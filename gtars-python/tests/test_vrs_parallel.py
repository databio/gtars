import gzip
import random
import struct
import zlib

import pytest

from gtars.refget import RefgetStore


def _bgzf_block(chunk: bytes) -> bytes:
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    cdata = c.compress(chunk) + c.flush()
    bsize = 18 + len(cdata) + 8 - 1
    header = struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, bsize)
    return header + cdata + struct.pack("<II", zlib.crc32(chunk), len(chunk))


def _write_bgzf(path, data: bytes, block_size: int = 997):
    """Small blocks so VCF lines straddle BGZF block boundaries."""
    with open(path, "wb") as f:
        for i in range(0, len(data), block_size):
            f.write(_bgzf_block(data[i : i + block_size]))
        f.write(_bgzf_block(b""))


@pytest.fixture(scope="module")
def vrs_inputs(tmp_path_factory):
    d = tmp_path_factory.mktemp("vrs_parallel")
    rng = random.Random(42)
    seqs = {c: "".join(rng.choice("ACGT") for _ in range(3000)) for c in ("chr1", "chr2")}
    fasta = d / "ref.fa"
    fasta.write_text("".join(f">{c}\n{s}\n" for c, s in seqs.items()))

    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for chrom, seq in seqs.items():
        for pos in range(10, 2900, 7):
            ref = seq[pos - 1 : pos - 1 + rng.randint(1, 4)]
            alts = [rng.choice("ACGT"), ref[0], ref[0] + "AC", "*"]
            alts = [a for a in rng.sample(alts, rng.randint(1, 3)) if a != ref] or ["N"]
            lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{','.join(alts)}\t.\t.\tAC=1")
    lines.append("chrUn\t5\t.\tA\tG\t.\t.\t.")
    text = ("\n".join(lines) + "\n").encode()

    plain = d / "in.vcf"
    plain.write_bytes(text)
    gz = d / "in.vcf.gz"
    gz.write_bytes(gzip.compress(text))
    bgz = d / "in.vcf.bgz"
    _write_bgzf(bgz, text)
    return d, fasta, {"plain": plain, "gzip": gz, "bgzf": bgz}


def _read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in f]
    return header, rows


@pytest.mark.parametrize("fmt", ["plain", "gzip", "bgzf"])
@pytest.mark.parametrize("threads", [1, 4, None])
def test_parallel_tsv_matches_serial(vrs_inputs, tmp_path, fmt, threads):
    _, fasta, vcfs = vrs_inputs
    store = RefgetStore.in_memory()
    meta, _ = store.add_sequence_collection_from_fasta(str(fasta))
    serial = store.compute_vrs_ids(meta.digest, str(vcfs["plain"]))
    assert len(serial) > 500

    out = tmp_path / "out.tsv"
    n = store.compute_vrs_ids_parallel(meta.digest, str(vcfs[fmt]), str(out), threads=threads)

    header, rows = _read_tsv(out)
    assert header == ["chrom", "pos", "ref", "alt", "vrs_id"]
    assert n == len(rows) == len(serial)
    expected = [[r["chrom"], str(r["pos"]), r["ref"], r["alt"], r["vrs_id"]] for r in serial]
    assert rows == expected


def test_parallel_tsv_lazy_disk_store(vrs_inputs, tmp_path):
    """A freshly opened local store has only stubs; sequences load on demand."""
    d, fasta, vcfs = vrs_inputs
    store_dir = tmp_path / "store"
    writer = RefgetStore.on_disk(str(store_dir))
    meta, _ = writer.add_sequence_collection_from_fasta(str(fasta))
    serial = writer.compute_vrs_ids(meta.digest, str(vcfs["plain"]))

    store = RefgetStore.open_local(str(store_dir))
    out = tmp_path / "out.tsv"
    n = store.compute_vrs_ids_parallel(meta.digest, str(vcfs["bgzf"]), str(out))
    _, rows = _read_tsv(out)
    assert n == len(serial)
    assert rows == [[r["chrom"], str(r["pos"]), r["ref"], r["alt"], r["vrs_id"]] for r in serial]
