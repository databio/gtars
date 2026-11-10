"""Binding tests for export_fasta's keyword signature, labels and errors."""

import pytest

from gtars.refget import RefgetStore


def _collection(store, path, body):
    path.write_text(body)
    metadata, _ = store.add_sequence_collection_from_fasta(str(path))
    return metadata.digest


def _headers(path):
    return [l[1:] for l in path.read_text().splitlines() if l.startswith(">")]


def test_export_fasta_collection_default_label_is_digest(tmp_path):
    store = RefgetStore.in_memory()
    digest = _collection(store, tmp_path / "a.fa", ">chr1\nACGT\n>chr2\nGGCC\n")

    out_path = tmp_path / "out_default.fa"
    summary = store.export_fasta(str(out_path), collections=digest)

    assert summary["n_records"] == 2
    assert summary["n_bases"] == 8
    assert all(h.startswith("SQ.") for h in _headers(out_path))


def test_export_fasta_names_and_label(tmp_path):
    store = RefgetStore.in_memory()
    digest = _collection(store, tmp_path / "a.fa", ">chr1\nACGT\n>chr2\nGGCC\n")

    out_path = tmp_path / "out_names.fa"
    store.export_fasta(out_path, collections=[digest], names=["chr1"], label="name", line_width=2)

    assert out_path.read_text() == ">chr1\nAC\nGT\n"


def test_export_fasta_seq_digests_and_all(tmp_path):
    store = RefgetStore.in_memory()
    digest = _collection(store, tmp_path / "a.fa", ">chr1\nACGT\n")
    _collection(store, tmp_path / "b.fa", ">chr1\nACGT\n>chrB\nTTTT\n")
    seq_digest = store.get_collection(digest).sequences[0].metadata.sha512t24u

    out_path = tmp_path / "out_by_digest.fa"
    store.export_fasta(str(out_path), seq_digests=[seq_digest])
    assert _headers(out_path) == [f"SQ.{seq_digest}"]

    summary = store.export_fasta(str(out_path), all=True)
    assert summary["n_records"] == 2

    summary = store.export_fasta(str(out_path), all=True, label="{col}/{name}")
    assert summary["n_records"] == 3

    summary = store.export_fasta(str(out_path), all=True, label="name")
    assert summary["n_deduplicated"] == 1


def test_export_fasta_bad_inputs_raise_value_error(tmp_path):
    store = RefgetStore.in_memory()
    digest = _collection(store, tmp_path / "a.fa", ">chr1\nACGT\n>chr1\nGGCC\n")
    out = str(tmp_path / "x.fa")

    with pytest.raises(ValueError, match="exactly one"):
        store.export_fasta(out)
    with pytest.raises(ValueError, match="exactly one"):
        store.export_fasta(out, collections=digest, all=True)
    with pytest.raises(ValueError, match="unknown field"):
        store.export_fasta(out, collections=digest, label="{nope}")
    with pytest.raises(ValueError, match="namespace"):
        store.export_fasta(out, collections=digest, label="seq:nope")
    with pytest.raises(ValueError, match="molecule"):
        store.export_fasta(out, collections=digest, molecule="rna")
    with pytest.raises(ValueError, match="two different sequences"):
        store.export_fasta(out, collections=digest, label="name")


def test_export_fasta_from_regions_label(tmp_path):
    store = RefgetStore.in_memory()
    digest = _collection(store, tmp_path / "a.fa", ">chr1\nACGTACGT\n")
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t2\t6\n")
    out = tmp_path / "r.fa"

    summary = store.export_fasta_from_regions(digest, str(bed), str(out), label="name")

    assert out.read_text() == ">chr1:2-6\nGTAC\n"
    assert summary["n_records"] == 1
