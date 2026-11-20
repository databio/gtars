"""The packed RefgetStore layout (alpha) and the standalone file limit."""

import pytest

from gtars.refget import RefgetStore

SEQS = {
    "chr1": "ACGTACGTACGTACGTAAAA",
    "chr2": "GGGGCCCCAAAATTTTACGT",
}


def test_packed_store_roundtrip(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in SEQS.items()))
    store_dir = tmp_path / "store"

    store = RefgetStore.on_disk(store_dir, packed=True)
    store.set_quiet(True)
    assert store.layout == "packed"
    store.add_sequence_collection_from_fasta(fasta)
    del store

    assert not list(store_dir.rglob("*.seq")), "a packed store has no .seq files"
    assert list((store_dir / "packs").glob("*.pack")), "a packed store has pack files"

    reopened = RefgetStore.open_local(store_dir)
    assert reopened.layout == "packed"
    by_name = {m.name: m for m in reopened.list_sequences()}
    assert sorted(by_name) == sorted(SEQS)
    for name, bases in SEQS.items():
        meta = by_name[name]
        assert reopened.get_substring(meta.sha512t24u, 0, meta.length) == bases

    report = reopened.compact()
    assert set(report) == {"packs_before", "packs_after", "bytes_reclaimed", "packs_pruned"}
    assert report["packs_pruned"] == 0
    for name, bases in SEQS.items():
        meta = by_name[name]
        assert reopened.get_substring(meta.sha512t24u, 0, meta.length) == bases


def test_compact_prune_unreferenced(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in SEQS.items()))
    store_dir = tmp_path / "store"
    store = RefgetStore.on_disk(store_dir, packed=True)
    store.set_quiet(True)
    store.add_sequence_collection_from_fasta(fasta)
    del store

    packs = list((store_dir / "packs").glob("*.pack"))
    assert len(packs) == 1
    unlisted = store_dir / "packs" / "0123456789abcdef0123456789abcdef.pack"
    unlisted.write_bytes(packs[0].read_bytes())

    store = RefgetStore.open_local(store_dir)
    assert store.compact()["packs_pruned"] == 0
    assert unlisted.exists(), "compact() must keep a pack no index lists"

    report = store.compact(prune_unreferenced=True)
    assert report["packs_pruned"] == 1
    assert not unlisted.exists()
    by_name = {m.name: m for m in store.list_sequences()}
    for name, bases in SEQS.items():
        meta = by_name[name]
        assert store.get_substring(meta.sha512t24u, 0, meta.length) == bases


def test_set_pack_cache_reads_the_same_bases(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in SEQS.items()))
    store_dir = tmp_path / "store"
    store = RefgetStore.on_disk(store_dir, packed=True)
    store.set_quiet(True)
    store.add_sequence_collection_from_fasta(fasta)
    del store

    store = RefgetStore.open_local(store_dir)
    by_name = {m.name: m for m in store.list_sequences()}
    # On, off, then on again with a different size.
    for megabytes in [16, 0, 1]:
        store.set_pack_cache(megabytes)
        for name, bases in SEQS.items():
            meta = by_name[name]
            assert store.get_substring(meta.sha512t24u, 0, meta.length) == bases
            assert store.get_substring(meta.sha512t24u, 3, 9) == bases[3:9]


def test_standalone_store_reports_its_layout(tmp_path):
    assert RefgetStore.on_disk(tmp_path / "store").layout == "standalone"


def test_standalone_file_limit_is_honored(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in SEQS.items()))

    store = RefgetStore.on_disk(tmp_path / "store")
    store.set_quiet(True)
    # force does not skip the limit.
    with pytest.raises(IOError, match="standalone_file_limit"):
        store.add_sequence_collection_from_fasta(fasta, force=True, standalone_file_limit=1)
    with pytest.raises(IOError, match="standalone_file_limit"):
        store.add_sequence_collections_from_fastas([str(fasta)], standalone_file_limit=1)
    assert not list((tmp_path / "store").rglob("*.seq"))

    _, was_new = store.add_sequence_collection_from_fasta(fasta, standalone_file_limit=2)
    assert was_new

    # Packed stores ignore the limit.
    packed = RefgetStore.on_disk(tmp_path / "packed", packed=True)
    packed.set_quiet(True)
    packed.add_sequence_collection_from_fasta(fasta, standalone_file_limit=1)
