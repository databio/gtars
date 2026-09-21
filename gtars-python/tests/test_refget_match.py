"""Binding tests for content-based collection name matching."""

import pytest

from gtars.refget import RefgetStore


def _collection(store, path, body):
    path.write_text(body)
    metadata, _ = store.add_sequence_collection_from_fasta(str(path))
    return metadata.digest


def test_match_sequence_names_shape_order_and_readonly(tmp_path):
    store = RefgetStore.in_memory()
    digest_a = _collection(
        store,
        tmp_path / "a.fa",
        ">chr1\nACGT\n>copy\nACGT\n>chr2\nGGCC\n>chr4\nCCCC\n",
    )
    digest_b = _collection(
        store,
        tmp_path / "b.fa",
        ">four\nCCCC\n>one\nACGT\n>one_alt\nACGT\n>extra\nTTAA\n",
    )

    expected = store.match_sequence_names(digest_a, digest_b)
    assert expected["collection_a"] == digest_a
    assert expected["collection_b"] == digest_b
    assert [row["names_a"] for row in expected["matches"]] == [
        ["chr1", "copy"],
        ["chr4"],
    ]
    assert set(expected["matches"][0]) == {"digest", "length", "names_a", "names_b"}
    assert not expected["matches"][0]["digest"].startswith("SQ.")
    assert expected["matches"][0]["names_a"] == ["chr1", "copy"]
    assert expected["matches"][0]["names_b"] == ["one", "one_alt"]
    assert expected["a_only"][0]["names_a"] == ["chr2"]
    assert expected["b_only"][0]["names_b"] == ["extra"]

    readonly = store.into_readonly()
    assert readonly.match_sequence_names(digest_a, digest_b) == expected


def test_match_sequence_names_rejects_unknown_digest():
    store = RefgetStore.in_memory()
    with pytest.raises(KeyError, match="missing-a"):
        store.match_sequence_names("missing-a", "missing-b")

    readonly = store.into_readonly()
    with pytest.raises(KeyError, match="missing-a"):
        readonly.match_sequence_names("missing-a", "missing-b")
