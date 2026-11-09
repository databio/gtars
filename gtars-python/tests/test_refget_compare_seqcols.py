"""Binding tests for compare_seqcols, the store-free level-2 comparison."""

import json
from pathlib import Path

import pytest

from gtars.refget import compare_seqcols

DATA = Path(__file__).parents[2] / "gtars-refget" / "tests" / "data" / "comparison"


def _load(path):
    return json.loads(path.read_text())


@pytest.mark.parametrize(
    "fixture", sorted(p.name for p in (DATA / "expected").glob("compare_*.json"))
)
def test_compare_seqcols_matches_fixture(fixture):
    a_name, b_name = fixture[len("compare_") : -len(".json")].split(".fa_")
    a = _load(DATA / "collections" / f"{a_name}.fa.json")
    b = _load(DATA / "collections" / f"{b_name}.json")
    expected = _load(DATA / "expected" / fixture)
    result = compare_seqcols(a, b, expected["digests"]["a"], expected["digests"]["b"])
    assert result == expected


def test_compare_seqcols_without_digests_and_objects():
    a = {"name_length_pairs": [{"name": "chr1", "length": 10}, {"name": "chr2", "length": 20}]}
    b = {"name_length_pairs": [{"length": 10, "name": "chr1"}, {"length": 20, "name": "chr2"}]}
    result = compare_seqcols(a, b)
    assert result["digests"] == {"a": None, "b": None}
    assert result["array_elements"]["a_and_b_same_order"] == {"name_length_pairs": True}


def test_compare_seqcols_one_shared_element_is_none():
    result = compare_seqcols({"names": ["chr1", "chr2"]}, {"names": ["chr1", "chrM"]})
    assert result["array_elements"]["a_and_b_count"] == {"names": 1}
    assert result["array_elements"]["a_and_b_same_order"] == {"names": None}


def test_compare_seqcols_rejects_non_list():
    with pytest.raises(ValueError, match="names"):
        compare_seqcols({"names": "chr1"}, {})
