"""Tests for gtars.vrs.hgvs PyO3 bindings."""

import json

import pytest

from gtars.vrs.hgvs import (
    Datum,
    Edit,
    HgvsError,
    HgvsVariant,
    PosEdit,
    Position,
    ReferenceType,
    parse_hgvs,
)


def test_parse_substitution():
    v = parse_hgvs("NM_000546.6:c.215C>G")
    assert v.accession == "NM_000546.6"
    assert v.reference_type == ReferenceType.Coding
    assert v.pos_edit.edit.kind == "substitution"
    assert v.pos_edit.edit.ref == "C"
    assert v.pos_edit.edit.alt == "G"
    assert v.pos_edit.start.base == 215
    assert v.pos_edit.start.offset == 0
    assert v.pos_edit.end is None


def test_parse_genomic_deletion():
    v = parse_hgvs("NC_000017.11:g.43124027del")
    assert v.reference_type == ReferenceType.Genomic
    assert v.pos_edit.edit.kind == "deletion"


def test_parse_insertion():
    v = parse_hgvs("NM_000546.6:c.215_216insATG")
    assert v.pos_edit.edit.kind == "insertion"
    assert v.pos_edit.edit.alt == "ATG"
    assert v.pos_edit.end is not None


def test_parse_delins():
    v = parse_hgvs("NM_000546.6:c.215delinsATGC")
    assert v.pos_edit.edit.kind == "delins"
    assert v.pos_edit.edit.alt == "ATGC"


def test_parse_dup():
    v = parse_hgvs("NM_000546.6:c.215dup")
    assert v.pos_edit.edit.kind == "duplication"


def test_parse_invalid_raises_hgvs_error():
    with pytest.raises(HgvsError):
        parse_hgvs("not a valid hgvs string")


def test_to_dict_roundtrip():
    v = parse_hgvs("NM_000546.6:c.215C>G")
    d = v.to_dict()
    # Must be JSON-serializable.
    serialized = json.dumps(d)
    parsed = json.loads(serialized)
    assert parsed["accession"] == "NM_000546.6"
    assert parsed["reference_type"] == "c"
    assert parsed["pos_edit"]["edit"]["kind"] == "substitution"
    assert parsed["pos_edit"]["edit"]["ref"] == "C"
    assert parsed["pos_edit"]["edit"]["alt"] == "G"
    assert parsed["pos_edit"]["start"]["base"] == 215


def test_intronic_offset():
    v = parse_hgvs("NM_000546.6:c.215+5C>G")
    assert v.pos_edit.start.base == 215
    assert v.pos_edit.start.offset == 5


def test_three_prime_utr():
    v = parse_hgvs("NM_000546.6:c.*37C>G")
    assert v.pos_edit.start.datum == Datum.CdsStop
    assert v.pos_edit.start.base == 37


def test_five_prime_utr():
    v = parse_hgvs("NM_000546.6:c.-14C>G")
    assert v.pos_edit.start.base == -14
    assert v.pos_edit.start.datum == Datum.Cds
