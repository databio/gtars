"""Tests for gtars.reftx PyO3 bindings."""

import pytest

from gtars.refget import sha512t24u_digest
from gtars.reftx import (
    CoordinateMapper,
    Exon,
    ManeStatus,
    MappingError,
    ReadonlyTxStore,
    Strand,
    Transcript,
    TxStoreBuilder,
    TxStoreError,
)


def _chrom_accession(seq: bytes = b"A" * 100) -> str:
    """Create a refget SQ.<digest> accession from a sequence."""
    return "SQ." + sha512t24u_digest(seq)


def _build_simple_store(tmp_path, chrom_acc=None, mane=False, plus_chrom=True):
    """Build a tiny store with a single forward-strand transcript and return path."""
    chrom_acc = chrom_acc or _chrom_accession()
    builder = TxStoreBuilder()
    tx_dict = {
        "accession": "NM_TEST.1",
        "gene": "TEST",
        "chrom": chrom_acc,
        "strand": "+" if plus_chrom else "-",
        "cds_start": 100,
        "cds_end": 200,
        "exons": [(50, 250)],
        "mane": {"select": True, "plus_clinical": False} if mane else None,
    }
    builder.add_transcript(tx_dict)
    out = str(tmp_path / "test.reftx")
    builder.build(out)
    return out, chrom_acc


def test_builder_roundtrip(tmp_path):
    chrom_acc = _chrom_accession()
    builder = TxStoreBuilder()
    builder.add_transcript({
        "accession": "NM_A.1",
        "gene": "A",
        "chrom": chrom_acc,
        "strand": "+",
        "cds_start": 100,
        "cds_end": 200,
        "exons": [(50, 250)],
    })
    builder.add_transcript({
        "accession": "NM_B.1",
        "gene": "B",
        "chrom": chrom_acc,
        "strand": "-",
        "cds_start": 1000,
        "cds_end": 1100,
        "exons": [(900, 1200)],
    })
    out = str(tmp_path / "two.reftx")
    builder.build(out)

    store = ReadonlyTxStore.open(out)
    assert len(store) == 2
    a = store.lookup("NM_A.1")
    assert a is not None
    assert a.accession == "NM_A.1"
    assert a.gene == "A"
    assert a.strand == Strand.Plus
    b = store.lookup("NM_B.1")
    assert b.strand == Strand.Minus
    assert store.lookup("NM_NOPE.1") is None


def test_mane_index(tmp_path):
    out, _ = _build_simple_store(tmp_path, mane=True)
    store = ReadonlyTxStore.open(out)
    assert store.has_mane_index()
    tx = store.lookup_mane("TEST")
    assert tx is not None
    assert tx.accession == "NM_TEST.1"
    assert tx.mane is not None
    assert tx.mane.select


def test_mane_missing_returns_none(tmp_path):
    out, _ = _build_simple_store(tmp_path, mane=False)
    store = ReadonlyTxStore.open(out)
    assert not store.has_mane_index()
    assert store.lookup_mane("TEST") is None


def test_coordinate_mapper_c_to_g(tmp_path):
    out, _ = _build_simple_store(tmp_path)
    store = ReadonlyTxStore.open(out)
    mapper = CoordinateMapper(store)

    # CDS starts at genomic 100; c.1 → 100.
    assert mapper.c_to_g("NM_TEST.1", 1) == 100
    assert mapper.c_to_g("NM_TEST.1", 50) == 149


def test_coordinate_mapper_full_metadata(tmp_path):
    out, chrom_acc = _build_simple_store(tmp_path)
    store = ReadonlyTxStore.open(out)
    mapper = CoordinateMapper(store)

    full = mapper.c_to_g_full("NM_TEST.1", 1)
    assert set(full.keys()) >= {"chrom", "chrom_accession", "genomic_pos", "strand"}
    assert full["chrom_accession"] == chrom_acc
    assert full["genomic_pos"] == 100
    assert full["strand"] == Strand.Plus


def test_c_to_g_by_gene_uses_mane(tmp_path):
    out, _ = _build_simple_store(tmp_path, mane=True)
    store = ReadonlyTxStore.open(out)
    mapper = CoordinateMapper(store)
    full = mapper.c_to_g_by_gene("TEST", 1)
    assert full["accession"] == "NM_TEST.1"
    assert full["genomic_pos"] == 100


def test_mapping_error_propagates(tmp_path):
    out, _ = _build_simple_store(tmp_path)
    store = ReadonlyTxStore.open(out)
    mapper = CoordinateMapper(store)
    with pytest.raises(MappingError):
        # Non-existent transcript
        mapper.c_to_g("NM_NOPE.1", 1)


def test_strand_from_str():
    assert Strand.from_str("+") == Strand.Plus
    assert Strand.from_str("-") == Strand.Minus
    with pytest.raises(ValueError):
        Strand.from_str("?")


def test_transcript_to_dict(tmp_path):
    out, chrom_acc = _build_simple_store(tmp_path, mane=True)
    store = ReadonlyTxStore.open(out)
    tx = store.lookup("NM_TEST.1")
    d = tx.to_dict()
    assert d["accession"] == "NM_TEST.1"
    assert d["chrom"] == chrom_acc
    assert d["strand"] == "+"
    assert d["cds_start"] == 100
    assert isinstance(d["exons"], list)
    assert d["mane"] is not None


def test_exon_to_dict():
    e = Exon(start=10, end=20)
    d = e.to_dict()
    assert d == {"start": 10, "end": 20, "cds_start": None, "cds_end": None}


def test_mane_status_to_dict():
    m = ManeStatus(select=True, plus_clinical=False)
    assert m.to_dict() == {"select": True, "plus_clinical": False}


def test_invalid_chrom_accession_raises(tmp_path):
    builder = TxStoreBuilder()
    with pytest.raises(TxStoreError):
        builder.add_transcript({
            "accession": "NM_X.1",
            "chrom": "SQ.notbase64!!!",
            "strand": "+",
            "exons": [(0, 10)],
        })
