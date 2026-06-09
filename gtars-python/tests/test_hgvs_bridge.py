"""Tests for the gtars.vrs.hgvs.hgvs_to_vrs_id bridge.

End-to-end fixture setup (refget collection + transcript store sharing a
chromosome digest) is involved. The full positive-path test is tracked as a
follow-up; for now verify the binding is wired and that error paths surface
HgvsError.
"""

import pytest

from gtars.refget import RefgetStore, sha512t24u_digest
from gtars.reftx import ReadonlyTxStore, ReftxProvider, TxStoreBuilder
from gtars.vrs.hgvs import HgvsError, hgvs_to_vrs_id


@pytest.fixture
def tiny_provider(tmp_path):
    chrom_acc = "SQ." + sha512t24u_digest(b"A" * 200)
    builder = TxStoreBuilder()
    builder.add_transcript({
        "accession": "NM_TEST.1",
        "gene": "TEST",
        "chrom": chrom_acc,
        "strand": "+",
        "cds_start": 100,
        "cds_end": 150,
        "exons": [(50, 200)],
    })
    out = str(tmp_path / "tiny.reftx")
    builder.build(out)
    store = ReadonlyTxStore.open(out)
    return ReftxProvider(store), chrom_acc


def test_unknown_collection_raises(tiny_provider):
    provider, _ = tiny_provider
    refget = RefgetStore.in_memory()
    with pytest.raises(HgvsError):
        hgvs_to_vrs_id(
            "NM_TEST.1:c.5A>T",
            provider,
            refget,
            "no_such_collection_digest",
        )


def test_invalid_hgvs_raises(tiny_provider):
    provider, _ = tiny_provider
    refget = RefgetStore.in_memory()
    with pytest.raises(HgvsError):
        hgvs_to_vrs_id(
            "definitely not hgvs",
            provider,
            refget,
            "any",
        )
