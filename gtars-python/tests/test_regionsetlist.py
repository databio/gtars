"""Tests for RegionSetList Python bindings."""

import os
import tempfile

import pytest

from gtars.models import Region, RegionSet
from gtars.models import RegionSetList


def _rs(*specs):
    """Helper: _rs(("chr1", 100, 200), ("chr1", 300, 400))"""
    return RegionSet.from_regions([Region(s[0], s[1], s[2], None) for s in specs])


# =========================================================================
# Constructor
# =========================================================================


class TestRegionSetListConstructor:
    def test_create_empty(self):
        rsl = RegionSetList([])
        assert len(rsl) == 0

    def test_create_from_region_sets(self):
        rs1 = _rs(("chr1", 0, 100), ("chr1", 200, 300))
        rs2 = _rs(("chr1", 50, 150))
        rsl = RegionSetList([rs1, rs2])
        assert len(rsl) == 2


# =========================================================================
# Indexing
# =========================================================================


class TestRegionSetListIndexing:
    def test_getitem(self):
        rs1 = _rs(("chr1", 0, 100))
        rs2 = _rs(("chr2", 0, 50))
        rsl = RegionSetList([rs1, rs2])

        first = rsl[0]
        assert isinstance(first, RegionSet)
        assert len(first) == 1

    def test_getitem_negative_index(self):
        rs1 = _rs(("chr1", 0, 100))
        rs2 = _rs(("chr2", 0, 50))
        rsl = RegionSetList([rs1, rs2])

        last = rsl[-1]
        assert isinstance(last, RegionSet)
        assert len(last) == 1

    def test_getitem_out_of_bounds(self):
        rsl = RegionSetList([_rs(("chr1", 0, 100))])
        with pytest.raises(IndexError):
            rsl[5]


# =========================================================================
# Iteration
# =========================================================================


class TestRegionSetListIter:
    def test_iter(self):
        rs1 = _rs(("chr1", 0, 100))
        rs2 = _rs(("chr2", 0, 50))
        rsl = RegionSetList([rs1, rs2])

        sets = list(rsl)
        assert len(sets) == 2
        assert all(isinstance(s, RegionSet) for s in sets)


# =========================================================================
# concat
# =========================================================================


class TestRegionSetListConcat:
    def test_concat_flattens(self):
        rs1 = _rs(("chr1", 0, 100), ("chr1", 200, 300))
        rs2 = _rs(("chr1", 50, 150))
        rsl = RegionSetList([rs1, rs2])

        combined = rsl.concat()
        assert isinstance(combined, RegionSet)
        assert len(combined) == 3  # no merging

    def test_concat_empty(self):
        rsl = RegionSetList([])
        combined = rsl.concat()
        assert len(combined) == 0

    def test_concat_preserves_duplicates(self):
        rs1 = _rs(("chr1", 0, 100))
        rs2 = _rs(("chr1", 0, 100))
        rsl = RegionSetList([rs1, rs2])

        combined = rsl.concat()
        assert len(combined) == 2  # raw union, no dedup


# =========================================================================
# names
# =========================================================================


class TestRegionSetListNames:
    def test_names_none_by_default(self):
        rsl = RegionSetList([_rs(("chr1", 0, 100))])
        assert rsl.names is None


# =========================================================================
# RegionDB integration
# =========================================================================


@pytest.fixture
def simple_db():
    """Create a simple RegionDB from temp BED files."""
    from gtars.lola import RegionDB

    tmpdir = tempfile.mkdtemp()
    bed1 = os.path.join(tmpdir, "a.bed")
    bed2 = os.path.join(tmpdir, "b.bed")
    bed3 = os.path.join(tmpdir, "c.bed")

    with open(bed1, "w") as f:
        f.write("chr1\t100\t200\nchr1\t500\t600\n")
    with open(bed2, "w") as f:
        f.write("chr1\t150\t250\nchr1\t700\t800\n")
    with open(bed3, "w") as f:
        f.write("chr2\t0\t1000\n")

    return RegionDB.from_bed_files([bed1, bed2, bed3])


class TestRegionSetListFromDB:
    def test_get_region_sets_by_index(self, simple_db):
        rsl = simple_db.get_region_sets([0, 1])
        assert isinstance(rsl, RegionSetList)
        assert len(rsl) == 2

    def test_get_region_sets_all(self, simple_db):
        rsl = simple_db.get_region_sets()
        assert len(rsl) == 3

    def test_universe_from_db(self, simple_db):
        """Build universe: get_region_sets -> concat."""
        rsl = simple_db.get_region_sets()
        combined = rsl.concat()
        assert isinstance(combined, RegionSet)
        # 2 + 2 + 1 = 5 regions total
        assert len(combined) == 5
