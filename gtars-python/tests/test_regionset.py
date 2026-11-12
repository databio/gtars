import os
from pathlib import Path

import pytest

from gtars.models import RegionSet

DATA_DIR = Path(__file__).resolve().parent.parent.parent / "tests" / "data" / "regionset"


class TestRegionSet:
    # TODO: re-enable when HTTP loading is reliable in CI
    # @pytest.mark.parametrize(
    #     "bed_file",
    #     [
    #         "https://raw.githubusercontent.com/databio/gtars/refs/heads/master/tests/data/regionset/dummy.narrowPeak",
    #     ],
    # )
    # def test_mean_region_width(self, bed_file):
    #     rs = RegionSet(bed_file)
    #
    #     assert rs.mean_region_width() == 4.22

    def test_from_regions(self):
        from gtars.models import Region

        regions = [
            Region(chr="chr1", start=14, end=514, rest=None),
            Region(chr="chr19", start=19, end=810, rest=None),
        ]
        rs = RegionSet.from_regions(regions)

        assert isinstance(rs, RegionSet)
        assert len(rs) == 2


def _rs(*specs):
    """Helper: _rs(("chr1", 100, 200), ("chr1", 300, 400))"""
    from gtars.models import Region

    return RegionSet.from_regions([Region(s[0], s[1], s[2], None) for s in specs])


class TestOverlapOps:
    def setup_method(self):
        self.a = _rs(("chr1", 100, 200), ("chr1", 300, 400), ("chr1", 500, 600))
        self.b = _rs(("chr1", 150, 250), ("chr1", 550, 650))

    def test_subset_by_overlaps(self):
        result = self.a.subset_by_overlaps(self.b)
        assert len(result) == 2

    def test_count_overlaps(self):
        assert self.a.count_overlaps(self.b) == [1, 0, 1]

    def test_any_overlaps(self):
        assert self.a.any_overlaps(self.b) == [True, False, True]

    def test_find_overlaps(self):
        result = self.a.find_overlaps(self.b)
        assert result == [[0], [], [1]]

    def test_intersect_all(self):
        result = self.a.intersect_all(self.b)
        assert len(result) == 2
        assert result[0].start == 150
        assert result[0].end == 200

    def test_setdiff(self):
        result = self.a.setdiff(self.b)
        assert len(result) > 0

    def test_closest(self):
        result = self.a.closest(self.b)
        assert len(result) == 3

    def test_cluster(self):
        assert self.a.cluster(0) == [0, 1, 2]
        assert self.a.cluster(200) == [0, 0, 0]

    def test_coverage(self):
        assert self.a.coverage(self.a) == 1.0
        empty = _rs()
        assert self.a.coverage(empty) == 0.0

    def test_overlap_coefficient(self):
        assert self.a.overlap_coefficient(self.a) == 1.0
        assert self.a.overlap_coefficient(self.b) == self.b.overlap_coefficient(self.a)


class TestDisjoin:
    def test_disjoin_breaks_overlaps(self):
        """Disjoin should split overlapping regions at every boundary."""
        rs = _rs(("chr1", 0, 100), ("chr1", 50, 150))
        result = rs.disjoin()
        # Two overlapping regions [0,100) and [50,150) produce three disjoint pieces:
        # [0,50), [50,100), [100,150)
        assert len(result) == 3

    def test_disjoin_nonoverlapping_unchanged(self):
        """Non-overlapping regions should pass through unchanged."""
        rs = _rs(("chr1", 0, 100), ("chr1", 200, 300))
        result = rs.disjoin()
        assert len(result) == 2

    def test_disjoin_empty(self):
        rs = _rs()
        result = rs.disjoin()
        assert len(result) == 0


class TestGaps:
    def test_gaps_basic(self):
        """Gaps should emit leading, inter-region, and trailing gaps."""
        rs = _rs(("chr1", 10, 20), ("chr1", 30, 40), ("chr1", 50, 60))
        result = rs.gaps({"chr1": 100})
        coords = [(r.start, r.end) for r in result]
        assert coords == [(0, 10), (20, 30), (40, 50), (60, 100)]

    def test_gaps_full_chrom_for_empty(self):
        """A chromosome in chrom_sizes with no regions yields a full-chrom gap."""
        rs = _rs(("chr1", 10, 20))
        result = rs.gaps({"chr1": 100, "chr2": 50})
        chr2 = [(r.start, r.end) for r in result if r.chr == "chr2"]
        assert chr2 == [(0, 50)]

    def test_gaps_empty(self):
        rs = _rs()
        result = rs.gaps({})
        assert len(result) == 0


class TestToBigBed:
    BED = DATA_DIR / "dummy.narrowPeak"
    CHROM_SIZES_FILE = DATA_DIR / "dummy_chrom_sizes"

    @staticmethod
    def _sizes_dict(path):
        sizes = {}
        with open(path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 2:
                    sizes[parts[0]] = int(parts[1])
        return sizes

    def test_to_bigbed_from_path(self, tmp_path):
        out = tmp_path / "from_path.bigbed"
        RegionSet(str(self.BED)).to_bigbed(str(out), str(self.CHROM_SIZES_FILE))
        assert out.stat().st_size > 0

    def test_to_bigbed_from_dict(self, tmp_path):
        out = tmp_path / "from_dict.bigbed"
        RegionSet(str(self.BED)).to_bigbed(str(out), self._sizes_dict(self.CHROM_SIZES_FILE))
        assert out.stat().st_size > 0

    def test_to_bigbed_from_pathlib(self, tmp_path):
        """A pathlib.Path (os.PathLike) is accepted as the chrom sizes argument."""
        out = tmp_path / "from_pathlib.bigbed"
        RegionSet(str(self.BED)).to_bigbed(str(out), self.CHROM_SIZES_FILE)
        assert out.stat().st_size > 0

    def test_to_bigbed_path_and_dict_match(self, tmp_path):
        """File-path and dict inputs for the same sizes produce identical output."""
        rs = RegionSet(str(self.BED))
        p_path = tmp_path / "path.bigbed"
        p_dict = tmp_path / "dict.bigbed"
        rs.to_bigbed(str(p_path), str(self.CHROM_SIZES_FILE))
        rs.to_bigbed(str(p_dict), self._sizes_dict(self.CHROM_SIZES_FILE))
        assert p_path.read_bytes() == p_dict.read_bytes()

    def test_to_bigbed_missing_file_raises(self, tmp_path):
        """A bad chrom sizes path raises a clean exception, not a panic."""
        with pytest.raises(Exception):
            RegionSet(str(self.BED)).to_bigbed(
                str(tmp_path / "x.bigbed"), str(tmp_path / "does_not_exist.sizes")
            )

    def test_to_bigbed_bad_dict_value_raises(self, tmp_path):
        with pytest.raises(Exception):
            RegionSet(str(self.BED)).to_bigbed(
                str(tmp_path / "y.bigbed"), {"chr1": "not-an-int"}
            )


