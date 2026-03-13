import os
from pathlib import Path

import pytest

from gtars.models import RegionSet


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

    def test_subtract(self):
        result = self.a.subtract(self.b)
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
