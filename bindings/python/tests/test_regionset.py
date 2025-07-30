import os
from pathlib import Path

import pytest

from gtars.models import RegionSet


class TestRegionSet:
    @pytest.mark.parametrize(
        "bed_file",
        [
            "https://raw.githubusercontent.com/databio/gtars/refs/heads/master/gtars/tests/data/regionset/dummy.narrowPeak",
        ],
    )
    def test_mean_region_width(self, bed_file):
        rs = RegionSet(bed_file)

        assert rs.mean_region_width() == 4.22

    @pytest
    def test_from_regions(self):
        from gtars.models import Region

        regions = [
            Region(chrom="chr1", start=14, end=514),
            Region(chrom="chr19", start=19, end=810),
        ]
        rs = RegionSet.from_regions(regions)

        assert isinstance(rs, RegionSet)
        assert len(rs) == 2
