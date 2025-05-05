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