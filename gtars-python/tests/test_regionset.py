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


class TestSpatialStats:
    def test_inter_peak_spacing_regular_array(self):
        """Evenly-spaced peaks produce zero std / zero IQR."""
        rs = _rs(
            ("chr1", 0, 10),
            ("chr1", 20, 30),
            ("chr1", 40, 50),
            ("chr1", 60, 70),
        )
        s = rs.inter_peak_spacing()
        assert s.n_gaps == 3
        assert s.mean == 10.0
        assert s.median == 10.0
        assert s.std == 0.0
        assert s.iqr == 0.0

    def test_inter_peak_spacing_empty(self):
        s = _rs().inter_peak_spacing()
        assert s.n_gaps == 0
        import math
        assert math.isnan(s.mean)
        assert math.isnan(s.std)

    def test_peak_clusters_mixed_default(self):
        """Two clusters + one singleton at radius 5. Default
        min_cluster_size=2 excludes the singleton uniformly: n_clusters,
        n_clustered_peaks, mean_cluster_size, and fraction_clustered
        all describe multi-peak clusters only."""
        rs = _rs(
            ("chr1", 0, 10),
            ("chr1", 13, 20),
            ("chr1", 100, 110),
            ("chr1", 113, 120),
            ("chr1", 122, 130),
            ("chr1", 500, 510),
        )
        c = rs.peak_clusters(5)
        assert c.n_clusters == 2  # multi-peak only, singleton excluded
        assert c.n_clustered_peaks == 5  # peaks in the two multi-peak clusters
        assert c.max_cluster_size == 3  # unfiltered, still the biggest cluster
        assert c.mean_cluster_size == 2.5  # (2 + 3) / 2
        # Arithmetic identity holds: 2 * 2.5 == 5
        assert c.n_clusters * c.mean_cluster_size == c.n_clustered_peaks
        # fraction_clustered uses raw total (6), not filtered count.
        assert abs(c.fraction_clustered - 5 / 6) < 1e-10

    def test_peak_clusters_mixed_min_1_simple_average(self):
        """Same fixture, min_cluster_size=1 includes every connected
        component including the singleton. Mean degenerates to
        total_peaks / n_clusters = 6 / 3 = 2.0."""
        rs = _rs(
            ("chr1", 0, 10),
            ("chr1", 13, 20),
            ("chr1", 100, 110),
            ("chr1", 113, 120),
            ("chr1", 122, 130),
            ("chr1", 500, 510),
        )
        c = rs.peak_clusters(5, min_cluster_size=1)
        assert c.n_clusters == 3  # includes the singleton
        assert c.n_clustered_peaks == 6  # == total_peaks under min=1
        assert c.max_cluster_size == 3
        assert c.mean_cluster_size == 2.0
        # Identity still holds: 3 * 2.0 == 6
        assert c.n_clusters * c.mean_cluster_size == c.n_clustered_peaks
        assert c.fraction_clustered == 1.0  # tautological under min=1

    def test_peak_clusters_empty(self):
        c = _rs().peak_clusters(100)
        assert c.n_clusters == 0
        assert c.max_cluster_size == 0
        assert c.radius_bp == 100

    def test_density_vector_single_chrom(self):
        """5 peaks in 5 consecutive bins → counts [1,1,1,1,1]."""
        rs = _rs(
            ("chr1", 5, 15),
            ("chr1", 25, 35),
            ("chr1", 45, 55),
            ("chr1", 65, 75),
            ("chr1", 85, 95),
        )
        dv = rs.density_vector({"chr1": 100}, 5)
        assert list(dv.counts) == [1, 1, 1, 1, 1]
        assert dv.bin_width == 20
        assert len(dv) == 5
        assert dv.chrom_offsets == [("chr1", 0)]

    def test_density_vector_zero_padding(self):
        """Empty RegionSet produces a zero-padded vector."""
        dv = _rs().density_vector({"chr1": 100}, 5)
        assert list(dv.counts) == [0, 0, 0, 0, 0]

    def test_density_homogeneity_even_distribution(self):
        """Perfectly even distribution → CV = 0, Gini = 0."""
        rs = _rs(
            ("chr1", 5, 15),
            ("chr1", 25, 35),
            ("chr1", 45, 55),
            ("chr1", 65, 75),
            ("chr1", 85, 95),
        )
        h = rs.density_homogeneity({"chr1": 100}, 5)
        assert h.n_windows == 5
        assert h.n_nonzero_windows == 5
        assert h.mean_count == 1.0
        assert h.cv == 0.0
        assert h.gini == 0.0

    def test_density_homogeneity_empty_regionset(self):
        """Empty RegionSet over populated chrom_sizes: Gini=0, CV=NaN."""
        import math

        h = _rs().density_homogeneity({"chr1": 100}, 5)
        assert h.n_windows == 5
        assert h.n_nonzero_windows == 0
        assert h.mean_count == 0.0
        assert h.gini == 0.0
        assert math.isnan(h.cv)
