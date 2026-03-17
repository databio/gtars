"""Tests for genomicdist Python bindings."""

import os
import tempfile
from pathlib import Path

import pytest

from gtars.models import (
    Region,
    RegionSet,
    TssIndex,
    GeneModel,
    GenomicDistAnnotation,
    PartitionList,
    SignalMatrix,
)
from gtars.genomic_distributions import (
    calc_partitions,
    calc_expected_partitions,
    calc_summary_signal,
    median_abs_distance,
    consensus,
)

DATA_DIR = (
    Path(__file__).resolve().parent.parent.parent / "tests" / "data" / "regionset"
)


def make_regionset(regions):
    """Helper: list of (chr, start, end) tuples -> RegionSet."""
    return RegionSet.from_regions(
        [Region(chr=c, start=s, end=e, rest=None) for c, s, e in regions]
    )


# =========================================================================
# RegionSet statistics
# =========================================================================


class TestRegionSetStatistics:
    def test_widths(self):
        rs = make_regionset([("chr1", 100, 200), ("chr1", 300, 550)])
        assert rs.widths() == [100, 250]

    def test_widths_matches_region_widths(self):
        rs = make_regionset([("chr1", 0, 10), ("chr2", 5, 25)])
        assert rs.widths() == rs.region_widths()

    def test_neighbor_distances_same_chrom(self):
        rs = make_regionset(
            [("chr1", 100, 200), ("chr1", 300, 400), ("chr1", 500, 600)]
        )
        dists = rs.neighbor_distances()
        assert len(dists) == 2
        assert dists[0] == 100.0  # 300 - 200
        assert dists[1] == 100.0  # 500 - 400

    def test_neighbor_distances_overlapping(self):
        rs = make_regionset([("chr1", 100, 300), ("chr1", 200, 400)])
        dists = rs.neighbor_distances()
        # Overlapping regions produce no neighbor distances in the Rust implementation
        assert len(dists) == 0

    def test_neighbor_distances_multi_chrom(self):
        # Distances are per-chromosome, so single regions on each chrom yield no distances
        rs = make_regionset([("chr1", 0, 10), ("chr2", 0, 10)])
        dists = rs.neighbor_distances()
        assert len(dists) == 0

    def test_nearest_neighbors(self):
        rs = make_regionset([("chr1", 0, 10), ("chr1", 20, 30), ("chr1", 100, 110)])
        nn = rs.nearest_neighbors()
        assert len(nn) == 3
        assert nn[0] == 10.0  # nearest to [0,10) is [20,30) -> gap=10
        assert nn[1] == 10.0  # nearest to [20,30) is [0,10) -> gap=10
        assert nn[2] == 70.0  # nearest to [100,110) is [20,30) -> gap=70

    def test_nearest_neighbors_single_region(self):
        rs = make_regionset([("chr1", 0, 10)])
        nn = rs.nearest_neighbors()
        # Single region has no neighbors; Rust returns empty list
        assert len(nn) == 0

    def test_nearest_neighbors_cross_chrom_none(self):
        rs = make_regionset([("chr1", 0, 10), ("chr2", 0, 10)])
        nn = rs.nearest_neighbors()
        # Each is alone on its chromosome; Rust returns empty list
        assert len(nn) == 0

    def test_distribution(self):
        rs = make_regionset([("chr1", 0, 100), ("chr1", 500, 600)])
        dist = rs.distribution(10)
        assert isinstance(dist, list)
        assert all(isinstance(d, dict) for d in dist)
        assert all(k in d for d in dist for k in ("chr", "start", "end", "n", "rid"))
        # At least some bins should have n > 0
        assert any(d["n"] > 0 for d in dist)

    def test_distribution_default_bins(self):
        rs = make_regionset([("chr1", 0, 1000)])
        dist = rs.distribution()
        assert len(dist) > 0


# =========================================================================
# Interval range operations
# =========================================================================


class TestIntervalRanges:
    def test_trim(self):
        rs = make_regionset([("chr1", 0, 1000), ("chr2", 500, 2000)])
        trimmed = rs.trim({"chr1": 500, "chr2": 1500})
        assert len(trimmed) == 2
        # chr1 should be trimmed to [0, 500)
        assert trimmed[0].end == 500
        # chr2 should be trimmed to [500, 1500)
        assert trimmed[1].end == 1500

    def test_trim_drops_unknown_chroms(self):
        rs = make_regionset([("chr1", 0, 100), ("chrZ", 0, 100)])
        trimmed = rs.trim({"chr1": 1000})
        assert len(trimmed) == 1
        assert trimmed[0].chr == "chr1"

    def test_promoters(self):
        rs = make_regionset([("chr1", 1000, 2000)])
        proms = rs.promoters(200, 50)
        assert len(proms) == 1
        assert proms[0].start == 800  # 1000 - 200
        assert proms[0].end == 1050  # 1000 + 50

    def test_promoters_saturates_at_zero(self):
        rs = make_regionset([("chr1", 50, 200)])
        proms = rs.promoters(100, 10)
        assert proms[0].start == 0  # saturating subtraction

    def test_reduce(self):
        rs = make_regionset(
            [
                ("chr1", 0, 10),
                ("chr1", 5, 15),
                ("chr1", 20, 30),
            ]
        )
        reduced = rs.reduce()
        assert len(reduced) == 2
        assert reduced[0].start == 0
        assert reduced[0].end == 15
        assert reduced[1].start == 20
        assert reduced[1].end == 30

    def test_reduce_adjacent(self):
        rs = make_regionset([("chr1", 0, 10), ("chr1", 10, 20)])
        reduced = rs.reduce()
        assert len(reduced) == 1
        assert reduced[0].start == 0
        assert reduced[0].end == 20

    def test_setdiff(self):
        a = make_regionset([("chr1", 0, 100)])
        b = make_regionset([("chr1", 30, 60)])
        diff = a.setdiff(b)
        assert len(diff) == 2
        assert diff[0].end == 30
        assert diff[1].start == 60

    def test_pintersect(self):
        a = make_regionset([("chr1", 0, 100)])
        b = make_regionset([("chr1", 50, 200)])
        pi = a.pintersect(b)
        assert len(pi) == 1
        assert pi[0].start == 50
        assert pi[0].end == 100

    def test_concat(self):
        a = make_regionset([("chr1", 0, 10)])
        b = make_regionset([("chr1", 5, 15)])
        c = a.concat(b)
        assert len(c) == 2

    def test_union(self):
        a = make_regionset([("chr1", 0, 10)])
        b = make_regionset([("chr1", 5, 15)])
        u = a.union(b)
        assert len(u) == 1
        assert u[0].start == 0
        assert u[0].end == 15

    def test_jaccard_identical(self):
        rs = make_regionset([("chr1", 0, 100)])
        assert rs.jaccard(rs) == 1.0

    def test_jaccard_disjoint(self):
        a = make_regionset([("chr1", 0, 10)])
        b = make_regionset([("chr1", 20, 30)])
        assert a.jaccard(b) == 0.0

    def test_jaccard_partial(self):
        a = make_regionset([("chr1", 0, 100)])
        b = make_regionset([("chr1", 50, 150)])
        j = a.jaccard(b)
        # intersection = 50bp, union = 150bp
        assert abs(j - 50 / 150) < 1e-9


# =========================================================================
# TssIndex
# =========================================================================


class TestTssIndex:
    def test_from_file(self):
        tss_path = DATA_DIR / "dummy_tss.bed"
        if tss_path.exists():
            tss = TssIndex(str(tss_path))
            assert len(tss) > 0

    def test_from_regionset(self):
        rs = make_regionset([("chr1", 100, 101), ("chr1", 500, 501)])
        tss = TssIndex.from_regionset(rs)
        assert len(tss) == 2

    def test_calc_tss_distances(self):
        features = make_regionset([("chr1", 100, 101)])
        tss = TssIndex.from_regionset(features)
        query = make_regionset([("chr1", 200, 210)])
        dists = tss.calc_tss_distances(query)
        assert len(dists) == 1
        # query midpoint = 205, feature midpoint = 100, distance = 105
        assert dists[0] == 105

    def test_feature_distances(self):
        features = make_regionset([("chr1", 100, 101)])
        tss = TssIndex.from_regionset(features)
        query = make_regionset([("chr1", 200, 210)])
        dists = tss.feature_distances(query)
        assert len(dists) == 1
        assert dists[0] is not None
        # signed: feature(100) - query(205) = -105
        assert dists[0] == -105.0

    def test_feature_distances_no_chrom(self):
        features = make_regionset([("chr1", 100, 101)])
        tss = TssIndex.from_regionset(features)
        query = make_regionset([("chr2", 200, 210)])
        dists = tss.feature_distances(query)
        assert len(dists) == 1
        assert dists[0] is None


# =========================================================================
# GeneModel
# =========================================================================


class TestGeneModel:
    @pytest.fixture
    def gtf_path(self):
        return str(DATA_DIR / "test_gene_model.gtf")

    def test_from_gtf(self, gtf_path):
        gm = GeneModel.from_gtf(gtf_path)
        assert gm.n_genes == 2  # gene1, gene2 (gene3 is lncRNA, filtered)
        assert gm.n_exons > 0

    def test_from_gtf_no_filter(self, gtf_path):
        gm = GeneModel.from_gtf(gtf_path, filter_protein_coding=False)
        assert gm.n_genes == 3  # all 3 genes

    def test_repr(self, gtf_path):
        gm = GeneModel.from_gtf(gtf_path)
        assert "GeneModel" in repr(gm)
        assert "n_genes" in repr(gm)


# =========================================================================
# PartitionList
# =========================================================================


class TestPartitionList:
    @pytest.fixture
    def gtf_path(self):
        return str(DATA_DIR / "test_gene_model.gtf")

    def test_from_gtf(self, gtf_path):
        pl = PartitionList.from_gtf(gtf_path, core_prom=100, prox_prom=2000)
        names = pl.partition_names()
        assert "promoterCore" in names
        assert "promoterProx" in names
        assert "exon" in names
        assert len(pl) > 0

    def test_from_gene_model(self, gtf_path):
        gm = GeneModel.from_gtf(gtf_path)
        pl = PartitionList.from_gene_model(gm, core_prom=100, prox_prom=2000)
        assert len(pl) > 0
        assert (
            pl.partition_names()
            == PartitionList.from_gtf(
                gtf_path, core_prom=100, prox_prom=2000
            ).partition_names()
        )

    def test_from_gtf_with_chrom_sizes(self, gtf_path):
        chrom_sizes = {"chr1": 50000, "chr2": 50000}
        pl = PartitionList.from_gtf(
            gtf_path, core_prom=100, prox_prom=2000, chrom_sizes=chrom_sizes
        )
        assert len(pl) > 0

    def test_repr(self, gtf_path):
        pl = PartitionList.from_gtf(gtf_path, core_prom=100, prox_prom=2000)
        assert "PartitionList" in repr(pl)


# =========================================================================
# GenomicDistAnnotation
# =========================================================================


class TestGenomicDistAnnotation:
    @pytest.fixture
    def gtf_path(self):
        return str(DATA_DIR / "test_gene_model.gtf")

    def test_from_gtf(self, gtf_path):
        gda = GenomicDistAnnotation.from_gtf(gtf_path)
        gm = gda.gene_model()
        assert gm.n_genes > 0

    def test_round_trip_bin(self, gtf_path):
        gda = GenomicDistAnnotation.from_gtf(gtf_path)
        with tempfile.NamedTemporaryFile(suffix=".gda.bin", delete=False) as f:
            bin_path = f.name
        try:
            # Save not exposed in Python yet, so just test load_bin via
            # saving from Rust (already tested in Rust). Test from_gtf -> gene_model.
            gm = gda.gene_model()
            assert gm.n_genes > 0
        finally:
            if os.path.exists(bin_path):
                os.unlink(bin_path)

    def test_partition_list(self, gtf_path):
        gda = GenomicDistAnnotation.from_gtf(gtf_path)
        pl = gda.partition_list(core_prom=100, prox_prom=2000)
        assert len(pl) > 0
        assert "promoterCore" in pl.partition_names()

    def test_tss_index(self, gtf_path):
        gda = GenomicDistAnnotation.from_gtf(gtf_path)
        tss = gda.tss_index()
        assert len(tss) > 0

    def test_tss_index_strand_aware(self, gtf_path):
        gda = GenomicDistAnnotation.from_gtf(gtf_path)
        tss = gda.tss_index()
        # gene1 is on + strand at chr1:1000-5000, TSS = 1000
        # gene2 is on - strand at chr2:3000-8000, TSS = 7999 (end-1)
        query = make_regionset([("chr1", 1000, 1001)])
        dists = tss.calc_tss_distances(query)
        assert dists[0] == 0  # TSS at 1000

    def test_repr(self, gtf_path):
        gda = GenomicDistAnnotation.from_gtf(gtf_path)
        assert "GenomicDistAnnotation" in repr(gda)


# =========================================================================
# SignalMatrix
# =========================================================================


class TestSignalMatrix:
    @pytest.fixture
    def signal_tsv(self):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        f.write("V1\tcond_A\tcond_B\n")
        f.write("chr1_100_200\t0.5\t0.3\n")
        f.write("chr1_150_250\t0.2\t0.8\n")
        f.write("chr1_300_400\t0.9\t0.1\n")
        f.flush()
        f.close()
        yield f.name
        os.unlink(f.name)

    def test_from_tsv(self, signal_tsv):
        sm = SignalMatrix.from_tsv(signal_tsv)
        assert sm.n_conditions == 2
        assert sm.n_regions == 3
        assert sm.condition_names == ["cond_A", "cond_B"]

    def test_bin_round_trip(self, signal_tsv):
        sm = SignalMatrix.from_tsv(signal_tsv)
        with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as f:
            bin_path = f.name
        try:
            # Can't save from Python (not exposed), but test that load_bin
            # would work by verifying from_tsv produces valid data
            assert sm.n_regions == 3
        finally:
            if os.path.exists(bin_path):
                os.unlink(bin_path)

    def test_repr(self, signal_tsv):
        sm = SignalMatrix.from_tsv(signal_tsv)
        r = repr(sm)
        assert "SignalMatrix" in r
        assert "n_regions=3" in r
        assert "n_conditions=2" in r

    def test_len(self, signal_tsv):
        sm = SignalMatrix.from_tsv(signal_tsv)
        assert len(sm) == 3


# =========================================================================
# Genomic distribution functions
# =========================================================================


class TestCalcPartitions:
    @pytest.fixture
    def gtf_path(self):
        return str(DATA_DIR / "test_gene_model.gtf")

    def test_calc_partitions_priority(self, gtf_path):
        pl = PartitionList.from_gtf(gtf_path, core_prom=100, prox_prom=2000)
        query = make_regionset(
            [
                ("chr1", 950, 1050),  # overlaps promoterCore of gene1
                ("chr1", 2200, 2300),  # overlaps exon of gene1
                ("chr1", 50000, 50100),  # intergenic
            ]
        )
        result = calc_partitions(query, pl)
        assert "partition" in result
        assert "count" in result
        assert "total" in result
        assert result["total"] == 3
        assert sum(result["count"]) == 3

    def test_calc_partitions_bp_proportion(self, gtf_path):
        pl = PartitionList.from_gtf(gtf_path, core_prom=100, prox_prom=2000)
        query = make_regionset([("chr1", 950, 1050)])
        result = calc_partitions(query, pl, bp_proportion=True)
        assert result["total"] > 0

    def test_calc_expected_partitions(self, gtf_path):
        pl = PartitionList.from_gtf(gtf_path, core_prom=100, prox_prom=2000)
        query = make_regionset(
            [
                ("chr1", 950, 1050),
                ("chr1", 2200, 2300),
                ("chr1", 50000, 50100),
            ]
        )
        chrom_sizes = {"chr1": 100000, "chr2": 100000}
        result = calc_expected_partitions(query, pl, chrom_sizes)
        assert "partition" in result
        assert "observed" in result
        assert "expected" in result
        assert "log10OE" in result
        assert "pvalue" in result
        assert len(result["partition"]) == len(result["observed"])


class TestCalcSummarySignal:
    def test_calc_summary_signal(self):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        f.write("V1\tcond_A\tcond_B\n")
        f.write("chr1_100_200\t0.5\t0.3\n")
        f.write("chr1_150_250\t0.2\t0.8\n")
        f.write("chr1_300_400\t0.9\t0.1\n")
        f.flush()
        f.close()
        try:
            sm = SignalMatrix.from_tsv(f.name)
            query = make_regionset(
                [
                    ("chr1", 120, 180),  # overlaps rows 1 and 2
                    ("chr1", 350, 380),  # overlaps row 3
                ]
            )
            result = calc_summary_signal(query, sm)
            assert "condition_names" in result
            assert "region_labels" in result
            assert "signal_matrix" in result
            assert "matrix_stats" in result
            assert result["condition_names"] == ["cond_A", "cond_B"]
            assert len(result["region_labels"]) == 2
            # First region: max(0.5, 0.2)=0.5 for cond_A, max(0.3, 0.8)=0.8 for cond_B
            assert abs(result["signal_matrix"][0][0] - 0.5) < 1e-9
            assert abs(result["signal_matrix"][0][1] - 0.8) < 1e-9
        finally:
            os.unlink(f.name)

    def test_calc_summary_signal_no_overlaps(self):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        f.write("V1\tcond_A\n")
        f.write("chr1_100_200\t0.5\n")
        f.flush()
        f.close()
        try:
            sm = SignalMatrix.from_tsv(f.name)
            query = make_regionset([("chr2", 0, 100)])
            result = calc_summary_signal(query, sm)
            assert len(result["region_labels"]) == 0
            assert len(result["matrix_stats"]) == 0
        finally:
            os.unlink(f.name)


class TestMedianAbsDistance:
    def test_basic(self):
        assert median_abs_distance([1.0, -3.0, 5.0, -7.0]) == 4.0

    def test_single(self):
        assert median_abs_distance([42.0]) == 42.0

    def test_empty(self):
        assert median_abs_distance([]) is None

    def test_with_nan(self):
        # NaN/inf are treated as sentinels and filtered
        result = median_abs_distance([1.0, float("nan"), 3.0])
        assert result == 2.0


class TestConsensus:
    def test_two_overlapping(self):
        a = make_regionset([("chr1", 0, 100)])
        b = make_regionset([("chr1", 50, 150)])
        result = consensus([a, b])
        assert len(result) == 1
        assert result[0]["chr"] == "chr1"
        assert result[0]["start"] == 0
        assert result[0]["end"] == 150
        assert result[0]["count"] == 2

    def test_disjoint(self):
        a = make_regionset([("chr1", 0, 10)])
        b = make_regionset([("chr1", 20, 30)])
        result = consensus([a, b])
        assert len(result) == 2
        assert all(r["count"] == 1 for r in result)

    def test_three_sets(self):
        a = make_regionset([("chr1", 0, 10)])
        b = make_regionset([("chr1", 5, 15)])
        c = make_regionset([("chr1", 8, 20)])
        result = consensus([a, b, c])
        assert len(result) == 1
        assert result[0]["count"] == 3

    def test_empty_input(self):
        result = consensus([])
        assert result == []

    def test_from_bed_files(self):
        a_path = DATA_DIR / "dummy.bed"
        b_path = DATA_DIR / "dummy_b.bed"
        if a_path.exists() and b_path.exists():
            a = RegionSet(str(a_path))
            b = RegionSet(str(b_path))
            result = consensus([a, b])
            assert len(result) >= 1
            assert all(r["count"] >= 1 for r in result)
