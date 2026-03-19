"""Tests for the LOLA (Locus Overlap Analysis) Python bindings."""

import pytest

from gtars.lola import RegionDB, run_lola, check_universe, redefine_user_sets


def make_regions(specs):
    """Convert list of (chr, start, end) tuples to the format expected by LOLA."""
    return [(chr, start, end) for chr, start, end in specs]


@pytest.fixture
def simple_db():
    """Create a simple RegionDB from in-memory BED-like data."""
    import tempfile
    import os

    tmpdir = tempfile.mkdtemp()
    # Create two BED files for the database
    bed1_path = os.path.join(tmpdir, "db1.bed")
    bed2_path = os.path.join(tmpdir, "db2.bed")

    with open(bed1_path, "w") as f:
        f.write("chr1\t100\t200\n")
        f.write("chr1\t500\t600\n")

    with open(bed2_path, "w") as f:
        f.write("chr1\t150\t250\n")
        f.write("chr1\t700\t800\n")

    db = RegionDB.from_bed_files([bed1_path, bed2_path])
    return db


@pytest.fixture
def universe():
    """A universe covering a broad range."""
    return make_regions(
        [
            ("chr1", 50, 250),
            ("chr1", 450, 650),
            ("chr1", 650, 850),
            ("chr1", 900, 1000),
            ("chr1", 1100, 1200),
        ]
    )


@pytest.fixture
def lolacore_db():
    """Load a RegionDB from the LOLA-format test folder (has collection.txt)."""
    import os

    db_path = os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "tests",
        "data",
        "lola_multi_db",
    )
    return RegionDB.from_folder(db_path)


class TestCollectionAnno:
    def test_collection_anno_returns_list(self, lolacore_db):
        annos = lolacore_db.collection_anno
        assert isinstance(annos, list)
        assert len(annos) >= 1

    def test_collection_anno_keys(self, lolacore_db):
        annos = lolacore_db.collection_anno
        expected_keys = {"collectionname", "collector", "date", "source", "description"}
        assert set(annos[0].keys()) == expected_keys

    def test_collection_anno_values(self, lolacore_db):
        annos = lolacore_db.collection_anno
        # collection1 has: Nathan, 2015-04-01, UCSC Genome Browser, ...
        names = [a["collectionname"] for a in annos]
        assert "collection1" in names

    def test_collection_anno_empty_for_bed_files(self, simple_db):
        """RegionDB built from BED files should have no collection annotations."""
        assert simple_db.collection_anno == []


class TestRunLola:
    def test_run_lola_basic(self, simple_db, universe):
        user_set = make_regions([("chr1", 120, 180)])
        result = run_lola([user_set], universe, simple_db)

        assert "pValueLog" in result
        assert "oddsRatio" in result
        assert "support" in result
        assert "qValue" in result
        assert "filename" in result
        assert len(result["pValueLog"]) == 2  # 1 user set x 2 db sets

    def test_run_lola_has_q_values(self, simple_db, universe):
        user_set = make_regions([("chr1", 120, 180)])
        result = run_lola([user_set], universe, simple_db)

        # q_value should be populated (not all None)
        q_values = result["qValue"]
        assert len(q_values) == 2
        # At least some q_values should be numeric (not None)
        assert any(q is not None for q in q_values)

    def test_run_lola_depletion(self, simple_db, universe):
        user_set = make_regions([("chr1", 120, 180)])
        result_enrich = run_lola(
            [user_set], universe, simple_db, direction="enrichment"
        )
        result_deplete = run_lola(
            [user_set], universe, simple_db, direction="depletion"
        )

        # p-values should differ between enrichment and depletion
        assert result_enrich["pValueLog"] != result_deplete["pValueLog"]


class TestCheckUniverse:
    def test_check_universe(self, universe):
        user_set = make_regions(
            [
                ("chr1", 120, 180),
                ("chr1", 500, 550),
            ]
        )
        result = check_universe([user_set], universe)

        assert "totalRegions" in result
        assert "coverage" in result
        assert result["totalRegions"][0] == 2


class TestRedefineUserSets:
    def test_redefine_user_sets(self, universe):
        user_set = make_regions([("chr1", 120, 180)])
        redefined = redefine_user_sets([user_set], universe)

        assert len(redefined) == 1
        # The redefined set should contain universe regions that overlap the user set
        assert len(redefined[0]) >= 1
