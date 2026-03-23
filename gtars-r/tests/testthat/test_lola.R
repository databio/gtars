library(testthat)
library(gtars)

# Small LOLA database shipped in tests/data (copied from LOLA R package extdata/multi)
db_path <- normalizePath(
  file.path(test_path(), "..", "..", "..", "tests", "data", "lola_multi_db"),
  mustWork = FALSE
)

skip_if_no_db <- function() {
  skip_if_not(dir.exists(db_path), "lola_multi_db test data not found")
}

# Helper: grab BED file paths from a collection
bed_files <- function(collection = "collection1") {
  list.files(
    file.path(db_path, collection, "regions"),
    pattern = "\\.bed$", full.names = TRUE
  )
}

# =========================================================================
# Database loading
# =========================================================================

test_that("loadRegionDB loads from folder", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path)
  expect_type(regionDB, "externalptr")
})

test_that("loadRegionDB filters by collection", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path, collections = "collection1")
  expect_type(regionDB, "externalptr")
  files <- listRegionSets(regionDB)
  expect_true(length(files) > 0)
})

test_that("loadRegionDBFromBeds works with BED file vector", {
  skip_if_no_db()
  beds <- bed_files("collection1")
  regionDB <- loadRegionDBFromBeds(beds)
  expect_type(regionDB, "externalptr")
})

# =========================================================================
# Database accessors
# =========================================================================

test_that("regionDBAnno returns data.frame with filename/collection columns", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path, collections = "collection1")
  anno <- regionDBAnno(regionDB)
  expect_s3_class(anno, "data.frame")
  expect_true("filename" %in% names(anno))
  expect_true("collection" %in% names(anno))
})

test_that("regionDBCollectionAnno returns data.frame", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path, collections = "collection1")
  coll_anno <- regionDBCollectionAnno(regionDB)
  expect_s3_class(coll_anno, "data.frame")
})

test_that("regionDBRegionSet returns a RegionSet", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path, collections = "collection1")
  rs <- regionDBRegionSet(regionDB, 1L)
  expect_s4_class(rs, "RegionSet")
  expect_true(length(rs) > 0)
})

test_that("listRegionSets returns character vector matching file count", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path, collections = "collection1")
  files <- listRegionSets(regionDB)
  expect_type(files, "character")
  expect_equal(length(files), length(bed_files("collection1")))
})

# =========================================================================
# Enrichment analysis
# =========================================================================

test_that("runLOLA returns data.frame with expected columns", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path, collections = "collection1")

  userSet <- RegionSet(bed_files("collection1")[1])
  beds <- bed_files("collection1")
  universe <- buildRestrictedUniverse(lapply(beds, RegionSet))

  result <- runLOLA(userSet, universe, regionDB)
  expect_s3_class(result, "data.frame")
  expected_cols <- c("userSet", "dbSet", "pValueLog", "oddsRatio",
                     "support", "b", "c", "d", "filename")
  for (col in expected_cols) {
    expect_true(col %in% names(result), info = paste("Missing column:", col))
  }
})

test_that("runLOLA support values are non-negative integers", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path, collections = "collection1")

  userSet <- RegionSet(bed_files("collection1")[1])
  beds <- bed_files("collection1")
  universe <- buildRestrictedUniverse(lapply(beds, RegionSet))

  result <- runLOLA(userSet, universe, regionDB)
  expect_true(all(result$support >= 0))
  expect_true(all(result$support == as.integer(result$support)))
})

test_that("runLOLA with direction='depletion' runs without error", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path, collections = "collection1")

  userSet <- RegionSet(bed_files("collection1")[1])
  beds <- bed_files("collection1")
  universe <- buildRestrictedUniverse(lapply(beds, RegionSet))

  expect_no_error(runLOLA(userSet, universe, regionDB, direction = "depletion"))
})

# =========================================================================
# Universe functions
# =========================================================================

test_that("checkUniverseAppropriateness runs without error", {
  skip_if_no_db()

  userSet <- RegionSet(bed_files("collection1")[1])
  beds <- bed_files("collection1")
  universe <- buildRestrictedUniverse(lapply(beds, RegionSet))

  expect_no_error(checkUniverseAppropriateness(userSet, universe))
})

test_that("redefineUserSets returns list of RegionSets", {
  skip_if_no_db()

  userSet <- RegionSet(bed_files("collection1")[1])
  beds <- bed_files("collection1")
  universe <- buildRestrictedUniverse(lapply(beds, RegionSet))

  result <- redefineUserSets(userSet, universe)
  expect_type(result, "list")
  expect_true(length(result) >= 1)
  expect_s4_class(result[[1]], "RegionSet")
})

test_that("buildRestrictedUniverse returns a RegionSet", {
  skip_if_no_db()

  beds <- bed_files("collection1")
  sets <- lapply(beds, RegionSet)

  result <- buildRestrictedUniverse(sets)
  expect_s4_class(result, "RegionSet")
  expect_true(length(result) > 0)
})
