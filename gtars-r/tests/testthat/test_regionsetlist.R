library(testthat)
library(gtars)

# Small LOLA database shipped in tests/data
db_path <- normalizePath(
  file.path(test_path(), "..", "..", "..", "tests", "data", "lola_multi_db"),
  mustWork = FALSE
)

skip_if_no_db <- function() {
  skip_if_not(dir.exists(db_path), "lola_multi_db test data not found")
}

# =========================================================================
# Constructor
# =========================================================================

test_that("RegionSetList can be created from RegionSet objects", {
  rs1 <- RegionSet(data.frame(
    chr = c("chr1", "chr1"), start = c(0L, 200L), end = c(100L, 300L)
  ))
  rs2 <- RegionSet(data.frame(
    chr = "chr1", start = 50L, end = 150L
  ))

  rsl <- RegionSetList(rs1, rs2)
  expect_s4_class(rsl, "RegionSetList")
  expect_equal(length(rsl), 2)
})

test_that("RegionSetList can be created from a list of RegionSets", {
  rs1 <- RegionSet(data.frame(chr = "chr1", start = 0L, end = 100L))
  rs2 <- RegionSet(data.frame(chr = "chr2", start = 0L, end = 50L))

  rsl <- RegionSetList(list(rs1, rs2))
  expect_s4_class(rsl, "RegionSetList")
  expect_equal(length(rsl), 2)
})

test_that("RegionSetList can be created from file paths", {
  skip_if_no_db()
  beds <- list.files(
    file.path(db_path, "collection1", "regions"),
    pattern = "\\.bed$", full.names = TRUE
  )
  rsl <- RegionSetList(lapply(beds, RegionSet))
  expect_s4_class(rsl, "RegionSetList")
  expect_equal(length(rsl), length(beds))
})

# =========================================================================
# [[ accessor
# =========================================================================

test_that("[[ returns a RegionSet", {
  rs1 <- RegionSet(data.frame(chr = "chr1", start = 0L, end = 100L))
  rs2 <- RegionSet(data.frame(chr = "chr2", start = 0L, end = 50L))

  rsl <- RegionSetList(rs1, rs2)
  first <- rsl[[1]]
  expect_s4_class(first, "RegionSet")
  expect_equal(length(first), 1)
})

test_that("[[ returns correct element by index", {
  rs1 <- RegionSet(data.frame(
    chr = c("chr1", "chr1"), start = c(0L, 200L), end = c(100L, 300L)
  ))
  rs2 <- RegionSet(data.frame(chr = "chr1", start = 50L, end = 150L))

  rsl <- RegionSetList(rs1, rs2)
  expect_equal(length(rsl[[1]]), 2)
  expect_equal(length(rsl[[2]]), 1)
})

# =========================================================================
# concat
# =========================================================================

test_that("concat flattens RegionSetList into single RegionSet", {
  rs1 <- RegionSet(data.frame(
    chr = c("chr1", "chr1"), start = c(0L, 200L), end = c(100L, 300L)
  ))
  rs2 <- RegionSet(data.frame(chr = "chr1", start = 50L, end = 150L))

  rsl <- RegionSetList(rs1, rs2)
  combined <- concat(rsl)
  expect_s4_class(combined, "RegionSet")
  expect_equal(length(combined), 3)  # 2 + 1, no merging
})

test_that("concat preserves all regions without deduplication", {
  rs1 <- RegionSet(data.frame(chr = "chr1", start = 0L, end = 100L))
  rs2 <- RegionSet(data.frame(chr = "chr1", start = 0L, end = 100L))

  rsl <- RegionSetList(rs1, rs2)
  combined <- concat(rsl)
  # Identical regions are NOT merged — concat is raw union
  expect_equal(length(combined), 2)
})

# =========================================================================
# getRegionSets from DB
# =========================================================================

test_that("getRegionSets extracts from DB by index", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path)

  rsl <- getRegionSets(regionDB, 1:3)
  expect_s4_class(rsl, "RegionSetList")
  expect_equal(length(rsl), 3)

  for (i in seq_len(length(rsl))) {
    expect_s4_class(rsl[[i]], "RegionSet")
    expect_true(length(rsl[[i]]) > 0)
  }
})

test_that("getRegionSets with no indices returns all", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path)
  anno <- regionDBAnno(regionDB)

  rsl <- getRegionSets(regionDB)
  expect_equal(length(rsl), nrow(anno))
})

# =========================================================================
# Universe construction
# =========================================================================

test_that("universe construction: DB -> getRegionSets -> concat -> reduce", {
  skip_if_no_db()
  regionDB <- loadRegionDB(db_path)

  rsl <- getRegionSets(regionDB)
  combined <- concat(rsl)
  expect_s4_class(combined, "RegionSet")
  expect_true(length(combined) > 0)

  universe <- reduce(combined)
  expect_s4_class(universe, "RegionSet")
  expect_true(length(universe) > 0)
  # reduce should not increase the number of regions
  expect_true(length(universe) <= length(combined))
})
