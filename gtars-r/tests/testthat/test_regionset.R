library(testthat)
library(gtars)

data_dir <- normalizePath(
  file.path(test_path(), "..", "..", "..", "tests", "data", "regionset"),
  mustWork = FALSE
)

skip_if_no_data <- function() {
  skip_if_not(dir.exists(data_dir), "regionset test data not found")
}

# =========================================================================
# Constructor
# =========================================================================

test_that("RegionSet from file path works", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  expect_s4_class(rs, "RegionSet")
  expect_true(length(rs) > 0)
})

test_that("RegionSet from data.frame works", {
  df <- data.frame(chr = "chr1", start = c(0L, 20L), end = c(10L, 30L))
  rs <- RegionSet(df)
  expect_s4_class(rs, "RegionSet")
  expect_equal(length(rs), 2L)
})

test_that("RegionSet from RegionSet is identity", {
  df <- data.frame(chr = "chr1", start = 0L, end = 10L)
  rs <- RegionSet(df)
  expect_identical(RegionSet(rs), rs)
})

test_that("as_regionset is an alias for RegionSet", {
  df <- data.frame(chr = "chr1", start = 0L, end = 10L)
  expect_s4_class(as_regionset(df), "RegionSet")
})

# =========================================================================
# Basic S4 methods
# =========================================================================

test_that("length returns correct count", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  expect_equal(length(rs), 4L)
})

test_that("as.data.frame round-trips correctly", {
  df <- data.frame(chr = "chr1", start = c(0L, 20L), end = c(10L, 30L))
  rs <- RegionSet(df)
  out <- as.data.frame(rs)
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 2L)
  expect_true(all(c("chr", "start", "end") %in% names(out)))
  expect_equal(out$start, c(0L, 20L))
})

test_that("show does not error", {
  df <- data.frame(chr = "chr1", start = 0L, end = 10L)
  rs <- RegionSet(df)
  expect_no_error(capture.output(show(rs)))
})

test_that("subsetting with [ works", {
  df <- data.frame(chr = "chr1", start = c(0L, 10L, 20L), end = c(5L, 15L, 25L))
  rs <- RegionSet(df)
  rs_sub <- rs[1:2]
  expect_s4_class(rs_sub, "RegionSet")
  expect_equal(length(rs_sub), 2L)
})

# =========================================================================
# Unary statistics
# =========================================================================

test_that("widths returns correct values", {
  skip_if_no_data()
  # dummy.bed: chr1:2-6, chr1:4-7, chr1:5-9, chr1:7-12
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  w <- widths(rs)
  # doubles instead of integers: u32 widths can exceed i32::MAX (2.1 Gbp)
  expect_type(w, "double")
  expect_equal(length(w), 4L)
  expect_equal(w[1], 4)
})

test_that("neighborDistances returns numeric vector", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  nd <- neighborDistances(rs)
  expect_true(is.numeric(nd))
})

test_that("nearestNeighbors returns numeric vector", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  nn <- nearestNeighbors(rs)
  # doubles instead of integers: distances can be chromosome-scale
  expect_type(nn, "double")
})

test_that("chromosomeStatistics returns data.frame", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  cs <- chromosomeStatistics(rs)
  expect_true(is.data.frame(cs))
})

# =========================================================================
# GenomicDistributions drop-in wrappers
# =========================================================================

test_that("calcWidth delegates to widths", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  expect_equal(calcWidth(rs), widths(rs))
})

test_that("calcNeighborDist delegates to neighborDistances", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  expect_equal(calcNeighborDist(rs), neighborDistances(rs))
})

test_that("calcNearestNeighbors delegates to nearestNeighbors", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  expect_equal(calcNearestNeighbors(rs), nearestNeighbors(rs))
})

# =========================================================================
# Unary interval operations
# =========================================================================

test_that("reduce merges overlapping regions", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  rs_red <- reduce(rs)
  expect_s4_class(rs_red, "RegionSet")
  expect_equal(length(rs_red), 1L)
  df <- as.data.frame(rs_red)
  expect_equal(df$start, 2L)
  expect_equal(df$end, 12L)
})

test_that("promoters returns same region count", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  rs_prom <- promoters(rs, upstream = 100L, downstream = 50L)
  expect_s4_class(rs_prom, "RegionSet")
  expect_equal(length(rs_prom), length(rs))
})

test_that("trim clips to chromosome boundary", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  rs_trim <- trim(rs, chromSizes = c(chr1 = 8L))
  expect_s4_class(rs_trim, "RegionSet")
  df <- as.data.frame(rs_trim)
  expect_true(all(df$end <= 8))
})

test_that("shift moves regions", {
  df <- data.frame(chr = "chr1", start = c(10L, 20L), end = c(15L, 25L))
  rs <- RegionSet(df)
  rs_shifted <- shift(rs, shift = 5L)
  out <- as.data.frame(rs_shifted)
  expect_equal(out$start, c(15L, 25L))
  expect_equal(out$end, c(20L, 30L))
})

test_that("disjoin returns non-overlapping pieces", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  rs_dj <- disjoin(rs)
  expect_s4_class(rs_dj, "RegionSet")
  expect_true(length(rs_dj) >= 1L)
})

test_that("gaps returns inter-region gaps bounded by chrom_sizes", {
  skip_if_no_data()
  # Three non-overlapping peaks → leading, 2 inter-region, trailing = 4 gaps.
  df <- data.frame(
    chr = "chr1",
    start = c(10L, 30L, 50L),
    end   = c(20L, 40L, 60L)
  )
  rs <- RegionSet(df)
  chrom_sizes <- c(chr1 = 100L)
  rs_gaps <- gaps(rs, chrom_sizes = chrom_sizes)
  expect_s4_class(rs_gaps, "RegionSet")
  expect_equal(length(rs_gaps), 4L)
})

test_that("gaps errors when chrom_sizes is missing", {
  df <- data.frame(chr = "chr1", start = 10L, end = 20L)
  rs <- RegionSet(df)
  expect_error(gaps(rs), "chrom_sizes")
})

test_that("gaps emits full-chromosome gap for empty chromosomes", {
  df <- data.frame(chr = "chr1", start = 10L, end = 20L)
  rs <- RegionSet(df)
  chrom_sizes <- c(chr1 = 100L, chr2 = 50L)
  rs_gaps <- gaps(rs, chrom_sizes = chrom_sizes)
  gaps_df <- as.data.frame(rs_gaps)
  chr2 <- gaps_df[gaps_df$chr == "chr2", ]
  expect_equal(nrow(chr2), 1L)
  expect_equal(chr2$start, 0L)
  expect_equal(chr2$end, 50L)
})

# =========================================================================
# Spatial-arrangement summary statistics
# =========================================================================

test_that("clusterRegions returns cluster IDs", {
  df <- data.frame(
    chr = "chr1",
    start = c(0L, 13L, 100L),
    end   = c(10L, 20L, 110L)
  )
  rs <- RegionSet(df)
  ids <- clusterRegions(rs, maxGap = 5L)
  expect_length(ids, 3L)
  expect_equal(ids[1], ids[2])   # first two are within 5bp → same cluster
  expect_false(ids[1] == ids[3]) # third is far away
})

# =========================================================================
# Binary interval operations
# =========================================================================

test_that("concat combines region counts", {
  skip_if_no_data()
  rs_a <- RegionSet(file.path(data_dir, "dummy.bed"))
  rs_b <- RegionSet(file.path(data_dir, "dummy_b.bed"))
  rs_cat <- concat(rs_a, rs_b)
  expect_s4_class(rs_cat, "RegionSet")
  expect_equal(length(rs_cat), length(rs_a) + length(rs_b))
})

test_that("union merges overlapping regions", {
  skip_if_no_data()
  rs_a <- RegionSet(file.path(data_dir, "dummy.bed"))
  rs_b <- RegionSet(file.path(data_dir, "dummy_b.bed"))
  rs_u <- union(rs_a, rs_b)
  expect_s4_class(rs_u, "RegionSet")
  df <- as.data.frame(rs_u)
  expect_equal(nrow(df), 1L)
  expect_equal(df$start, 2L)
  expect_equal(df$end, 12L)
})

test_that("setdiff returns RegionSet", {
  skip_if_no_data()
  rs_a <- RegionSet(file.path(data_dir, "dummy.bed"))
  rs_b <- RegionSet(file.path(data_dir, "dummy_b.bed"))
  rs_sd <- setdiff(rs_a, rs_b)
  expect_s4_class(rs_sd, "RegionSet")
})

test_that("pintersect returns RegionSet", {
  df <- data.frame(chr = "chr1", start = c(0L, 20L), end = c(10L, 30L))
  rs <- RegionSet(df)
  rs_pi <- pintersect(rs, rs)
  expect_s4_class(rs_pi, "RegionSet")
})

test_that("jaccard returns correct values", {
  skip_if_no_data()
  rs_a <- RegionSet(file.path(data_dir, "dummy.bed"))
  rs_b <- RegionSet(file.path(data_dir, "dummy_b.bed"))
  j <- jaccard(rs_a, rs_b)
  expect_type(j, "double")
  expect_true(j >= 0 && j <= 1)

  j_self <- jaccard(rs_a, rs_a)
  expect_equal(j_self, 1.0, tolerance = 1e-10)
})

test_that("findOverlaps returns data.frame with hits", {
  df_a <- data.frame(chr = "chr1", start = c(0L, 20L), end = c(10L, 30L))
  df_b <- data.frame(chr = "chr1", start = c(5L, 25L), end = c(15L, 35L))
  rs_a <- RegionSet(df_a)
  rs_b <- RegionSet(df_b)
  hits <- findOverlaps(rs_a, rs_b)
  expect_true(is.data.frame(hits))
  expect_true(all(c("queryHits", "subjectHits") %in% names(hits)))
  expect_true(nrow(hits) >= 2)
})

test_that("countOverlaps returns integer vector", {
  df_a <- data.frame(chr = "chr1", start = c(0L, 20L), end = c(10L, 30L))
  df_b <- data.frame(chr = "chr1", start = c(5L, 25L), end = c(15L, 35L))
  rs_a <- RegionSet(df_a)
  rs_b <- RegionSet(df_b)
  counts <- countOverlaps(rs_a, rs_b)
  expect_type(counts, "integer")
  expect_equal(length(counts), 2L)
  expect_true(all(counts >= 1))
})

# =========================================================================
# Consensus
# =========================================================================

test_that("consensus returns data.frame with count column", {
  skip_if_no_data()
  rs_a <- RegionSet(file.path(data_dir, "dummy.bed"))
  rs_b <- RegionSet(file.path(data_dir, "dummy_b.bed"))
  cons <- consensus(list(rs_a, rs_b))
  expect_true(is.data.frame(cons))
  expect_true("count" %in% names(cons))
})

test_that("consensus with disjoint sets has multiple regions", {
  skip_if_no_data()
  rs_a <- RegionSet(file.path(data_dir, "dummy.bed"))
  rs_far <- RegionSet(data.frame(chr = "chr2", start = 100L, end = 200L))
  cons <- consensus(list(rs_a, rs_far))
  expect_equal(nrow(cons), 2L)
})
