#!/usr/bin/env Rscript
# Smoke test for RegionSet S4 class
# Run: Rscript tests/test_regionset.R

library(gtars)

setwd('/Users/sam/Documents/Work/ai-sandbox/workspaces/sam/bedbase/repos/gtars/gtars-r')

pass <- 0
fail <- 0
errors <- list()

check <- function(label, expr) {
  result <- tryCatch(
    {
      val <- eval(expr)
      if (isTRUE(val) || (!is.logical(val) && !is.null(val))) {
        pass <<- pass + 1
        cat(sprintf("  PASS  %s\n", label))
        TRUE
      } else {
        fail <<- fail + 1
        errors[[length(errors) + 1]] <<- label
        cat(sprintf("  FAIL  %s\n", label))
        FALSE
      }
    },
    error = function(e) {
      fail <<- fail + 1
      errors[[length(errors) + 1]] <<- paste0(label, ": ", e$message)
      cat(sprintf("  ERR   %s: %s\n", label, e$message))
      FALSE
    }
  )
  result
}

data_dir <- normalizePath(file.path("..", "tests", "data", "regionset"))

# =========================================================================
cat("\n=== 1. Constructor ===\n")
# =========================================================================

rs <- RegionSet(file.path(data_dir, "dummy.bed"))
check("RegionSet from file path", is(rs, "RegionSet"))
check("RegionSet has ptr slot", is(rs@ptr, "externalptr"))

df <- data.frame(chr = "chr1", start = c(0L, 20L), end = c(10L, 30L))
rs_df <- RegionSet(df)
check("RegionSet from data.frame", is(rs_df, "RegionSet"))

check("as_regionset alias works", is(as_regionset(df), "RegionSet"))

rs_same <- RegionSet(rs)
check("RegionSet from RegionSet is identity", identical(rs, rs_same))

# =========================================================================
cat("\n=== 2. Basic S4 Methods ===\n")
# =========================================================================

check("length(rs) returns integer", length(rs) == 4L)
check("length(rs_df) returns integer", length(rs_df) == 2L)

rs_df_out <- as.data.frame(rs_df)
check("as.data.frame returns data.frame", is.data.frame(rs_df_out))
check("as.data.frame has 2 rows", nrow(rs_df_out) == 2L)
check("as.data.frame has chr/start/end", all(c("chr", "start", "end") %in% names(rs_df_out)))
check("as.data.frame coords correct", rs_df_out$start[1] == 0L && rs_df_out$end[1] == 10L)

# show() should not error
tryCatch(
  { capture.output(show(rs)); check("show(RegionSet) works", TRUE) },
  error = function(e) check("show(RegionSet) works", FALSE)
)

# =========================================================================
cat("\n=== 3. Unary Statistics ===\n")
# =========================================================================

# dummy.bed: chr1:2-6, chr1:4-7, chr1:5-9, chr1:7-12
w <- widths(rs)
check("widths returns numeric vector", is.numeric(w))
check("widths has 4 values", length(w) == 4L)
check("widths first region is 4", w[1] == 4)

nd <- neighborDistances(rs)
check("neighborDistances returns numeric", is.numeric(nd))

nn <- nearestNeighbors(rs)
check("nearestNeighbors returns numeric", is.numeric(nn))

cs <- chromosomeStatistics(rs)
check("chromosomeStatistics returns data.frame", is.data.frame(cs))

# =========================================================================
cat("\n=== 4. GenomicDistributions Drop-in Wrappers ===\n")
# =========================================================================

w2 <- calcWidth(rs)
check("calcWidth delegates to widths", all(w == w2))

nd2 <- calcNeighborDist(rs)
check("calcNeighborDist delegates", all(nd == nd2))

nn2 <- calcNearestNeighbors(rs)
check("calcNearestNeighbors delegates", all(nn == nn2))

# =========================================================================
cat("\n=== 5. Unary Interval Operations ===\n")
# =========================================================================

# reduce: dummy.bed regions overlap, should merge to single region
rs_red <- reduce(rs)
check("reduce returns RegionSet", is(rs_red, "RegionSet"))
check("reduce merges overlapping regions", length(rs_red) == 1L)
df_red <- as.data.frame(rs_red)
check("reduced region spans 2-12", df_red$start == 2 && df_red$end == 12)

# promoters
rs_prom <- promoters(rs, upstream = 100L, downstream = 50L)
check("promoters returns RegionSet", is(rs_prom, "RegionSet"))
check("promoters same region count", length(rs_prom) == length(rs))

# trim
cs_vec <- c(chr1 = 8L)
rs_trim <- trim(rs, chromSizes = cs_vec)
check("trim returns RegionSet", is(rs_trim, "RegionSet"))
df_trim <- as.data.frame(rs_trim)
check("trim clips to chrom boundary", all(df_trim$end <= 8))

# =========================================================================
cat("\n=== 6. Binary Interval Operations ===\n")
# =========================================================================

rs_b <- RegionSet(file.path(data_dir, "dummy_b.bed"))

# concat
rs_cat <- concat(rs, rs_b)
check("concat returns RegionSet", is(rs_cat, "RegionSet"))
check("concat count is sum", length(rs_cat) == length(rs) + length(rs_b))

# union
rs_u <- union(rs, rs_b)
check("union returns RegionSet", is(rs_u, "RegionSet"))
df_u <- as.data.frame(rs_u)
check("union merges to 1 region", nrow(df_u) == 1L)
check("union spans 2-12", df_u$start == 2 && df_u$end == 12)

# setdiff
rs_sd <- setdiff(rs, rs_b)
check("setdiff returns RegionSet", is(rs_sd, "RegionSet"))

# pintersect
rs_pi <- pintersect(rs_df, rs_df)
check("pintersect returns RegionSet", is(rs_pi, "RegionSet"))

# jaccard
j <- jaccard(rs, rs_b)
check("jaccard returns numeric", is.numeric(j))
check("jaccard is 0.4", abs(j - 0.4) < 1e-10)

j_self <- jaccard(rs, rs)
check("jaccard self is 1.0", abs(j_self - 1.0) < 1e-10)

# =========================================================================
cat("\n=== 7. Consensus ===\n")
# =========================================================================

cons <- consensus(list(rs, rs_b))
check("consensus returns data.frame", is.data.frame(cons))
check("consensus has count column", "count" %in% names(cons))
check("consensus 1 region", nrow(cons) == 1L)
check("consensus count is 2", cons$count == 2)

# disjoint set
rs_far <- RegionSet(data.frame(chr = "chr2", start = 100L, end = 200L))
cons3 <- consensus(list(rs, rs_b, rs_far))
check("consensus with disjoint: 2 regions", nrow(cons3) == 2L)

# =========================================================================
cat("\n=== 8. Mixed Input Types ===\n")
# =========================================================================

# data.frame input to binary ops
j_df <- jaccard(
  data.frame(chr = "chr1", start = c(0L, 20L), end = c(10L, 30L)),
  data.frame(chr = "chr1", start = 5L, end = 25L)
)
check("jaccard with data.frames", is.numeric(j_df))

# union(df, df) dispatches to base::union, not gtars — need at least one RegionSet
rs_from_df <- RegionSet(data.frame(chr = "chr1", start = c(0L, 20L), end = c(10L, 30L)))
u_df <- union(rs_from_df, data.frame(chr = "chr1", start = 5L, end = 25L))
check("union RegionSet + data.frame", is(u_df, "RegionSet"))
df_u2 <- as.data.frame(u_df)
check("union df result: [0,30)", df_u2$start == 0 && df_u2$end == 30)

# mixed: RegionSet + data.frame
u_mixed <- union(rs, data.frame(chr = "chr1", start = 100L, end = 200L))
check("union mixed RegionSet + data.frame", is(u_mixed, "RegionSet"))

# =========================================================================
cat("\n=== 9. GRanges Integration ===\n")
# =========================================================================

if (requireNamespace("GenomicRanges", quietly = TRUE) &&
    requireNamespace("IRanges", quietly = TRUE)) {

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(3, 21), end = c(10, 30))
  )

  # Constructor from GRanges
  rs_gr <- RegionSet(gr)
  check("RegionSet from GRanges", is(rs_gr, "RegionSet"))
  check("GRanges coord conversion", length(rs_gr) == 2L)
  df_gr <- as.data.frame(rs_gr)
  check("GRanges 1-based to 0-based", df_gr$start[1] == 2L)

  # Methods on GRanges directly
  w_gr <- widths(gr)
  check("widths(GRanges) works", is.numeric(w_gr))
  check("widths match RegionSet", all(w_gr == widths(rs_gr)))

  # calcWidth drop-in with GRanges
  w_calc <- calcWidth(gr)
  check("calcWidth(GRanges) drop-in", all(w_calc == w_gr))

  # as_granges round-trip
  gr_rt <- as_granges(rs_gr)
  check("as_granges returns GRanges", inherits(gr_rt, "GRanges"))
  check("round-trip preserves regions", length(gr_rt) == 2L)

  # Binary ops with GRanges
  j_gr <- jaccard(gr, rs)
  check("jaccard(GRanges, RegionSet) works", is.numeric(j_gr))

} else {
  cat("  SKIP  GRanges tests (GenomicRanges not installed)\n")
}

# =========================================================================
cat("\n============================\n")
cat(sprintf("RESULTS: %d passed, %d failed\n", pass, fail))
if (fail > 0) {
  cat("\nFailures:\n")
  for (e in errors) cat(sprintf("  - %s\n", e))
  quit(status = 1)
} else {
  cat("All RegionSet tests passed!\n")
}
cat("============================\n")
