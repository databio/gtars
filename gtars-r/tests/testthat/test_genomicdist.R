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
# Width, distance, and distribution wrappers
# =========================================================================

test_that("calcWidth returns correct widths", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  w <- calcWidth(rs)
  # doubles instead of integers: u32 widths can exceed i32::MAX (2.1 Gbp)
  expect_type(w, "double")
  expect_equal(length(w), 4L)
  expect_equal(w[1], 4)
})

test_that("calcNeighborDist returns distances", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  nd <- calcNeighborDist(rs)
  expect_true(is.numeric(nd))
})

test_that("calcNearestNeighbors returns distances", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  nn <- calcNearestNeighbors(rs)
  # doubles instead of integers: distances can be chromosome-scale,
  # approaching or exceeding i32::MAX for concatenated/genome-scale regions
  expect_type(nn, "double")
  expect_equal(length(nn), length(rs))
})

test_that("regionDistribution returns data.table", {
  skip_if_no_data()
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  dist <- regionDistribution(rs, nBins = 10L)
  expect_s3_class(dist, "data.table")
  expect_true("N" %in% names(dist))
  expect_true("chr" %in% names(dist))
})

# =========================================================================
# GC content and dinucleotide frequency
# =========================================================================

test_that("loadGenomeAssembly and calcGCContent work", {
  fasta_path <- file.path(test_path(), "..", "..", "..", "tests", "data", "fasta", "base.fa")
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")

  assembly <- loadGenomeAssembly(fasta_path)
  expect_type(assembly, "externalptr")

  # base.fa: chr1 = "GGAA" (4bp), chrX = "TTGGGGAA" (8bp)
  df <- data.frame(chr = "chrX", start = 0L, end = 4L)
  gc <- calcGCContent(df, assembly, ignoreUnkChroms = TRUE)
  expect_type(gc, "double")
  expect_equal(length(gc), 1L)
  expect_true(gc[1] >= 0 && gc[1] <= 1)
})

test_that("calcDinuclFreq returns data.frame with 16 dinucleotide columns", {
  fasta_path <- file.path(test_path(), "..", "..", "..", "tests", "data", "fasta", "base.fa")
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")

  assembly <- loadGenomeAssembly(fasta_path)
  # base.fa: chrX = "TTGGGGAA" (8bp) — need at least 2bp for a dinucleotide
  df <- data.frame(chr = "chrX", start = 0L, end = 8L)
  freq <- calcDinuclFreq(df, assembly)
  expect_true(is.data.frame(freq))
  # 16 dinucleotides + region column
  expect_true(ncol(freq) >= 16)
  # Default rawCounts=FALSE: each row sums to 100
  num_cols <- freq[, sapply(freq, is.numeric)]
  expect_equal(sum(num_cols[1, ]), 100, tolerance = 1e-6)
})

test_that("calcDinuclFreq rawCounts=TRUE returns integer counts", {
  fasta_path <- file.path(test_path(), "..", "..", "..", "tests", "data", "fasta", "base.fa")
  skip_if_not(file.exists(fasta_path), "Test FASTA file not found")

  assembly <- loadGenomeAssembly(fasta_path)
  # base.fa: chrX = "TTGGGGAA" (8bp) → 7 dinucleotide windows
  df <- data.frame(chr = "chrX", start = 0L, end = 8L)
  freq <- calcDinuclFreq(df, assembly, rawCounts = TRUE)
  expect_true(is.data.frame(freq))
  num_cols <- freq[, sapply(freq, is.numeric)]
  # 8 bases → 7 dinucleotides
  expect_equal(sum(num_cols[1, ]), 7)
})

# =========================================================================
# Partitions
# =========================================================================

test_that("genomePartitionList from BED files works", {
  skip_if_no_data()
  genes_path <- file.path(data_dir, "test_genes.bed")
  exons_path <- file.path(data_dir, "test_exons.bed")
  skip_if_not(file.exists(genes_path), "test_genes.bed not found")
  skip_if_not(file.exists(exons_path), "test_exons.bed not found")

  genes <- RegionSet(genes_path)
  exons <- RegionSet(exons_path)
  pl <- genomePartitionList(genes, exons)
  expect_type(pl, "externalptr")
})

test_that("partitionListFromGTF works", {
  skip_if_no_data()
  gtf_path <- file.path(data_dir, "test_gene_model.gtf")
  skip_if_not(file.exists(gtf_path), "test_gene_model.gtf not found")

  pl <- partitionListFromGTF(gtf_path, filterProteinCoding = FALSE)
  expect_type(pl, "externalptr")
})

test_that("calcPartitions returns partition table", {
  skip_if_no_data()
  gtf_path <- file.path(data_dir, "test_gene_model.gtf")
  skip_if_not(file.exists(gtf_path), "test_gene_model.gtf not found")

  pl <- partitionListFromGTF(gtf_path, filterProteinCoding = FALSE)
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  result <- calcPartitions(rs, pl)
  expect_true(is.data.frame(result))
  expect_true("partition" %in% names(result))
  expect_true("Freq" %in% names(result))
})

test_that("calcExpectedPartitions returns enrichment table", {
  skip_if_no_data()
  gtf_path <- file.path(data_dir, "test_gene_model.gtf")
  chrom_sizes_path <- file.path(data_dir, "dummy_chrom_sizes")
  skip_if_not(file.exists(gtf_path), "test_gene_model.gtf not found")
  skip_if_not(file.exists(chrom_sizes_path), "dummy_chrom_sizes not found")

  # Read chrom sizes
  cs <- read.table(chrom_sizes_path, header = FALSE, stringsAsFactors = FALSE)
  chrom_sizes <- setNames(as.integer(cs$V2), cs$V1)

  pl <- partitionListFromGTF(gtf_path, filterProteinCoding = FALSE,
                              chromSizes = chrom_sizes)
  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  result <- calcExpectedPartitions(rs, pl, chrom_sizes)
  expect_true(is.data.frame(result))
  expect_true(all(c("partition", "observed", "expected", "log10OE", "pvalue") %in% names(result)))
})

# =========================================================================
# TSS / Feature distances
# =========================================================================

test_that("calcFeatureDist returns signed distances", {
  skip_if_no_data()
  tss_path <- file.path(data_dir, "dummy_tss.bed")
  skip_if_not(file.exists(tss_path), "dummy_tss.bed not found")

  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  tss <- RegionSet(tss_path)
  dists <- calcFeatureDist(rs, tss)
  expect_true(is.numeric(dists))
  expect_equal(length(dists), length(rs))
})

test_that("calcTSSDist returns absolute distances", {
  skip_if_no_data()
  tss_path <- file.path(data_dir, "dummy_tss.bed")
  skip_if_not(file.exists(tss_path), "dummy_tss.bed not found")

  rs <- RegionSet(file.path(data_dir, "dummy.bed"))
  tss <- RegionSet(tss_path)
  dists <- calcTSSDist(rs, tss)
  expect_true(is.numeric(dists))
  expect_equal(length(dists), length(rs))
  expect_true(all(dists >= 0))
})
