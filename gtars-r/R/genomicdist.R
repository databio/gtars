#' @useDynLib gtars, .registration = TRUE
NULL

# =========================================================================
# RegionSet conversion utilities
# =========================================================================

#' Convert to RegionSet
#'
#' @description Converts a GRanges object, file path, or data.frame to an
#'   ExternalPtr<RegionSet> for use with gtars genomicdist functions.
#'
#' @param x One of:
#'   - A GRanges object (requires GenomicRanges package)
#'   - A character string file path to a BED/narrowPeak/gzip file
#'   - A data.frame with chr, start, end columns (0-based half-open, BED convention)
#'   - An existing RegionSet external pointer (returned as-is)
#' @return An external pointer to a RegionSet
#'
#' @details GRanges coordinates (1-based closed) are automatically converted to
#'   BED coordinates (0-based half-open). Data.frame inputs are assumed to
#'   already be in 0-based half-open (BED) format.
#'
#' @export
as_regionset <- function(x) {
  if (inherits(x, "GRanges")) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
      stop("GenomicRanges package required for GRanges conversion")
    }
    chrs <- as.character(GenomicRanges::seqnames(x))
    starts <- BiocGenerics::start(x) - 1L  # 1-based -> 0-based
    ends <- as.integer(BiocGenerics::end(x))  # closed -> half-open (same value)
    .Call(wrap__r_regionset_from_vectors, chrs, starts, ends)
  } else if (is.character(x) && length(x) == 1) {
    .Call(wrap__r_load_regionset, x)
  } else if (is.data.frame(x)) {
    chrs <- as.character(x$chr)
    starts <- as.integer(x$start)
    ends <- as.integer(x$end)
    .Call(wrap__r_regionset_from_vectors, chrs, starts, ends)
  } else if (is(x, "externalptr")) {
    x
  } else {
    stop("Cannot convert to RegionSet: expected GRanges, file path, data.frame, or externalptr")
  }
}

#' Convert RegionSet to GRanges
#'
#' @description Converts an ExternalPtr<RegionSet> to a GRanges object.
#'
#' @param rs An external pointer to a RegionSet
#' @return A GRanges object
#'
#' @details Coordinates are converted from 0-based half-open (BED) to 1-based
#'   closed (GRanges).
#'
#' @export
as_granges <- function(rs) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges package required for GRanges conversion")
  }
  vecs <- .Call(wrap__r_regionset_to_vectors, rs)
  GenomicRanges::GRanges(
    seqnames = vecs$chr,
    ranges = IRanges::IRanges(
      start = vecs$start + 1L,  # 0-based -> 1-based
      end = vecs$end            # half-open -> closed (same value)
    )
  )
}

#' Convert RegionSet to data.frame
#'
#' @description Converts an ExternalPtr<RegionSet> to a data.frame with
#'   chr, start, end columns in 0-based half-open (BED) coordinates.
#'
#' @param rs An external pointer to a RegionSet
#' @return A data.frame with chr, start, end columns
#'
#' @export
regionset_to_df <- function(rs) {
  vecs <- .Call(wrap__r_regionset_to_vectors, rs)
  data.frame(chr = vecs$chr, start = vecs$start, end = vecs$end,
             stringsAsFactors = FALSE)
}

#' Get number of regions in a RegionSet
#'
#' @param rs An external pointer to a RegionSet
#' @return Integer count of regions
#'
#' @export
regionset_length <- function(rs) {
  .Call(wrap__r_regionset_length, rs)
}

# =========================================================================
# Helper: ensure query is a RegionSet pointer
# =========================================================================

.ensure_regionset <- function(x) {
  if (is(x, "externalptr")) return(x)
  as_regionset(x)
}

# =========================================================================
# Statistics
# =========================================================================

#' Calculate region widths
#'
#' @description Calculates the width (end - start) of each region.
#'   Port of GenomicDistributions::calcWidth.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @return Numeric vector of widths
#'
#' @export
calcWidth <- function(query) {
  query <- .ensure_regionset(query)
  .Call(wrap__r_calc_widths, query)
}

#' Calculate neighbor distances
#'
#' @description Computes distances between consecutive regions on each chromosome.
#'   Port of GenomicDistributions::calcNeighborDist.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @return Numeric vector of distances
#'
#' @export
calcNeighborDist <- function(query) {
  query <- .ensure_regionset(query)
  .Call(wrap__r_calc_neighbor_distances, query)
}

#' Calculate nearest neighbor distances
#'
#' @description Computes the distance from each region to its nearest neighbor.
#'   Port of GenomicDistributions::calcNearestNeighbors.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @return Numeric vector of nearest neighbor distances
#'
#' @export
calcNearestNeighbors <- function(query) {
  query <- .ensure_regionset(query)
  .Call(wrap__r_calc_nearest_neighbors, query)
}

#' Calculate per-chromosome statistics
#'
#' @description Computes summary statistics (region count, bounds, length
#'   min/max/mean/median) for each chromosome.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @return A data.frame with one row per chromosome
#'
#' @export
chromosomeStatistics <- function(query) {
  query <- .ensure_regionset(query)
  result <- .Call(wrap__r_chromosome_statistics, query)
  as.data.frame(result, stringsAsFactors = FALSE)
}

#' Calculate region distribution across bins
#'
#' @description Partitions the genome into fixed-size windows and counts
#'   region overlaps per bin. Output is compatible with
#'   GenomicDistributions::plotChromBins.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @param nBins Number of bins (default 250)
#' @return A data.table with chr (factor, karyotype-ordered), start, end,
#'   N, regionID, and withinGroupID columns — compatible with plotChromBins.
#'
#' @export
regionDistribution <- function(query, nBins = 250L) {
  query <- .ensure_regionset(query)
  result <- .Call(wrap__r_region_distribution, query, as.integer(nBins))
  dt <- data.table::as.data.table(result)
  data.table::setnames(dt, c("n", "rid"), c("N", "regionID"))
  data.table::set(dt, j = "withinGroupID", value = dt[["regionID"]])
  data.table::set(dt, j = "chr", value = factor(dt[["chr"]], levels = unique(dt[["chr"]])))
  dt
}

# =========================================================================
# GC Content & Dinucleotide Frequency
# =========================================================================

#' Load a genome assembly from FASTA
#'
#' @description Loads a FASTA file into memory as a GenomeAssembly pointer.
#'   This replaces the BSgenome object used in the original GenomicDistributions.
#'
#' @param fasta_path Path to a FASTA file
#' @return An external pointer to a GenomeAssembly
#'
#' @export
loadGenomeAssembly <- function(fasta_path) {
  if (!is.character(fasta_path) || length(fasta_path) != 1) {
    stop("fasta_path must be a single character string")
  }
  .Call(wrap__r_load_genome_assembly, fasta_path)
}

#' Calculate GC content
#'
#' @description Calculates GC content (proportion of G and C bases) for each
#'   region. Port of GenomicDistributions::calcGCContent.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @param ref A GenomeAssembly pointer (from loadGenomeAssembly)
#' @param ignoreUnkChroms If TRUE, skip regions on chromosomes not in the
#'   assembly (default FALSE)
#' @return Numeric vector of GC content values (0 to 1)
#'
#' @details Unlike the original GenomicDistributions which accepts a BSgenome
#'   object, this function requires a GenomeAssembly loaded from a FASTA file
#'   via \code{\link{loadGenomeAssembly}}.
#'
#' @export
calcGCContent <- function(query, ref, ignoreUnkChroms = FALSE) {
  query <- .ensure_regionset(query)
  .Call(wrap__r_calc_gc_content, query, ref, ignoreUnkChroms)
}

#' Calculate dinucleotide frequency
#'
#' @description Calculates the frequency of all 16 dinucleotides for each
#'   region as percentages. Port of GenomicDistributions::calcDinuclFreq.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @param ref A GenomeAssembly pointer (from loadGenomeAssembly)
#' @return A data.frame with a \code{region} column and 16 dinucleotide
#'   columns (AA, AC, ..., TT) containing frequency percentages.
#'   Compatible with \code{GenomicDistributions::plotDinuclFreq}.
#'
#' @export
calcDinuclFreq <- function(query, ref) {
  query <- .ensure_regionset(query)
  result <- .Call(wrap__r_calc_dinucl_freq, query, ref)
  as.data.frame(result, stringsAsFactors = FALSE)
}

# =========================================================================
# Interval Ranges
# =========================================================================

#' Trim regions to chromosome boundaries
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @param chromSizes A named integer vector of chromosome sizes
#'   (e.g., from GenomeInfoDb::seqlengths)
#' @return A RegionSet pointer with trimmed regions
#'
#' @export
gtars_trim <- function(query, chromSizes) {
  query <- .ensure_regionset(query)
  .Call(wrap__r_trim, query, names(chromSizes), as.integer(chromSizes))
}

#' Generate promoter regions
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @param upstream Number of bases upstream of start (default 2000)
#' @param downstream Number of bases downstream of start (default 200)
#' @return A RegionSet pointer with promoter regions
#'
#' @export
gtars_promoters <- function(query, upstream = 2000L, downstream = 200L) {
  query <- .ensure_regionset(query)
  .Call(wrap__r_promoters, query, as.integer(upstream), as.integer(downstream))
}

#' Merge overlapping intervals
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @return A RegionSet pointer with merged regions
#'
#' @export
gtars_reduce <- function(query) {
  query <- .ensure_regionset(query)
  .Call(wrap__r_reduce, query)
}

#' Set difference of two region sets
#'
#' @param a A GRanges object, file path, or RegionSet pointer
#' @param b A GRanges object, file path, or RegionSet pointer
#' @return A RegionSet pointer with the difference (a minus b)
#'
#' @export
gtars_setdiff <- function(a, b) {
  a <- .ensure_regionset(a)
  b <- .ensure_regionset(b)
  .Call(wrap__r_setdiff, a, b)
}

#' Pairwise intersection by position
#'
#' @param a A GRanges object, file path, or RegionSet pointer
#' @param b A GRanges object, file path, or RegionSet pointer
#' @return A RegionSet pointer with pairwise intersections
#'
#' @export
gtars_pintersect <- function(a, b) {
  a <- .ensure_regionset(a)
  b <- .ensure_regionset(b)
  .Call(wrap__r_pintersect, a, b)
}

# =========================================================================
# Partitions
# =========================================================================

#' Create a genome partition list
#'
#' @description Builds genomic partitions (promoter core, promoter proximal,
#'   3'UTR, 5'UTR, exon, intron) from gene model components.
#'   Port of GenomicDistributions::genomePartitionList.
#'
#' @param genesGR GRanges, file path, or RegionSet pointer for gene boundaries
#' @param exonsGR GRanges, file path, or RegionSet pointer for exon boundaries
#' @param threeUTRGR GRanges, file path, or RegionSet pointer for 3'UTR regions
#'   (NULL to omit)
#' @param fiveUTRGR GRanges, file path, or RegionSet pointer for 5'UTR regions
#'   (NULL to omit)
#' @param corePromSize Core promoter size in bp (default 100)
#' @param proxPromSize Proximal promoter size in bp (default 2000)
#' @param chromSizes Optional named numeric vector of chromosome sizes
#'   (e.g., from GenomeInfoDb::seqlengths). When provided, promoter regions
#'   are trimmed at chromosome boundaries before merging.
#' @return An external pointer to a PartitionList
#'
#' @export
genomePartitionList <- function(genesGR, exonsGR, threeUTRGR = NULL,
                                fiveUTRGR = NULL, corePromSize = 100L,
                                proxPromSize = 2000L, chromSizes = NULL) {
  # If any input is GRanges, extract strand and use the strand-aware path
  has_granges <- inherits(genesGR, "GRanges") || inherits(exonsGR, "GRanges") ||
    (!is.null(threeUTRGR) && inherits(threeUTRGR, "GRanges")) ||
    (!is.null(fiveUTRGR) && inherits(fiveUTRGR, "GRanges"))

  if (has_granges) {
    .extract_stranded <- function(gr) {
      if (is.null(gr)) {
        return(list(chrs = character(0), starts = integer(0),
                    ends = integer(0), strands = character(0)))
      }
      if (!inherits(gr, "GRanges")) {
        stop("When using strand-aware partitions, all inputs must be GRanges")
      }
      list(
        chrs = as.character(GenomicRanges::seqnames(gr)),
        starts = as.integer(BiocGenerics::start(gr) - 1L),
        ends = as.integer(BiocGenerics::end(gr)),
        strands = as.character(BiocGenerics::strand(gr))
      )
    }
    g <- .extract_stranded(genesGR)
    e <- .extract_stranded(exonsGR)
    t3 <- .extract_stranded(threeUTRGR)
    t5 <- .extract_stranded(fiveUTRGR)
    cs_names <- if (!is.null(chromSizes)) names(chromSizes) else character(0)
    cs_sizes <- if (!is.null(chromSizes)) as.integer(chromSizes) else integer(0)
    .Call(wrap__r_partition_list_from_regions_stranded,
          g$chrs, g$starts, g$ends, g$strands,
          e$chrs, e$starts, e$ends, e$strands,
          t3$chrs, t3$starts, t3$ends, t3$strands,
          t5$chrs, t5$starts, t5$ends, t5$strands,
          as.integer(corePromSize), as.integer(proxPromSize),
          cs_names, cs_sizes)
  } else {
    # Unstranded path: file paths, data.frames, or externalptrs
    genes <- .ensure_regionset(genesGR)
    exons <- .ensure_regionset(exonsGR)
    three_utr <- if (!is.null(threeUTRGR)) .ensure_regionset(threeUTRGR) else NULL
    five_utr <- if (!is.null(fiveUTRGR)) .ensure_regionset(fiveUTRGR) else NULL
    cs_names <- if (!is.null(chromSizes)) names(chromSizes) else character(0)
    cs_sizes <- if (!is.null(chromSizes)) as.integer(chromSizes) else integer(0)
    .Call(wrap__r_partition_list_from_regions, genes, exons, three_utr, five_utr,
          as.integer(corePromSize), as.integer(proxPromSize),
          cs_names, cs_sizes)
  }
}

#' Create a partition list from a GTF file
#'
#' @description Convenience function that loads a gene model from GTF and
#'   builds the partition list in one step.
#'
#' @param gtfPath Path to a GTF or GTF.gz file
#' @param filterProteinCoding Keep only protein-coding genes (default TRUE)
#' @param convertEnsemblUCSC Prepend "chr" to bare chromosome names (default FALSE)
#' @param corePromSize Core promoter size in bp (default 100)
#' @param proxPromSize Proximal promoter size in bp (default 2000)
#' @param chromSizes Optional named numeric vector of chromosome sizes
#'   (e.g., from GenomeInfoDb::seqlengths). When provided, promoter regions
#'   are trimmed at chromosome boundaries before merging.
#' @return An external pointer to a PartitionList
#'
#' @export
partitionListFromGTF <- function(gtfPath, filterProteinCoding = TRUE,
                                  convertEnsemblUCSC = FALSE,
                                  corePromSize = 100L, proxPromSize = 2000L,
                                  chromSizes = NULL) {
  cs_names <- if (!is.null(chromSizes)) names(chromSizes) else character(0)
  cs_sizes <- if (!is.null(chromSizes)) as.integer(chromSizes) else integer(0)
  .Call(wrap__r_partition_list_from_gtf, gtfPath, filterProteinCoding,
        convertEnsemblUCSC, as.integer(corePromSize), as.integer(proxPromSize),
        cs_names, cs_sizes)
}

#' Assign regions to genomic partitions
#'
#' @description Classifies query regions into genomic partitions.
#'   Port of GenomicDistributions::calcPartitions.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @param partitionList A PartitionList pointer (from genomePartitionList or
#'   partitionListFromGTF)
#' @param bpProportion If TRUE, count overlapping base pairs; if FALSE, count
#'   regions (default FALSE)
#' @return A data.frame with partition and Freq columns, plus a total attribute
#'
#' @export
calcPartitions <- function(query, partitionList, bpProportion = FALSE) {
  query <- .ensure_regionset(query)
  result <- .Call(wrap__r_calc_partitions, query, partitionList, bpProportion)
  df <- data.frame(partition = result$partition, Freq = result$count,
                   stringsAsFactors = FALSE)
  attr(df, "total") <- result$total
  df
}

#' Calculate expected partition enrichment
#'
#' @description Computes observed vs expected partition overlaps with
#'   chi-square p-values. Port of GenomicDistributions::calcExpectedPartitions.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @param partitionList A PartitionList pointer
#' @param genomeSize A named numeric vector of chromosome sizes
#'   (e.g., from GenomeInfoDb::seqlengths). Required.
#' @param bpProportion If TRUE, count base pairs; if FALSE, count regions
#'   (default FALSE)
#' @return A data.frame with partition, observed, expected, log10OE, pvalue columns
#'
#' @export
calcExpectedPartitions <- function(query, partitionList, genomeSize,
                                    bpProportion = FALSE) {
  query <- .ensure_regionset(query)
  result <- .Call(wrap__r_calc_expected_partitions, query, partitionList,
                  names(genomeSize), as.integer(genomeSize), bpProportion)
  data.frame(
    partition = result$partition,
    observed = result$observed,
    expected = result$expected,
    log10OE = result$log10OE,
    pvalue = result$pvalue,
    stringsAsFactors = FALSE
  )
}

# =========================================================================
# Signal
# =========================================================================

#' Calculate summary signal
#'
#' @description Overlaps query regions with a signal matrix and computes
#'   per-condition summary statistics. Port of GenomicDistributions::calcSummarySignal.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @param signalMatrix A data.frame or data.table where the first column contains
#'   region IDs in "chr_start_end" format and remaining columns contain signal
#'   values for each condition
#' @return A list compatible with GenomicDistributions' plotSummarySignal:
#'   \item{signalSummaryMatrix}{data.table with queryPeak column + one column per condition}
#'   \item{matrixStats}{data.frame with stats as rownames and conditions as columns}
#'
#' @export
calcSummarySignal <- function(query, signalMatrix) {
  query <- .ensure_regionset(query)

  signalMatrix <- as.data.frame(signalMatrix)
  region_ids <- as.character(signalMatrix[[1]])
  condition_names <- colnames(signalMatrix)[-1]
  values <- as.matrix(signalMatrix[, -1, drop = FALSE])

  n_regions <- nrow(values)
  n_conditions <- ncol(values)

  # Flatten to row-major vector for Rust
  values_flat <- as.numeric(t(values))

  result <- .Call(wrap__r_calc_summary_signal, query, region_ids,
                  condition_names, values_flat,
                  as.integer(n_regions), as.integer(n_conditions))

  # Reshape to match GD's plotSummarySignal expected format:
  # signalSummaryMatrix: wide data.table with queryPeak + condition columns
  signal_mat <- do.call(cbind, result$signal_matrix)
  colnames(signal_mat) <- result$condition_names
  ssm <- data.table::as.data.table(
    cbind(data.frame(queryPeak = result$region_labels,
                     stringsAsFactors = FALSE),
          as.data.frame(signal_mat))
  )

  # matrixStats: wide data.frame with stats as rownames, conditions as columns
  ms_long <- as.data.frame(result$matrixStats, stringsAsFactors = FALSE)
  stat_names <- c("lowerWhisker", "lowerHinge", "median", "upperHinge", "upperWhisker")
  stat_cols <- c("lower_whisker", "lower_hinge", "median", "upper_hinge", "upper_whisker")
  ms_wide <- as.data.frame(
    t(as.matrix(ms_long[, stat_cols])),
    stringsAsFactors = FALSE
  )
  colnames(ms_wide) <- ms_long$condition
  rownames(ms_wide) <- stat_names

  list(signalSummaryMatrix = ssm, matrixStats = ms_wide)
}

# =========================================================================
# TSS / Feature Distances
# =========================================================================

#' Calculate feature distances (signed)
#'
#' @description Calculates the signed distance from each query region to its
#'   nearest feature. Positive = feature is downstream, negative = upstream.
#'   Port of GenomicDistributions::calcFeatureDist.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @param features A GRanges object, file path, or RegionSet pointer of
#'   feature locations (e.g., TSS sites)
#' @return Numeric vector of signed distances
#'
#' @export
calcFeatureDist <- function(query, features) {
  query <- .ensure_regionset(query)
  features <- .ensure_regionset(features)
  .Call(wrap__r_calc_feature_distances, query, features)
}

#' Calculate TSS distances (absolute)
#'
#' @description Calculates the absolute distance from each query region midpoint
#'   to the nearest feature midpoint.
#'
#' @param query A GRanges object, file path, or RegionSet pointer
#' @param features A GRanges object, file path, or RegionSet pointer of
#'   feature locations (e.g., TSS sites)
#' @return Integer vector of absolute distances
#'
#' @export
calcTSSDist <- function(query, features) {
  query <- .ensure_regionset(query)
  features <- .ensure_regionset(features)
  .Call(wrap__r_calc_tss_distances, query, features)
}
