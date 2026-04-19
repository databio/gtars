#' @include RegionSet-methods.R
#' @useDynLib gtars, .registration = TRUE
NULL

# =========================================================================
# GenomicDistributions drop-in wrappers
#
# These delegate to S4 methods on RegionSet, accepting both GRanges and
# RegionSet inputs for backwards compatibility.
# =========================================================================

#' Calculate region widths
#'
#' @description Drop-in replacement for GenomicDistributions::calcWidth.
#'   Delegates to \code{\link{widths}}.
#'
#' @param query A GRanges object, file path, data.frame, or RegionSet
#' @return Numeric vector of widths
#'
#' @export
calcWidth <- function(query) widths(query)

#' Calculate neighbor distances
#'
#' @description Drop-in replacement for GenomicDistributions::calcNeighborDist.
#'   Delegates to \code{\link{neighborDistances}}.
#'
#' @param query A GRanges object, file path, data.frame, or RegionSet
#' @return Numeric vector of distances
#'
#' @export
calcNeighborDist <- function(query) neighborDistances(query)

#' Calculate nearest neighbor distances
#'
#' @description Drop-in replacement for GenomicDistributions::calcNearestNeighbors.
#'   Delegates to \code{\link{nearestNeighbors}}.
#'
#' @param query A GRanges object, file path, data.frame, or RegionSet
#' @return Numeric vector of nearest neighbor distances
#'
#' @export
calcNearestNeighbors <- function(query) nearestNeighbors(query)

#' Calculate region distribution across bins
#'
#' @description Drop-in replacement for GenomicDistributions::calcChromBins.
#'   Delegates to \code{\link{distribution}}.
#'
#' @param query A GRanges object, file path, data.frame, or RegionSet
#' @param nBins Number of bins (default 250)
#' @return A data.table compatible with plotChromBins
#'
#' @export
regionDistribution <- function(query, nBins = 250L, chromSizes = NULL) {
  distribution(query, nBins = nBins, chromSizes = chromSizes)
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
#'   region. Drop-in replacement for GenomicDistributions::calcGCContent.
#'
#' @param query A GRanges object, file path, data.frame, or RegionSet
#' @param ref A GenomeAssembly pointer (from loadGenomeAssembly)
#' @param ignoreUnkChroms If TRUE, skip regions on chromosomes not in the
#'   assembly (default FALSE)
#' @return Numeric vector of GC content values (0 to 1)
#'
#' @export
calcGCContent <- function(query, ref, ignoreUnkChroms = FALSE) {
  query <- .ensure_regionset(query)
  .Call(wrap__r_calc_gc_content, query, ref, ignoreUnkChroms)
}

#' Calculate dinucleotide frequency
#'
#' @description Calculates the frequency of all 16 dinucleotides for each
#'   region. Drop-in replacement for GenomicDistributions::calcDinuclFreq,
#'   including the \code{rawCounts} parameter.
#'
#' @param query A GRanges object, file path, data.frame, or RegionSet
#' @param ref A GenomeAssembly pointer (from loadGenomeAssembly)
#' @param rawCounts If TRUE, return raw integer counts; if FALSE (default),
#'   return percentages per row (each row sums to 100)
#' @param ignoreUnknownChroms If TRUE, skip regions on chromosomes not in
#'   the assembly. If FALSE (default), error on unknown chromosomes.
#' @return A data.frame with a \code{region} column and 16 dinucleotide
#'   columns. For pooled global counts, call \code{colSums(result[, -1])}
#'   on a \code{rawCounts=TRUE} result.
#'
#' @export
calcDinuclFreq <- function(query, ref, rawCounts = FALSE, ignoreUnknownChroms = FALSE) {
  query <- .ensure_regionset(query)
  result <- .Call(
    wrap__r_calc_dinucl_freq,
    query,
    ref,
    rawCounts,
    ignoreUnknownChroms
  )
  as.data.frame(result, stringsAsFactors = FALSE)
}

# =========================================================================
# Partitions
# =========================================================================

#' Create a genome partition list
#'
#' @description Builds genomic partitions (promoter core, promoter proximal,
#'   3'UTR, 5'UTR, exon, intron) from gene model components.
#'   Drop-in replacement for GenomicDistributions::genomePartitionList.
#'
#' @param genesGR GRanges, file path, data.frame, or RegionSet for gene boundaries
#' @param exonsGR GRanges, file path, data.frame, or RegionSet for exon boundaries
#' @param threeUTRGR GRanges, file path, data.frame, or RegionSet for 3'UTR regions
#'   (NULL to omit)
#' @param fiveUTRGR GRanges, file path, data.frame, or RegionSet for 5'UTR regions
#'   (NULL to omit)
#' @param corePromSize Core promoter size in bp (default 100, matching R
#'   GenomicDistributions)
#' @param proxPromSize Proximal promoter size in bp (default 2000, matching R
#'   GenomicDistributions)
#' @param chromSizes Optional named numeric vector of chromosome sizes.
#'   When provided, promoter regions are trimmed at chromosome boundaries.
#' @return An external pointer to a PartitionList
#'
#' @export
genomePartitionList <- function(genesGR, exonsGR, threeUTRGR = NULL,
                                fiveUTRGR = NULL, corePromSize = 100L,
                                proxPromSize = 2000L, chromSizes = NULL) {
  # Extract chr/start/end/strand from any supported input type
  .extract <- function(x) {
    if (is.null(x)) {
      return(list(chrs = character(0), starts = integer(0),
                  ends = integer(0), strands = character(0)))
    }
    rs <- RegionSet(x)
    vecs <- .Call(wrap__r_regionset_to_vectors, .ptr(rs))
    list(chrs = vecs$chr, starts = vecs$start, ends = vecs$end,
         strands = rs@strand)
  }

  g <- .extract(genesGR)
  e <- .extract(exonsGR)
  t3 <- .extract(threeUTRGR)
  t5 <- .extract(fiveUTRGR)

  # Check if any input has real strand info
  has_strand <- any(g$strands != "*") || any(e$strands != "*") ||
    any(t3$strands != "*") || any(t5$strands != "*")

  cs_names <- if (!is.null(chromSizes)) names(chromSizes) else character(0)
  cs_sizes <- if (!is.null(chromSizes)) as.integer(chromSizes) else integer(0)

  if (has_strand) {
    .Call(wrap__r_partition_list_from_regions_stranded,
          g$chrs, as.integer(g$starts), as.integer(g$ends), g$strands,
          e$chrs, as.integer(e$starts), as.integer(e$ends), e$strands,
          t3$chrs, as.integer(t3$starts), as.integer(t3$ends), t3$strands,
          t5$chrs, as.integer(t5$starts), as.integer(t5$ends), t5$strands,
          as.integer(corePromSize), as.integer(proxPromSize),
          cs_names, cs_sizes)
  } else {
    .Call(wrap__r_partition_list_from_regions,
          .ptr(RegionSet(genesGR)), .ptr(RegionSet(exonsGR)),
          if (!is.null(threeUTRGR)) .ptr(RegionSet(threeUTRGR)) else NULL,
          if (!is.null(fiveUTRGR)) .ptr(RegionSet(fiveUTRGR)) else NULL,
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
#' @param chromSizes Optional named numeric vector of chromosome sizes.
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
#'   Drop-in replacement for GenomicDistributions::calcPartitions.
#'
#' @param query A GRanges object, file path, data.frame, or RegionSet
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
#'   chi-square p-values. Drop-in replacement for
#'   GenomicDistributions::calcExpectedPartitions.
#'
#' @param query A GRanges object, file path, data.frame, or RegionSet
#' @param partitionList A PartitionList pointer
#' @param genomeSize A named numeric vector of chromosome sizes
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
#'   per-condition summary statistics. Drop-in replacement for
#'   GenomicDistributions::calcSummarySignal.
#'
#' @param query A GRanges object, file path, data.frame, or RegionSet
#' @param signalMatrix A data.frame or data.table where the first column contains
#'   region IDs in "chr_start_end" format and remaining columns contain signal
#'   values for each condition
#' @return A list compatible with GenomicDistributions' plotSummarySignal
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
  values_flat <- as.numeric(t(values))

  result <- .Call(wrap__r_calc_summary_signal, query, region_ids,
                  condition_names, values_flat,
                  as.integer(n_regions), as.integer(n_conditions))

  signal_mat <- do.call(cbind, result$signal_matrix)
  colnames(signal_mat) <- result$condition_names
  ssm <- data.table::as.data.table(
    cbind(data.frame(queryPeak = result$region_labels,
                     stringsAsFactors = FALSE),
          as.data.frame(signal_mat))
  )

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
#'   nearest feature. Drop-in replacement for
#'   GenomicDistributions::calcFeatureDist.
#'
#' @param query A GRanges object, file path, data.frame, or RegionSet
#' @param features A GRanges object, file path, data.frame, or RegionSet
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
#' @param query A GRanges object, file path, data.frame, or RegionSet
#' @param features A GRanges object, file path, data.frame, or RegionSet
#' @return Integer vector of absolute distances
#'
#' @export
calcTSSDist <- function(query, features) {
  query <- .ensure_regionset(query)
  features <- .ensure_regionset(features)
  .Call(wrap__r_calc_tss_distances, query, features)
}
