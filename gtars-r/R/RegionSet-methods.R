#' @include RegionSet-class.R
#' @importFrom BiocGenerics union setdiff intersect
#' @importFrom IRanges reduce promoters trim shift narrow resize flank disjoin gaps findOverlaps countOverlaps
NULL

# Helper: coerce-and-dispatch for Bioconductor generics.
# We register methods for "character" and "data.frame" rather than
# "ANY", to avoid hijacking dispatch for GRanges/IRanges objects
# (which have their own methods on these generics).
.coerce_classes <- c("character", "data.frame")

# =========================================================================
# Own generics: unary statistics
# =========================================================================

#' Region widths
#'
#' @description Calculates the width (end - start) of each region.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param ... ignored
#' @return Numeric vector of widths
#' @export
setGeneric("widths", function(x, ...) standardGeneric("widths"))

#' @rdname widths
#' @export
setMethod("widths", "RegionSet", function(x, ...) {
  .Call(wrap__r_calc_widths, .ptr(x))
})

#' @rdname widths
#' @export
setMethod("widths", "ANY", function(x, ...) {
  widths(RegionSet(x))
})

#' Neighbor distances
#'
#' @description Computes distances between consecutive regions on each
#'   chromosome.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param ... ignored
#' @return Numeric vector of distances
#' @export
setGeneric("neighborDistances", function(x, ...) standardGeneric("neighborDistances"))

#' @rdname neighborDistances
#' @export
setMethod("neighborDistances", "RegionSet", function(x, ...) {
  .Call(wrap__r_calc_neighbor_distances, .ptr(x))
})

#' @rdname neighborDistances
#' @export
setMethod("neighborDistances", "ANY", function(x, ...) {
  neighborDistances(RegionSet(x))
})

#' Nearest neighbor distances
#'
#' @description Computes the distance from each region to its nearest
#'   neighbor.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param ... ignored
#' @return Numeric vector of nearest neighbor distances
#' @export
setGeneric("nearestNeighbors", function(x, ...) standardGeneric("nearestNeighbors"))

#' @rdname nearestNeighbors
#' @export
setMethod("nearestNeighbors", "RegionSet", function(x, ...) {
  .Call(wrap__r_calc_nearest_neighbors, .ptr(x))
})

#' @rdname nearestNeighbors
#' @export
setMethod("nearestNeighbors", "ANY", function(x, ...) {
  nearestNeighbors(RegionSet(x))
})

#' Per-chromosome statistics
#'
#' @description Computes summary statistics (region count, bounds, length
#'   min/max/mean/median) for each chromosome.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param ... ignored
#' @return A data.frame with one row per chromosome
#' @export
setGeneric("chromosomeStatistics", function(x, ...) standardGeneric("chromosomeStatistics"))

#' @rdname chromosomeStatistics
#' @export
setMethod("chromosomeStatistics", "RegionSet", function(x, ...) {
  result <- .Call(wrap__r_chromosome_statistics, .ptr(x))
  as.data.frame(result, stringsAsFactors = FALSE)
})

#' @rdname chromosomeStatistics
#' @export
setMethod("chromosomeStatistics", "ANY", function(x, ...) {
  chromosomeStatistics(RegionSet(x))
})

#' Region distribution across bins
#'
#' @description Partitions the genome into fixed-size windows and counts
#'   region overlaps per bin. Output is compatible with
#'   GenomicDistributions::plotChromBins.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param nBins Number of bins (default 250)
#' @param ... ignored
#' @return A data.table compatible with plotChromBins
#' @export
setGeneric("distribution", function(x, nBins = 250L, chromSizes = NULL, ...) standardGeneric("distribution"))

#' @rdname distribution
#' @export
setMethod("distribution", "RegionSet", function(x, nBins = 250L, chromSizes = NULL, ...) {
  chrom_names <- NULL
  chrom_lengths <- NULL
  if (!is.null(chromSizes)) {
    chrom_names <- names(chromSizes)
    chrom_lengths <- as.numeric(chromSizes)
  }
  result <- .Call(wrap__r_region_distribution, .ptr(x), as.integer(nBins), chrom_names, chrom_lengths)
  dt <- data.table::as.data.table(result)
  data.table::setnames(dt, c("n", "rid"), c("N", "regionID"))
  data.table::set(dt, j = "withinGroupID", value = dt[["regionID"]])
  data.table::set(dt, j = "chr",
                  value = factor(dt[["chr"]], levels = unique(dt[["chr"]])))
  dt
})

#' @rdname distribution
#' @export
setMethod("distribution", "ANY", function(x, nBins = 250L, chromSizes = NULL, ...) {
  distribution(RegionSet(x), nBins = nBins, chromSizes = chromSizes)
})

#' Cluster regions via single-linkage with a stitching radius
#'
#' @description Assigns a cluster ID to each region. Regions within
#'   ``maxGap`` bp of another region on the same chromosome are linked
#'   into the same cluster; chromosome boundaries always break clusters.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param maxGap Maximum bp gap between regions to link (default 0)
#' @param ... ignored
#' @return Integer vector of 0-based cluster IDs in original region order
#' @export
setGeneric("clusterRegions", function(x, maxGap = 0L, ...) standardGeneric("clusterRegions"))

#' @rdname clusterRegions
#' @export
setMethod("clusterRegions", "RegionSet", function(x, maxGap = 0L, ...) {
  .Call(wrap__r_cluster, .ptr(x), as.integer(maxGap))
})

#' @rdname clusterRegions
#' @export
setMethod("clusterRegions", "ANY", function(x, maxGap = 0L, ...) {
  clusterRegions(RegionSet(x), maxGap = maxGap)
})

#' Inter-peak spacing summary statistics
#'
#' @description Scalar summary over the distribution of inter-region
#'   spacings: mean, median, standard deviation, IQR, and log-space
#'   mean / std. Overlapping and abutting neighbors are excluded, and
#'   cross-chromosome pairs are never counted.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param ... ignored
#' @return A named list with fields: n_gaps, mean, median, std, iqr,
#'   log_mean, log_std. Float fields are NaN when ``n_gaps`` is 0.
#' @export
setGeneric("interPeakSpacing", function(x, ...) standardGeneric("interPeakSpacing"))

#' @rdname interPeakSpacing
#' @export
setMethod("interPeakSpacing", "RegionSet", function(x, ...) {
  .Call(wrap__r_calc_inter_peak_spacing, .ptr(x))
})

#' @rdname interPeakSpacing
#' @export
setMethod("interPeakSpacing", "ANY", function(x, ...) {
  interPeakSpacing(RegionSet(x))
})

#' Peak cluster summary statistics
#'
#' @description Cluster-level summary at a given stitching radius:
#'   number of clusters, clustered-peak count, mean / max cluster size,
#'   and fraction of peaks belonging to a cluster.
#'
#' @details
#' \code{min_cluster_size} applies uniformly to every size-dependent
#' field except \code{max_cluster_size} (which is always the biggest
#' cluster in the input regardless of filter). \code{n_clusters},
#' \code{n_clustered_peaks}, \code{mean_cluster_size}, and
#' \code{fraction_clustered} all restrict to clusters with size at
#' least \code{min_cluster_size}. The arithmetic identity
#' \code{n_clusters * mean_cluster_size == n_clustered_peaks} holds at
#' any threshold.
#'
#' \strong{Default \code{min_cluster_size = 2L}}: every reported field
#' describes "clusters with at least 2 peaks". This matches the typical
#' enhancer-clustering / super-enhancer-stitching use case;
#' \code{fraction_clustered} is then the fraction of peaks with at
#' least one neighbor within \code{radius_bp}, and
#' \code{mean_cluster_size} is the average size of multi-peak clusters.
#'
#' \strong{Pass \code{min_cluster_size = 1L}} to include singletons. Under
#' this threshold \code{mean_cluster_size} becomes the simple
#' \code{total_peaks / n_clusters} average, but \code{n_clustered_peaks}
#' degenerates to \code{total_peaks} and \code{fraction_clustered} to
#' \code{1.0} (every peak is in a cluster of size >= 1).
#'
#' Higher values progressively restrict to larger clusters.
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param radius_bp Stitching radius in bp (non-negative)
#' @param min_cluster_size Minimum cluster size that qualifies clusters
#'   for inclusion in all size-dependent fields except
#'   \code{max_cluster_size}. Default \code{2L}.
#' @param ... ignored
#' @return A named list with fields: radius_bp, n_clusters,
#'   n_clustered_peaks, mean_cluster_size, max_cluster_size,
#'   fraction_clustered.
#' @export
setGeneric("peakClusters",
           function(x, radius_bp, min_cluster_size = 2L, ...) standardGeneric("peakClusters"))

#' @rdname peakClusters
#' @export
setMethod("peakClusters", "RegionSet", function(x, radius_bp, min_cluster_size = 2L, ...) {
  .Call(
    wrap__r_calc_peak_clusters,
    .ptr(x),
    as.integer(radius_bp),
    as.integer(min_cluster_size)
  )
})

#' @rdname peakClusters
#' @export
setMethod("peakClusters", "ANY", function(x, radius_bp, min_cluster_size = 2L, ...) {
  peakClusters(RegionSet(x), radius_bp = radius_bp, min_cluster_size = min_cluster_size)
})

#' Dense per-window peak count vector
#'
#' @description Returns the full zero-padded peak count vector across the
#'   genome, ordered by karyotypic chromosome order and bin index. Unlike
#'   \code{distribution}, empty bins are included. Suitable for ML feature
#'   extraction.
#'
#' @details
#' \strong{\code{nBins} is a target, not the total window count.} It is
#' the target bin count for the longest chromosome in \code{chrom_sizes},
#' not the length of the returned \code{counts} vector. Bin width is
#' derived as \code{max(chrom_sizes) / nBins} (floored, minimum 1 bp),
#' and every chromosome is tiled with windows of that width. The total
#' number of bins returned is \code{sum(ceil(size / bin_width))} across
#' \code{chrom_sizes}, which can substantially exceed \code{nBins} when
#' many chromosomes are present. To target a specific bin width in bp,
#' pass \code{nBins = floor(max(chrom_sizes) / desired_width_bp)}.
#'
#' \strong{Per-chromosome bin width.} The last bin on each chromosome is
#' narrower than \code{bin_width} whenever \code{chrom_size} is not an
#' exact multiple of \code{bin_width}. Chromosomes shorter than
#' \code{bin_width} (common with UCSC alt / random / unplaced contigs
#' such as \code{chrUn_*}, \code{*_random}, \code{*_alt}) reduce to a
#' single bin whose effective width equals the chromosome size rather
#' than \code{bin_width}. Entries in \code{counts} are therefore counts
#' per bin, not counts per \code{bin_width} bp — bins of different
#' effective widths are not directly comparable as densities when
#' \code{chrom_sizes} contains contigs significantly shorter than
#' \code{bin_width}.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param chrom_sizes Named integer vector of chromosome sizes (required)
#' @param nBins Target bin count for the longest chromosome — see Details.
#' @param ... ignored
#' @return A named list with fields: n_bins, bin_width, counts (numeric
#'   vector), chrom_offset_names, chrom_offset_indices.
#' @export
setGeneric("densityVector",
           function(x, chrom_sizes, nBins = 250L, ...) standardGeneric("densityVector"))

#' @rdname densityVector
#' @export
setMethod("densityVector", "RegionSet", function(x, chrom_sizes, nBins = 250L, ...) {
  if (missing(chrom_sizes) || is.null(chrom_sizes)) {
    stop("'chrom_sizes' is required: pass a named integer vector of chromosome sizes")
  }
  .Call(
    wrap__r_calc_density_vector,
    .ptr(x),
    as.integer(nBins),
    names(chrom_sizes),
    as.integer(chrom_sizes)
  )
})

#' @rdname densityVector
#' @export
setMethod("densityVector", "ANY", function(x, chrom_sizes, nBins = 250L, ...) {
  densityVector(RegionSet(x), chrom_sizes = chrom_sizes, nBins = nBins)
})

#' Density homogeneity summary statistics
#'
#' @description Summary of how evenly peaks are distributed across genome
#'   windows: mean, population variance, coefficient of variation, Gini
#'   coefficient, and count of nonzero windows.
#'
#' @details
#' See \code{\link{densityVector}} for the definition of \code{nBins}
#' (target bin count for the longest chromosome, not the total window
#' count) and for the treatment of chromosomes shorter than the derived
#' bin width. Both affect the interpretation of the statistics below —
#' short contigs in \code{chrom_sizes} each contribute a narrow
#' single-bin entry which dilutes \code{mean_count}, inflates
#' \code{n_windows}, and raises \code{gini}.
#'
#' The Gini coefficient is biased high for very sparse count
#' distributions (many zero-count windows). Check \code{n_nonzero_windows}
#' before interpreting Gini on sparse peak sets.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param chrom_sizes Named integer vector of chromosome sizes (required)
#' @param nBins Target bin count for the longest chromosome — see
#'   \code{\link{densityVector}} for the full semantic.
#' @param ... ignored
#' @return A named list with fields: bin_width, n_windows,
#'   n_nonzero_windows, mean_count, variance, cv, gini.
#' @export
setGeneric("densityHomogeneity",
           function(x, chrom_sizes, nBins = 250L, ...) standardGeneric("densityHomogeneity"))

#' @rdname densityHomogeneity
#' @export
setMethod("densityHomogeneity", "RegionSet", function(x, chrom_sizes, nBins = 250L, ...) {
  if (missing(chrom_sizes) || is.null(chrom_sizes)) {
    stop("'chrom_sizes' is required: pass a named integer vector of chromosome sizes")
  }
  .Call(
    wrap__r_calc_density_homogeneity,
    .ptr(x),
    as.integer(nBins),
    names(chrom_sizes),
    as.integer(chrom_sizes)
  )
})

#' @rdname densityHomogeneity
#' @export
setMethod("densityHomogeneity", "ANY", function(x, chrom_sizes, nBins = 250L, ...) {
  densityHomogeneity(RegionSet(x), chrom_sizes = chrom_sizes, nBins = nBins)
})

# =========================================================================
# IRanges generics: trim, reduce, promoters
#
# Imported from IRanges so that these generics dispatch to gtars for
# RegionSet and to IRanges for GRanges.
# =========================================================================

#' Trim regions to chromosome boundaries
#'
#' @description Clips regions to fit within chromosome boundaries. Regions
#'   on unknown chromosomes are dropped.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param ... Must include \code{chromSizes}: a named integer vector of
#'   chromosome sizes
#' @return A RegionSet with trimmed regions
#' @rdname trim
#' @export
setMethod("trim", "RegionSet", function(x, chromSizes, ...) {
  # Trim may drop regions, so strand mapping is lost
  .rs_from_ptr(.Call(wrap__r_trim, .ptr(x),
                     names(chromSizes), as.integer(chromSizes)))
})

#' @rdname trim
#' @export
setMethod("trim", "character", function(x, chromSizes, ...) {
  trim(RegionSet(x), chromSizes = chromSizes)
})

#' @rdname trim
#' @export
setMethod("trim", "data.frame", function(x, chromSizes, ...) {
  trim(RegionSet(x), chromSizes = chromSizes)
})

#' Merge overlapping and adjacent intervals
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param ... ignored
#' @return A RegionSet with merged regions
#' @rdname reduce
#' @export
setMethod("reduce", "RegionSet", function(x, ...) {
  .rs_from_ptr(.Call(wrap__r_reduce, .ptr(x)))
})

#' @rdname reduce
#' @export
setMethod("reduce", "character", function(x, ...) {
  reduce(RegionSet(x))
})

#' @rdname reduce
#' @export
setMethod("reduce", "data.frame", function(x, ...) {
  reduce(RegionSet(x))
})

#' Generate promoter regions
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param upstream Bases upstream of start (default 2000)
#' @param downstream Bases downstream of start (default 200)
#' @param ... ignored
#' @return A RegionSet with promoter regions
#' @rdname promoters
#' @export
setMethod("promoters", "RegionSet", function(x, upstream = 2000L,
                                              downstream = 200L, ...) {
  ptr <- .Call(wrap__r_promoters, .ptr(x), as.integer(upstream),
               as.integer(downstream))
  # Same regions in same order — preserve strand
  new("RegionSet", ptr = ptr, strand = x@strand)
})

#' @rdname promoters
#' @export
setMethod("promoters", "character", function(x, upstream = 2000L,
                                              downstream = 200L, ...) {
  promoters(RegionSet(x), upstream = upstream,
            downstream = downstream)
})

#' @rdname promoters
#' @export
setMethod("promoters", "data.frame", function(x, upstream = 2000L,
                                               downstream = 200L, ...) {
  promoters(RegionSet(x), upstream = upstream,
            downstream = downstream)
})

#' Shift regions by a fixed offset
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param shift Integer offset (positive = downstream, negative = upstream)
#' @param use.names ignored
#' @return A RegionSet with shifted regions
#' @rdname shift
#' @export
setMethod("shift", "RegionSet", function(x, shift = 0L, use.names = TRUE) {
  ptr <- .Call(wrap__r_shift, .ptr(x), as.integer(shift))
  new("RegionSet", ptr = ptr, strand = x@strand)
})

#' @rdname shift
#' @export
setMethod("shift", "character", function(x, shift = 0L,
                                          use.names = TRUE) {
  shift(RegionSet(x), shift = shift)
})

#' @rdname shift
#' @export
setMethod("shift", "data.frame", function(x, shift = 0L,
                                           use.names = TRUE) {
  shift(RegionSet(x), shift = shift)
})

#' Generate flanking regions
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param width Flank width in base pairs
#' @param start If TRUE (default), flank upstream of start; if FALSE,
#'   flank downstream of end
#' @param both If TRUE, flank on both sides
#' @param ... ignored
#' @return A RegionSet with flanking regions
#' @rdname flank
#' @export
setMethod("flank", "RegionSet", function(x, width, start = TRUE,
                                          both = FALSE, use.names = TRUE, ...) {
  ptr <- .Call(wrap__r_flank, .ptr(x), as.integer(width),
               as.logical(start), as.logical(both))
  new("RegionSet", ptr = ptr, strand = x@strand)
})

#' @rdname flank
#' @export
setMethod("flank", "character", function(x, width, start = TRUE,
                                          both = FALSE,
                                          use.names = TRUE, ...) {
  flank(RegionSet(x), width = width, start = start,
        both = both)
})

#' @rdname flank
#' @export
setMethod("flank", "data.frame", function(x, width,
                                           start = TRUE,
                                           both = FALSE,
                                           use.names = TRUE, ...) {
  flank(RegionSet(x), width = width, start = start,
        both = both)
})

#' Resize regions to a fixed width
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param width New width in base pairs
#' @param fix Anchor point: "start" (default), "end", or "center"
#' @param ... ignored
#' @return A RegionSet with resized regions
#' @rdname resize
#' @export
setMethod("resize", "RegionSet", function(x, width, fix = "start", ...) {
  ptr <- .Call(wrap__r_resize, .ptr(x), as.integer(width), as.character(fix))
  new("RegionSet", ptr = ptr, strand = x@strand)
})

#' @rdname resize
#' @export
setMethod("resize", "character", function(x, width,
                                           fix = "start", ...) {
  resize(RegionSet(x), width = width, fix = fix)
})

#' @rdname resize
#' @export
setMethod("resize", "data.frame", function(x, width,
                                            fix = "start", ...) {
  resize(RegionSet(x), width = width, fix = fix)
})

#' Narrow regions by specifying a relative sub-range
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param start 1-based relative start within each region (NA to omit)
#' @param end 1-based relative end within each region (NA to omit)
#' @param width Width of the sub-range (NA to omit)
#' @param use.names ignored
#' @return A RegionSet with narrowed regions
#' @rdname narrow
#' @export
setMethod("narrow", "RegionSet", function(x, start = NA, end = NA,
                                           width = NA, use.names = TRUE) {
  s <- if (is.na(start)) NA_integer_ else as.integer(start)
  e <- if (is.na(end)) NA_integer_ else as.integer(end)
  w <- if (is.na(width)) NA_integer_ else as.integer(width)
  ptr <- .Call(wrap__r_narrow, .ptr(x), s, e, w)
  new("RegionSet", ptr = ptr, strand = x@strand)
})

#' @rdname narrow
#' @export
setMethod("narrow", "character", function(x, start = NA,
                                           end = NA, width = NA,
                                           use.names = TRUE) {
  narrow(RegionSet(x), start = start, end = end,
         width = width)
})

#' @rdname narrow
#' @export
setMethod("narrow", "data.frame", function(x, start = NA,
                                            end = NA, width = NA,
                                            use.names = TRUE) {
  narrow(RegionSet(x), start = start, end = end,
         width = width)
})

#' Break regions into non-overlapping disjoint pieces
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param ... ignored
#' @return A RegionSet with disjoint regions
#' @rdname disjoin
#' @export
setMethod("disjoin", "RegionSet", function(x, ...) {
  .rs_from_ptr(.Call(wrap__r_disjoin, .ptr(x)))
})

#' @rdname disjoin
#' @export
setMethod("disjoin", "character", function(x, ...) {
  disjoin(RegionSet(x))
})

#' @rdname disjoin
#' @export
setMethod("disjoin", "data.frame", function(x, ...) {
  disjoin(RegionSet(x))
})

#' Return gaps between regions per chromosome, bounded by chromosome sizes
#'
#' Emits the peak-free intervals of each chromosome listed in ``chrom_sizes``:
#' leading gaps from 0, inter-region gaps, trailing gaps to the chromosome
#' end, and full-chromosome gaps for any chromosome with no regions.
#' Regions on chromosomes not listed in ``chrom_sizes`` are skipped.
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param chrom_sizes Named integer vector of chromosome sizes (name =
#'   chromosome, value = size in bp). Required.
#' @param start,end Ignored. Accepted only for compatibility with the
#'   \code{IRanges::gaps} generic signature; use ``chrom_sizes`` instead.
#' @param ... ignored
#' @return A RegionSet containing the gap intervals
#' @rdname gaps
#' @export
setMethod("gaps", "RegionSet", function(x, start = NA, end = NA, chrom_sizes = NULL, ...) {
  if (is.null(chrom_sizes)) {
    stop("'chrom_sizes' is required: pass a named integer vector of chromosome sizes")
  }
  .rs_from_ptr(.Call(
    wrap__r_gaps,
    .ptr(x),
    names(chrom_sizes),
    as.integer(chrom_sizes)
  ))
})

#' @rdname gaps
#' @export
setMethod("gaps", "character", function(x, start = NA, end = NA, chrom_sizes = NULL, ...) {
  gaps(RegionSet(x), chrom_sizes = chrom_sizes)
})

#' @rdname gaps
#' @export
setMethod("gaps", "data.frame", function(x, start = NA, end = NA, chrom_sizes = NULL, ...) {
  gaps(RegionSet(x), chrom_sizes = chrom_sizes)
})

# =========================================================================
# BiocGenerics methods: binary (merge ops — strand lost)
# =========================================================================

#' @rdname RegionSet-class
#' @export
setMethod("union", c("RegionSet", "RegionSet"), function(x, y) {
  .rs_from_ptr(.Call(wrap__r_union, .ptr(x), .ptr(y)))
})

#' @rdname RegionSet-class
#' @export
setMethod("union", c("RegionSet", "ANY"), function(x, y) {
  union(x, RegionSet(y))
})

#' @rdname RegionSet-class
#' @export
setMethod("union", c("ANY", "RegionSet"), function(x, y) {
  union(RegionSet(x), y)
})

#' @rdname RegionSet-class
#' @export
setMethod("intersect", c("RegionSet", "RegionSet"), function(x, y) {
  .rs_from_ptr(.Call(wrap__r_intersect, .ptr(x), .ptr(y)))
})

#' @rdname RegionSet-class
#' @export
setMethod("intersect", c("RegionSet", "ANY"), function(x, y) {
  intersect(x, RegionSet(y))
})

#' @rdname RegionSet-class
#' @export
setMethod("intersect", c("ANY", "RegionSet"), function(x, y) {
  intersect(RegionSet(x), y)
})

#' @rdname RegionSet-class
#' @export
setMethod("setdiff", c("RegionSet", "RegionSet"), function(x, y) {
  .rs_from_ptr(.Call(wrap__r_setdiff, .ptr(x), .ptr(y)))
})

#' @rdname RegionSet-class
#' @export
setMethod("setdiff", c("RegionSet", "ANY"), function(x, y) {
  setdiff(x, RegionSet(y))
})

#' @rdname RegionSet-class
#' @export
setMethod("setdiff", c("ANY", "RegionSet"), function(x, y) {
  setdiff(RegionSet(x), y)
})

# =========================================================================
# Own generics: binary interval operations
# =========================================================================

#' Pairwise intersection by position
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param y A RegionSet, GRanges, file path, or data.frame
#' @return A RegionSet with pairwise intersections
#' @export
setGeneric("pintersect", function(x, y, ...) standardGeneric("pintersect"))

#' @rdname pintersect
#' @export
setMethod("pintersect", c("ANY", "ANY"), function(x, y, ...) {
  x <- RegionSet(x)
  y <- RegionSet(y)
  # Same number of regions, same order — take strand from x
  ptr <- .Call(wrap__r_pintersect, .ptr(x), .ptr(y))
  new("RegionSet", ptr = ptr, strand = x@strand)
})

#' Combine two region sets without merging
#'
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param y A RegionSet, GRanges, file path, or data.frame
#' @return A RegionSet with all regions from both sets
#' @export
setGeneric("concat", function(x, y) standardGeneric("concat"))

#' @rdname concat
#' @export
setMethod("concat", c("ANY", "ANY"), function(x, y) {
  x <- RegionSet(x)
  y <- RegionSet(y)
  ptr <- .Call(wrap__r_concat, .ptr(x), .ptr(y))
  # Concatenate strand vectors
  new("RegionSet", ptr = ptr, strand = c(x@strand, y@strand))
})

#' Nucleotide-level Jaccard similarity
#'
#' @description Computes |intersection| / |union| in base pairs between two
#'   region sets. Both sets are reduced (merged) before comparison.
#' @param x A RegionSet, GRanges, file path, or data.frame
#' @param y A RegionSet, GRanges, file path, or data.frame
#' @return Numeric value between 0 and 1
#' @export
setGeneric("jaccard", function(x, y) standardGeneric("jaccard"))

#' @rdname jaccard
#' @export
setMethod("jaccard", c("ANY", "ANY"), function(x, y) {
  x <- RegionSet(x)
  y <- RegionSet(y)
  .Call(wrap__r_jaccard, .ptr(x), .ptr(y))
})

# =========================================================================
# Overlap queries (findOverlaps / countOverlaps)
# =========================================================================

#' Find overlapping pairs between two region sets
#'
#' @description Returns all (queryHits, subjectHits) pairs where the query
#'   region overlaps the subject region by at least \code{minoverlap} base
#'   pairs. Uses an IGD (Integrated Genome Database) index for fast queries.
#' @param query A RegionSet, GRanges, file path, or data.frame
#' @param subject A RegionSet, GRanges, file path, or data.frame
#' @param minoverlap Minimum overlap in base pairs (default 1)
#' @param type ignored (present for generic compatibility)
#' @param maxgap ignored
#' @param select ignored
#' @param ... ignored
#' @return A data.frame with columns queryHits and subjectHits (1-based)
#' @rdname findOverlaps
#' @export
setMethod("findOverlaps", c("RegionSet", "RegionSet"),
  function(query, subject, maxgap = -1L, minoverlap = 0L,
           type = c("any", "start", "end", "within", "equal"),
           select = c("all", "first", "last", "arbitrary"), ...) {
    mo <- max(as.integer(minoverlap), 1L)
    result <- .Call(wrap__r_find_overlaps, .ptr(query), .ptr(subject), mo)
    data.frame(queryHits = result$queryHits, subjectHits = result$subjectHits)
})

#' @rdname findOverlaps
#' @export
setMethod("findOverlaps", c("RegionSet", "ANY"),
  function(query, subject, maxgap = -1L, minoverlap = 0L,
           type = c("any", "start", "end", "within", "equal"),
           select = c("all", "first", "last", "arbitrary"), ...) {
    findOverlaps(query, RegionSet(subject),
                 minoverlap = minoverlap)
})

#' @rdname findOverlaps
#' @export
setMethod("findOverlaps", c("ANY", "RegionSet"),
  function(query, subject, maxgap = -1L, minoverlap = 0L,
           type = c("any", "start", "end", "within", "equal"),
           select = c("all", "first", "last", "arbitrary"), ...) {
    findOverlaps(RegionSet(query), subject,
                 minoverlap = minoverlap)
})

#' Count overlaps per query region
#'
#' @description For each query region, counts the number of subject regions
#'   that overlap by at least \code{minoverlap} base pairs.
#' @param query A RegionSet, GRanges, file path, or data.frame
#' @param subject A RegionSet, GRanges, file path, or data.frame
#' @param minoverlap Minimum overlap in base pairs (default 1)
#' @param type ignored
#' @param maxgap ignored
#' @param ... ignored
#' @return Integer vector of length equal to the number of query regions
#' @rdname countOverlaps
#' @export
setMethod("countOverlaps", c("RegionSet", "RegionSet"),
  function(query, subject, maxgap = -1L, minoverlap = 0L,
           type = c("any", "start", "end", "within", "equal"), ...) {
    mo <- max(as.integer(minoverlap), 1L)
    .Call(wrap__r_count_overlaps, .ptr(query), .ptr(subject), mo)
})

#' @rdname countOverlaps
#' @export
setMethod("countOverlaps", c("RegionSet", "ANY"),
  function(query, subject, maxgap = -1L, minoverlap = 0L,
           type = c("any", "start", "end", "within", "equal"), ...) {
    countOverlaps(query, RegionSet(subject),
                  minoverlap = minoverlap)
})

#' @rdname countOverlaps
#' @export
setMethod("countOverlaps", c("ANY", "RegionSet"),
  function(query, subject, maxgap = -1L, minoverlap = 0L,
           type = c("any", "start", "end", "within", "equal"), ...) {
    countOverlaps(RegionSet(query), subject,
                  minoverlap = minoverlap)
})

# =========================================================================
# Consensus (list-based, not a method on a single RegionSet)
# =========================================================================

#' Compute consensus regions from multiple region sets
#'
#' @description Given a list of region sets, computes the union of all regions
#'   and annotates each union region with the number of input sets that overlap
#'   it.
#'
#' @param sets A list of RegionSet objects, GRanges objects, file paths,
#'   or data.frames
#' @return A data.frame with chr, start, end, count columns (0-based half-open)
#'
#' @export
consensus <- function(sets) {
  sets <- lapply(sets, function(s) .ptr(RegionSet(s)))
  result <- .Call(wrap__r_consensus, sets)
  data.frame(chr = result$chr, start = result$start,
             end = result$end, count = result$count,
             stringsAsFactors = FALSE)
}
