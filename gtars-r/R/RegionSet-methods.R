#' @include RegionSet-class.R
#' @importFrom BiocGenerics union setdiff
#' @importFrom IRanges reduce promoters trim
NULL

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
setGeneric("distribution", function(x, nBins = 250L, ...) standardGeneric("distribution"))

#' @rdname distribution
#' @export
setMethod("distribution", "RegionSet", function(x, nBins = 250L, ...) {
  result <- .Call(wrap__r_region_distribution, .ptr(x), as.integer(nBins))
  dt <- data.table::as.data.table(result)
  data.table::setnames(dt, c("n", "rid"), c("N", "regionID"))
  data.table::set(dt, j = "withinGroupID", value = dt[["regionID"]])
  data.table::set(dt, j = "chr",
                  value = factor(dt[["chr"]], levels = unique(dt[["chr"]])))
  dt
})

#' @rdname distribution
#' @export
setMethod("distribution", "ANY", function(x, nBins = 250L, ...) {
  distribution(RegionSet(x), nBins = nBins)
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
setMethod("trim", "ANY", function(x, chromSizes, ...) {
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
setMethod("reduce", "ANY", function(x, ...) {
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
setMethod("promoters", "ANY", function(x, upstream = 2000L,
                                        downstream = 200L, ...) {
  promoters(RegionSet(x), upstream = upstream, downstream = downstream)
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
