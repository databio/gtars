# =========================================================================
# RegionSetList S4 class
# =========================================================================

#' RegionSetList class
#'
#' @description S4 class representing a collection of RegionSets, backed by
#'   a Rust RegionSetList via externalptr. The gtars equivalent of GRangesList.
#'
#' @slot ptr An externalptr to the Rust RegionSetList
#'
#' @exportClass RegionSetList
setClass("RegionSetList", slots = list(
  ptr = "externalptr"
))

#' Create a RegionSetList
#'
#' @description Constructs a RegionSetList from RegionSet objects, file paths,
#'   GRanges objects, or an existing pointer.
#'
#' @param ... RegionSet objects, file paths, GRanges objects, or a single list
#'   of such objects. If no arguments, creates an empty RegionSetList.
#' @return A \code{RegionSetList} object
#'
#' @export
RegionSetList <- function(...) {
  args <- list(...)

  # Handle single list argument
  if (length(args) == 1 && is.list(args[[1]]) && !is(args[[1]], "RegionSet")) {
    args <- args[[1]]
  }

  # Empty case
  if (length(args) == 0) {
    sets <- list()
  } else {
    sets <- lapply(args, function(x) {
      if (is(x, "RegionSet")) return(x)
      RegionSet(x)
    })
  }

  # Extract pointers and pass to Rust
  ptrs <- lapply(sets, function(s) s@ptr)
  ptr <- regionsetlist_from_sets(ptrs)
  new("RegionSetList", ptr = ptr)
}

# =========================================================================
# Basic S4 methods
# =========================================================================

#' @export
setMethod("show", "RegionSetList", function(object) {
  n <- regionsetlist_length(object@ptr)
  cat(sprintf("RegionSetList with %d region sets\n", n))
  if (n == 0L) return(invisible(NULL))
  nms <- regionsetlist_names(object@ptr)
  show_n <- min(n, 5L)
  for (i in seq_len(show_n)) {
    rs_ptr <- regionsetlist_get(object@ptr, as.integer(i))
    rs_len <- .Call(wrap__r_regionset_length, rs_ptr)
    label <- if (!is.null(nms)) nms[i] else paste0("[[", i, "]]")
    cat(sprintf("  %s: %d regions\n", label, rs_len))
  }
  if (n > 5L) cat(sprintf("  ... and %d more\n", n - 5L))
})

#' @export
setMethod("length", "RegionSetList", function(x) {
  regionsetlist_length(x@ptr)
})

#' @export
setMethod("[[", "RegionSetList", function(x, i, j, ...) {
  ptr <- regionsetlist_get(x@ptr, as.integer(i))
  .rs_from_ptr(ptr)
})

#' @export
setMethod("names", "RegionSetList", function(x) {
  regionsetlist_names(x@ptr)
})

# =========================================================================
# concat method for RegionSetList
# =========================================================================

#' Flatten a RegionSetList into a single RegionSet
#'
#' @description Flatten all region sets into a single RegionSet without
#'   merging or deduplication. Apply \code{disjoin()} or \code{reduce()}
#'   on the result as needed.
#'
#' When called as \code{concat(rsl)} (single argument), flattens the list.
#' The generic \code{concat(x, y)} is defined in RegionSet-methods.R.
#'
#' @param x A RegionSetList object
#' @param y missing (not used for RegionSetList)
#' @return A \code{RegionSet} object
#'
#' @export
setMethod("concat", c("RegionSetList", "missing"), function(x, y) {
  ptr <- regionsetlist_concat(x@ptr)
  .rs_from_ptr(ptr)
})
