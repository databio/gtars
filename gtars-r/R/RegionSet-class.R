#' @importFrom methods setClass setGeneric setMethod new is show
#' @importFrom BiocGenerics union setdiff start end
NULL

# =========================================================================
# RegionSet S4 class
# =========================================================================

#' RegionSet class
#'
#' @description S4 class representing a set of genomic regions, backed by
#'   a Rust RegionSet via externalptr. Carries an R-side strand vector.
#'
#' @slot ptr An externalptr to the Rust RegionSet
#' @slot strand A character vector of strand values ("+", "-", or "*")
#'
#' @exportClass RegionSet
setClass("RegionSet", slots = list(
  ptr = "externalptr",
  strand = "character"
))

#' Create a RegionSet
#'
#' @description Constructs a RegionSet from a file path, GRanges object,
#'   data.frame, or existing pointer.
#'
#' @param x One of:
#'   \itemize{
#'     \item A character string file path to a BED/narrowPeak/gzip file
#'     \item A GRanges object (requires GenomicRanges package; coordinates
#'       converted from 1-based closed to 0-based half-open)
#'     \item A data.frame with chr, start, end columns (0-based half-open);
#'       optionally a strand column
#'     \item An existing RegionSet (returned as-is)
#'     \item An externalptr (wrapped in RegionSet with strand = "*")
#'   }
#' @return A \code{RegionSet} object
#'
#' @export
RegionSet <- function(x) {
  if (is(x, "RegionSet")) return(x)

  if (is(x, "externalptr")) {
    n <- .Call(wrap__r_regionset_length, x)
    return(new("RegionSet", ptr = x, strand = rep("*", n)))
  }

  if (inherits(x, "GRanges")) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
      stop("GenomicRanges package required for GRanges conversion")
    }
    chrs <- as.character(GenomicRanges::seqnames(x))
    starts <- BiocGenerics::start(x) - 1L
    ends <- as.integer(BiocGenerics::end(x))
    strand_vec <- as.character(BiocGenerics::strand(x))
    ptr <- .Call(wrap__r_regionset_from_vectors, chrs, starts, ends)
    return(new("RegionSet", ptr = ptr, strand = strand_vec))
  }

  if (is.character(x) && length(x) == 1) {
    ptr <- .Call(wrap__r_load_regionset, x)
    n <- .Call(wrap__r_regionset_length, ptr)
    return(new("RegionSet", ptr = ptr, strand = rep("*", n)))
  }

  if (is.data.frame(x)) {
    strand_vec <- if ("strand" %in% names(x)) {
      as.character(x$strand)
    } else {
      rep("*", nrow(x))
    }
    ptr <- .Call(wrap__r_regionset_from_vectors,
                 as.character(x$chr), as.integer(x$start), as.integer(x$end))
    return(new("RegionSet", ptr = ptr, strand = strand_vec))
  }

  stop("Cannot create RegionSet: expected file path, GRanges, ",
       "data.frame, or externalptr")
}

# =========================================================================
# Basic S4 methods
# =========================================================================

#' @export
setMethod("show", "RegionSet", function(object) {
  n <- .Call(wrap__r_regionset_length, .ptr(object))
  cat(sprintf("RegionSet with %d regions (0-based half-open)\n", n))
  if (n == 0L) return(invisible(NULL))
  vecs <- .Call(wrap__r_regionset_to_vectors, .ptr(object))
  show_n <- min(n, 5L)
  has_strand <- any(object@strand != "*")
  df <- data.frame(
    chr    = vecs$chr[seq_len(show_n)],
    start  = vecs$start[seq_len(show_n)],
    end    = vecs$end[seq_len(show_n)],
    stringsAsFactors = FALSE
  )
  if (has_strand) {
    df$strand <- object@strand[seq_len(show_n)]
  }
  rownames(df) <- paste0("  [", seq_len(show_n), "]")
  print(df, right = FALSE)
  if (n > 5L) cat(sprintf("  ... and %d more regions\n", n - 5L))
})

#' @export
setMethod("length", "RegionSet", function(x) {
  .Call(wrap__r_regionset_length, .ptr(x))
})

#' @rdname RegionSet-class
#' @export
setMethod("as.data.frame", "RegionSet", function(x, ...) {
  vecs <- .Call(wrap__r_regionset_to_vectors, .ptr(x))
  data.frame(chr = vecs$chr, start = vecs$start, end = vecs$end,
             strand = x@strand, stringsAsFactors = FALSE)
})

#' Subset a RegionSet
#'
#' @param x A RegionSet
#' @param i Numeric, logical, or character index
#' @param j ignored
#' @param ... ignored
#' @param drop ignored
#' @return A RegionSet with selected regions
#' @rdname RegionSet-class
#' @export
setMethod("[", "RegionSet", function(x, i, j, ..., drop = TRUE) {
  df <- as.data.frame(x)
  df <- df[i, , drop = FALSE]
  RegionSet(df)
})

# =========================================================================
# Conversion utilities
# =========================================================================

#' Convert to RegionSet
#'
#' @description Converts a GRanges object, file path, data.frame, or
#'   externalptr to a RegionSet. Alias for \code{\link{RegionSet}}.
#'
#' @param x Input to convert (see \code{\link{RegionSet}} for accepted types)
#' @return A \code{RegionSet} object
#'
#' @export
as_regionset <- function(x) RegionSet(x)

#' Convert RegionSet to GRanges
#'
#' @description Converts a RegionSet to a GRanges object. Coordinates are
#'   converted from 0-based half-open (BED) to 1-based closed (GRanges).
#'   Strand information is preserved.
#'
#' @param rs A RegionSet object
#' @return A GRanges object
#'
#' @export
as_granges <- function(rs) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges package required for GRanges conversion")
  }
  rs <- RegionSet(rs)
  vecs <- .Call(wrap__r_regionset_to_vectors, .ptr(rs))
  GenomicRanges::GRanges(
    seqnames = vecs$chr,
    ranges = IRanges::IRanges(
      start = vecs$start + 1L,
      end = vecs$end
    ),
    strand = rs@strand
  )
}

# =========================================================================
# Internal helpers
# =========================================================================

# Extract the raw externalptr for passing to .Call()
.ptr <- function(x) {
  if (is(x, "RegionSet")) return(x@ptr)
  x
}

# Ensure input is a raw externalptr for .Call()
.ensure_regionset <- function(x) {
  if (is(x, "RegionSet")) return(x@ptr)
  if (is(x, "externalptr")) return(x)
  RegionSet(x)@ptr
}

# Create a RegionSet from a Rust pointer with default strand
.rs_from_ptr <- function(ptr) {
  n <- .Call(wrap__r_regionset_length, ptr)
  new("RegionSet", ptr = ptr, strand = rep("*", n))
}
