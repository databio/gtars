# LOLA (Locus Overlap Analysis) R interface
#
# Drop-in replacement for the LOLA R/Bioconductor package,
# powered by gtars-lola (Rust) for performance.

#' Load a LOLA region database from a folder
#'
#' @param dbLocation Path to the LOLA database folder
#' @param useCache Logical, currently ignored (reserved for future caching)
#' @param limit Optional integer limit on files per collection
#' @param collections Optional character vector of collection names to load
#' @return An external pointer to a RegionDB object
#' @export
loadRegionDB <- function(dbLocation, useCache = TRUE, limit = NULL,
                          collections = NULL) {
  if (!is.character(dbLocation) || length(dbLocation) != 1) {
    stop("dbLocation must be a single character string")
  }
  if (!dir.exists(dbLocation)) {
    stop("Database directory does not exist: ", dbLocation)
  }
  .Call(wrap__load_region_db, dbLocation, collections, limit)
}

#' Load a region database from BED file paths
#'
#' @param bedFiles Character vector of paths to BED files
#' @param filenames Optional character vector of display names
#' @return An external pointer to a RegionDB object
#' @export
loadRegionDBFromBeds <- function(bedFiles, filenames = NULL) {
  if (!is.character(bedFiles)) {
    stop("bedFiles must be a character vector")
  }
  .Call(wrap__load_region_db_from_beds, bedFiles, filenames)
}

#' Run LOLA enrichment analysis
#'
#' @param userSets A RegionSet object or list of RegionSet objects
#' @param userUniverse A RegionSet object representing the universe
#' @param regionDB An external pointer to a RegionDB (from loadRegionDB)
#' @param minOverlap Minimum base-pair overlap (default 1)
#' @param cores Number of cores (currently unused, reserved for future parallelism)
#' @param redefineUserSets Logical, whether to redefine user sets against universe
#' @param direction "enrichment" or "depletion"
#' @return A data.frame with LOLA results
#' @export
runLOLA <- function(userSets, userUniverse, regionDB,
                     minOverlap = 1, cores = 1,
                     redefineUserSets = FALSE,
                     direction = "enrichment") {
  # Normalize userSets to a list of RegionSet pointers
  if (is(userSets, "RegionSet")) {
    user_ptrs <- list(.ptr(userSets))
  } else if (is.list(userSets)) {
    user_ptrs <- lapply(userSets, function(x) {
      if (is(x, "RegionSet")) .ptr(x)
      else stop("All elements of userSets must be RegionSet objects")
    })
  } else {
    stop("userSets must be a RegionSet or list of RegionSets")
  }

  # Get universe pointer
  if (is(userUniverse, "RegionSet")) {
    universe_ptr <- .ptr(userUniverse)
  } else {
    stop("userUniverse must be a RegionSet object")
  }

  # Optionally redefine user sets
  if (redefineUserSets) {
    user_ptrs <- .Call(wrap__redefine_user_sets, user_ptrs, universe_ptr)
  }

  # Validate direction
  direction <- match.arg(direction, c("enrichment", "depletion"))

  # Run LOLA
  result <- .Call(wrap__run_lola, user_ptrs, universe_ptr, regionDB,
                  as.integer(minOverlap), direction)

  # Convert to data.frame
  as.data.frame(result, stringsAsFactors = FALSE)
}

#' List region set filenames in a database
#'
#' @param regionDB An external pointer to a RegionDB
#' @param collections Optional character vector filter
#' @return Character vector of filenames
#' @export
listRegionSets <- function(regionDB, collections = NULL) {
  .Call(wrap__list_region_sets, regionDB, collections)
}

#' Check universe appropriateness
#'
#' @param userSets A RegionSet or list of RegionSets
#' @param userUniverse A RegionSet
#' @param cores Number of cores (currently unused)
#' @param fast Logical (currently unused)
#' @return A data.frame with coverage diagnostics per user set
#' @export
checkUniverseAppropriateness <- function(userSets, userUniverse,
                                          cores = 1, fast = FALSE) {
  if (is(userSets, "RegionSet")) {
    user_ptrs <- list(.ptr(userSets))
  } else if (is.list(userSets)) {
    user_ptrs <- lapply(userSets, .ptr)
  } else {
    stop("userSets must be a RegionSet or list of RegionSets")
  }

  if (is(userUniverse, "RegionSet")) {
    universe_ptr <- .ptr(userUniverse)
  } else {
    stop("userUniverse must be a RegionSet object")
  }

  result <- .Call(wrap__check_universe, user_ptrs, universe_ptr)
  df <- as.data.frame(result[names(result) != "warnings"],
                      stringsAsFactors = FALSE)

  # Print warnings
  for (w in result$warnings) {
    warning(w, call. = FALSE)
  }

  df
}

#' Redefine user sets in terms of universe regions
#'
#' @param userSets A RegionSet or list of RegionSets
#' @param userUniverse A RegionSet
#' @param cores Number of cores (currently unused)
#' @return A list of RegionSet objects
#' @export
redefineUserSets <- function(userSets, userUniverse, cores = 1) {
  if (is(userSets, "RegionSet")) {
    user_ptrs <- list(.ptr(userSets))
  } else if (is.list(userSets)) {
    user_ptrs <- lapply(userSets, .ptr)
  } else {
    stop("userSets must be a RegionSet or list of RegionSets")
  }

  if (is(userUniverse, "RegionSet")) {
    universe_ptr <- .ptr(userUniverse)
  } else {
    stop("userUniverse must be a RegionSet object")
  }

  ptrs <- .Call(wrap__redefine_user_sets, user_ptrs, universe_ptr)
  lapply(ptrs, function(p) RegionSet(p))
}

#' Build a restricted universe from user sets
#'
#' @param userSets A list of RegionSets
#' @return A RegionSet representing the restricted universe
#' @export
buildRestrictedUniverse <- function(userSets) {
  if (is(userSets, "RegionSet")) {
    user_ptrs <- list(.ptr(userSets))
  } else if (is.list(userSets)) {
    user_ptrs <- lapply(userSets, .ptr)
  } else {
    stop("userSets must be a RegionSet or list of RegionSets")
  }

  ptr <- .Call(wrap__build_restricted_universe, user_ptrs)
  RegionSet(ptr)
}

# Helper to extract the external pointer from a RegionSet S4 object
.ptr <- function(x) {
  if (is(x, "RegionSet")) x@ptr
  else if (is(x, "externalptr")) x
  else stop("Expected a RegionSet object")
}
