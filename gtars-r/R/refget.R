#' @useDynLib gtars, .registration = TRUE
#' @importFrom utils capture.output
NULL

# =========================================================================
# S4 Class Definitions
# =========================================================================

#' SequenceMetadata class
#'
#' @description Metadata for one sequence, without its bases. Returned by
#'   \code{\link{get_sequence_metadata}} and \code{\link{list_sequences}}.
#'
#' @slot name Sequence name
#' @slot length Sequence length
#' @slot sha512t24u GA4GH sha512t24u digest
#' @slot md5 MD5 digest
#' @slot alphabet Sequence alphabet
#' @slot description Description from the FASTA header, or ""
#' @param object A SequenceMetadata object
setClass(
  'SequenceMetadata',
  slots = list(
    name = 'character',
    length = 'integer',
    sha512t24u = 'character',
    md5 = 'character',
    alphabet = 'character',
    description = 'character'
  )
)

#' SequenceRecord class
#'
#' @description A sequence's metadata plus its bases. Returned by
#'   \code{\link{get_sequence}}.
#'
#' @slot metadata A \code{SequenceMetadata} object
#' @slot data The sequence bases, or "" when they are not loaded
#' @param object A SequenceRecord object
setClass(
  'SequenceRecord',
  slots = list(
    metadata = 'SequenceMetadata',
    data = 'character'
  )
)

#' SeqColDigestLvl1 class
#'
#' @description The level 1 digests of a sequence collection: one digest per
#'   attribute array.
#'
#' @slot sequences_digest Digest of the sequences array
#' @slot names_digest Digest of the names array
#' @slot lengths_digest Digest of the lengths array
#' @param object A SeqColDigestLvl1 object
setClass(
  'SeqColDigestLvl1',
  slots = list(
    sequences_digest = 'character',
    names_digest = 'character',
    lengths_digest = 'character'
  )
)

#' SequenceCollectionMetadata class
#'
#' @description Metadata for one sequence collection, without its sequences.
#'   Returned by \code{\link{get_collection_metadata}} and
#'   \code{\link{list_collections}}.
#'
#' @slot digest Collection digest
#' @slot n_sequences Number of sequences in the collection
#' @slot names_digest Digest of the names array
#' @slot sequences_digest Digest of the sequences array
#' @slot lengths_digest Digest of the lengths array
#' @slot name_length_pairs_digest Digest of the name_length_pairs array
#' @slot sorted_name_length_pairs_digest Digest of the
#'   sorted_name_length_pairs array
#' @slot sorted_sequences_digest Digest of the sorted_sequences array
#' @param object A SequenceCollectionMetadata object
setClass(
  'SequenceCollectionMetadata',
  slots = list(
    digest = 'character',
    n_sequences = 'integer',
    names_digest = 'character',
    sequences_digest = 'character',
    lengths_digest = 'character',
    name_length_pairs_digest = 'character',
    sorted_name_length_pairs_digest = 'character',
    sorted_sequences_digest = 'character'
  )
)

#' SequenceCollection class
#'
#' @description A sequence collection: its digests and a list of
#'   \code{SequenceRecord} objects. Returned by \code{\link{digest_fasta}} and
#'   \code{\link{get_collection}}. Use \code{length()} for the number of
#'   sequences and \code{x[i]} for the i-th \code{SequenceRecord}.
#'
#' @slot sequences A list of \code{SequenceRecord} objects
#' @slot digest Collection digest
#' @slot lvl1 A \code{SeqColDigestLvl1} object
#' @slot file_path Source FASTA path, or ""
#' @param object,x A SequenceCollection object
#' @param i Index of the sequence to return
setClass(
  'SequenceCollection',
  slots = list(
    sequences = 'list',
    digest = 'character',
    lvl1 = 'SeqColDigestLvl1',
    file_path = 'character'
  )
)

#' RefgetStore class
#'
#' @description A store of sequences and sequence collections, backed by a
#'   Rust RefgetStore. Create one with \code{\link{refget_store}},
#'   \code{\link{refget_store_on_disk}}, \code{\link{refget_store_open_local}}
#'   or \code{\link{refget_store_open_remote}}.
#'
#' @slot ptr An externalptr to the Rust RefgetStore
#' @param object,x A RefgetStore object
setClass(
  'RefgetStore',
  slots = list(
    ptr = 'externalptr'
  )
)

#' RetrievedSequence class
#'
#' @description The bases of one region, as returned by
#'   \code{\link{get_seqs_bed_file_to_vec}}.
#'
#' @slot sequence The region's bases
#' @slot chrom_name Sequence (chromosome) name
#' @slot start Start position (0-based, inclusive)
#' @slot end End position (0-based, exclusive)
#' @param object A RetrievedSequence object
setClass(
  'RetrievedSequence',
  slots = list(
    sequence = 'character',
    chrom_name = 'character',
    start = 'integer',
    end = 'integer'
  )
)

# =========================================================================
# Show Methods
# =========================================================================

#' @rdname SequenceMetadata-class
setMethod('show', 'SequenceMetadata', function(object) {
  cat(sprintf('SequenceMetadata for sequence %s\n', object@name))
  cat(sprintf('  name: %s\n', object@name))
  cat(sprintf('  length: %d\n', object@length))
  cat(sprintf('  sha512t24u: %s\n', object@sha512t24u))
  cat(sprintf('  md5: %s\n', object@md5))
  cat(sprintf('  alphabet: %s\n', object@alphabet))
  if (nzchar(object@description)) {
    cat(sprintf('  description: %s\n', object@description))
  }
})

#' @rdname SequenceRecord-class
setMethod('show', 'SequenceRecord', function(object) {
  cat(sprintf('SequenceRecord for %s\n', object@metadata@name))
})

#' @rdname SeqColDigestLvl1-class
setMethod('show', 'SeqColDigestLvl1', function(object) {
  cat('SeqColDigestLvl1: \n')
  cat(sprintf('  sequences digest: %s\n', object@sequences_digest))
  cat(sprintf('  names digest: %s\n', object@names_digest))
  cat(sprintf('  lengths digest: %s\n', object@lengths_digest))
})

#' @rdname SequenceCollectionMetadata-class
setMethod('show', 'SequenceCollectionMetadata', function(object) {
  cat(sprintf('SequenceCollectionMetadata: %s (%d sequences)\n',
              object@digest, object@n_sequences))
})

#' @rdname SequenceCollection-class
setMethod('show', 'SequenceCollection', function(object) {
  cat(sprintf('SequenceCollection with %d sequences, digest = %s\n',
              length(object@sequences), object@digest))
})

#' @rdname SequenceCollection-class
setMethod('length', 'SequenceCollection', function(x) {
  length(x@sequences)
})

#' @rdname RetrievedSequence-class
setMethod('show', 'RetrievedSequence', function(object) {
  cat(sprintf('RetrievedSequence: %s|%d-%d: %s\n',
              object@chrom_name, object@start, object@end, object@sequence))
})

setMethod('print', 'RetrievedSequence', function(x) {
  cat(sprintf('%s|%d-%d: %s\n', x@chrom_name, x@start, x@end, x@sequence))
})

#' @rdname SequenceCollection-class
setMethod('[', c('SequenceCollection', 'numeric', 'missing'),
          function(x, i) {
            if (i < 1 || i > length(x@sequences)) {
              stop('Index out of range')
            }
            return(x@sequences[[i]])
          })

#' @rdname RefgetStore-class
setMethod('show', 'RefgetStore', function(object) {
  stats <- tryCatch(
    .Call(wrap__stats_store, object@ptr),
    error = function(e) NULL
  )
  if (!is.null(stats)) {
    cat('RefgetStore:\n')
    cat(sprintf('  sequences: %d\n', stats$n_sequences))
    cat(sprintf('  collections: %d\n', stats$n_collections))
    cat(sprintf('  storage_mode: %s\n', stats$storage_mode))
  } else {
    ptr_address <- sub('<pointer: (.*)>', '\\1', capture.output(object@ptr))
    cat('RefgetStore: \n')
    cat(sprintf('  ptr: %s\n', ptr_address))
  }
})

#' @rdname RefgetStore-class
setMethod('as.character', 'RefgetStore', function(x) {
  ptr_address <- sub('<pointer: (.*)>', '\\1', capture.output(x@ptr))
  sprintf('RefgetStore (ptr = %s)', ptr_address)
})

# =========================================================================
# Helper Functions
# =========================================================================

convert_to_sequence_metadata <- function(raw_result) {
  if (is.null(raw_result)) return(NULL)

  new(
    'SequenceMetadata',
    name = raw_result$name %||% '',
    length = as.integer(raw_result$length %||% 0L),
    sha512t24u = raw_result$sha512t24u %||% '',
    md5 = raw_result$md5 %||% '',
    alphabet = raw_result$alphabet %||% '',
    description = raw_result$description %||% ''
  )
}

convert_to_sequence_record <- function(raw_result) {
  if (is.null(raw_result)) return(NULL)

  metadata <- convert_to_sequence_metadata(raw_result$metadata)

  new(
    'SequenceRecord',
    metadata = metadata,
    data = raw_result$data %||% ''
  )
}

convert_to_collection_metadata <- function(raw_result) {
  if (is.null(raw_result)) return(NULL)

  new(
    'SequenceCollectionMetadata',
    digest = raw_result$digest %||% '',
    n_sequences = as.integer(raw_result$n_sequences %||% 0L),
    names_digest = raw_result$names_digest %||% '',
    sequences_digest = raw_result$sequences_digest %||% '',
    lengths_digest = raw_result$lengths_digest %||% '',
    name_length_pairs_digest = raw_result$name_length_pairs_digest %||% '',
    sorted_name_length_pairs_digest = raw_result$sorted_name_length_pairs_digest %||% '',
    sorted_sequences_digest = raw_result$sorted_sequences_digest %||% ''
  )
}

convert_to_sequence_collection <- function(raw_result) {
  if (is.null(raw_result)) return(NULL)

  lvl1 <- new(
    'SeqColDigestLvl1',
    sequences_digest = raw_result$lvl1$sequences_digest %||% '',
    names_digest = raw_result$lvl1$names_digest %||% '',
    lengths_digest = raw_result$lvl1$lengths_digest %||% ''
  )

  sequence_records <- lapply(raw_result$sequences, function(seq) {
    convert_to_sequence_record(seq)
  })

  new(
    'SequenceCollection',
    sequences = sequence_records,
    digest = raw_result$digest %||% '',
    lvl1 = lvl1,
    file_path = raw_result$file_path %||% ''
  )
}

convert_to_retrieved_sequence <- function(raw_result) {
  if (is.null(raw_result)) return(NULL)

  new(
    'RetrievedSequence',
    sequence = raw_result$sequence,
    chrom_name = raw_result$chrom_name,
    start = as.integer(raw_result$start),
    end = as.integer(raw_result$end)
  )
}

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# =========================================================================
# Digest Functions
# =========================================================================

#' @title digest fasta
#'
#' @description digest a fasta file given a filepath
#'
#' @param fasta a filepath string to a fasta file
#'
#' @return a SequenceCollection object
#'
#' @examples
#' \dontrun{
#' fasta_path <- 'tests/data/fasta/base.fa'
#' result <- gtars::digest_fasta(fasta_path)
#' }
#'
#' @export
digest_fasta <- function(fasta) {
  result <- .Call(wrap__digest_fasta_raw, fasta)
  convert_to_sequence_collection(result)
}

# =========================================================================
# Store Constructors
# =========================================================================

#' @title Create an in-memory RefgetStore
#'
#' @description Create a new in-memory RefgetStore for storing sequences.
#'
#' @param mode Storage mode: "raw" or "encoded" (default: "encoded")
#'
#' @return A RefgetStore object
#'
#' @examples
#' \dontrun{
#' store <- refget_store()
#' store <- refget_store("raw")
#' }
#'
#' @export
refget_store <- function(mode = "encoded") {
  result <- .Call(wrap__refget_store_raw, mode)
  invisible(new('RefgetStore', ptr = result))
}

#' @title Create a disk-backed RefgetStore
#'
#' @description Create a RefgetStore that persists to disk, or open the one
#' already at \code{path}.
#'
#' @param path Directory path for storing sequences and metadata
#' @param packed Create the store with the packed layout (alpha): sequences are
#'   written into a few large pack files under \code{packs/} instead of one
#'   \code{.seq} file each. A store's layout is fixed when it is created;
#'   opening an existing standalone store with \code{packed = TRUE} is an
#'   error. The packed format may change before it is declared stable.
#'
#' @return A RefgetStore object
#'
#' @examples
#' \dontrun{
#' store <- refget_store_on_disk("/path/to/store")
#' packed <- refget_store_on_disk("/path/to/packed_store", packed = TRUE)
#' }
#'
#' @export
refget_store_on_disk <- function(path, packed = FALSE) {
  result <- .Call(wrap__on_disk_store, path, packed)
  invisible(new('RefgetStore', ptr = result))
}

#' @title Open a local RefgetStore
#'
#' @description Load an existing RefgetStore from a directory.
#'
#' @param path Path to the RefgetStore directory
#'
#' @return A RefgetStore object
#'
#' @examples
#' \dontrun{
#' store <- refget_store_open_local("/path/to/store")
#' }
#'
#' @export
refget_store_open_local <- function(path) {
  result <- .Call(wrap__open_local_store, path)
  invisible(new('RefgetStore', ptr = result))
}

#' @title Open a remote RefgetStore
#'
#' @description Open a remote RefgetStore. Metadata downloads in parallel at
#' open; sequences download on demand over parallel requests, each checked
#' against its digest.
#' \code{cache_path} is the shared cache root (see
#' \code{refget_default_cache_dir()}): metadata goes to \code{remotes/<key>/}
#' and downloaded sequences to \code{sequences/}, shared by every remote store
#' opened with the same root.
#'
#' @param cache_path Cache root shared by all remote stores
#' @param remote_url URL of the remote RefgetStore
#' @param fetch_jobs Parallel HTTP requests
#'
#' @return A RefgetStore object
#'
#' @examples
#' \dontrun{
#' store <- refget_store_open_remote(refget_default_cache_dir(), "https://example.com/store")
#' }
#'
#' @export
refget_store_open_remote <- function(cache_path, remote_url, fetch_jobs = 32L) {
  result <- .Call(wrap__open_remote_store, cache_path, remote_url, as.integer(fetch_jobs))
  invisible(new('RefgetStore', ptr = result))
}

#' @title Default cache root for remote stores
#' @description \code{REFGET_CACHE_DIR}, else \code{$XDG_CACHE_HOME/refget}, else \code{~/.cache/refget}.
#' @return A path
#' @export
refget_default_cache_dir <- function() {
  .Call(wrap__default_cache_dir_raw)
}

#' @title Summarize a cache
#' @param root Cache root
#' @return A named list: root, sequence_files, sequence_bytes, and remotes
#'   (a data frame of url and bytes)
#' @export
refget_cache_info <- function(root) {
  .Call(wrap__cache_info_raw, root)
}

#' @title Delete a cache
#' @description Delete everything in a cache: sequences and remote metadata.
#' @param root Cache root
#' @return Bytes freed
#' @export
refget_clear_cache <- function(root) {
  .Call(wrap__clear_cache_raw, root)
}

# Keep old function name for backwards compatibility
#' @rdname refget_store_open_local
#' @export
load_from_directory <- refget_store_open_local

# =========================================================================
# Store Configuration Methods
# =========================================================================

#' RefgetStore settings
#'
#' @description Turn store features on or off, and read store settings.
#'
#' @details
#' * `enable_encoding()` / `disable_encoding()`: store sequences 2-bit
#'   encoded or as raw bytes. `storage_mode()` returns the current mode
#'   (`"raw"`, `"encoded"` or `"zstd"`).
#' * `set_quiet()` / `is_quiet()`: silence progress messages.
#' * `enable_persistence()` / `disable_persistence()` / `is_persisting()`:
#'   write sequences and metadata to a directory as they are added.
#' * `cache_path()`: the store's local directory, or NULL.
#' * `remote_url()`: the URL of a remote store, or NULL.
#' * `enable_ancillary_digests()` / `disable_ancillary_digests()` /
#'   `has_ancillary_digests()`: compute the extra seqcol digests
#'   (name_length_pairs, sorted_name_length_pairs, sorted_sequences).
#' * `enable_attribute_index()` / `disable_attribute_index()` /
#'   `has_attribute_index()`: keep the index used by
#'   [find_collections_by_attribute()].
#'
#' @param store A RefgetStore object
#' @param quiet TRUE to suppress progress messages
#' @param path Directory for storing sequences and metadata
#' @return The `enable_*`, `disable_*` and `set_quiet` functions return the
#'   store, invisibly. The others return a logical, a string, or NULL.
#' @name store_settings
NULL

#' @rdname store_settings
#' @export
setGeneric('enable_encoding', function(store) standardGeneric('enable_encoding'))
#' @rdname store_settings
setMethod('enable_encoding', 'RefgetStore', function(store) {
  .Call(wrap__enable_encoding_store, store@ptr)
  invisible(store)
})

#' @rdname store_settings
#' @export
setGeneric('disable_encoding', function(store) standardGeneric('disable_encoding'))
#' @rdname store_settings
setMethod('disable_encoding', 'RefgetStore', function(store) {
  .Call(wrap__disable_encoding_store, store@ptr)
  invisible(store)
})

#' @rdname store_settings
#' @export
setGeneric('storage_mode', function(store) standardGeneric('storage_mode'))
#' @rdname store_settings
setMethod('storage_mode', 'RefgetStore', function(store) {
  .Call(wrap__get_storage_mode_store, store@ptr)
})

#' @rdname store_settings
#' @export
setGeneric('set_quiet', function(store, quiet) standardGeneric('set_quiet'))
#' @rdname store_settings
setMethod('set_quiet', 'RefgetStore', function(store, quiet) {
  .Call(wrap__set_quiet_store, store@ptr, quiet)
  invisible(store)
})

#' @rdname store_settings
#' @export
setGeneric('is_quiet', function(store) standardGeneric('is_quiet'))
#' @rdname store_settings
setMethod('is_quiet', 'RefgetStore', function(store) {
  .Call(wrap__get_quiet_store, store@ptr)
})

#' @rdname store_settings
#' @export
setGeneric('is_persisting', function(store) standardGeneric('is_persisting'))
#' @rdname store_settings
setMethod('is_persisting', 'RefgetStore', function(store) {
  .Call(wrap__is_persisting_store, store@ptr)
})

#' @rdname store_settings
#' @export
setGeneric('enable_persistence', function(store, path) standardGeneric('enable_persistence'))
#' @rdname store_settings
setMethod('enable_persistence', 'RefgetStore', function(store, path) {
  .Call(wrap__enable_persistence_store, store@ptr, path)
  invisible(store)
})

#' @rdname store_settings
#' @export
setGeneric('disable_persistence', function(store) standardGeneric('disable_persistence'))
#' @rdname store_settings
setMethod('disable_persistence', 'RefgetStore', function(store) {
  .Call(wrap__disable_persistence_store, store@ptr)
  invisible(store)
})

#' @rdname store_settings
#' @export
setGeneric('cache_path', function(store) standardGeneric('cache_path'))
#' @rdname store_settings
setMethod('cache_path', 'RefgetStore', function(store) {
  .Call(wrap__get_cache_path_store, store@ptr)
})

#' @rdname store_settings
#' @export
setGeneric('remote_url', function(store) standardGeneric('remote_url'))
#' @rdname store_settings
setMethod('remote_url', 'RefgetStore', function(store) {
  .Call(wrap__get_remote_url_store, store@ptr)
})

# =========================================================================
# Adding Sequences/Collections
# =========================================================================

#' Import a FASTA file into a store
#'
#' @description Add the sequences of a FASTA file to the store as one
#'   collection. See \code{\link{add_fasta}} for more options.
#'
#' @param store A RefgetStore object
#' @param file_path Path to the FASTA file
#' @return The store, invisibly
#' @export
setGeneric('import_fasta', function(store, file_path) standardGeneric('import_fasta'))
#' @rdname import_fasta
setMethod('import_fasta', 'RefgetStore', function(store, file_path) {
  .Call(wrap__import_fasta_store, store@ptr, file_path)
  invisible(store)
})

#' Add a FASTA file to a store.
#'
#' @param store A RefgetStore object
#' @param file_path Path to the FASTA file
#' @param force Replace collections and aliases that already exist
#' @param standalone_file_limit Most \code{.seq} files a standard on-disk
#'   store may hold after the import. \code{NULL} (the default) means 500,000.
#'   The import stops before writing anything if it would go over. 0 is not
#'   special; pass a large number to lift the limit. \code{force} does not
#'   skip it. Packed stores ignore it.
#' @return A named list describing the collection, with \code{was_new}
#' @export
setGeneric('add_fasta', function(store, file_path, force = FALSE, standalone_file_limit = NULL) standardGeneric('add_fasta'))
#' @rdname add_fasta
setMethod('add_fasta', 'RefgetStore', function(store, file_path, force = FALSE, standalone_file_limit = NULL) {
  if (!is.null(standalone_file_limit)) {
    standalone_file_limit <- as.numeric(standalone_file_limit)
  }
  .Call(wrap__add_fasta_store, store@ptr, file_path, force, standalone_file_limit)
})

# =========================================================================
# Sequence Retrieval
# =========================================================================

#' Get sequences from a store
#'
#' @description
#' * `get_sequence()`: a sequence by digest (sha512t24u or MD5).
#' * `get_sequence_by_name()`: a sequence by its name in a collection.
#' * `get_sequence_metadata()`: a sequence's metadata, without its bases.
#' * `get_substring()`: part of a sequence, as a string.
#' * `get_sequence_by_id()` and `get_sequence_by_collection_and_name()` are
#'   old names for `get_sequence()` and `get_sequence_by_name()`.
#'
#' @param store A RefgetStore object
#' @param digest Sequence digest (sha512t24u or MD5)
#' @param collection_digest Sequence collection digest
#' @param sequence_name Sequence name in the collection
#' @param seq_digest Sequence digest
#' @param start Start position (0-based, inclusive)
#' @param end End position (0-based, exclusive)
#' @return A \code{SequenceRecord} object, a \code{SequenceMetadata} object
#'   (\code{get_sequence_metadata}), or a string or NULL
#'   (\code{get_substring}).
#' @export
setGeneric('get_sequence', function(store, digest) standardGeneric('get_sequence'))
#' @rdname get_sequence
setMethod('get_sequence', 'RefgetStore', function(store, digest) {
  result <- .Call(wrap__get_sequence_store, store@ptr, digest)
  convert_to_sequence_record(result)
})

#' @rdname get_sequence
#' @export
setGeneric('get_sequence_by_name', function(store, collection_digest, sequence_name) standardGeneric('get_sequence_by_name'))
#' @rdname get_sequence
setMethod('get_sequence_by_name', 'RefgetStore', function(store, collection_digest, sequence_name) {
  result <- .Call(wrap__get_sequence_by_name_store, store@ptr, collection_digest, sequence_name)
  convert_to_sequence_record(result)
})

#' @rdname get_sequence
#' @export
setGeneric('get_sequence_metadata', function(store, digest) standardGeneric('get_sequence_metadata'))
#' @rdname get_sequence
setMethod('get_sequence_metadata', 'RefgetStore', function(store, digest) {
  result <- .Call(wrap__get_sequence_metadata_store, store@ptr, digest)
  convert_to_sequence_metadata(result)
})

#' @rdname get_sequence
#' @export
setGeneric('get_substring', function(store, seq_digest, start, end) standardGeneric('get_substring'))
#' @rdname get_sequence
setMethod('get_substring', 'RefgetStore', function(store, seq_digest, start, end) {
  result <- .Call(wrap__get_substring_store, store@ptr, seq_digest, as.integer(start), as.integer(end))
  if (is.null(result)) NULL else as.character(result)
})

#' Stream a (sub)sequence's decoded bases (flow 2).
#'
#' Reads only the bytes covering \code{[start, end)} from resident memory, a
#' local \code{.seq} file, or a remote HTTP byte-range, and returns the decoded
#' region as a character string. Pass \code{NULL} for \code{start}/\code{end} to
#' stream the whole sequence. Works on remote-backed stores with no preload.
#' @param store A RefgetStore.
#' @param seq_digest Sequence digest (optionally "SQ."-prefixed).
#' @param start Start position (0-based, inclusive), or NULL for the start.
#' @param end End position (0-based, exclusive), or NULL for the sequence length.
#' @export
setGeneric('stream_sequence', function(store, seq_digest, start = NULL, end = NULL) standardGeneric('stream_sequence'))
#' @rdname stream_sequence
setMethod('stream_sequence', 'RefgetStore', function(store, seq_digest, start = NULL, end = NULL) {
  if (!is.null(start)) start <- as.integer(start)
  if (!is.null(end)) end <- as.integer(end)
  as.character(.Call(wrap__stream_sequence_store, store@ptr, seq_digest, start, end))
})

#' Download and cache a single sequence's whole byte data (flow 3).
#'
#' Opt-in for repeat-heavy workloads: after this the sequence is resident (and
#' persisted to the local cache dir when the store persists), so subsequent
#' \code{get_substring}/\code{stream_sequence} calls are served from memory.
#' @param store A RefgetStore.
#' @param digest Sequence digest (optionally "SQ."-prefixed).
#' @export
setGeneric('load_sequence', function(store, digest) standardGeneric('load_sequence'))
#' @rdname load_sequence
setMethod('load_sequence', 'RefgetStore', function(store, digest) {
  invisible(.Call(wrap__load_sequence_store, store@ptr, digest))
})

#' Download and cache every sequence in the store (flow 3, in bulk).
#' @param store A RefgetStore.
#' @export
setGeneric('load_all_sequences', function(store) standardGeneric('load_all_sequences'))
#' @rdname load_all_sequences
setMethod('load_all_sequences', 'RefgetStore', function(store) {
  invisible(.Call(wrap__load_all_sequences_store, store@ptr))
})

#' Load several collections; a remote store downloads the missing files in parallel.
#' @param store A RefgetStore.
#' @param digests Collection digests.
#' @export
setGeneric('load_collections', function(store, digests) standardGeneric('load_collections'))
#' @rdname load_collections
setMethod('load_collections', 'RefgetStore', function(store, digests) {
  invisible(.Call(wrap__load_collections_store, store@ptr, as.character(digests)))
})

#' Download sequences into the shared cache, in parallel, without loading them.
#'
#' Already-cached sequences are skipped; each download is checked against its
#' digest before it is cached.
#' @param store A RefgetStore.
#' @param digests Sequence digests.
#' @param progress TRUE to print a progress line to stderr every 5 s.
#' @return A named list of counters: requests, retries, bytes_downloaded,
#'   sequences_fetched, cache_hits, cache_hit_bytes, local_reads,
#'   local_read_bytes, missing, elapsed_secs.
#' @export
setGeneric('cache_sequences', function(store, digests, progress = TRUE) standardGeneric('cache_sequences'))
#' @rdname cache_sequences
setMethod('cache_sequences', 'RefgetStore', function(store, digests, progress = TRUE) {
  .Call(wrap__cache_sequences_store, store@ptr, as.character(digests), isTRUE(progress))
})

#' Network and cache counters since the store opened (zeros for a local store).
#' @param store A RefgetStore.
#' @return A named list of counters (see \code{cache_sequences}).
#' @export
setGeneric('fetch_stats', function(store) standardGeneric('fetch_stats'))
#' @rdname fetch_stats
setMethod('fetch_stats', 'RefgetStore', function(store) {
  .Call(wrap__fetch_stats_store, store@ptr)
})

# Keep old function name for backwards compatibility
#' @rdname get_sequence
#' @export
setGeneric('get_sequence_by_id', function(store, digest) standardGeneric('get_sequence_by_id'))
#' @rdname get_sequence
setMethod('get_sequence_by_id', 'RefgetStore', function(store, digest) {
  get_sequence(store, digest)
})

#' @rdname get_sequence
#' @export
setGeneric('get_sequence_by_collection_and_name', function(store, collection_digest, sequence_name) standardGeneric('get_sequence_by_collection_and_name'))
#' @rdname get_sequence
setMethod('get_sequence_by_collection_and_name', 'RefgetStore', function(store, collection_digest, sequence_name) {
  get_sequence_by_name(store, collection_digest, sequence_name)
})

# =========================================================================
# Collection Operations
# =========================================================================

#' Browse the contents of a store
#'
#' @description
#' * `list_collections()` / `list_sequences()`: metadata for every
#'   collection or sequence, without sequence data.
#' * `get_collection()`: a collection with its sequences.
#'   `get_collection_metadata()`: a collection's metadata only.
#' * `is_collection_loaded()`: whether a collection's sequences are loaded.
#' * `iter_collections()` / `iter_sequences()`: every collection or sequence,
#'   with data loaded.
#' * `stats()`: counts of sequences and collections, and the storage mode.
#'
#' @param store A RefgetStore object
#' @param digest Collection digest
#' @return A list of \code{SequenceCollectionMetadata}, \code{SequenceMetadata},
#'   \code{SequenceCollection} or \code{SequenceRecord} objects; a single such
#'   object; a logical (\code{is_collection_loaded}); or a named list
#'   (\code{stats}).
#' @name store_contents
NULL

#' @rdname store_contents
#' @export
setGeneric('list_collections', function(store) standardGeneric('list_collections'))
#' @rdname store_contents
setMethod('list_collections', 'RefgetStore', function(store) {
  result <- .Call(wrap__list_collections_store, store@ptr)
  lapply(result, convert_to_collection_metadata)
})

#' @rdname store_contents
#' @export
setGeneric('list_sequences', function(store) standardGeneric('list_sequences'))
#' @rdname store_contents
setMethod('list_sequences', 'RefgetStore', function(store) {
  result <- .Call(wrap__list_sequences_store, store@ptr)
  lapply(result, convert_to_sequence_metadata)
})

#' @rdname store_contents
#' @export
setGeneric('get_collection', function(store, digest) standardGeneric('get_collection'))
#' @rdname store_contents
setMethod('get_collection', 'RefgetStore', function(store, digest) {
  result <- .Call(wrap__get_collection_store, store@ptr, digest)
  convert_to_sequence_collection(result)
})

#' @rdname store_contents
#' @export
setGeneric('get_collection_metadata', function(store, digest) standardGeneric('get_collection_metadata'))
#' @rdname store_contents
setMethod('get_collection_metadata', 'RefgetStore', function(store, digest) {
  result <- .Call(wrap__get_collection_metadata_store, store@ptr, digest)
  convert_to_collection_metadata(result)
})

#' @rdname store_contents
#' @export
setGeneric('is_collection_loaded', function(store, digest) standardGeneric('is_collection_loaded'))
#' @rdname store_contents
setMethod('is_collection_loaded', 'RefgetStore', function(store, digest) {
  .Call(wrap__is_collection_loaded_store, store@ptr, digest)
})

#' @rdname store_contents
#' @export
setGeneric('iter_collections', function(store) standardGeneric('iter_collections'))
#' @rdname store_contents
setMethod('iter_collections', 'RefgetStore', function(store) {
  result <- .Call(wrap__iter_collections_store, store@ptr)
  lapply(result, convert_to_sequence_collection)
})

#' @rdname store_contents
#' @export
setGeneric('iter_sequences', function(store) standardGeneric('iter_sequences'))
#' @rdname store_contents
setMethod('iter_sequences', 'RefgetStore', function(store) {
  result <- .Call(wrap__iter_sequences_store, store@ptr)
  lapply(result, convert_to_sequence_record)
})

#' @rdname store_contents
#' @export
setGeneric('stats', function(store) standardGeneric('stats'))
#' @rdname store_contents
setMethod('stats', 'RefgetStore', function(store) {
  .Call(wrap__stats_store, store@ptr)
})

#' Verify sequences against their stored digests.
#'
#' Streams each sequence's decoded bytes through the normal decode path and
#' recomputes its sha512t24u (and md5, when one is stored) to check that the
#' stored bytes still hash to the digest they are stored under.
#'
#' Named \code{verify_store} (not \code{verify}) to avoid clashing with other
#' packages' \code{verify} generic.
#'
#' @param store A RefgetStore object
#' @param digests Optional character vector of sequence digests (sha512t24u
#'   or md5) to verify. NULL (the default) verifies every sequence in the
#'   store.
#' @param jobs Worker threads. 0 (default) = auto, 1 = serial.
#' @param check_md5 Also recompute and compare md5 when the store has one
#'   recorded (default TRUE).
#' @return A list with n_checked, n_ok, n_failed, and failures (a data.frame
#'   with columns digest, name, alphabet, length, kind,
#'   computed_sha512t24u, stored_md5, computed_md5, error -- one row per
#'   failing sequence).
#' @export
setGeneric('verify_store', function(store, digests = NULL, jobs = 0L, check_md5 = TRUE) standardGeneric('verify_store'))
#' @rdname verify_store
setMethod('verify_store', 'RefgetStore', function(store, digests = NULL, jobs = 0L, check_md5 = TRUE) {
  result <- .Call(wrap__verify_store, store@ptr, digests, as.integer(jobs), isTRUE(check_md5))
  result$failures <- data.frame(
    digest = as.character(result$failures$digest),
    name = as.character(result$failures$name),
    alphabet = as.character(result$failures$alphabet),
    length = as.integer(result$failures$length),
    kind = as.character(result$failures$kind),
    computed_sha512t24u = as.character(result$failures$computed_sha512t24u),
    stored_md5 = as.character(result$failures$stored_md5),
    computed_md5 = as.character(result$failures$computed_md5),
    error = as.character(result$failures$error),
    stringsAsFactors = FALSE
  )
  result
})

# =========================================================================
# Seqcol Spec Operations
# =========================================================================

#' Sequence collection (seqcol) operations
#'
#' @description
#' * `get_level1()`: a collection's attribute digests (seqcol level 1).
#' * `get_level2()`: a collection's full attribute arrays (seqcol level 2).
#' * `compare()`: compare two collections, as in the seqcol comparison
#'   function.
#' * `find_collections_by_attribute()`: digests of the collections whose
#'   attribute has the given digest.
#' * `get_attribute()`: the array for an attribute digest.
#'
#' @param store A RefgetStore object
#' @param digest Collection digest
#' @param digest_a,digest_b Digests of the two collections to compare
#' @param attr_name Attribute name (names, lengths, sequences, ...)
#' @param attr_digest Attribute digest
#' @return A named list (\code{get_level1}, \code{get_level2},
#'   \code{compare}), a character vector of collection digests
#'   (\code{find_collections_by_attribute}), or the attribute array
#'   (\code{get_attribute}).
#' @name seqcol_operations
NULL

#' @rdname seqcol_operations
#' @export
setGeneric('get_level1', function(store, digest) standardGeneric('get_level1'))
#' @rdname seqcol_operations
setMethod('get_level1', 'RefgetStore', function(store, digest) {
  .Call(wrap__get_collection_level1_store, store@ptr, digest)
})

#' @rdname seqcol_operations
#' @export
setGeneric('get_level2', function(store, digest) standardGeneric('get_level2'))
#' @rdname seqcol_operations
setMethod('get_level2', 'RefgetStore', function(store, digest) {
  .Call(wrap__get_collection_level2_store, store@ptr, digest)
})

#' @rdname seqcol_operations
#' @export
setGeneric('compare', function(store, digest_a, digest_b) standardGeneric('compare'))
#' @rdname seqcol_operations
setMethod('compare', 'RefgetStore', function(store, digest_a, digest_b) {
  .Call(wrap__compare_store, store@ptr, digest_a, digest_b)
})

#' @rdname seqcol_operations
#' @export
setGeneric('find_collections_by_attribute', function(store, attr_name, attr_digest) standardGeneric('find_collections_by_attribute'))
#' @rdname seqcol_operations
setMethod('find_collections_by_attribute', 'RefgetStore', function(store, attr_name, attr_digest) {
  .Call(wrap__find_collections_by_attribute_store, store@ptr, attr_name, attr_digest)
})

#' @rdname seqcol_operations
#' @export
setGeneric('get_attribute', function(store, attr_name, attr_digest) standardGeneric('get_attribute'))
#' @rdname seqcol_operations
setMethod('get_attribute', 'RefgetStore', function(store, attr_name, attr_digest) {
  .Call(wrap__get_attribute_store, store@ptr, attr_name, attr_digest)
})

#' @rdname store_settings
#' @export
setGeneric('enable_ancillary_digests', function(store) standardGeneric('enable_ancillary_digests'))
#' @rdname store_settings
setMethod('enable_ancillary_digests', 'RefgetStore', function(store) {
  .Call(wrap__enable_ancillary_digests_store, store@ptr)
  invisible(store)
})

#' @rdname store_settings
#' @export
setGeneric('disable_ancillary_digests', function(store) standardGeneric('disable_ancillary_digests'))
#' @rdname store_settings
setMethod('disable_ancillary_digests', 'RefgetStore', function(store) {
  .Call(wrap__disable_ancillary_digests_store, store@ptr)
  invisible(store)
})

#' @rdname store_settings
#' @export
setGeneric('has_ancillary_digests', function(store) standardGeneric('has_ancillary_digests'))
#' @rdname store_settings
setMethod('has_ancillary_digests', 'RefgetStore', function(store) {
  .Call(wrap__has_ancillary_digests_store, store@ptr)
})

#' @rdname store_settings
#' @export
setGeneric('enable_attribute_index', function(store) standardGeneric('enable_attribute_index'))
#' @rdname store_settings
setMethod('enable_attribute_index', 'RefgetStore', function(store) {
  .Call(wrap__enable_attribute_index_store, store@ptr)
  invisible(store)
})

#' @rdname store_settings
#' @export
setGeneric('disable_attribute_index', function(store) standardGeneric('disable_attribute_index'))
#' @rdname store_settings
setMethod('disable_attribute_index', 'RefgetStore', function(store) {
  .Call(wrap__disable_attribute_index_store, store@ptr)
  invisible(store)
})

#' @rdname store_settings
#' @export
setGeneric('has_attribute_index', function(store) standardGeneric('has_attribute_index'))
#' @rdname store_settings
setMethod('has_attribute_index', 'RefgetStore', function(store) {
  .Call(wrap__has_attribute_index_store, store@ptr)
})

# =========================================================================
# Write/Export Operations
# =========================================================================

#' Write a store to disk
#'
#' @description `write_store_to_directory()` writes the store to a directory.
#'   `write_store()` writes it to the path it was created or opened with.
#'
#' @param store A RefgetStore object
#' @param root_path Directory to write the store to
#' @param seqdata_path_template Template for sequence file paths; "" for the
#'   default
#' @return The store, invisibly
#' @rdname write_store
#' @export
setGeneric('write_store_to_directory', function(store, root_path, seqdata_path_template = '') standardGeneric('write_store_to_directory'))
#' @rdname write_store
setMethod('write_store_to_directory', 'RefgetStore', function(store, root_path, seqdata_path_template = '') {
  .Call(wrap__write_store_to_directory_store, store@ptr, root_path, seqdata_path_template)
  invisible(store)
})

#' @rdname write_store
#' @export
setGeneric('write_store', function(store) standardGeneric('write_store'))
#' @rdname write_store
setMethod('write_store', 'RefgetStore', function(store) {
  .Call(wrap__write_store, store@ptr)
  invisible(store)
})

#' Reclaim dead bytes in a packed store (alpha).
#'
#' Removing collections leaves their bytes inside pack files. This rewrites
#' the live sequences into fresh packs and deletes the old ones. Reads return
#' the same bases before and after. Errors on a store that is not packed.
#'
#' @param store A RefgetStore object created with \code{packed = TRUE}
#' @param prune_unreferenced Also delete pack files that no index entry points
#'   to, such as packs left by a crashed writer. Only use this when no other
#'   process is writing to the store: a pack another writer has sealed but not
#'   committed yet is deleted too.
#' @return A named list: packs_before, packs_after, bytes_reclaimed,
#'   packs_pruned (pack files deleted by \code{prune_unreferenced})
#' @export
setGeneric('compact', function(store, prune_unreferenced = FALSE) standardGeneric('compact'))
#' @rdname compact
setMethod('compact', 'RefgetStore', function(store, prune_unreferenced = FALSE) {
  .Call(wrap__compact_store, store@ptr, isTRUE(prune_unreferenced))
})

#' Turn the pack block cache on or off (alpha).
#'
#' Off by default. When on, reads from a packed store's local pack files go
#' through a cache of 1 MiB blocks. It helps clustered reads (for example a
#' sorted BED) on a cold network file system, and hurts scattered reads on
#' stores larger than the cache.
#'
#' @param store A RefgetStore object
#' @param megabytes Cache size in megabytes; 0 turns the cache off
#' @return The store, invisibly
#' @export
setGeneric('set_pack_cache', function(store, megabytes) standardGeneric('set_pack_cache'))
#' @rdname set_pack_cache
setMethod('set_pack_cache', 'RefgetStore', function(store, megabytes) {
  .Call(wrap__set_pack_cache_store, store@ptr, as.integer(megabytes))
  invisible(store)
})

#' Export sequences as FASTA, to a file or to stdout.
#'
#' Choose records with exactly one of \code{collections}, \code{seq_digests}
#' or \code{all}; choose how they are labeled with \code{label}. Exports
#' read sequences ahead of the writer in parallel, within a fixed memory
#' budget, and leave the store's loaded state unchanged.
#'
#' Label fields: \code{{seq}} sequence digest (\code{SQ.}...), the default;
#' \code{{seq:NS}} the sequence's alias in namespace NS; \code{{col}}
#' collection digest; \code{{col:NS}} the collection's alias in namespace NS;
#' \code{{name}} the sequence's name in its collection. A missing alias falls
#' back to the digest. The
#' same label for the same sequence is written once; the same label for
#' different sequences is an error, raised before any output.
#'
#' @param store A RefgetStore object
#' @param output_path Output FASTA path ('.gz' compresses), or "-" for stdout
#' @param collections Collection digest(s) to export
#' @param names Only these sequence names in each selected collection
#' @param seq_digests Sequence digests to export (no collection context)
#' @param all Export every collection plus sequences in no collection
#' @param label Label template (default `"{seq}"`)
#' @param molecule NULL, "dna" or "protein"
#' @param line_width Bases per line (default 80; 0 = one line per sequence)
#' @param write_cache Write sequences downloaded from a remote store to the
#'   shared cache, so a rerun downloads only what is missing
#' @return A named list of counts: n_records, n_bases, n_fallbacks,
#'   n_deduplicated, and \code{fetch}, the network and cache counters
#' @export
setGeneric('export_fasta', function(store, output_path, collections = NULL, names = NULL, seq_digests = NULL, all = FALSE, label = "{seq}", molecule = NULL, line_width = 80L, write_cache = FALSE) standardGeneric('export_fasta'))
#' @rdname export_fasta
setMethod('export_fasta', 'RefgetStore', function(store, output_path, collections = NULL, names = NULL, seq_digests = NULL, all = FALSE, label = "{seq}", molecule = NULL, line_width = 80L, write_cache = FALSE) {
  invisible(.Call(wrap__export_fasta_store, store@ptr, output_path, collections, names, seq_digests, all, label, molecule, as.integer(line_width), as.logical(write_cache)))
})

#' Export BED regions from a collection as FASTA.
#'
#' Writes one record per region, labeled by \code{label} (see
#' \code{export_fasta}) with \code{:start-end} appended. The default label is
#' the sequence name, so headers look like \code{>chr1:0-10}.
#'
#' @param store A RefgetStore object
#' @param collection_digest Sequence collection digest
#' @param bed_file_path Path to BED file
#' @param output_file_path Output FASTA path, or "-" for stdout
#' @param label Label template (default `"{name}"`); `"{seq}"` labels regions by
#'   sequence digest
#' @return A named list of counts
#' @export
setGeneric('get_seqs_bed_file', function(store, collection_digest, bed_file_path, output_file_path, label = "{name}") standardGeneric('get_seqs_bed_file'))
#' @rdname get_seqs_bed_file
setMethod('get_seqs_bed_file', 'RefgetStore', function(store, collection_digest, bed_file_path, output_file_path, label = "{name}") {
  invisible(.Call(wrap__get_seqs_bed_file_store, store@ptr, collection_digest, bed_file_path, output_file_path, label))
})

#' Extract BED file sequences as a list
#'
#' @description Returns one \code{RetrievedSequence} object per BED region.
#'   For many regions, \code{\link{get_seqs_bed_file_to_df}} is faster.
#'
#' @param store A RefgetStore object
#' @param collection_digest Sequence collection digest
#' @param bed_file_path Path to BED file
#' @return A list of \code{RetrievedSequence} objects
#' @export
setGeneric('get_seqs_bed_file_to_vec', function(store, collection_digest, bed_file_path) standardGeneric('get_seqs_bed_file_to_vec'))
#' @rdname get_seqs_bed_file_to_vec
setMethod('get_seqs_bed_file_to_vec', 'RefgetStore', function(store, collection_digest, bed_file_path) {
  result <- .Call(wrap__get_seqs_bed_file_to_vec_store, store@ptr, collection_digest, bed_file_path)
  if (!is.null(result) && length(result) > 0) {
    lapply(result, convert_to_retrieved_sequence)
  } else {
    list()
  }
})

#' Extract BED file sequences as a data.frame (flat structure-of-arrays).
#'
#' Returns one data.frame with columns sequence (character), chrom_name
#' (character), start (integer), end (integer) -- one row per region. This avoids
#' the per-region S4 construction of \code{get_seqs_bed_file_to_vec}.
#'
#' @param store A RefgetStore object
#' @param collection_digest Sequence collection digest
#' @param bed_file_path Path to BED file
#' @return A data.frame with columns sequence, chrom_name, start, end
#' @export
setGeneric('get_seqs_bed_file_to_df', function(store, collection_digest, bed_file_path) standardGeneric('get_seqs_bed_file_to_df'))
#' @rdname get_seqs_bed_file_to_df
setMethod('get_seqs_bed_file_to_df', 'RefgetStore', function(store, collection_digest, bed_file_path) {
  result <- .Call(wrap__get_seqs_bed_file_to_df_store, store@ptr, collection_digest, bed_file_path)
  data.frame(
    sequence = as.character(result$sequence),
    chrom_name = as.character(result$chrom_name),
    start = as.integer(result$start),
    end = as.integer(result$end),
    stringsAsFactors = FALSE
  )
})

#' Extract sequences for regions given as parallel vectors (flat structure-of-arrays).
#'
#' Pass parallel \code{chroms}/\code{starts}/\code{ends} vectors directly into the
#' store with no temp BED file written to disk. Returns one data.frame with columns
#' sequence (character), chrom_name (character), start (integer), end (integer) --
#' one row per region. Coordinates are 0-based, half-open (BED convention): a region
#' with start=0, end=4 yields the same bytes as a BED line \code{chrom\\t0\\t4}.
#'
#' @param store A RefgetStore object
#' @param collection_digest Sequence collection digest
#' @param chroms Character vector of chromosome / sequence names
#' @param starts Integer vector of 0-based start coordinates (BED convention)
#' @param ends Integer vector of 0-based half-open end coordinates (BED convention)
#' @return A data.frame with columns sequence, chrom_name, start, end
#' @export
setGeneric('get_seqs_from_regions', function(store, collection_digest, chroms, starts, ends) standardGeneric('get_seqs_from_regions'))
#' @rdname get_seqs_from_regions
setMethod('get_seqs_from_regions', 'RefgetStore', function(store, collection_digest, chroms, starts, ends) {
  stopifnot(length(chroms) == length(starts), length(chroms) == length(ends))
  result <- .Call(
    wrap__get_seqs_from_region_vectors_store,
    store@ptr, collection_digest,
    as.character(chroms), as.integer(starts), as.integer(ends)
  )
  data.frame(
    sequence = as.character(result$sequence),
    chrom_name = as.character(result$chrom_name),
    start = as.integer(result$start),
    end = as.integer(result$end),
    stringsAsFactors = FALSE
  )
})

# =========================================================================
# Sequence Alias Operations
# =========================================================================

#' Sequence and collection aliases
#'
#' @description Aliases are human-readable names for sequence or collection
#'   digests, grouped into namespaces (for example \code{"refseq"} or
#'   \code{"ucsc"}). The \code{*_sequence_*} functions work on sequence
#'   aliases and the \code{*_collection_*} functions on collection aliases.
#'
#' @details
#' * `add_*_alias()`: add an alias for a digest.
#' * `get_*_metadata_by_alias()`: look up metadata by alias.
#' * `get_aliases_for_*()`: all aliases of a digest.
#' * `list_*_alias_namespaces()` / `list_*_aliases()`: list namespaces, or
#'   the aliases in one namespace.
#' * `remove_*_alias()`: remove an alias.
#' * `load_*_aliases()`: load aliases into a namespace from a TSV file.
#'
#' @param store A RefgetStore object
#' @param namespace Alias namespace
#' @param alias Alias name
#' @param digest Sequence or collection digest
#' @param path Path to a TSV file of aliases
#' @return \code{add_*_alias} returns the store, invisibly. Lookups return
#'   a metadata object or NULL; listings return character vectors;
#'   \code{remove_*_alias} returns TRUE if an alias was removed; and
#'   \code{load_*_aliases} returns the number of aliases loaded.
#' @name refget_aliases
NULL

#' @rdname refget_aliases
#' @export
setGeneric('add_sequence_alias', function(store, namespace, alias, digest) standardGeneric('add_sequence_alias'))
#' @rdname refget_aliases
setMethod('add_sequence_alias', 'RefgetStore', function(store, namespace, alias, digest) {
  .Call(wrap__add_sequence_alias_store, store@ptr, namespace, alias, digest)
  invisible(store)
})

#' @rdname refget_aliases
#' @export
setGeneric('get_sequence_metadata_by_alias', function(store, namespace, alias) standardGeneric('get_sequence_metadata_by_alias'))
#' @rdname refget_aliases
setMethod('get_sequence_metadata_by_alias', 'RefgetStore', function(store, namespace, alias) {
  result <- .Call(wrap__get_sequence_metadata_by_alias_store, store@ptr, namespace, alias)
  convert_to_sequence_metadata(result)
})

#' @rdname refget_aliases
#' @export
setGeneric('get_aliases_for_sequence', function(store, digest) standardGeneric('get_aliases_for_sequence'))
#' @rdname refget_aliases
setMethod('get_aliases_for_sequence', 'RefgetStore', function(store, digest) {
  .Call(wrap__get_aliases_for_sequence_store, store@ptr, digest)
})

#' @rdname refget_aliases
#' @export
setGeneric('list_sequence_alias_namespaces', function(store) standardGeneric('list_sequence_alias_namespaces'))
#' @rdname refget_aliases
setMethod('list_sequence_alias_namespaces', 'RefgetStore', function(store) {
  .Call(wrap__list_sequence_alias_namespaces_store, store@ptr)
})

#' @rdname refget_aliases
#' @export
setGeneric('list_sequence_aliases', function(store, namespace) standardGeneric('list_sequence_aliases'))
#' @rdname refget_aliases
setMethod('list_sequence_aliases', 'RefgetStore', function(store, namespace) {
  .Call(wrap__list_sequence_aliases_store, store@ptr, namespace)
})

#' @rdname refget_aliases
#' @export
setGeneric('remove_sequence_alias', function(store, namespace, alias) standardGeneric('remove_sequence_alias'))
#' @rdname refget_aliases
setMethod('remove_sequence_alias', 'RefgetStore', function(store, namespace, alias) {
  .Call(wrap__remove_sequence_alias_store, store@ptr, namespace, alias)
})

#' @rdname refget_aliases
#' @export
setGeneric('load_sequence_aliases', function(store, namespace, path) standardGeneric('load_sequence_aliases'))
#' @rdname refget_aliases
setMethod('load_sequence_aliases', 'RefgetStore', function(store, namespace, path) {
  .Call(wrap__load_sequence_aliases_store, store@ptr, namespace, path)
})

# =========================================================================
# Collection Alias Operations
# =========================================================================

#' @rdname refget_aliases
#' @export
setGeneric('add_collection_alias', function(store, namespace, alias, digest) standardGeneric('add_collection_alias'))
#' @rdname refget_aliases
setMethod('add_collection_alias', 'RefgetStore', function(store, namespace, alias, digest) {
  .Call(wrap__add_collection_alias_store, store@ptr, namespace, alias, digest)
  invisible(store)
})

#' @rdname refget_aliases
#' @export
setGeneric('get_collection_metadata_by_alias', function(store, namespace, alias) standardGeneric('get_collection_metadata_by_alias'))
#' @rdname refget_aliases
setMethod('get_collection_metadata_by_alias', 'RefgetStore', function(store, namespace, alias) {
  result <- .Call(wrap__get_collection_metadata_by_alias_store, store@ptr, namespace, alias)
  convert_to_collection_metadata(result)
})

#' @rdname refget_aliases
#' @export
setGeneric('get_aliases_for_collection', function(store, digest) standardGeneric('get_aliases_for_collection'))
#' @rdname refget_aliases
setMethod('get_aliases_for_collection', 'RefgetStore', function(store, digest) {
  .Call(wrap__get_aliases_for_collection_store, store@ptr, digest)
})

#' @rdname refget_aliases
#' @export
setGeneric('list_collection_alias_namespaces', function(store) standardGeneric('list_collection_alias_namespaces'))
#' @rdname refget_aliases
setMethod('list_collection_alias_namespaces', 'RefgetStore', function(store) {
  .Call(wrap__list_collection_alias_namespaces_store, store@ptr)
})

#' @rdname refget_aliases
#' @export
setGeneric('list_collection_aliases', function(store, namespace) standardGeneric('list_collection_aliases'))
#' @rdname refget_aliases
setMethod('list_collection_aliases', 'RefgetStore', function(store, namespace) {
  .Call(wrap__list_collection_aliases_store, store@ptr, namespace)
})

#' @rdname refget_aliases
#' @export
setGeneric('remove_collection_alias', function(store, namespace, alias) standardGeneric('remove_collection_alias'))
#' @rdname refget_aliases
setMethod('remove_collection_alias', 'RefgetStore', function(store, namespace, alias) {
  .Call(wrap__remove_collection_alias_store, store@ptr, namespace, alias)
})

#' @rdname refget_aliases
#' @export
setGeneric('load_collection_aliases', function(store, namespace, path) standardGeneric('load_collection_aliases'))
#' @rdname refget_aliases
setMethod('load_collection_aliases', 'RefgetStore', function(store, namespace, path) {
  .Call(wrap__load_collection_aliases_store, store@ptr, namespace, path)
})

# =========================================================================
# FHR Metadata Operations
# =========================================================================

#' FHR metadata for collections
#'
#' @description Attach, read and remove FAIR Headers Reference genome (FHR)
#'   metadata for a sequence collection.
#'
#' @param store A RefgetStore object
#' @param collection_digest Collection digest
#' @param metadata_json FHR metadata as a JSON string, or an R list that is
#'   converted with \code{jsonlite::toJSON}
#' @param path Path to an FHR JSON file
#' @return \code{set_fhr_metadata} and \code{load_fhr_metadata} return the
#'   store, invisibly. \code{get_fhr_metadata} returns the parsed metadata or
#'   NULL; \code{remove_fhr_metadata} returns TRUE if metadata was removed;
#'   \code{list_fhr_metadata} returns the digests of collections with
#'   metadata.
#' @name fhr_metadata
NULL

#' @rdname fhr_metadata
#' @export
setGeneric('set_fhr_metadata', function(store, collection_digest, metadata_json) standardGeneric('set_fhr_metadata'))
#' @rdname fhr_metadata
setMethod('set_fhr_metadata', 'RefgetStore', function(store, collection_digest, metadata_json) {
  json_str <- if (is.character(metadata_json)) {
    metadata_json
  } else {
    jsonlite::toJSON(metadata_json, auto_unbox = TRUE)
  }
  .Call(wrap__set_fhr_metadata_store, store@ptr, collection_digest, json_str)
  invisible(store)
})

#' @rdname fhr_metadata
#' @export
setGeneric('get_fhr_metadata', function(store, collection_digest) standardGeneric('get_fhr_metadata'))
#' @rdname fhr_metadata
setMethod('get_fhr_metadata', 'RefgetStore', function(store, collection_digest) {
  result <- .Call(wrap__get_fhr_metadata_store, store@ptr, collection_digest)
  if (is.null(result) || result == '') {
    NULL
  } else {
    jsonlite::fromJSON(result)
  }
})

#' @rdname fhr_metadata
#' @export
setGeneric('remove_fhr_metadata', function(store, collection_digest) standardGeneric('remove_fhr_metadata'))
#' @rdname fhr_metadata
setMethod('remove_fhr_metadata', 'RefgetStore', function(store, collection_digest) {
  .Call(wrap__remove_fhr_metadata_store, store@ptr, collection_digest)
})

#' @rdname fhr_metadata
#' @export
setGeneric('list_fhr_metadata', function(store) standardGeneric('list_fhr_metadata'))
#' @rdname fhr_metadata
setMethod('list_fhr_metadata', 'RefgetStore', function(store) {
  .Call(wrap__list_fhr_metadata_store, store@ptr)
})

#' @rdname fhr_metadata
#' @export
setGeneric('load_fhr_metadata', function(store, collection_digest, path) standardGeneric('load_fhr_metadata'))
#' @rdname fhr_metadata
setMethod('load_fhr_metadata', 'RefgetStore', function(store, collection_digest, path) {
  .Call(wrap__load_fhr_metadata_store, store@ptr, collection_digest, path)
  invisible(store)
})
