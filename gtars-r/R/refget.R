#' @useDynLib gtars, .registration = TRUE
NULL

# =========================================================================
# S4 Class Definitions
# =========================================================================

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

setClass(
  'SequenceRecord',
  slots = list(
    metadata = 'SequenceMetadata',
    data = 'character'
  )
)

setClass(
  'SeqColDigestLvl1',
  slots = list(
    sequences_digest = 'character',
    names_digest = 'character',
    lengths_digest = 'character'
  )
)

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

setClass(
  'SequenceCollection',
  slots = list(
    sequences = 'list',
    digest = 'character',
    lvl1 = 'SeqColDigestLvl1',
    file_path = 'character'
  )
)

setClass(
  'RefgetStore',
  slots = list(
    ptr = 'externalptr'
  )
)

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

setMethod('show', 'SequenceRecord', function(object) {
  cat(sprintf('SequenceRecord for %s\n', object@metadata@name))
})

setMethod('show', 'SeqColDigestLvl1', function(object) {
  cat('SeqColDigestLvl1: \n')
  cat(sprintf('  sequences digest: %s\n', object@sequences_digest))
  cat(sprintf('  names digest: %s\n', object@names_digest))
  cat(sprintf('  lengths digest: %s\n', object@lengths_digest))
})

setMethod('show', 'SequenceCollectionMetadata', function(object) {
  cat(sprintf('SequenceCollectionMetadata: %s (%d sequences)\n',
              object@digest, object@n_sequences))
})

setMethod('show', 'SequenceCollection', function(object) {
  cat(sprintf('SequenceCollection with %d sequences, digest = %s\n',
              length(object@sequences), object@digest))
})

setMethod('length', 'SequenceCollection', function(x) {
  length(x@sequences)
})

setMethod('show', 'RetrievedSequence', function(object) {
  cat(sprintf('RetrievedSequence: %s|%d-%d: %s\n',
              object@chrom_name, object@start, object@end, object@sequence))
})

setMethod('print', 'RetrievedSequence', function(x) {
  cat(sprintf('%s|%d-%d: %s\n', x@chrom_name, x@start, x@end, x@sequence))
})

setMethod('[', c('SequenceCollection', 'numeric', 'missing'),
          function(x, i) {
            if (i < 1 || i > length(x@sequences)) {
              stop('Index out of range')
            }
            return(x@sequences[[i]])
          })

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
#' @description Create a RefgetStore that persists to disk.
#'
#' @param path Directory path for storing sequences and metadata
#'
#' @return A RefgetStore object
#'
#' @examples
#' \dontrun{
#' store <- refget_store_on_disk("/path/to/store")
#' }
#'
#' @export
refget_store_on_disk <- function(path) {
  result <- .Call(wrap__on_disk_store, path)
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
#' @description Open a remote RefgetStore with local caching.
#'
#' @param cache_path Local directory for caching downloaded data
#' @param remote_url URL of the remote RefgetStore
#'
#' @return A RefgetStore object
#'
#' @examples
#' \dontrun{
#' store <- refget_store_open_remote("/local/cache", "https://example.com/store")
#' }
#'
#' @export
refget_store_open_remote <- function(cache_path, remote_url) {
  result <- .Call(wrap__open_remote_store, cache_path, remote_url)
  invisible(new('RefgetStore', ptr = result))
}

# Keep old function name for backwards compatibility
#' @export
load_from_directory <- refget_store_open_local

# =========================================================================
# Store Configuration Methods
# =========================================================================

#' @export
setGeneric('enable_encoding', function(store) standardGeneric('enable_encoding'))
setMethod('enable_encoding', 'RefgetStore', function(store) {
  .Call(wrap__enable_encoding_store, store@ptr)
  invisible(store)
})

#' @export
setGeneric('disable_encoding', function(store) standardGeneric('disable_encoding'))
setMethod('disable_encoding', 'RefgetStore', function(store) {
  .Call(wrap__disable_encoding_store, store@ptr)
  invisible(store)
})

#' @export
setGeneric('storage_mode', function(store) standardGeneric('storage_mode'))
setMethod('storage_mode', 'RefgetStore', function(store) {
  .Call(wrap__get_storage_mode_store, store@ptr)
})

#' @export
setGeneric('set_quiet', function(store, quiet) standardGeneric('set_quiet'))
setMethod('set_quiet', 'RefgetStore', function(store, quiet) {
  .Call(wrap__set_quiet_store, store@ptr, quiet)
  invisible(store)
})

#' @export
setGeneric('is_quiet', function(store) standardGeneric('is_quiet'))
setMethod('is_quiet', 'RefgetStore', function(store) {
  .Call(wrap__get_quiet_store, store@ptr)
})

#' @export
setGeneric('is_persisting', function(store) standardGeneric('is_persisting'))
setMethod('is_persisting', 'RefgetStore', function(store) {
  .Call(wrap__is_persisting_store, store@ptr)
})

#' @export
setGeneric('enable_persistence', function(store, path) standardGeneric('enable_persistence'))
setMethod('enable_persistence', 'RefgetStore', function(store, path) {
  .Call(wrap__enable_persistence_store, store@ptr, path)
  invisible(store)
})

#' @export
setGeneric('disable_persistence', function(store) standardGeneric('disable_persistence'))
setMethod('disable_persistence', 'RefgetStore', function(store) {
  .Call(wrap__disable_persistence_store, store@ptr)
  invisible(store)
})

#' @export
setGeneric('cache_path', function(store) standardGeneric('cache_path'))
setMethod('cache_path', 'RefgetStore', function(store) {
  .Call(wrap__get_cache_path_store, store@ptr)
})

#' @export
setGeneric('remote_url', function(store) standardGeneric('remote_url'))
setMethod('remote_url', 'RefgetStore', function(store) {
  .Call(wrap__get_remote_url_store, store@ptr)
})

# =========================================================================
# Adding Sequences/Collections
# =========================================================================

#' @export
setGeneric('import_fasta', function(store, file_path) standardGeneric('import_fasta'))
setMethod('import_fasta', 'RefgetStore', function(store, file_path) {
  .Call(wrap__import_fasta_store, store@ptr, file_path)
  invisible(store)
})

#' @export
setGeneric('add_fasta', function(store, file_path, force = FALSE) standardGeneric('add_fasta'))
setMethod('add_fasta', 'RefgetStore', function(store, file_path, force = FALSE) {
  result <- .Call(wrap__add_fasta_store, store@ptr, file_path, force)
  result
})

# =========================================================================
# Sequence Retrieval
# =========================================================================

#' @export
setGeneric('get_sequence', function(store, digest) standardGeneric('get_sequence'))
setMethod('get_sequence', 'RefgetStore', function(store, digest) {
  result <- .Call(wrap__get_sequence_store, store@ptr, digest)
  convert_to_sequence_record(result)
})

#' @export
setGeneric('get_sequence_by_name', function(store, collection_digest, sequence_name) standardGeneric('get_sequence_by_name'))
setMethod('get_sequence_by_name', 'RefgetStore', function(store, collection_digest, sequence_name) {
  result <- .Call(wrap__get_sequence_by_name_store, store@ptr, collection_digest, sequence_name)
  convert_to_sequence_record(result)
})

#' @export
setGeneric('get_sequence_metadata', function(store, digest) standardGeneric('get_sequence_metadata'))
setMethod('get_sequence_metadata', 'RefgetStore', function(store, digest) {
  result <- .Call(wrap__get_sequence_metadata_store, store@ptr, digest)
  convert_to_sequence_metadata(result)
})

#' @export
setGeneric('get_substring', function(store, seq_digest, start, end) standardGeneric('get_substring'))
setMethod('get_substring', 'RefgetStore', function(store, seq_digest, start, end) {
  result <- .Call(wrap__get_substring_store, store@ptr, seq_digest, as.integer(start), as.integer(end))
  if (is.null(result)) NULL else as.character(result)
})

# Keep old function name for backwards compatibility
#' @export
setGeneric('get_sequence_by_id', function(store, digest) standardGeneric('get_sequence_by_id'))
setMethod('get_sequence_by_id', 'RefgetStore', function(store, digest) {
  get_sequence(store, digest)
})

#' @export
setGeneric('get_sequence_by_collection_and_name', function(store, collection_digest, sequence_name) standardGeneric('get_sequence_by_collection_and_name'))
setMethod('get_sequence_by_collection_and_name', 'RefgetStore', function(store, collection_digest, sequence_name) {
  get_sequence_by_name(store, collection_digest, sequence_name)
})

# =========================================================================
# Collection Operations
# =========================================================================

#' @export
setGeneric('list_collections', function(store) standardGeneric('list_collections'))
setMethod('list_collections', 'RefgetStore', function(store) {
  result <- .Call(wrap__list_collections_store, store@ptr)
  lapply(result, convert_to_collection_metadata)
})

#' @export
setGeneric('list_sequences', function(store) standardGeneric('list_sequences'))
setMethod('list_sequences', 'RefgetStore', function(store) {
  result <- .Call(wrap__list_sequences_store, store@ptr)
  lapply(result, convert_to_sequence_metadata)
})

#' @export
setGeneric('get_collection', function(store, digest) standardGeneric('get_collection'))
setMethod('get_collection', 'RefgetStore', function(store, digest) {
  result <- .Call(wrap__get_collection_store, store@ptr, digest)
  convert_to_sequence_collection(result)
})

#' @export
setGeneric('get_collection_metadata', function(store, digest) standardGeneric('get_collection_metadata'))
setMethod('get_collection_metadata', 'RefgetStore', function(store, digest) {
  result <- .Call(wrap__get_collection_metadata_store, store@ptr, digest)
  convert_to_collection_metadata(result)
})

#' @export
setGeneric('is_collection_loaded', function(store, digest) standardGeneric('is_collection_loaded'))
setMethod('is_collection_loaded', 'RefgetStore', function(store, digest) {
  .Call(wrap__is_collection_loaded_store, store@ptr, digest)
})

#' @export
setGeneric('iter_collections', function(store) standardGeneric('iter_collections'))
setMethod('iter_collections', 'RefgetStore', function(store) {
  result <- .Call(wrap__iter_collections_store, store@ptr)
  lapply(result, convert_to_sequence_collection)
})

#' @export
setGeneric('iter_sequences', function(store) standardGeneric('iter_sequences'))
setMethod('iter_sequences', 'RefgetStore', function(store) {
  result <- .Call(wrap__iter_sequences_store, store@ptr)
  lapply(result, convert_to_sequence_record)
})

#' @export
setGeneric('stats', function(store) standardGeneric('stats'))
setMethod('stats', 'RefgetStore', function(store) {
  .Call(wrap__stats_store, store@ptr)
})

# =========================================================================
# Seqcol Spec Operations
# =========================================================================

#' @export
setGeneric('get_level1', function(store, digest) standardGeneric('get_level1'))
setMethod('get_level1', 'RefgetStore', function(store, digest) {
  .Call(wrap__get_collection_level1_store, store@ptr, digest)
})

#' @export
setGeneric('get_level2', function(store, digest) standardGeneric('get_level2'))
setMethod('get_level2', 'RefgetStore', function(store, digest) {
  .Call(wrap__get_collection_level2_store, store@ptr, digest)
})

#' @export
setGeneric('compare', function(store, digest_a, digest_b) standardGeneric('compare'))
setMethod('compare', 'RefgetStore', function(store, digest_a, digest_b) {
  .Call(wrap__compare_store, store@ptr, digest_a, digest_b)
})

#' @export
setGeneric('find_collections_by_attribute', function(store, attr_name, attr_digest) standardGeneric('find_collections_by_attribute'))
setMethod('find_collections_by_attribute', 'RefgetStore', function(store, attr_name, attr_digest) {
  .Call(wrap__find_collections_by_attribute_store, store@ptr, attr_name, attr_digest)
})

#' @export
setGeneric('get_attribute', function(store, attr_name, attr_digest) standardGeneric('get_attribute'))
setMethod('get_attribute', 'RefgetStore', function(store, attr_name, attr_digest) {
  .Call(wrap__get_attribute_store, store@ptr, attr_name, attr_digest)
})

#' @export
setGeneric('enable_ancillary_digests', function(store) standardGeneric('enable_ancillary_digests'))
setMethod('enable_ancillary_digests', 'RefgetStore', function(store) {
  .Call(wrap__enable_ancillary_digests_store, store@ptr)
  invisible(store)
})

#' @export
setGeneric('disable_ancillary_digests', function(store) standardGeneric('disable_ancillary_digests'))
setMethod('disable_ancillary_digests', 'RefgetStore', function(store) {
  .Call(wrap__disable_ancillary_digests_store, store@ptr)
  invisible(store)
})

#' @export
setGeneric('has_ancillary_digests', function(store) standardGeneric('has_ancillary_digests'))
setMethod('has_ancillary_digests', 'RefgetStore', function(store) {
  .Call(wrap__has_ancillary_digests_store, store@ptr)
})

#' @export
setGeneric('enable_attribute_index', function(store) standardGeneric('enable_attribute_index'))
setMethod('enable_attribute_index', 'RefgetStore', function(store) {
  .Call(wrap__enable_attribute_index_store, store@ptr)
  invisible(store)
})

#' @export
setGeneric('disable_attribute_index', function(store) standardGeneric('disable_attribute_index'))
setMethod('disable_attribute_index', 'RefgetStore', function(store) {
  .Call(wrap__disable_attribute_index_store, store@ptr)
  invisible(store)
})

#' @export
setGeneric('has_attribute_index', function(store) standardGeneric('has_attribute_index'))
setMethod('has_attribute_index', 'RefgetStore', function(store) {
  .Call(wrap__has_attribute_index_store, store@ptr)
})

# =========================================================================
# Write/Export Operations
# =========================================================================

#' @export
setGeneric('write_store_to_directory', function(store, root_path, seqdata_path_template = '') standardGeneric('write_store_to_directory'))
setMethod('write_store_to_directory', 'RefgetStore', function(store, root_path, seqdata_path_template = '') {
  .Call(wrap__write_store_to_directory_store, store@ptr, root_path, seqdata_path_template)
  invisible(store)
})

#' @export
setGeneric('write_store', function(store) standardGeneric('write_store'))
setMethod('write_store', 'RefgetStore', function(store) {
  .Call(wrap__write_store, store@ptr)
  invisible(store)
})

#' @export
setGeneric('export_fasta', function(store, collection_digest, output_path, sequence_names = NULL, line_width = 80L) standardGeneric('export_fasta'))
setMethod('export_fasta', 'RefgetStore', function(store, collection_digest, output_path, sequence_names = NULL, line_width = 80L) {
  .Call(wrap__export_fasta_store, store@ptr, collection_digest, output_path, sequence_names, as.integer(line_width))
  invisible(store)
})

#' @export
setGeneric('export_fasta_by_digests', function(store, seq_digests, output_path, line_width = 80L) standardGeneric('export_fasta_by_digests'))
setMethod('export_fasta_by_digests', 'RefgetStore', function(store, seq_digests, output_path, line_width = 80L) {
  .Call(wrap__export_fasta_by_digests_store, store@ptr, seq_digests, output_path, as.integer(line_width))
  invisible(store)
})

#' @export
setGeneric('get_seqs_bed_file', function(store, collection_digest, bed_file_path, output_file_path) standardGeneric('get_seqs_bed_file'))
setMethod('get_seqs_bed_file', 'RefgetStore', function(store, collection_digest, bed_file_path, output_file_path) {
  .Call(wrap__get_seqs_bed_file_store, store@ptr, collection_digest, bed_file_path, output_file_path)
  invisible(store)
})

#' @export
setGeneric('get_seqs_bed_file_to_vec', function(store, collection_digest, bed_file_path) standardGeneric('get_seqs_bed_file_to_vec'))
setMethod('get_seqs_bed_file_to_vec', 'RefgetStore', function(store, collection_digest, bed_file_path) {
  result <- .Call(wrap__get_seqs_bed_file_to_vec_store, store@ptr, collection_digest, bed_file_path)
  if (!is.null(result) && length(result) > 0) {
    lapply(result, convert_to_retrieved_sequence)
  } else {
    list()
  }
})

# =========================================================================
# Sequence Alias Operations
# =========================================================================

#' @export
setGeneric('add_sequence_alias', function(store, namespace, alias, digest) standardGeneric('add_sequence_alias'))
setMethod('add_sequence_alias', 'RefgetStore', function(store, namespace, alias, digest) {
  .Call(wrap__add_sequence_alias_store, store@ptr, namespace, alias, digest)
  invisible(store)
})

#' @export
setGeneric('get_sequence_metadata_by_alias', function(store, namespace, alias) standardGeneric('get_sequence_metadata_by_alias'))
setMethod('get_sequence_metadata_by_alias', 'RefgetStore', function(store, namespace, alias) {
  result <- .Call(wrap__get_sequence_metadata_by_alias_store, store@ptr, namespace, alias)
  convert_to_sequence_metadata(result)
})

#' @export
setGeneric('get_aliases_for_sequence', function(store, digest) standardGeneric('get_aliases_for_sequence'))
setMethod('get_aliases_for_sequence', 'RefgetStore', function(store, digest) {
  .Call(wrap__get_aliases_for_sequence_store, store@ptr, digest)
})

#' @export
setGeneric('list_sequence_alias_namespaces', function(store) standardGeneric('list_sequence_alias_namespaces'))
setMethod('list_sequence_alias_namespaces', 'RefgetStore', function(store) {
  .Call(wrap__list_sequence_alias_namespaces_store, store@ptr)
})

#' @export
setGeneric('list_sequence_aliases', function(store, namespace) standardGeneric('list_sequence_aliases'))
setMethod('list_sequence_aliases', 'RefgetStore', function(store, namespace) {
  .Call(wrap__list_sequence_aliases_store, store@ptr, namespace)
})

#' @export
setGeneric('remove_sequence_alias', function(store, namespace, alias) standardGeneric('remove_sequence_alias'))
setMethod('remove_sequence_alias', 'RefgetStore', function(store, namespace, alias) {
  .Call(wrap__remove_sequence_alias_store, store@ptr, namespace, alias)
})

#' @export
setGeneric('load_sequence_aliases', function(store, namespace, path) standardGeneric('load_sequence_aliases'))
setMethod('load_sequence_aliases', 'RefgetStore', function(store, namespace, path) {
  .Call(wrap__load_sequence_aliases_store, store@ptr, namespace, path)
})

# =========================================================================
# Collection Alias Operations
# =========================================================================

#' @export
setGeneric('add_collection_alias', function(store, namespace, alias, digest) standardGeneric('add_collection_alias'))
setMethod('add_collection_alias', 'RefgetStore', function(store, namespace, alias, digest) {
  .Call(wrap__add_collection_alias_store, store@ptr, namespace, alias, digest)
  invisible(store)
})

#' @export
setGeneric('get_collection_metadata_by_alias', function(store, namespace, alias) standardGeneric('get_collection_metadata_by_alias'))
setMethod('get_collection_metadata_by_alias', 'RefgetStore', function(store, namespace, alias) {
  result <- .Call(wrap__get_collection_metadata_by_alias_store, store@ptr, namespace, alias)
  convert_to_collection_metadata(result)
})

#' @export
setGeneric('get_aliases_for_collection', function(store, digest) standardGeneric('get_aliases_for_collection'))
setMethod('get_aliases_for_collection', 'RefgetStore', function(store, digest) {
  .Call(wrap__get_aliases_for_collection_store, store@ptr, digest)
})

#' @export
setGeneric('list_collection_alias_namespaces', function(store) standardGeneric('list_collection_alias_namespaces'))
setMethod('list_collection_alias_namespaces', 'RefgetStore', function(store) {
  .Call(wrap__list_collection_alias_namespaces_store, store@ptr)
})

#' @export
setGeneric('list_collection_aliases', function(store, namespace) standardGeneric('list_collection_aliases'))
setMethod('list_collection_aliases', 'RefgetStore', function(store, namespace) {
  .Call(wrap__list_collection_aliases_store, store@ptr, namespace)
})

#' @export
setGeneric('remove_collection_alias', function(store, namespace, alias) standardGeneric('remove_collection_alias'))
setMethod('remove_collection_alias', 'RefgetStore', function(store, namespace, alias) {
  .Call(wrap__remove_collection_alias_store, store@ptr, namespace, alias)
})

#' @export
setGeneric('load_collection_aliases', function(store, namespace, path) standardGeneric('load_collection_aliases'))
setMethod('load_collection_aliases', 'RefgetStore', function(store, namespace, path) {
  .Call(wrap__load_collection_aliases_store, store@ptr, namespace, path)
})

# =========================================================================
# FHR Metadata Operations
# =========================================================================

#' @export
setGeneric('set_fhr_metadata', function(store, collection_digest, metadata_json) standardGeneric('set_fhr_metadata'))
setMethod('set_fhr_metadata', 'RefgetStore', function(store, collection_digest, metadata_json) {
  json_str <- if (is.character(metadata_json)) {
    metadata_json
  } else {
    jsonlite::toJSON(metadata_json, auto_unbox = TRUE)
  }
  .Call(wrap__set_fhr_metadata_store, store@ptr, collection_digest, json_str)
  invisible(store)
})

#' @export
setGeneric('get_fhr_metadata', function(store, collection_digest) standardGeneric('get_fhr_metadata'))
setMethod('get_fhr_metadata', 'RefgetStore', function(store, collection_digest) {
  result <- .Call(wrap__get_fhr_metadata_store, store@ptr, collection_digest)
  if (is.null(result) || result == '') {
    NULL
  } else {
    jsonlite::fromJSON(result)
  }
})

#' @export
setGeneric('remove_fhr_metadata', function(store, collection_digest) standardGeneric('remove_fhr_metadata'))
setMethod('remove_fhr_metadata', 'RefgetStore', function(store, collection_digest) {
  .Call(wrap__remove_fhr_metadata_store, store@ptr, collection_digest)
})

#' @export
setGeneric('list_fhr_metadata', function(store) standardGeneric('list_fhr_metadata'))
setMethod('list_fhr_metadata', 'RefgetStore', function(store) {
  .Call(wrap__list_fhr_metadata_store, store@ptr)
})

#' @export
setGeneric('load_fhr_metadata', function(store, collection_digest, path) standardGeneric('load_fhr_metadata'))
setMethod('load_fhr_metadata', 'RefgetStore', function(store, collection_digest, path) {
  .Call(wrap__load_fhr_metadata_store, store@ptr, collection_digest, path)
  invisible(store)
})
