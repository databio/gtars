#' @useDynLib gtars, .registration = TRUE
NULL

setClass(
  'SequenceMetadata', 
  slots = list(
    name = 'character',
    length = 'integer',
    sha512t24u = 'character',
    md5 = 'character',
    alphabet = 'character'
  )
)

setClass(
  'SequenceRecord', 
  slots = list(
    metadata = 'SequenceMetadata',
    data = 'raw'
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
  'SequenceCollection', 
  slots = list(
    sequences = 'list',
    digest = 'character',
    lvl1 = 'SeqColDigestLvl1',
    file_path = 'character',
    has_data = 'logical'
  )
)

setMethod('show', 'SequenceMetadata', function(object) {
  cat(sprintf('SequenceMetadata for sequence %s\n', object@name))
  cat(sprintf('  name: %s\n', object@name))
  cat(sprintf('  length: %d\n', object@length))
  cat(sprintf('  sha512t24u: %s\n', object@sha512t24u))
  cat(sprintf('  md5: %s\n', object@md5))
  cat(sprintf('  alphabet: %s\n', object@alphabet))
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

setMethod('show', 'SequenceCollection', function(object) {
  cat(sprintf('SequenceCollection with %d sequences, digest = %s\n', length(object@sequences), object@digest))
})

setMethod('length', 'SequenceCollection', function(x) {
  length(x@sequences)
})

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
#' r_md5_digest('ATCG')
#' }
#' 
#' @export
digest_fasta <- function(fasta) {
  result <- .Call(wrap__digest_fasta_raw, fasta)
  
  lvl1 <- new(
    'SeqColDigestLvl1', 
    sequences_digest = result$lvl1$sequences_digest, 
    names_digest = result$lvl1$names_digest,
    lengths_digest = result$lvl1$lengths_digest
  )
  
  sequence_records <- lapply(result$sequences, function(seq) {
    metadata <- new(
      'SequenceMetadata', 
      name = seq$metadata$name,
      length = seq$metadata$length,
      sha512t24u = seq$metadata$sha512t24u,
      md5 = seq$metadata$md5,
      alphabet = seq$metadata$alphabet
    )
    
    new(
      'SequenceRecord',
      metadata = metadata,
      data = seq$data
    )
  })
  
  invisible(new(
    'SequenceCollection',
    sequences = sequence_records,
    digest = result$digest,
    lvl1 = lvl1,
    file_path = result$file_path,
    has_data = result$has_data
  ))
}


setClass(
  'GlobalRefgetStore',
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


convert_to_sequence_record <- function(raw_result) {
  if (is.null(raw_result)) return(NULL)
  
  metadata <- new(
    'SequenceMetadata',
    name = raw_result$metadata$name,
    length = as.integer(raw_result$metadata$length),
    sha512t24u = raw_result$metadata$sha512t24u,
    md5 = raw_result$metadata$md5,
    alphabet = raw_result$metadata$alphabet
  )
  
  new(
    'SequenceRecord',
    metadata = metadata,
    data = as.raw(raw_result$data)
  )
}

setGeneric('import_fasta', function(store, file_path) standardGeneric('import_fasta'))
setMethod('import_fasta', 'GlobalRefgetStore', function(store, file_path) {
  .Call(wrap__import_fasta_store, store@ptr, file_path)
  invisible(store)
})

setGeneric('get_sequence_by_id', function(store, digest) standardGeneric('get_sequence_by_id'))
setMethod('get_sequence_by_id', 'GlobalRefgetStore', function(store, digest) {
  result <- .Call(wrap__get_sequence_by_id_store, store@ptr, digest)
  convert_to_sequence_record(result)
})

setGeneric('get_sequence_by_collection_and_name', function(store, collection_digest, sequence_name) standardGeneric('get_sequence_by_collection_and_name'))
setMethod('get_sequence_by_collection_and_name', 'GlobalRefgetStore', function(store, collection_digest, sequence_name) {
  result <- .Call(wrap__get_sequence_by_collection_and_name_store, store@ptr, collection_digest, sequence_name)
  convert_to_sequence_record(result)
})

setGeneric('get_substring', function(store, seq_digest, start, end) standardGeneric('get_substring'))
setMethod('get_substring', 'GlobalRefgetStore', function(store, seq_digest, start, end) {
  result <- .Call(wrap__get_substring_store, store@ptr, seq_digest, as.integer(start), as.integer(end))
  if (is.null(result)) NULL else as.character(result)
})

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

# Write store to directory
setGeneric('write_store_to_directory', function(store, root_path, seqdata_path_template) standardGeneric('write_store_to_directory'))
setMethod('write_store_to_directory', 'GlobalRefgetStore', function(store, root_path, seqdata_path_template) {
  .Call(wrap__write_store_to_directory_store, 
        store@ptr, root_path, seqdata_path_template)
  invisible(store)  # Return for chaining
})

# Get sequences from BED file (writes to file)
setGeneric('get_seqs_bed_file', function(store, collection_digest, bed_file_path, output_file_path) standardGeneric('get_seqs_bed_file'))
setMethod('get_seqs_bed_file', 'GlobalRefgetStore', function(store, collection_digest, bed_file_path, output_file_path) {
  .Call(wrap__get_seqs_bed_file_store, 
        store@ptr, collection_digest, bed_file_path, output_file_path)
  invisible(store)  # Return for chaining
})

# Get sequences from BED file (returns vector)
setGeneric('get_seqs_bed_file_to_vec', function(store, collection_digest, bed_file_path) standardGeneric('get_seqs_bed_file_to_vec'))
setMethod('get_seqs_bed_file_to_vec', 'GlobalRefgetStore', function(store, collection_digest, bed_file_path) {
  result <- .Call(wrap__get_seqs_bed_file_to_vec_store, 
                  store@ptr, collection_digest, bed_file_path)
  if (!is.null(result) && length(result) > 0) {
    # Convert each element to RetrievedSequence S4 object
    lapply(result, convert_to_retrieved_sequence)
  } else {
    list()
  }
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

setMethod('show', 'GlobalRefgetStore', function(object) {
  cat('GlobalRefgetStore object\n')
  cat('  Use methods: import_fasta(), get_sequence_by_id(), get_substring(), etc.\n')
})


#' @title global refget store
#' 
#' @description create a pointer to a new global refget store rust object
#' 
#' @param mode either 'raw' or 'encoded'
#' 
#' @return a GlobalRefgetStore object
#' 
#' @examples
#' \dontrun{
#' r_md5_digest('ATCG')
#' }
#' 
#' @export
global_refget_store <- function(mode) {
  result <- .Call(wrap__global_refget_store_raw, mode)
  
  invisible(new('GlobalRefgetStore', ptr = result))
}

#' @title load refget store from directory
#' 
#' @description load refget store from directory
#' 
#' @param root_path path to refget store
#' 
#' @return a GlobalRefgetStore object
#' 
#' @examples
#' \dontrun{
#' r_md5_digest('ATCG')
#' }
#' 
#' @export
load_from_directory <- function(root_path) {
  result <- .Call(wrap__load_from_directory_store, root_path)
  
  invisible(new('GlobalRefgetStore', ptr = result))
}



