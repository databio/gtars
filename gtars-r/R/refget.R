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

setClass(
  'RetrievedSequence',
  slots = list(
    sequence = 'character',
    chrom_name = 'character',
    start = 'integer', 
    end = 'integer'
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
  cat(sprintf('SequenceRecord for %s:\n', object@metadata@name))
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


