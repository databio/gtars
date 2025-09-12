#' @useDynLib gtars, .registration = TRUE
NULL

#' @title create sha512t24u_digest
#' 
#' @description create a sha512t24u digest from a readable string
#' 
#' @param readable a readable string
#' 
#' @return a sha512t24u digest on success
#' 
#' @examples
#' \dontrun{
#' r_sha512t24u_digest('ATCG')
#' }
#' 
#' @export
r_sha512t24u_digest <- function(readable) {
  result <- .Call(wrap__sha512t24u_digest, readable)
  invisible(result)
}


#' @title create md5_digest
#' 
#' @description create an md5 digest from a readable string
#' 
#' @param readable a readable string
#' 
#' @return an md5 digest on success
#' 
#' @examples
#' \dontrun{
#' r_md5_digest('ATCG')
#' }
#' 
#' @export
r_md5_digest <- function(readable) {
  result <- .Call(wrap__md5_digest, readable)
  invisible(result)
}


setClass(
  'SequenceCollection', 
  slots =
    list(sequences = 'list',
         digest = 'character',
         lvl1 = 'character',
         file_path = 'character',
         has_data = 'boolean'
    )
)

setClass(
  'SequenceMetadata', 
  slots =
    list(name = 'character',
         length = 'double',
         sha512t24u = 'character',
         md5 = 'character',
         alphabet = 'character'
    )
)

setClass(
  'SeqColDigestLvl1', 
  slots =
    list(sequences_digest = 'character',
         names_digest = 'character',
         lengths_digest = 'character',
    )
)

setClass(
  'SequenceRecord', 
  slots =
    list(metadata = 'SequenceMetadata',
         data = 'list'
    )
)

# setClass(
#   'GlobalRefgetStore', 
#   slots =
#     list(sequences = 'list',
#          digest = 'character',
#          lvl1 = 'character',
#          file_path = 'character',
#          has_data = 'boolean'
#     )
# )
