#' @useDynLib gtars, .registration = TRUE
NULL

#' @title Create IGD Database
#' 
#' @description Creates an IGD (Indexed Genomic Data) database from a collection of BED files.
#' 
#' @param output_path Character string specifying the directory where the IGD database will be saved
#' @param filelist Character string specifying either:
#'   - Path to a text file containing paths to BED files (one per line)
#'   - Path to a directory containing BED files
#'   - "-" or "stdin" to read paths from standard input
#' @param db_name Character string specifying the name for the database (will be used in output filenames).
#'   Defaults to "igd_database"
#' 
#' @return NULL invisibly on success
#' 
#' @examples
#' \dontrun{
#' # Create database with default name
#' igd_create("path/to/output", "path/to/bed/files")
#' }
#' 
#' @export
igd_create <- function(output_path, filelist, db_name = "igd_database") {
  # Input validation
  if (!is.character(output_path) || length(output_path) != 1) {
    stop("output_path must be a single character string")
  }
  if (!is.character(filelist) || length(filelist) != 1) {
    stop("filelist must be a single character string")
  }

  # Call Rust function
  .Call(wrap__igd_create, output_path, filelist, db_name)

  invisible(NULL)
}


#' @title Search IGD Database
#' 
#' @description Searches an IGD database for region overlaps with an input BED file
#' 
#' @param database_path path to .igd database
#' @param query_path path to .bed file
#' 
#' @return dataframe of overlap hits
#' 
#' @examples
#' \dontrun{
#' # Search database with default name
#' igd_search("path/to/database", "path/to/query/file")
#' }
#' 
#' @export
igd_search <- function(database_path, query_path) {
  
  # Input validation
  if (!is.character(database_path) || length(database_path) != 1) {
    stop("database_path must be a single character string")
  }
  if (!is.character(query_path) || length(query_path) != 1) {
    stop("query_path must be a single character string")
  }
  
  # Call Rust function
  chr_vector <- .Call(wrap__igd_search, database_path, query_path)
  
  split_result <- strsplit(chr_vector, split = '\t')
  df <- data.frame(matrix(unlist(split_result[-1]), nrow = length(chr_vector)-1, byrow = TRUE))
  colnames(df) <- split_result[[1]]
  
  invisible(df)
}
