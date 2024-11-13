#' @useDynLib gtars, .registration = TRUE
#' @importFrom methods new
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
#' 
#' # Create database with custom name
#' igd_create("path/to/output", "path/to/bed/files", "my_database")
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
  .Call(wrap__r_igd_create, output_path, filelist, db_name)
  
  invisible(NULL)
}
