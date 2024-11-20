use extendr_api::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;
use gtars::igd::search::{igd_search, get_igd_info, get_file_info_tsv};
use gtars::igd::create::create_igd_f;



/// RUST WRAPPERS SHOULD BE MINIMAL. HANDLE DATA STRUCTURES IN IGD

/// Search igd with a bed file
/// @param database_path A string representing the path to the database igd file.
/// @param query_path A string representing the path to the query bed file.
/// @export
#[extendr]
pub fn r_igd_search(database_path: &str, query_path: &str) -> std::result::Result<Vec<String>, extendr_api::Error> {

    let dbpath = String::from(database_path);
    let qpath = String::from(query_path);

    let result = igd_search(&dbpath, &qpath);

    match result {
        Ok(vector_strings) => return Ok(vector_strings),
        Err(e) => return Err(Error::from(e)),
    }

}
//
/// Create an IGD database from a directory of bed files
/// @param output_path String path where the IGD database will be saved
/// @param filelist String path to either a text file containing paths to bed files, or a directory containing bed files
/// @param db_name String name for the database (will be used in output filenames)
/// @export
#[extendr]
fn r_igd_create(output_path: &str, filelist: &str, db_name: &str) -> std::result::Result<(), extendr_api::Error> {

    // Call the underlying create function
    create_igd_f(&output_path.to_string(), &filelist.to_string(), &db_name.to_string());

    Ok(())
}

extendr_module! {
    mod igd;
    fn r_igd_search;
    fn r_igd_create;
}
