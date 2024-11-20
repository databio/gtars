use extendr_api::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;
use gtars::igd::search::{get_igd_info, get_file_info_tsv};
use gtars::igd::create::create_igd_f;



/// RUST WRAPPERS SHOULD BE MINIMAL. HANDLE DATA STRUCTURES IN IGD

/// Search igd with a bed file
/// @param database_path A string representing the path to the database igd file.
/// @param query_path A string representing the path to the query bed file.
/// @export
#[extendr(r_name = "igd_search")]
pub fn r_igd_search(database_path: &str, query_path: &str) -> extendr_api::Result<Robj> {
    
    // Create data structures
    let mut hash_table: HashMap<String, i32> = HashMap::new();
    
    // Get IGD info
    let mut igd = get_igd_info(&database_path.to_string(), &mut hash_table)
        .map_err(|e| Error::Other(format!("Failed to open IGD: {}", e)))?;
    
    // Get TSV info
    let tsv_path = {
        let path = std::path::Path::new(database_path);
        let stem = path.file_stem()
            .ok_or_else(|| Error::Other("Invalid database path".into()))?;
        let mut tsv_path = path.with_file_name(stem);
        tsv_path.set_extension("tsv");
        tsv_path
    };
    
    get_file_info_tsv(tsv_path, &mut igd)
        .map_err(|e| Error::Other(format!("Failed to get file info: {}", e)))?;

    // Initialize hits vector
    let mut hits = vec![0; igd.nFiles as usize];
    
    // Process the search THIS IS MOST IMPORTANT IN WRAPPER
    gtars::igd::search::getOverlaps(
        &igd,
        &database_path.to_string(),
        &query_path.to_string(),
        &mut hits,
        &mut hash_table,
    );
    
    // Prepare the data
    let mut file_names = vec![];
    let mut region_counts = vec![];
    let mut hit_counts = vec![];
    
    for (i, hit) in hits.iter().enumerate() {
        if *hit > 0 {
            file_names.push(&igd.file_info[i].fileName);
            region_counts.push(igd.file_info[i].nr);
            hit_counts.push(*hit);
        }
    }
    
    // Create R list using the named list function
    let result = call!("list", 
        file_name = file_names,
        n_regions = region_counts,
        n_hits = hit_counts
    )?;
    
    Ok(result)
}

/// Create an IGD database from a directory of bed files
/// @param output_path String path where the IGD database will be saved
/// @param filelist String path to either a text file containing paths to bed files, or a directory containing bed files
/// @param db_name String name for the database (will be used in output filenames)
#[extendr]
fn r_igd_create(output_path: &str, filelist: &str, db_name: &str) -> std::result::Result<(), extendr_api::Error> {
    // Validate inputs
    if output_path.is_empty() {
        return Err(Error::from("output_path cannot be empty"));
    }
    if filelist.is_empty() {
        return Err(Error::from("filelist cannot be empty"));
    }
    if db_name.is_empty() {
        return Err("db_name cannot be empty".into());
    }

    // Ensure output path exists
    let output_pathbuf = PathBuf::from(output_path);
    if !output_pathbuf.exists() {
        if let Err(e) = std::fs::create_dir_all(&output_pathbuf) {
            return Err(Error::from(format!("Failed to create output directory: {}", e)));
        }
    }

    // Call the underlying create function
    create_igd_f(&output_path.to_string(), &filelist.to_string(), &db_name.to_string());
    
    Ok(())
}

extendr_module! {
    mod igd;
    fn r_igd_search;
    fn r_igd_create;
}
