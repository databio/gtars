// bindings/r/src/rust/src/lib.rs
use extendr_api::prelude::*;
pub mod io;

/// @export
#[extendr]
fn igd_search(database_path: &str, query_path: &str) -> extendr_api::Result<Robj> {
    use std::collections::HashMap;
    use gtars::igd::search::{
        get_igd_info, 
        get_file_info_tsv
    };

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
    
    // Process the search
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

#[extendr]
fn __init__() {}

extendr_module! {
    mod gtars;
    fn __init__;
    fn igd_search;
}
