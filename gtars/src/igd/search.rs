use clap::ArgMatches;
use std::path::Path;
use crate::common::consts::{BED_FILE_EXTENSION, IGD_FILE_EXTENSION};

/// Searches IGD database
pub fn igd_get_search_matches(matches: &ArgMatches) {

    let database_path = matches
        .get_one::<String>("database")
        .expect("Database path is required");

    let query = matches
        .get_one::<String>("query")
        .expect("Query bed file path is required");


    igd_search(database_path, query).expect("Error:");
}

pub fn igd_search(database_path: &String, query_file_path: &String) -> Result<(), String> {

    // First check that BOTH the igd database and the query are the proper file types
    // else raise error

    let mode = 1;

    match check_file_extension(database_path, IGD_FILE_EXTENSION) {
        Ok(_) => (),
        Err(e) => return Err(e),
    }

    match check_file_extension(query_file_path, BED_FILE_EXTENSION) {
        Ok(_) => (),
        Err(e) => {;
            return Err(e);
        
        }
        ,
    }

    println!("\n {} \n {}", database_path,query_file_path);


    //Get file info from the associated TSV



    // If query "-q" used set to mode 1

    match mode {

        1 => {



        },
        _ => {
            println!("Invalid mode selected, exiting");
            return Ok(());
        },


    }


    println!("FINISHED");
    
    Ok(())

}

fn check_file_extension(path: &str, expected_extension: &str) -> Result<(), String> {
    let path = Path::new(path);
    let actual_extension = path
        .extension()
        .and_then(|ext| ext.to_str())
        .ok_or_else(|| format!("Invalid file path: {}", path.display()))?;

    if actual_extension != expected_extension {
        return Err(format!("Incorrect file extension. Expected: {}, got: {}", expected_extension, actual_extension));
    }

    Ok(())
}
