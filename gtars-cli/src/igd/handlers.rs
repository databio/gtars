use clap::ArgMatches;

use igd::create::create_igd_f;
use igd::search::igd_search;

pub fn igd_get_create_matches(matches: &ArgMatches) {
    //println!("HELLO FROM IGD CREATE SUBMODULE!");

    let output_path = matches
        .get_one::<String>("output")
        .expect("Output path is required");

    let filelist = matches
        .get_one::<String>("filelist")
        .expect("File list path is required");

    let db_output_name = matches
        .get_one::<String>("dbname")
        .expect("File list path is required");

    let _igd = create_igd_f(output_path, filelist, db_output_name);
}

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