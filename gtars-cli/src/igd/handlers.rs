use std::fs;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

use clap::ArgMatches;

use gtars_core::consts::{BED_FILE_EXTENSION, GZ_FILE_EXTENSION};
use gtars_core::models::RegionSet;
use gtars_igd::igd::Igd;

/// Resolve a filelist argument into a list of BED file paths.
///
/// Supports three input modes (matching legacy `create_igd_f` behavior):
/// - A `.txt` file: each line is a path to a BED file
/// - `"-"` or `"stdin"`: read paths from stdin, one per line
/// - A directory: enumerate BED/gz files in it
fn resolve_bed_paths(filelist: &str) -> anyhow::Result<Vec<PathBuf>> {
    if filelist.ends_with(".txt") {
        let file = fs::File::open(filelist)?;
        let reader = BufReader::new(file);
        Ok(reader
            .lines()
            .filter_map(|l| l.ok())
            .map(|l| PathBuf::from(l.trim()))
            .filter(|p| !p.as_os_str().is_empty())
            .collect())
    } else if filelist == "-" || filelist == "stdin" {
        let stdin = io::stdin();
        let reader = BufReader::new(stdin.lock());
        Ok(reader
            .lines()
            .filter_map(|l| l.ok())
            .map(|l| PathBuf::from(l.trim()))
            .filter(|p| !p.as_os_str().is_empty())
            .collect())
    } else {
        let mut paths: Vec<PathBuf> = Vec::new();
        for entry in fs::read_dir(filelist)? {
            let p = entry?.path();
            if let Some(ext) = p.extension().and_then(|e| e.to_str()) {
                if (ext == BED_FILE_EXTENSION.trim_start_matches('.')
                    || ext == GZ_FILE_EXTENSION.trim_start_matches('.'))
                    && p.is_file()
                {
                    paths.push(p);
                }
            }
        }
        paths.sort();
        Ok(paths)
    }
}

pub fn igd_get_create_matches(matches: &ArgMatches) {
    let output_path = matches
        .get_one::<String>("output")
        .expect("Output path is required");

    let filelist = matches
        .get_one::<String>("filelist")
        .expect("File list path is required");

    let db_output_name = matches
        .get_one::<String>("dbname")
        .expect("Database name is required");

    let paths = resolve_bed_paths(filelist).expect("Failed to resolve BED file paths");
    let igd = Igd::from_bed_files(paths).expect("Failed to create IGD from BED files");

    let save_path = Path::new(output_path).join(format!("{}.igd", db_output_name));
    igd.save(&save_path).expect("Failed to save IGD database");
}

/// Searches IGD database
pub fn igd_get_search_matches(matches: &ArgMatches) {
    let database_path = matches
        .get_one::<String>("database")
        .expect("Database path is required");

    let query = matches
        .get_one::<String>("query")
        .expect("Query bed file path is required");

    let igd = Igd::from_igd_file(Path::new(database_path)).expect("Failed to load IGD database");
    let query_rs =
        RegionSet::try_from(Path::new(query).to_path_buf()).expect("Failed to parse query BED file");
    let hits = igd.count_set_overlaps(&query_rs, 1);

    // Format output matching legacy TSV: "index\tnRegions\tnHits\tfilename"
    println!("index\t number of regions\t number of hits\t File_name");
    for (i, fi) in igd.file_info.iter().enumerate() {
        if hits[i] > 0 {
            println!("{}\t{}\t{}\t{}", i, fi.num_regions, hits[i], fi.filename);
        }
    }
    let total: u64 = hits.iter().sum();
    println!("Total: {}", total);
}
