use std::fs;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

use extendr_api::prelude::*;
use gtars_core::consts::{BED_FILE_EXTENSION, GZ_FILE_EXTENSION};
use gtars_core::models::RegionSet;
use gtars_igd::igd::Igd;

/// Resolve a filelist argument into a list of BED file paths.
///
/// Supports three input modes (matching legacy `create_igd_f` behavior):
/// - A `.txt` file: each line is a path to a BED file
/// - `"-"` or `"stdin"`: read paths from stdin, one per line
/// - A directory: enumerate BED/gz files in it
fn resolve_bed_paths(filelist: &str) -> Result<Vec<PathBuf>, String> {
    if filelist.ends_with(".txt") {
        let file = fs::File::open(filelist).map_err(|e| e.to_string())?;
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
        for entry in fs::read_dir(filelist).map_err(|e| e.to_string())? {
            let p = entry.map_err(|e| e.to_string())?.path();
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

/// Create an IGD database from a directory of bed files
/// @param output_path String path where the IGD database will be saved
/// @param filelist String path to either a text file containing paths to bed files, or a directory containing bed files
/// @param db_name String name for the database (will be used in output filenames)
#[extendr]
fn r_igd_create(
    output_path: &str,
    filelist: &str,
    db_name: &str,
) -> std::result::Result<(), extendr_api::Error> {
    let paths = resolve_bed_paths(filelist).map_err(Error::from)?;
    let igd = Igd::from_bed_files(paths).map_err(|e| Error::from(e.to_string()))?;

    let save_path = Path::new(output_path).join(format!("{}.igd", db_name));
    igd.save(&save_path).map_err(|e| Error::from(e.to_string()))?;

    Ok(())
}

/// Search igd with a bed file
/// @param database_path A string representing the path to the database igd file.
/// @param query_path A string representing the path to the query bed file.
/// @return A named list with columns: filename, numRegions, hits (convertible to data.frame)
#[extendr]
pub fn r_igd_search(
    database_path: &str,
    query_path: &str,
) -> std::result::Result<List, extendr_api::Error> {
    let igd =
        Igd::from_igd_file(Path::new(database_path)).map_err(|e| Error::from(e.to_string()))?;
    let query_rs =
        RegionSet::try_from(Path::new(query_path)).map_err(|e| Error::from(e.to_string()))?;
    let hits = igd.count_set_overlaps(&query_rs, 1);

    // Build column vectors (only files with hits > 0)
    let mut filenames: Vec<String> = Vec::new();
    let mut num_regions: Vec<i32> = Vec::new();
    let mut hit_counts: Vec<i32> = Vec::new();

    for (i, fi) in igd.file_info.iter().enumerate() {
        if hits[i] > 0 {
            filenames.push(fi.filename.clone());
            num_regions.push(fi.num_regions as i32);
            hit_counts.push(hits[i] as i32);
        }
    }

    Ok(list!(
        filename = filenames,
        numRegions = num_regions,
        hits = hit_counts
    ))
}

extendr_module! {
    mod igd;
    fn r_igd_create;
    fn r_igd_search;
}
