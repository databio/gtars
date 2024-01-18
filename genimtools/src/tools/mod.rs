use std::fs::File;
use std::io::{BufReader, Write, BufRead};

use walkdir::WalkDir;

pub mod cli;

pub mod consts {
    pub const DATA_DIR_STAT_CMD: &str = "dir-stat";
    pub const DEFAULT_OUTPUT: &str = "output.tsv";
}

///
/// Collect statistics of a given data directory. This
/// will look through a potential data directory and
/// report some statistics about the data contained
/// within.
/// 
/// It creates a report with three columns:
/// - `file`: The name of the file.
/// - `size`: The size of the file in bytes.
/// - `lines`: The number of lines in the file.
/// 
/// # Arguments
/// - `path`: The path to the data directory.
/// - `out`: The path to the output file.
///
pub fn data_dir_stat(path: &str, out: &str) {
    // create the writer
    let mut writer = File::create(out).unwrap();

    // iterate over the files in the directory
    for entry in WalkDir::new(path).into_iter().filter_map(|e| e.ok()) {
        if entry.file_type().is_file() {

            // get the file name, size, and number of lines
            let fname = entry.file_name().to_str().unwrap();
            let size = entry.metadata().unwrap().len();
            let num_lines = BufReader::new(File::open(entry.path()).unwrap()).lines().count();

            // create the line
            let line = format!("{}\t{}\t{}\n", fname, size, num_lines);

            // write the line
            writer.write_all(line.as_bytes()).unwrap();
        }
    }
}