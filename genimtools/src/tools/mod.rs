use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use indicatif::{ProgressBar, ProgressStyle};
use walkdir::WalkDir;

pub mod cli;

pub mod consts {
    pub const TOOLS_CMD: &str = "tools";
    pub const DATA_DIR_STAT_CMD: &str = "dir-stat";
    pub const DEFAULT_DATA_DIR_STAT_OUTPUT: &str = "output.tsv";
    pub const DEFAULT_PRETOKENIZE_OUT: &str = "gtoks";
    pub const PRE_TOKENIZATION_CMD: &str = "pretokenize";
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

    // get num files in dir
    println!("Counting files in directory...");
    let num_files = WalkDir::new(path)
        .into_iter()
        .filter_map(|e| e.ok())
        .count();

    // create the progress bar
    let pb = ProgressBar::new(num_files as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-"),
    );

    // iterate over the files in the directory
    println!("Iterating over files in directory...");
    for entry in WalkDir::new(path).into_iter().filter_map(|e| e.ok()) {
        if entry.file_type().is_file() {
            // get the file name, size, and number of lines
            let fname = entry.file_name().to_str().unwrap();
            let size = entry.metadata().unwrap().len();
            let num_lines = BufReader::new(File::open(entry.path()).unwrap())
                .lines()
                .count();

            // create the line
            let line = format!("{}\t{}\t{}\n", fname, size, num_lines);

            // write the line
            writer.write_all(line.as_bytes()).unwrap();
        }

        pb.inc(1);
    }
}
