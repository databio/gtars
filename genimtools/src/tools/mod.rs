use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use walkdir::WalkDir;

use crate::common::models::RegionSet;
use crate::io::write_tokens_to_gtok;
use crate::tokenizers::{Tokenizer, TreeTokenizer};

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

pub fn pre_tokenize_data(
    path_to_data: &Path,
    outdir: &str,
    tokenizer: &TreeTokenizer,
) -> Result<()> {
    if path_to_data.is_file() {
        pre_tokenize_file(path_to_data, outdir, tokenizer)?;
    } else {
        let bed_file_path = Path::new(path_to_data).join("**/*.bed");
        let zipped_bed_file_path = Path::new(path_to_data).join("**/*.bed.gz");

        let num_beds = glob::glob(bed_file_path.to_str().unwrap())
            .unwrap()
            .filter_map(|e| e.ok())
            .count();

        let num_zipped = glob::glob(zipped_bed_file_path.to_str().unwrap())
            .unwrap()
            .filter_map(|e| e.ok())
            .count();

        let bed_file_matches = glob::glob(bed_file_path.to_str().unwrap())
            .unwrap()
            .filter_map(|e| e.ok());

        let zipped_bed_file_matches = glob::glob(zipped_bed_file_path.to_str().unwrap())
            .unwrap()
            .filter_map(|e| e.ok());

        let num_files = num_beds + num_zipped;

        let pb = ProgressBar::new(num_files as u64);
        pb.set_style(
            ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
            )
            .unwrap()
            .progress_chars("##-"),
        );

        for entry in bed_file_matches {
            let path = entry;
            pre_tokenize_file(&path, outdir, tokenizer)?;
            pb.inc(1);
        }

        for entry in zipped_bed_file_matches {
            let path = entry;
            pre_tokenize_file(&path, outdir, tokenizer)?;
            pb.inc(1);
        }
    }

    Ok(())
}

fn pre_tokenize_file(
    path_to_bedfile: &Path,
    outdir: &str,
    tokenizer: &TreeTokenizer,
) -> Result<()> {
    // make sure the file ends in .bed or .bed.gz
    let ext = path_to_bedfile.extension().unwrap();
    if ext != OsStr::new("bed") && ext != OsStr::new("gz") {
        println!("Skipping file: {}", path_to_bedfile.display());
        println!(
            "File must end in .bed or .bed.gz, ends with .{}",
            ext.to_str().unwrap()
        );
        return Ok(());
    }

    let out_file = Path::new(outdir).join(format!(
        "{}.{}",
        path_to_bedfile.file_stem().unwrap().to_str().unwrap(),
        crate::common::consts::GTOK_EXT
    ));

    let out_file = out_file.to_str().unwrap();

    let regions = RegionSet::try_from(path_to_bedfile).expect("Failed to read bed file");

    let tokens = tokenizer.tokenize_region_set(&regions);

    write_tokens_to_gtok(out_file, &tokens.ids)?;

    Ok(())
}
