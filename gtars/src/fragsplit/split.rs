use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;
use std::{collections::HashMap, fs};

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use indicatif::{ProgressBar, ProgressStyle};

use crate::common::utils::get_dynamic_reader;
use crate::fragsplit::map::BarcodeToClusterMap;
use crate::fragsplit::utils::remove_all_extensions;

use super::map::ClusterLookup;

///
/// Psuedobulks fragment files accoring to a specified mapping.
///
/// Given a folder of fragment files, this function will read them in
/// and split off the reads into new files according to a mapping
/// specified by the user:
///
/// | barcode1 | A |
/// |----------|---|
/// | barcode2 | B |
/// | barcode3 | A |
/// | barcode4 | A |
/// | barcode5 | B |
///
/// # Arguments:
/// - files: path to fragment files
/// - mapping: path to mapping (a tsv file)
/// - output: path to the output folder where new files should go
///
pub fn pseudobulk_fragment_files(
    files: &Path,
    mapping: &BarcodeToClusterMap,
    output: &Path,
) -> Result<()> {
    let files = fs::read_dir(files).with_context(|| {
        format!(
            "There was an error reading the specifed fragment file directory: {:?}",
            files
        )
    })?;

    // convert files to Path -- consume iterator
    let files: Vec<Result<PathBuf>> = files
        .map(|f| {
            let f = f?;
            Ok(f.path())
        })
        .collect();

    // create actual output directory
    fs::create_dir_all(output).with_context(|| {
        format!(
            "There was an error creating the output directory: {:?}",
            output
        )
    })?;

    let mut handle_map: HashMap<String, BufWriter<GzEncoder<File>>> = HashMap::new();
    for cluster_id in mapping.get_cluster_labels() {
        let file_name = format!("cluster_{cluster_id}.bed.gz");
        let file_path = output.join(file_name);
        let file_path = Path::new(&file_path);
        let file = File::create(file_path)?;

        let buf_writer = BufWriter::new(GzEncoder::new(file, Compression::default()));

        handle_map.insert(cluster_id, buf_writer);
    }

    let total_files = files.len();

    let pb = ProgressBar::new(total_files as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} files ({eta})")?
            .progress_chars("##-"),
    );

    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed}] {msg} ({per_sec})")
            .unwrap()
            .tick_strings(&["-", "\\", "|", "/"]),
    );

    spinner.set_message("Processing fragment files...");

    let _start_time = Instant::now();
    let mut processed_reads: u64 = 0;

    for file in files {
        let file = file?;
        let reader = get_dynamic_reader(file.as_path())?;

        // strip out any *.*.gz
        let file_path = file.as_path();
        let file_stem = remove_all_extensions(file_path);

        for (index, line) in reader.lines().enumerate() {
            let line = line?;
            let mut parts = line.split_whitespace();

            let chr = parts.next();
            let start = parts.next();
            let end = parts.next();
            let barcode = parts.next();
            let read_support = parts.next();

            if let (Some(chr), Some(start), Some(end), Some(barcode), Some(read_support)) =
                (chr, start, end, barcode, read_support)
            {
                // merge file stem + barcode to get lookup values
                let lookup_value = format!("{}+{}", file_stem, barcode);
                let cluster = mapping.get_cluster_from_barcode(&lookup_value);
                if let Some(cluster) = cluster {
                    let cluster_file = handle_map.get_mut(&cluster).unwrap();
                    cluster_file.write_all(
                        format!("{chr}\t{start}\t{end}\t{barcode}\t{read_support}\n").as_bytes(),
                    )?;
                }
                // pass on else, since the barcode was most likely a cell tossed in QC/processing
            } else {
                anyhow::bail!(format!(
                    "Failed to parse fragments file at line {index}: {}",
                    line
                ))
            }

            // let elapsed = start_time.elapsed().as_secs();
            processed_reads += 1;
            if processed_reads % 10_000 == 0 {
                spinner.set_message(format!("Processed {} reads", processed_reads));
            }

            spinner.inc(1);
        }

        pb.inc(1);
    }

    spinner.finish_with_message("Done!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn barcode_cluster_map_file() -> &'static str {
        "tests/data/barcode_cluster_map.tsv"
    }

    #[fixture]
    fn path_to_fragment_files() -> &'static str {
        "tests/data/fragments/fragsplit"
    }

    #[fixture]
    fn path_to_output() -> &'static str {
        "tests/data/out"
    }

    #[fixture]
    fn filtered_out_barcode() -> &'static str {
        "AAACGCAAGCAAAGGATCGGCT"
    }

    #[rstest]
    fn test_fragment_file_splitter(
        barcode_cluster_map_file: &str,
        path_to_fragment_files: &str,
        path_to_output: &str,
    ) {
        let barcode_cluster_map_file = Path::new(barcode_cluster_map_file);
        let mapping = BarcodeToClusterMap::from_file(barcode_cluster_map_file).unwrap();

        let path_to_fragment_files = Path::new(path_to_fragment_files);
        let path_to_output = Path::new(path_to_output);

        let res = pseudobulk_fragment_files(path_to_fragment_files, &mapping, path_to_output);

        assert_eq!(res.is_ok(), true);
    }
}
