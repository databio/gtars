use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;
use std::{collections::HashMap, fs};

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::common::utils::get_dynamic_reader;
use crate::fragsplit::map::BarcodeToClusterMap;

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

    // create actual output directory
    fs::create_dir_all(output).with_context(|| {
        format!(
            "There was an error creating the output directory: {:?}",
            output
        )
    })?;

    let mut handle_map: HashMap<char, BufWriter<GzEncoder<File>>> = HashMap::new();
    for cluster_id in mapping.get_cluster_labels() {
        let file_name = format!("cluster_{cluster_id}.bed.gz");
        let file_path = output.join(file_name);
        let file_path = Path::new(&file_path);
        let file = File::create(file_path)?;

        let buf_writer = BufWriter::new(GzEncoder::new(file, Compression::default()));

        handle_map.insert(cluster_id, buf_writer);
    }

    for file in files {
        let file = file?;
        let reader = get_dynamic_reader(&file.path())?;

        // strip out any *.*.gz
        let file_path = file.path();
        let file_stem = file_path.file_stem().unwrap();
        let file_stem = file_stem.to_string_lossy();

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
                let cluster_id = mapping.get_cluster_from_barcode(&lookup_value);
                if let Some(cluster_id) = cluster_id {
                    let cluster_file = handle_map.get_mut(&cluster_id).unwrap();
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
        }
    }

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
        "tests/data/fragments"
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
