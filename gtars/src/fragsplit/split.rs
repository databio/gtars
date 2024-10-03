use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::{collections::HashMap, fs};
use std::path::Path;


use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::fragsplit::map::BarcodeToClusterMap;
use crate::common::utils::get_dynamic_reader;

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
/// - mapping: path to mapping (a csv file)
/// - output: path to the output folder where new files should go
/// 
pub fn pseudobulk_fragment_files(files: &Path, mapping: &BarcodeToClusterMap, output: &Path) -> Result<()> {
    let files = fs::read_dir(files)
        .with_context(|| {
            format!("There was an error reading the specifed fragment file directory: {:?}", files)
        })?;
    
    fs::create_dir_all(output)
        .with_context(|| {
            format!("There was an error creating the output directory: {:?}", output)
        })?;

    let mut handle_map: HashMap<char, BufWriter<GzEncoder<File>>> = HashMap::new();
    for cluster_id in mapping.get_cluster_labels() {

        let file_name = format!("cluster_{cluster_id}.bed.gz");
        let file_path = output.join(file_name);
        let file_path = Path::new(&file_path);
        let file = File::create(file_path)?;

        let buf_writer = BufWriter::new(
            GzEncoder::new(file, Compression::default())
        );

        handle_map.insert(cluster_id, buf_writer);
    }

    for file in files {
        let file = file?;
        let reader = get_dynamic_reader(&file.path())?;
        for (index, line) in reader.lines().enumerate() {
            let line = line?;

            let mut parts = line.split('\t');
            
            let chr = parts.next();
            let start = parts.next();
            let end = parts.next();
            let barcode = parts.next();
            let read_support = parts.next();

            if let (
                Some(chr),
                Some(start),
                Some(end),
                Some(barcode),
                Some(read_support)
            ) = (chr, start, end, barcode, read_support) {

                let cluster_id = mapping.get_cluster_from_barcode(barcode);
                if let Some(cluster_id) = cluster_id {
                    let cluster_file = handle_map.get_mut(&cluster_id).unwrap();
                    cluster_file.write_all(
                        format!("{chr}\t{start}\t{end}\t{barcode}\t{read_support}\n").as_bytes()
                    )?;
                } else {
                    anyhow::bail!(
                        format!("No cluster assignment found for barcode: {barcode}")
                    )
                }
            } else {
                anyhow::bail!(
                    format!("Failed to parse fragments file at line {index}: {}", line)
                )
            }
        }
    }


    Ok(())

}      

