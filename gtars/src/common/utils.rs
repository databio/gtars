use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};

use flate2::read::MultiGzDecoder;

use crate::common::models::region::Region;

///
/// Get a reader for either a gzip'd or non-gzip'd file.
///
/// # Arguments
///
/// - path: path to the file to read
///
pub fn get_dynamic_reader(path: &Path) -> Result<BufReader<Box<dyn Read>>> {
    let is_gzipped = path.extension() == Some(OsStr::new("gz"));
    let file = File::open(path).with_context(|| format!("Failed to open file: {:?}", path))?;

    let file: Box<dyn Read> = match is_gzipped {
        true => Box::new(MultiGzDecoder::new(file)),
        false => Box::new(file),
    };

    let reader = BufReader::new(file);

    Ok(reader)
}

/// Get a reader for either a gzipped, non-gzipped file, or stdin
///
/// # Arguments
///
/// - file_path: path to the file to read, or '-' for stdin
///
/// # Returns
///
/// A `BufReader` object for a given file path or stdin.
pub fn get_dynamic_reader_w_stdin(file_path_str: &str) -> Result<BufReader<Box<dyn Read>>> {
    if file_path_str == "-" {
        Ok(BufReader::new(Box::new(std::io::stdin()) as Box<dyn Read>))
    } else {
        let file_path = Path::new(file_path_str);
        get_dynamic_reader(file_path)
    }
}

///
/// Create a region-to-id hash-map from a list of regions
///
/// # Arguments:
/// - regions: vec![] of [Region] structs
pub fn generate_region_to_id_map(regions: &[Region]) -> HashMap<Region, u32> {
    let mut current_id = 0;
    let mut region_to_id: HashMap<Region, u32> = HashMap::new();
    for region in regions.iter() {
        region_to_id.entry(region.to_owned()).or_insert_with(|| {
            let old_id = current_id;
            current_id += 1;
            old_id
        });
    }

    region_to_id
}

///
/// Generate an id-to-region hash-map from a list of regions
///
/// # Arguments:
/// - regions: vec![] of [Region] structs
pub fn generate_id_to_region_map(regions: &[Region]) -> HashMap<u32, Region> {
    let mut current_id = 0;
    let mut id_to_region: HashMap<u32, Region> = HashMap::new();

    for region in regions.iter() {
        id_to_region.entry(current_id).or_insert_with(|| {
            current_id += 1;
            region.clone()
        });
    }

    id_to_region
}

///
/// Create a region-to-id hash-map from a list of region strings
///
/// # Arguments:
/// - regions: vec![] of region strings in the form `chr:start-end`
pub fn generate_region_string_to_id_map(regions: &[String]) -> HashMap<String, u32> {
    let mut current_id = 0;
    let mut region_to_id: HashMap<String, u32> = HashMap::new();
    for region in regions.iter() {
        region_to_id.entry(region.to_owned()).or_insert_with(|| {
            let old_id = current_id;
            current_id += 1;
            old_id
        });
    }

    region_to_id
}

///
/// Generate an id-to-region string hash-map from a list of region strings
///
/// # Arguments:
/// - regions: vec![] of region strings in the form `chr:start-end`
pub fn generate_id_to_region_string_map(regions: &[String]) -> HashMap<u32, String> {
    let mut current_id = 0;
    let mut id_to_region: HashMap<u32, String> = HashMap::new();

    for region in regions.iter() {
        id_to_region.entry(current_id).or_insert_with(|| {
            current_id += 1;
            region.clone()
        });
    }

    id_to_region
}

pub fn get_chrom_sizes<T: AsRef<Path>>(path: T) -> HashMap<String, u32> {
    let chrom_sizes_file = File::open(path.as_ref())
        .with_context(|| "Failed to open chrom sizes file.")
        .unwrap();

    let mut chrom_sizes: HashMap<String, u32> = HashMap::new();

    let file_buf = BufReader::new(chrom_sizes_file);

    for line in file_buf.lines() {
        let line_string: String = match line {
            Ok(value) => value,
            Err(_) => panic!("Error while reading chrom sizes file"),
        };

        let line_parts: Vec<String> = line_string
            .split_whitespace()
            .map(|s| s.to_string())
            .collect();

        chrom_sizes.insert(line_parts[0].clone(), line_parts[1].parse::<u32>().unwrap());
    }

    chrom_sizes
}

///
/// Gen
pub fn generate_ordering_map_for_universe_regions<T: AsRef<Path>>(
    path: T,
) -> Result<HashMap<Region, f64>> {
    let mut map = HashMap::new();

    let reader = get_dynamic_reader(path.as_ref())?;

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() < 5 {
            anyhow::bail!("BED file line does not have at least 5 fields: {}. It needs to have chr, start, end, name, and score.", line);
        }

        // parse the fields
        let chr = parts[0];
        let start = parts[1].parse::<u32>().with_context(|| {
            format!("Failed to parse start position in BED file line: {}", line)
        })?;

        let end = parts[2]
            .parse::<u32>()
            .with_context(|| format!("Failed to parse end position in BED file line: {}", line))?;

        let score = parts[4]
            .parse::<f64>()
            .with_context(|| format!("Failed to parse score in BED file line: {}", line))?;

        let rest = Some(parts[3..].join("\t")).filter(|s| !s.is_empty());

        let region = Region {
            chr: chr.to_owned(),
            start,
            end,
            rest,
        };

        map.insert(region, score);
    }

    Ok(map)
}
