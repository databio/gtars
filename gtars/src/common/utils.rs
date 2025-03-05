use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::hash::Hash;
use std::io::prelude::*;
// use std::io::{BufRead, BufReader, Cursor};
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
// use flate2::read::{GzDecoder, MultiGzDecoder};
use flate2::read::MultiGzDecoder;
// use reqwest::blocking::Client;
use rust_lapper::{Interval, Lapper};
// use std::error::Error;

use crate::common::models::region::Region;
use crate::common::models::universe::Universe;

///
/// Get a reader for either a gzip'd or non-gzip'd file.
///
/// # Arguments
///
/// - path: path to the file to read
///
pub fn get_dynamic_reader(path: &Path) -> Result<BufReader<Box<dyn Read>>> {
    let is_gzipped = path.extension() == Some(OsStr::new("gz"));
    let file = File::open(path).with_context(|| "Failed to open bed file.")?;

    let file: Box<dyn Read> = match is_gzipped {
        true => Box::new(MultiGzDecoder::new(file)),
        false => Box::new(file),
    };

    let reader = BufReader::new(file);

    Ok(reader)
}

///
/// Get a reader for url ling. Either for gzip'd or non-gzip'd file
///
/// # Arguments
///
/// - path: path to the file to read
///
// pub fn get_dynamic_reader_from_url(
//     url: &Path,
// ) -> Result<BufReader<Box<dyn std::io::Read>>, Box<dyn Error>> {
//     // Create an HTTP client and fetch the content
//     let mut url: String = url.to_str().unwrap().to_string();

//     let is_ftp: bool = url.starts_with("ftp");

//     if is_ftp {
//         println!("ftp is not fully implemented. Bugs could appear");
//         url = url.replacen("ftp://", "http://", 1);
//     }

//     let response = Client::new()
//         .get(&url)
//         .send()
//         .with_context(|| format!("Failed to fetch content from URL: {}", &url))?
//         .error_for_status()?
//         .bytes()?;

//     // Convert the response into a cursor for reading
//     let cursor = Cursor::new(response);

//     let is_gzipped = url.ends_with(".gz");

//     let reader: Box<dyn std::io::Read> = match is_gzipped {
//         true => Box::new(GzDecoder::new(cursor)),
//         false => Box::new(cursor),
//     };

//     Ok(BufReader::new(reader))
// }

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
        return get_dynamic_reader(&file_path);
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
/// Simple wrapper function that will create a [Lapper] object (an interval tree)
/// from a [Universe] struct.
///
/// # Arguments:
/// - universe: the universe to create the interval tree for.
pub fn create_interval_tree_from_universe(
    universe: &Universe,
) -> HashMap<String, Lapper<u32, u32>> {
    // instantiate the tree and list of intervals
    let mut tree: HashMap<String, Lapper<u32, u32>> = HashMap::new();
    let mut intervals: HashMap<String, Vec<Interval<u32, u32>>> = HashMap::new();

    for region in universe.regions.iter() {
        // create interval
        let interval = Interval {
            start: region.start,
            stop: region.end,
            val: universe.convert_region_to_id(region).unwrap(),
        };

        // use chr to get the vector of intervals
        let chr_intervals = intervals.entry(region.chr.clone()).or_default();

        // push interval to vector
        chr_intervals.push(interval);
    }

    // build the tree
    for (chr, chr_intervals) in intervals.into_iter() {
        let lapper: Lapper<u32, u32> = Lapper::new(chr_intervals);
        tree.insert(chr.to_string(), lapper);
    }

    tree
}

pub fn get_chrom_sizes<T: AsRef<Path>>(path: T) -> HashMap<String, u32> {
    let chrom_sizes_file = File::open(path.as_ref())
        .with_context(|| format!("Failed to open chrom sizes file."))
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

        let region = Region {
            chr: chr.to_owned(),
            start,
            end,
            rest: None,
        };

        map.insert(region, score);
    }

    Ok(map)
}
