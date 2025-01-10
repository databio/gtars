use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use rust_lapper::{Interval, Lapper};

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
/// Read in a bed file into a vector of [Region] structs. It handles detecting
/// the file-type, verifying each line, and error handling.
///
/// # Arguments:
/// - path: path to the bed file to read in.
pub fn extract_regions_from_bed_file(path: &Path) -> Result<Vec<Region>> {
    let reader = get_dynamic_reader(path)?;

    let mut regions = Vec::new();

    for line in reader.lines() {
        let line = line.with_context(|| "Failed parsing line in BED file")?;
        let fields: Vec<&str> = line.split('\t').collect();

        // check length of fields
        if fields.len() < 3 {
            anyhow::bail!("BED file line does not have at least 3 fields: {}", line);
        }

        let chr = fields[0];
        let start = fields[1].parse::<u32>().with_context(|| {
            format!("Failed to parse start position in BED file line: {}", line)
        })?;
        let end = fields[2]
            .parse::<u32>()
            .with_context(|| format!("Failed to parse end position in BED file line: {}", line))?;

        let region = Region {
            chr: chr.to_string(),
            start,
            end,
        };

        regions.push(region);
    }

    Ok(regions)
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
