use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;

use crate::common::models::region::Region;

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

pub fn extract_regions_from_bed_file(path: &Path) -> Result<Vec<Region>> {
    let mut regions = Vec::new();

    let is_gzipped = path.extension() == Some(OsStr::new("gz"));
    let file = File::open(path).with_context(|| "Failed to open bed file.")?;

    let file: Box<dyn Read> = match is_gzipped {
        true => Box::new(GzDecoder::new(file)),
        false => Box::new(file),
    };

    let reader = BufReader::new(file);

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
