use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::models::Region;

pub fn extract_regions_from_bed_file(
    path: &Path,
) -> Result<Vec<Region>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut regions = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        let chr = fields[0];
        let start = fields[1].parse::<u32>()?;
        let end = fields[2].parse::<u32>()?;

        let region = Region {
            chr: chr.to_string(),
            start,
            end,
        };

        regions.push(region);
    }
    Ok(regions)
}
