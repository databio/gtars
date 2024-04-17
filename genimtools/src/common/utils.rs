use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::io::prelude::*;
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use polars::datatypes::DataType;
use polars::prelude::*;

use crate::common::consts::{CHR_COL_NAME, DELIMITER, END_COL_NAME, START_COL_NAME};
use crate::common::models::region::Region;

pub fn bed_file_to_df(path: &Path) -> Result<DataFrame> {
    let schema = Schema::from_iter(vec![
        Field::new(CHR_COL_NAME, DataType::Utf8),
        Field::new(START_COL_NAME, DataType::UInt32),
        Field::new(END_COL_NAME, DataType::UInt32),
    ]);

    let df = CsvReader::from_path(path)?
        .has_header(false)
        .with_schema(Some(Arc::new(schema)))
        .with_separator(DELIMITER as u8)
        .truncate_ragged_lines(true)
        .finish()?;

    Ok(df)
}

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
