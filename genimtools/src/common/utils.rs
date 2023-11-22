use std::path::Path;
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::collections::HashMap;

use polars::prelude::*;
use polars::datatypes::DataType;

use crate::common::models::region::Region;
use crate::common::consts::{
    CHR_COL_NAME,
    START_COL_NAME,
    END_COL_NAME,
    DELIMITER
};

pub fn bed_file_to_df(
    path: &Path,
) -> Result<DataFrame, Box<dyn std::error::Error>> {

    let schema = Schema::from_iter(vec![
        Field::new(CHR_COL_NAME, DataType::Utf8),
        Field::new(START_COL_NAME, DataType::UInt32),
        Field::new(END_COL_NAME, DataType::UInt32),
    ]);

    let df = CsvReader::from_path(path)?
        .has_header(false)
        .with_schema(Some(Arc::new(schema)))
        .with_separator(DELIMITER as u8) // the [0] is needed to convert from 
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

pub fn extract_regions_from_bed_file(path: &Path) -> Result<Vec<Region>, Box<dyn Error>> {
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