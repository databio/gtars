use std::path::Path;

use polars::prelude::*;
use polars::datatypes::DataType;
use super::consts::{
    CHR_COL_NAME,
    START_COL_NAME,
    END_COL_NAME,
    DELIMITER
};

pub fn extract_regions_from_bed_file(
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
