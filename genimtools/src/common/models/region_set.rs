use crate::common::utils::bed_file_to_df;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::common::models::Region;
use crate::common::consts::{CHR_COL_NAME, END_COL_NAME, START_COL_NAME};

use polars::prelude::*;

pub struct RegionSet {
    pub regions: DataFrame,
}

impl TryFrom<&Path> for RegionSet {
    type Error = Box<dyn std::error::Error>;

    fn try_from(value: &Path) -> Result<Self, Self::Error> {
        let regions = bed_file_to_df(value)?;
        Ok(RegionSet { regions })
    }
}

impl From<Vec<Region>> for RegionSet {
    fn from(regions: Vec<Region>) -> Self {
        let mut chrs = Vec::new();
        let mut starts = Vec::new();
        let mut ends = Vec::new();

        for region in regions {
            chrs.push(region.chr);
            starts.push(region.start);
            ends.push(region.end);
        }

        let regions = DataFrame::new(vec![
            Series::new(CHR_COL_NAME, chrs),
            Series::new(START_COL_NAME, starts),
            Series::new(END_COL_NAME, ends),
        ])
        .unwrap();

        RegionSet { regions }
    }
}

impl RegionSet {
    pub fn to_bed(&self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(path)?;

        let chrs = self.regions.column(CHR_COL_NAME).unwrap();
        let starts = self.regions.column(START_COL_NAME).unwrap();
        let ends = self.regions.column(END_COL_NAME).unwrap();

        // is there a better way to do this?
        for i in 0..self.regions.height() {
            let chr: String;
            let start: u32;
            let end: u32;

            // extract out the chr, start, and end values
            if let AnyValue::Utf8(v) = chrs.get(i)? {
                chr = v.to_string();
            } else {
                panic!("chr column must be of type Utf8");
            }

            if let AnyValue::UInt32(v) = starts.get(i)? {
                start = v;
            } else {
                panic!("start column must be of type UInt32");
            }

            if let AnyValue::UInt32(v) = ends.get(i)? {
                end = v;
            } else {
                panic!("end column must be of type UInt32");
            }

            let line = format!("{}\t{}\t{}\n", chr, start, end);
            file.write_all(line.as_bytes())?;
        }

        Ok(())
    }

    pub fn len(&self) -> usize {
        self.regions.height()
    }

    pub fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }

    pub fn chrs(&self) -> &Series {
        self.regions.column(CHR_COL_NAME).unwrap()
    }

    pub fn starts(&self) -> &Series {
        self.regions.column(START_COL_NAME).unwrap()
    }

    pub fn ends(&self) -> &Series {
        self.regions.column(END_COL_NAME).unwrap()
    }
}
