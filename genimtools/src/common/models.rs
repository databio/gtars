use super::utils::extract_regions_from_bed_file;
use core::panic;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use super::consts::CHR_COL_NAME;
use polars::prelude::*;

pub struct Region {
    pub chr: String,
    pub start: u32,
    pub end: u32,
}

pub struct RegionSet {
    pub regions: DataFrame,
}

impl TryFrom<&Path> for RegionSet {
    type Error = Box<dyn std::error::Error>;

    fn try_from(value: &Path) -> Result<Self, Self::Error> {
        let regions = extract_regions_from_bed_file(value)?;
        Ok(RegionSet { regions })
    }
}

impl RegionSet {
    pub fn to_bed(&self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(path)?;

        let chrs = self.regions.column(CHR_COL_NAME).unwrap();
        let starts = self.regions.column("start").unwrap();
        let ends = self.regions.column("end").unwrap();
        
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
}


pub struct BedSet {
    pub region_sets: Vec<RegionSet>,
}
