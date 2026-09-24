#![allow(dead_code)]

use gtars_core::models::{Region, RegionSet};
use std::collections::HashMap;
use std::path::PathBuf;

pub type Iv<'a> = (&'a str, u32, u32);

pub fn make_region(chr: &str, start: u32, end: u32) -> Region {
    Region {
        chr: chr.to_string(),
        start,
        end,
        rest: None,
    }
}

pub fn make_regionset(regions: Vec<Iv>) -> RegionSet {
    RegionSet::from(
        regions
            .into_iter()
            .map(|(chr, start, end)| make_region(chr, start, end))
            .collect::<Vec<Region>>(),
    )
}

pub fn coords(rs: &RegionSet) -> Vec<Iv<'_>> {
    rs.regions
        .iter()
        .map(|r| (r.chr.as_str(), r.start, r.end))
        .collect()
}

pub fn chrom_sizes(sizes: &[(&str, u32)]) -> HashMap<String, u32> {
    sizes.iter().map(|(c, s)| (c.to_string(), *s)).collect()
}

pub fn load_bed(file_name: &str) -> RegionSet {
    let path = PathBuf::from("../tests/data/regionset").join(file_name);
    RegionSet::try_from(path.as_path()).unwrap()
}

pub fn assert_approx(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-10,
        "expected {expected}, got {actual}"
    );
}
