use super::utils::extract_regions_from_bed_file;
use std::path::Path;

pub struct Region {
    pub chr: String,
    pub start: u32,
    pub end: u32,
}

pub struct RegionSet {
    pub regions: Vec<Region>,
}

pub struct BedSet {
    pub region_sets: Vec<RegionSet>,
}

impl TryFrom<&Path> for RegionSet {
    type Error = Box<dyn std::error::Error>;

    fn try_from(value: &Path) -> Result<Self, Self::Error> {
        let regions = extract_regions_from_bed_file(value)?;
        Ok(RegionSet { regions })
    }
}

impl From<Vec<Region>> for RegionSet {
    fn from(value: Vec<Region>) -> Self {
        RegionSet { regions: value }
    }
}
