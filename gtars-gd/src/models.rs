/// Models for used in calculation statistics for regionset
///
#[derive(Debug, Clone)]
pub struct ChromosomeStatistics {
    pub chromosome: String,
    pub number_of_regions: u32,
    pub start_nucleotide_position: u32,
    pub end_nucleotide_position: u32,
    pub minimum_region_length: u32,
    pub maximum_region_length: u32,
    pub mean_region_length: f64,
    pub median_region_length: f64,
}

#[derive(Debug, Clone)]
pub struct RegionBin {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub n: u32,
    pub rid: u32,
}
