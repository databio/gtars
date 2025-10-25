/// Models for used in calculation statistics for regionset
///
#[derive(Debug, Clone)]
pub struct ChromosomeStats {
    pub chromosome: String,
    pub count: u32,   // number of regions
    pub start: u32,   // first nucleotide
    pub end: u32,     // end nucleotide
    pub minimum: u32, // smallest region
    pub maximum: u32, // largest region
    pub mean: f64,    // mean region width
    pub median: f64,  // median region width
}

#[derive(Debug, Clone)]
pub struct RegionBin {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub n: u32,
    pub rid: u32,
}
