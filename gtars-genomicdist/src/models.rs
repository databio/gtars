/// Statistics summary for regions on a single chromosome.
///
/// Contains counts, bounds, and descriptive statistics for region lengths.
#[derive(Debug, Clone)]
pub struct ChromosomeStatistics {
    /// Chromosome name
    pub chromosome: String,
    /// Total number of regions on this chromosome
    pub number_of_regions: u32,
    /// Leftmost start position across all regions
    pub start_nucleotide_position: u32,
    /// Rightmost end position across all regions
    pub end_nucleotide_position: u32,
    /// Length of the shortest region
    pub minimum_region_length: u32,
    /// Length of the longest region
    pub maximum_region_length: u32,
    /// Average region length
    pub mean_region_length: f64,
    /// Median region length
    pub median_region_length: f64,
}

/// A genomic bin with a count of overlapping regions.
///
/// Used to represent distribution of regions across fixed-size windows.
#[derive(Debug, Clone)]
pub struct RegionBin {
    /// Chromosome name
    pub chr: String,
    /// Start position of the bin
    pub start: u32,
    /// End position of the bin
    pub end: u32,
    /// Number of regions overlapping this bin
    pub n: u32
}
