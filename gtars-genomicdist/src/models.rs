use crate::errors::GtarsGenomicDistError;
use bio::io::fasta;
use gtars_core::models::{Region, RegionSet};
// use gtars_overlaprs::Overlapper;
use std::{collections::HashMap, fmt::Debug};
// use num_traits::{PrimInt, Unsigned};
use std::fs::File;
use std::path::{Path, PathBuf};

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
    pub n: u32,
    /// Rid: needed for plot to have correct order
    pub rid: u32,
}

pub struct GenomeAssembly {
    seq_map: HashMap<String, Vec<u8>>,
}

impl TryFrom<&str> for GenomeAssembly {
    type Error = GtarsGenomicDistError;

    fn try_from(value: &str) -> Result<Self, GtarsGenomicDistError> {
        GenomeAssembly::try_from(Path::new(value))
    }
}

impl TryFrom<String> for GenomeAssembly {
    type Error = GtarsGenomicDistError;

    fn try_from(value: String) -> Result<Self, GtarsGenomicDistError> {
        // println!("Converting String to Path: {}", value);
        GenomeAssembly::try_from(Path::new(&value))
    }
}

impl TryFrom<&Path> for GenomeAssembly {
    type Error = GtarsGenomicDistError;

    ///
    /// Create a new [GenomeAssembly] from fasta file
    ///
    fn try_from(value: &Path) -> Result<GenomeAssembly, GtarsGenomicDistError> {
        let file = File::open(value)?;
        let genome = fasta::Reader::new(file);

        let records = genome.records();

        // store the genome in a hashmap
        let mut seq_map: HashMap<String, Vec<u8>> = HashMap::new();
        for record in records {
            match record {
                Ok(record) => {
                    seq_map.insert(record.id().to_string(), record.seq().to_owned());
                }
                Err(e) => {
                    return Err(GtarsGenomicDistError::CustomError(format!(
                        "Error reading genome file: {}",
                        e
                    )));
                }
            }
        }

        Ok(GenomeAssembly { seq_map })
    }
}

impl GenomeAssembly {
    pub fn seq_from_region<'a>(&self, coords: &Region) -> Result<&[u8], GtarsGenomicDistError> {
        let chr = &coords.chr;
        let start = coords.start as usize;
        let end = coords.end as usize;

        if let Some(seq) = self.seq_map.get(chr) {
            if end <= seq.len() && start <= end {
                Ok(&seq[start..end])
            } else {
                Err(GtarsGenomicDistError::CustomError(format!(
                    "Invalid range: start={}, end={} for chromosome {} with length {}",
                    start,
                    end,
                    chr,
                    seq.len()
                )))
            }
        } else {
            Err(GtarsGenomicDistError::CustomError(format!(
                "Unknown chromosome found in region set: {}",
                chr
            )))
        }
    }

    pub fn contains_chr(&self, chr: &str) -> bool {
        self.seq_map.contains_key(chr)
    }
}

#[derive(PartialEq, Eq, Hash)]
pub enum Dinucleotide {
    Aa,
    Ac,
    Ag,
    At,
    Ca,
    Cc,
    Cg,
    Ct,
    Ga,
    Gc,
    Gg,
    Gt,
    Ta,
    Tc,
    Tg,
    Tt,
}

impl Dinucleotide {
    pub fn from_bytes(bytes: &[u8]) -> Option<Dinucleotide> {
        if bytes.len() != 2 {
            return None;
        }
        // Normalize to uppercase for case-insensitive matching
        let normalized = [bytes[0].to_ascii_uppercase(), bytes[1].to_ascii_uppercase()];
        match &normalized {
            b"AA" => Some(Dinucleotide::Aa),
            b"AC" => Some(Dinucleotide::Ac),
            b"AG" => Some(Dinucleotide::Ag),
            b"AT" => Some(Dinucleotide::At),
            b"CA" => Some(Dinucleotide::Ca),
            b"CC" => Some(Dinucleotide::Cc),
            b"CG" => Some(Dinucleotide::Cg),
            b"CT" => Some(Dinucleotide::Ct),
            b"GA" => Some(Dinucleotide::Ga),
            b"GC" => Some(Dinucleotide::Gc),
            b"GG" => Some(Dinucleotide::Gg),
            b"GT" => Some(Dinucleotide::Gt),
            b"TA" => Some(Dinucleotide::Ta),
            b"TC" => Some(Dinucleotide::Tc),
            b"TG" => Some(Dinucleotide::Tg),
            b"TT" => Some(Dinucleotide::Tt),
            _ => None,
        }
    }

    pub fn to_string(&self) -> Result<String, GtarsGenomicDistError> {
        match self {
            Dinucleotide::Aa => Ok("Aa".to_string()),
            Dinucleotide::Ac => Ok("Ac".to_string()),
            Dinucleotide::Ag => Ok("Ag".to_string()),
            Dinucleotide::At => Ok("At".to_string()),
            Dinucleotide::Ca => Ok("Ca".to_string()),
            Dinucleotide::Cc => Ok("Cc".to_string()),
            Dinucleotide::Cg => Ok("Cg".to_string()),
            Dinucleotide::Ct => Ok("Ct".to_string()),
            Dinucleotide::Ga => Ok("Ga".to_string()),
            Dinucleotide::Gc => Ok("Gc".to_string()),
            Dinucleotide::Gg => Ok("Gg".to_string()),
            Dinucleotide::Gt => Ok("Gt".to_string()),
            Dinucleotide::Ta => Ok("Ta".to_string()),
            Dinucleotide::Tc => Ok("Tc".to_string()),
            Dinucleotide::Tg => Ok("Tg".to_string()),
            Dinucleotide::Tt => Ok("Tt".to_string()),
        }
    }
}

///
/// Struct to hold Tss information (RegionSet with additionally indexing) that is initialized from
/// RegionSet or BED file that holds tss regions
///
pub struct TssIndex {
    pub region_set: RegionSet,
    pub mid_points: HashMap<String, Vec<u32>>,
}

impl TryFrom<RegionSet> for TssIndex {
    type Error = GtarsGenomicDistError;
    fn try_from(value: RegionSet) -> Result<Self, GtarsGenomicDistError> {
        let mut mid_points = value.calc_mid_points();

        for points in mid_points.values_mut() {
            points.sort_unstable();
        }

        Ok(TssIndex {
            region_set: value,
            mid_points,
        })
    }
}

impl TryFrom<&Path> for TssIndex {
    type Error = GtarsGenomicDistError;
    fn try_from(value: &Path) -> Result<Self, GtarsGenomicDistError> {
        let region_set = match RegionSet::try_from(value) {
            Ok(region_set) => region_set,
            Err(_e) => {
                return Err(GtarsGenomicDistError::TSSContentError(String::from(
                    "Unable to open Tss file",
                )));
            }
        };
        TssIndex::try_from(region_set)
    }
}

impl TryFrom<&str> for TssIndex {
    type Error = GtarsGenomicDistError;
    fn try_from(value: &str) -> Result<Self, GtarsGenomicDistError> {
        let region_set = match RegionSet::try_from(value) {
            Ok(region_set) => region_set,
            Err(_e) => {
                return Err(GtarsGenomicDistError::TSSContentError(String::from(
                    "Unable to open Tss file",
                )));
            }
        };
        TssIndex::try_from(region_set)
    }
}

impl TryFrom<String> for TssIndex {
    type Error = GtarsGenomicDistError;
    fn try_from(value: String) -> Result<Self, GtarsGenomicDistError> {
        let region_set = match RegionSet::try_from(value) {
            Ok(region_set) => region_set,
            Err(_e) => {
                return Err(GtarsGenomicDistError::TSSContentError(String::from(
                    "Unable to open Tss file",
                )));
            }
        };
        TssIndex::try_from(region_set)
    }
}

impl TssIndex {
    ///
    /// Calculate the distance from each region to the nearest TSS mid-point.
    ///
    /// Uses binary search for O(R * log M) complexity instead of O(R * M),
    /// where R is the number of regions and M is the number of TSS midpoints.
    ///
    pub fn calc_tss_distances(&self, rs: &RegionSet) -> Result<Vec<u32>, GtarsGenomicDistError> {
        let mut distances: Vec<u32> = Vec::with_capacity(rs.len());

        for chromosome in rs.iter_chroms() {
            if let Some(chr_midpoints) = self.mid_points.get(chromosome.as_str()) {
                for region in rs.iter_chr_regions(chromosome.as_str()) {
                    let target = region.mid_point();

                    let min_distance = match chr_midpoints.binary_search(&target) {
                        Ok(_) => 0,
                        Err(idx) => {
                            let left = idx
                                .checked_sub(1)
                                .map(|i| target.abs_diff(chr_midpoints[i]));
                            let right = chr_midpoints.get(idx).map(|&v| target.abs_diff(v));

                            match (left, right) {
                                (Some(l), Some(r)) => l.min(r),
                                (Some(l), None) => l,
                                (None, Some(r)) => r,
                                (None, None) => continue,
                            }
                        }
                    };
                    distances.push(min_distance);
                }
            }
        }
        Ok(distances)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Error;
    use std::path::PathBuf;

    use pretty_assertions::assert_eq;
    use rstest::*;

    fn get_test_path(file_name: &str) -> Result<PathBuf, Error> {
        let file_path: PathBuf = std::env::current_dir()
            .unwrap()
            .join("../tests/data/regionset")
            .join(file_name);
        Ok(file_path)
    }

    #[rstest]
    fn test_calc_tss_distances() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let tss_path = get_test_path("dummy_tss.bed").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
        let tss_index = TssIndex::try_from(tss_path.to_str().unwrap()).unwrap();

        let distances = tss_index.calc_tss_distances(&region_set).unwrap();

        assert_eq!(distances.len(), 9);
        assert_eq!(distances.iter().min(), Some(&2));
    }
}
