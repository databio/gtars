use crate::errors::GtarsGenomicDistError;
use bio::io::fasta;
use gtars_core::models::{CoordinateMode, Region, RegionSet};
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, fmt::Debug};
use std::fs::File;
use std::path::Path;

/// A `RegionSet` that is guaranteed to be sorted by (chr, start).
///
/// Created by moving a `RegionSet` into `SortedRegionSet::new()`, which
/// sorts in place (no clone). Functions that require sorted input can
/// accept this type instead of re-sorting every time.
pub struct SortedRegionSet(pub RegionSet);

impl SortedRegionSet {
    /// Sort a RegionSet in place and wrap it. This is a move, not a clone.
    pub fn new(mut rs: RegionSet) -> Self {
        rs.sort();
        Self(rs)
    }
}

/// Genomic strand orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Strand {
    Plus,
    Minus,
    Unstranded,
}

impl Strand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => Strand::Plus,
            '-' => Strand::Minus,
            _ => Strand::Unstranded,
        }
    }
}

/// A `RegionSet` paired with a parallel `Vec<Strand>`.
///
/// Follows the `SortedRegionSet` wrapper pattern: wraps a `RegionSet` and
/// adds strand information without modifying `Region` itself (which lives
/// in gtars-core). Strand-aware operations (promoters, reduce, setdiff)
/// are implemented as methods on this type.
#[derive(Clone, Serialize, Deserialize)]
pub struct StrandedRegionSet {
    pub inner: RegionSet,
    pub strands: Vec<Strand>,
}

impl StrandedRegionSet {
    /// Create a new StrandedRegionSet from a RegionSet and parallel strand vector.
    ///
    /// # Panics
    /// Panics if `strands.len() != rs.regions.len()`.
    pub fn new(rs: RegionSet, strands: Vec<Strand>) -> Self {
        assert_eq!(
            rs.regions.len(),
            strands.len(),
            "StrandedRegionSet: regions and strands must have the same length"
        );
        StrandedRegionSet {
            inner: rs,
            strands,
        }
    }

    /// Wrap a RegionSet with all-Unstranded. Preserves existing behavior.
    pub fn unstranded(rs: RegionSet) -> Self {
        let n = rs.regions.len();
        StrandedRegionSet {
            strands: vec![Strand::Unstranded; n],
            inner: rs,
        }
    }

    /// Consume into the inner RegionSet, dropping strand information.
    pub fn into_regionset(self) -> RegionSet {
        self.inner
    }

    pub fn len(&self) -> usize {
        self.inner.regions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.regions.is_empty()
    }
}

/// Statistics summary for regions on a single chromosome.
///
/// Contains counts, bounds, and descriptive statistics for region lengths.
#[derive(Debug, Clone, Serialize, Deserialize)]
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
#[derive(Debug, Clone, Serialize, Deserialize)]
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

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
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
        TssIndex::from_region_set(value, CoordinateMode::Bed)
    }
}

impl TssIndex {
    /// Create a TssIndex from a RegionSet using the specified coordinate mode for midpoints.
    pub fn from_region_set(
        value: RegionSet,
        mode: CoordinateMode,
    ) -> Result<Self, GtarsGenomicDistError> {
        let mut mid_points = value.calc_mid_points_with_mode(mode);

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
    pub fn calc_tss_distances(
        &self,
        rs: &RegionSet,
        mode: CoordinateMode,
    ) -> Result<Vec<u32>, GtarsGenomicDistError> {
        let mut distances: Vec<u32> = Vec::with_capacity(rs.len());

        for chromosome in rs.iter_chroms() {
            if let Some(chr_midpoints) = self.mid_points.get(chromosome.as_str()) {
                for region in rs.iter_chr_regions(chromosome.as_str()) {
                    let target = region.mid_point_with_mode(mode);

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
            } else {
                // No features on this chromosome — push u32::MAX for each region
                for _ in rs.iter_chr_regions(chromosome.as_str()) {
                    distances.push(u32::MAX);
                }
            }
        }
        Ok(distances)
    }

    /// Calculate signed distances from each region to its nearest feature.
    ///
    /// Like `calc_tss_distances` but returns signed distances where:
    /// - Positive: nearest feature is downstream (right) of the query
    /// - Negative: nearest feature is upstream (left) of the query
    ///
    /// Sign convention: `nearest_feature_midpoint - query_midpoint`
    /// (matches R GenomicDistributions `calcFeatureDist()`).
    pub fn calc_feature_distances(
        &self,
        rs: &RegionSet,
        mode: CoordinateMode,
    ) -> Result<Vec<i64>, GtarsGenomicDistError> {
        let mut distances: Vec<i64> = Vec::with_capacity(rs.len());

        for chromosome in rs.iter_chroms() {
            if let Some(chr_midpoints) = self.mid_points.get(chromosome.as_str()) {
                for region in rs.iter_chr_regions(chromosome.as_str()) {
                    let target = region.mid_point_with_mode(mode) as i64;

                    let distance = match chr_midpoints.binary_search(&(target as u32)) {
                        Ok(_) => 0i64,
                        Err(idx) => {
                            // distance = feature_mid - query_mid
                            let left = idx
                                .checked_sub(1)
                                .map(|i| chr_midpoints[i] as i64 - target);
                            let right =
                                chr_midpoints.get(idx).map(|&v| v as i64 - target);

                            match (left, right) {
                                (Some(l), Some(r)) => {
                                    if l.unsigned_abs() <= r.unsigned_abs() {
                                        l
                                    } else {
                                        r
                                    }
                                }
                                (Some(l), None) => l,
                                (None, Some(r)) => r,
                                (None, None) => continue,
                            }
                        }
                    };
                    distances.push(distance);
                }
            } else {
                // No features on this chromosome — push i64::MAX for each region
                for _ in rs.iter_chr_regions(chromosome.as_str()) {
                    distances.push(i64::MAX);
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

    fn get_fasta_path(file_name: &str) -> PathBuf {
        std::env::current_dir()
            .unwrap()
            .join("../tests/data/fasta")
            .join(file_name)
    }

    // --- Strand ---

    #[test]
    fn test_strand_from_char() {
        assert_eq!(Strand::from_char('+'), Strand::Plus);
        assert_eq!(Strand::from_char('-'), Strand::Minus);
        assert_eq!(Strand::from_char('.'), Strand::Unstranded);
        assert_eq!(Strand::from_char('?'), Strand::Unstranded);
    }

    // --- Dinucleotide ---

    #[test]
    fn test_dinucleotide_from_bytes_all_variants() {
        let pairs = [
            (b"AA", Dinucleotide::Aa), (b"AC", Dinucleotide::Ac),
            (b"AG", Dinucleotide::Ag), (b"AT", Dinucleotide::At),
            (b"CA", Dinucleotide::Ca), (b"CC", Dinucleotide::Cc),
            (b"CG", Dinucleotide::Cg), (b"CT", Dinucleotide::Ct),
            (b"GA", Dinucleotide::Ga), (b"GC", Dinucleotide::Gc),
            (b"GG", Dinucleotide::Gg), (b"GT", Dinucleotide::Gt),
            (b"TA", Dinucleotide::Ta), (b"TC", Dinucleotide::Tc),
            (b"TG", Dinucleotide::Tg), (b"TT", Dinucleotide::Tt),
        ];
        for (bytes, expected) in &pairs {
            assert_eq!(Dinucleotide::from_bytes(&bytes[..]), Some(*expected));
        }
    }

    #[test]
    fn test_dinucleotide_case_insensitive() {
        assert_eq!(Dinucleotide::from_bytes(b"aa"), Some(Dinucleotide::Aa));
        assert_eq!(Dinucleotide::from_bytes(b"cG"), Some(Dinucleotide::Cg));
        assert_eq!(Dinucleotide::from_bytes(b"Tc"), Some(Dinucleotide::Tc));
    }

    #[test]
    fn test_dinucleotide_invalid() {
        assert_eq!(Dinucleotide::from_bytes(b"AN"), None);
        assert_eq!(Dinucleotide::from_bytes(b"A"), None);  // too short
        assert_eq!(Dinucleotide::from_bytes(b"ACG"), None); // too long
    }

    #[test]
    fn test_dinucleotide_to_string_round_trip() {
        let all = [
            Dinucleotide::Aa, Dinucleotide::Ac, Dinucleotide::Ag, Dinucleotide::At,
            Dinucleotide::Ca, Dinucleotide::Cc, Dinucleotide::Cg, Dinucleotide::Ct,
            Dinucleotide::Ga, Dinucleotide::Gc, Dinucleotide::Gg, Dinucleotide::Gt,
            Dinucleotide::Ta, Dinucleotide::Tc, Dinucleotide::Tg, Dinucleotide::Tt,
        ];
        for d in &all {
            let s = d.to_string().unwrap();
            assert_eq!(s.len(), 2);
            let round_tripped = Dinucleotide::from_bytes(s.as_bytes()).unwrap();
            assert_eq!(*d, round_tripped);
        }
    }

    // --- SortedRegionSet ---

    #[test]
    fn test_sorted_regionset_sorts_in_place() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 100, end: 200, rest: None },
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
            Region { chr: "chr2".into(), start: 5, end: 15, rest: None },
        ];
        let sorted = SortedRegionSet::new(RegionSet::from(regions));
        let starts: Vec<u32> = sorted.0.regions.iter().map(|r| r.start).collect();
        // chr1 regions should come first (sorted by chr then start)
        assert_eq!(starts, vec![10, 100, 5]);
    }

    // --- StrandedRegionSet ---

    #[test]
    fn test_stranded_regionset_new() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
            Region { chr: "chr1".into(), start: 30, end: 40, rest: None },
        ];
        let strands = vec![Strand::Plus, Strand::Minus];
        let srs = StrandedRegionSet::new(RegionSet::from(regions), strands);
        assert_eq!(srs.len(), 2);
        assert!(!srs.is_empty());
        assert_eq!(srs.strands[0], Strand::Plus);
        assert_eq!(srs.strands[1], Strand::Minus);
    }

    #[test]
    fn test_stranded_regionset_unstranded() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
        ];
        let srs = StrandedRegionSet::unstranded(RegionSet::from(regions));
        assert_eq!(srs.strands, vec![Strand::Unstranded]);
    }

    #[test]
    #[should_panic(expected = "regions and strands must have the same length")]
    fn test_stranded_regionset_mismatched_lengths() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
        ];
        StrandedRegionSet::new(RegionSet::from(regions), vec![]);
    }

    #[test]
    fn test_stranded_regionset_into_regionset() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
        ];
        let srs = StrandedRegionSet::unstranded(RegionSet::from(regions));
        let rs = srs.into_regionset();
        assert_eq!(rs.regions.len(), 1);
    }

    // --- GenomeAssembly ---

    #[test]
    fn test_genome_assembly_from_fasta() {
        // base.fa: chrX=TTGGGGAA, chr1=GGAA, chr2=GCGC
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.as_path()).unwrap();
        assert!(ga.contains_chr("chr1"));
        assert!(ga.contains_chr("chr2"));
        assert!(ga.contains_chr("chrX"));
        assert!(!ga.contains_chr("chr3"));
    }

    #[test]
    fn test_genome_assembly_seq_from_region() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.as_path()).unwrap();

        let region = Region { chr: "chr1".into(), start: 0, end: 4, rest: None };
        let seq = ga.seq_from_region(&region).unwrap();
        assert_eq!(seq, b"GGAA");

        let region2 = Region { chr: "chrX".into(), start: 2, end: 6, rest: None };
        let seq2 = ga.seq_from_region(&region2).unwrap();
        assert_eq!(seq2, b"GGGG");
    }

    #[test]
    fn test_genome_assembly_seq_unknown_chrom() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.as_path()).unwrap();

        let region = Region { chr: "chr99".into(), start: 0, end: 1, rest: None };
        assert!(ga.seq_from_region(&region).is_err());
    }

    #[test]
    fn test_genome_assembly_seq_out_of_bounds() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.as_path()).unwrap();

        // chr1 is only 4bp, request 0-100
        let region = Region { chr: "chr1".into(), start: 0, end: 100, rest: None };
        assert!(ga.seq_from_region(&region).is_err());
    }

    #[test]
    fn test_genome_assembly_try_from_str() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap());
        assert!(ga.is_ok());
    }

    #[test]
    fn test_genome_assembly_try_from_string() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap().to_string());
        assert!(ga.is_ok());
        assert!(ga.unwrap().contains_chr("chr1"));
    }

    #[test]
    fn test_tss_index_try_from_path() {
        let path = get_test_path("dummy_tss.bed").unwrap();
        let tss = TssIndex::try_from(path.as_path());
        assert!(tss.is_ok());
    }

    #[test]
    fn test_tss_index_try_from_path_invalid() {
        let path = PathBuf::from("/nonexistent/file.bed");
        let tss = TssIndex::try_from(path.as_path());
        assert!(tss.is_err());
    }

    #[test]
    fn test_tss_index_try_from_string() {
        let path = get_test_path("dummy_tss.bed").unwrap();
        let tss = TssIndex::try_from(path.to_str().unwrap().to_string());
        assert!(tss.is_ok());
    }

    // --- TssIndex sentinel behavior ---

    #[test]
    fn test_tss_distances_sentinel_for_missing_chrom() {
        // TSS features only on chr1
        let tss_regions = vec![
            Region { chr: "chr1".into(), start: 50, end: 51, rest: None },
        ];
        let tss_index = TssIndex::try_from(RegionSet::from(tss_regions)).unwrap();

        // Query has regions on chr1 and chr2 (no TSS on chr2)
        let query = RegionSet::from(vec![
            Region { chr: "chr1".into(), start: 40, end: 45, rest: None },
            Region { chr: "chr2".into(), start: 10, end: 20, rest: None },
        ]);

        let distances = tss_index.calc_tss_distances(&query, CoordinateMode::Bed).unwrap();
        assert_eq!(distances.len(), 2); // one per input region
        // One should be a real distance, the other should be u32::MAX sentinel
        // (order depends on HashSet iteration of iter_chroms)
        assert_eq!(distances.iter().filter(|&&d| d == u32::MAX).count(), 1);
        assert_eq!(distances.iter().filter(|&&d| d < u32::MAX).count(), 1);
    }

    #[test]
    fn test_feature_distances_sentinel_for_missing_chrom() {
        let tss_regions = vec![
            Region { chr: "chr1".into(), start: 50, end: 51, rest: None },
        ];
        let tss_index = TssIndex::try_from(RegionSet::from(tss_regions)).unwrap();

        let query = RegionSet::from(vec![
            Region { chr: "chr1".into(), start: 40, end: 45, rest: None },
            Region { chr: "chr2".into(), start: 10, end: 20, rest: None },
        ]);

        let distances = tss_index.calc_feature_distances(&query, CoordinateMode::Bed).unwrap();
        assert_eq!(distances.len(), 2);
        // One real distance, one i64::MAX sentinel
        assert_eq!(distances.iter().filter(|&&d| d == i64::MAX).count(), 1);
        assert_eq!(distances.iter().filter(|&&d| d != i64::MAX).count(), 1);
    }

    // --- Existing tests ---

    #[rstest]
    fn test_calc_tss_distances() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let tss_path = get_test_path("dummy_tss.bed").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
        let tss_index = TssIndex::try_from(tss_path.to_str().unwrap()).unwrap();

        let distances = tss_index.calc_tss_distances(&region_set, CoordinateMode::Bed).unwrap();

        assert_eq!(distances.len(), 9);
        assert_eq!(distances.iter().min(), Some(&2));
    }

    #[rstest]
    fn test_calc_feature_distances() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let tss_path = get_test_path("dummy_tss.bed").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
        let tss_index = TssIndex::try_from(tss_path.to_str().unwrap()).unwrap();

        let signed_distances = tss_index.calc_feature_distances(&region_set, CoordinateMode::Bed).unwrap();
        let abs_distances = tss_index.calc_tss_distances(&region_set, CoordinateMode::Bed).unwrap();

        // same number of results
        assert_eq!(signed_distances.len(), abs_distances.len());
        // absolute values should match calc_tss_distances
        for (signed, abs) in signed_distances.iter().zip(abs_distances.iter()) {
            assert_eq!(signed.unsigned_abs() as u32, *abs);
        }
        // should have both positive and negative distances
        assert!(signed_distances.iter().any(|d| *d > 0));
        assert!(signed_distances.iter().any(|d| *d < 0));
    }
}
