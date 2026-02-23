//! Consensus region analysis for multiple region sets.
//!
//! Given N region sets, computes the union of all regions and annotates each
//! union region with the number of input sets that overlap it. This enables
//! filtering by support threshold (e.g., "regions present in at least 2/3
//! replicates").

use gtars_core::models::{Region, RegionSet};
use gtars_overlaprs::{OverlapperType, multi_chrom_overlapper::IntoMultiChromOverlapper};

use crate::interval_ranges::IntervalRanges;

/// A region annotated with the number of input sets overlapping it.
#[derive(Debug, Clone, PartialEq)]
pub struct ConsensusRegion {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    /// Number of input region sets that overlap this region.
    pub count: u32,
}

/// Compute consensus regions from multiple region sets.
///
/// 1. Concatenates all input sets and reduces to a non-overlapping union.
/// 2. For each union region, counts how many input sets have at least one
///    overlapping interval (using AIList for fast queries).
/// 3. Returns `ConsensusRegion`s sorted by (chr, start).
///
/// Returns an empty vec if `sets` is empty.
pub fn consensus(sets: &[RegionSet]) -> Vec<ConsensusRegion> {
    if sets.is_empty() {
        return Vec::new();
    }

    // Build union of all input sets
    let mut all_regions: Vec<Region> = Vec::new();
    for set in sets {
        all_regions.extend(set.regions.iter().cloned());
    }
    let union = RegionSet::from(all_regions).reduce();

    // Build an AIList overlapper per input set
    let overlappers: Vec<_> = sets
        .iter()
        .map(|s| s.clone().into_multi_chrom_overlapper(OverlapperType::AIList))
        .collect();

    // For each union region, count how many input sets overlap it
    union
        .regions
        .iter()
        .map(|r| {
            let query = RegionSet::from(vec![Region {
                chr: r.chr.clone(),
                start: r.start,
                end: r.end,
                rest: None,
            }]);
            let count = overlappers
                .iter()
                .filter(|ov| !ov.find_overlaps(&query).is_empty())
                .count() as u32;
            ConsensusRegion {
                chr: r.chr.clone(),
                start: r.start,
                end: r.end,
                count,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use std::path::PathBuf;

    fn get_test_path(file_name: &str) -> PathBuf {
        std::env::current_dir()
            .unwrap()
            .join("../tests/data/regionset")
            .join(file_name)
    }

    fn make_regionset(regions: Vec<(&str, u32, u32)>) -> RegionSet {
        let regions: Vec<Region> = regions
            .into_iter()
            .map(|(chr, start, end)| Region {
                chr: chr.to_string(),
                start,
                end,
                rest: None,
            })
            .collect();
        RegionSet::from(regions)
    }

    #[rstest]
    fn test_consensus_two_sets() {
        // dummy.bed: chr1:2-6, chr1:4-7, chr1:5-9, chr1:7-12 -> reduces to chr1:2-12
        // dummy_b.bed: chr1:3-5, chr1:8-10
        // Union = chr1:2-12 (single region)
        // Both sets overlap chr1:2-12, so count = 2
        let path_a = get_test_path("dummy.bed");
        let path_b = get_test_path("dummy_b.bed");
        let a = RegionSet::try_from(path_a.to_str().unwrap()).unwrap();
        let b = RegionSet::try_from(path_b.to_str().unwrap()).unwrap();
        let result = consensus(&[a, b]);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].chr, "chr1");
        assert_eq!(result[0].start, 2);
        assert_eq!(result[0].end, 12);
        assert_eq!(result[0].count, 2);
    }

    #[rstest]
    fn test_consensus_identical() {
        let a = make_regionset(vec![("chr1", 10, 20), ("chr1", 30, 40)]);
        let result = consensus(&[a.clone(), a]);
        assert_eq!(result.len(), 2);
        assert!(result.iter().all(|r| r.count == 2));
    }

    #[rstest]
    fn test_consensus_single_set() {
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]);
        let result = consensus(&[a]);
        assert_eq!(result.len(), 2);
        assert!(result.iter().all(|r| r.count == 1));
    }

    #[rstest]
    fn test_consensus_empty() {
        let result = consensus(&[]);
        assert!(result.is_empty());
    }

    #[rstest]
    fn test_consensus_partial_overlap() {
        // Set A: chr1:0-10, chr1:20-30
        // Set B: chr1:5-15
        // Union: chr1:0-15, chr1:20-30
        // chr1:0-15 overlapped by both, chr1:20-30 overlapped by A only
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let result = consensus(&[a, b]);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].count, 2); // chr1:0-15
        assert_eq!(result[1].count, 1); // chr1:20-30
    }
}
