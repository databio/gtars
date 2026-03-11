//! Overlap query operations on RegionSet.
//!
//! These methods build a MultiChromOverlapper index from `other` and query
//! each region in `self` against it.

use gtars_core::models::RegionSet;

use crate::OverlapperType;
use crate::multi_chrom_overlapper::build_indexed_overlapper;

/// Overlap query operations on RegionSet.
///
/// These methods build a MultiChromOverlapper index from `other` and query
/// each region in `self` against it.
pub trait RegionSetOverlaps {
    /// Return a new RegionSet containing only regions from self that overlap
    /// at least one region in other.
    fn subset_by_overlaps(&self, other: &RegionSet) -> RegionSet;

    /// Return a Vec<usize> with one entry per region in self, counting
    /// how many regions in other overlap it.
    fn count_overlaps(&self, other: &RegionSet) -> Vec<usize>;

    /// Return a Vec<bool> with one entry per region in self, true if
    /// any region in other overlaps it.
    fn any_overlaps(&self, other: &RegionSet) -> Vec<bool>;

    /// Return a Vec<Vec<usize>> with one entry per region in self,
    /// containing indices into other of all overlapping regions.
    fn find_overlaps(&self, other: &RegionSet) -> Vec<Vec<usize>>;
}

impl RegionSetOverlaps for RegionSet {
    fn subset_by_overlaps(&self, other: &RegionSet) -> RegionSet {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        let flags = index.any_query_overlaps(self);
        let kept: Vec<_> = self
            .regions
            .iter()
            .zip(flags)
            .filter_map(|(r, hit)| if hit { Some(r.clone()) } else { None })
            .collect();
        RegionSet::from(kept)
    }

    fn count_overlaps(&self, other: &RegionSet) -> Vec<usize> {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        index.count_query_overlaps(self)
    }

    fn any_overlaps(&self, other: &RegionSet) -> Vec<bool> {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        index.any_query_overlaps(self)
    }

    fn find_overlaps(&self, other: &RegionSet) -> Vec<Vec<usize>> {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        index.find_query_overlaps(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_core::models::{Region, RegionSet};

    fn make_region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    #[test]
    fn test_subset_by_overlaps() {
        let a = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
            make_region("chr1", 500, 600),
        ]);
        let b = RegionSet::from(vec![
            make_region("chr1", 150, 250), // overlaps first
            make_region("chr1", 550, 650), // overlaps third
        ]);
        let result = a.subset_by_overlaps(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0].start, 100);
        assert_eq!(result.regions[1].start, 500);
    }

    #[test]
    fn test_count_overlaps() {
        let a = RegionSet::from(vec![make_region("chr1", 100, 300)]);
        let b = RegionSet::from(vec![
            make_region("chr1", 150, 200),
            make_region("chr1", 250, 350),
            make_region("chr1", 500, 600),
        ]);
        let counts = a.count_overlaps(&b);
        assert_eq!(counts, vec![2]);
    }

    #[test]
    fn test_any_overlaps() {
        let a = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
        ]);
        let b = RegionSet::from(vec![make_region("chr1", 150, 250)]);
        let any = a.any_overlaps(&b);
        assert_eq!(any, vec![true, false]);
    }

    #[test]
    fn test_find_overlaps_returns_indices() {
        let a = RegionSet::from(vec![make_region("chr1", 100, 300)]);
        let b = RegionSet::from(vec![
            make_region("chr1", 50, 150),  // index 0, overlaps
            make_region("chr1", 200, 250), // index 1, overlaps
            make_region("chr1", 400, 500), // index 2, no overlap
        ]);
        let indices = a.find_overlaps(&b);
        assert_eq!(indices.len(), 1);
        let mut sorted = indices[0].clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1]);
    }

    #[test]
    fn test_multi_chrom() {
        let a = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr2", 100, 200),
        ]);
        let b = RegionSet::from(vec![make_region("chr1", 150, 250)]);
        let any = a.any_overlaps(&b);
        assert_eq!(any, vec![true, false]);
    }

    #[test]
    fn test_empty_other() {
        let a = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let b = RegionSet::from(vec![]);
        assert_eq!(a.count_overlaps(&b), vec![0]);
        assert_eq!(a.any_overlaps(&b), vec![false]);
        assert_eq!(a.find_overlaps(&b), vec![vec![] as Vec<usize>]);
        assert_eq!(a.subset_by_overlaps(&b).regions.len(), 0);
    }
}
