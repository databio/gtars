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
///
/// All methods accept an optional `min_overlap` parameter specifying the minimum
/// number of overlapping base pairs required. `None` or `Some(0)` means any
/// overlap at all (matching R's `minoverlap = 0L` default).
pub trait RegionSetOverlaps {
    /// Return a new RegionSet containing only regions from self that overlap
    /// at least one region in other.
    fn subset_by_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> RegionSet;

    /// Return a Vec<usize> with one entry per region in self, counting
    /// how many regions in other overlap it.
    fn count_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> Vec<usize>;

    /// Return a Vec<bool> with one entry per region in self, true if
    /// any region in other overlaps it.
    fn any_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> Vec<bool>;

    /// Return a Vec<Vec<usize>> with one entry per region in self,
    /// containing indices into other of all overlapping regions.
    fn find_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> Vec<Vec<usize>>;
}

/// Check whether two regions overlap by at least `min_bp` base pairs.
#[inline]
fn meets_min_overlap(
    a_start: u32,
    a_end: u32,
    b_start: u32,
    b_end: u32,
    min_bp: i32,
) -> bool {
    if min_bp <= 0 {
        return true; // any overlap suffices
    }
    let overlap = a_end.min(b_end) as i64 - a_start.max(b_start) as i64;
    overlap >= min_bp as i64
}

impl RegionSetOverlaps for RegionSet {
    fn subset_by_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> RegionSet {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        let min_bp = min_overlap.unwrap_or(0);
        let mut kept = Vec::new();
        for region in &self.regions {
            if let Some(lapper) = index.get_chr_overlapper(&region.chr) {
                let has_hit = lapper
                    .find_iter(region.start, region.end)
                    .any(|iv| {
                        meets_min_overlap(
                            region.start,
                            region.end,
                            other.regions[iv.val].start,
                            other.regions[iv.val].end,
                            min_bp,
                        )
                    });
                if has_hit {
                    kept.push(region.clone());
                }
            }
        }
        RegionSet::from(kept)
    }

    fn count_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> Vec<usize> {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        let min_bp = min_overlap.unwrap_or(0);
        self.regions
            .iter()
            .map(|region| match index.get_chr_overlapper(&region.chr) {
                Some(lapper) => {
                    if min_bp <= 0 {
                        lapper.find_iter(region.start, region.end).count()
                    } else {
                        lapper
                            .find_iter(region.start, region.end)
                            .filter(|iv| {
                                meets_min_overlap(
                                    region.start,
                                    region.end,
                                    other.regions[iv.val].start,
                                    other.regions[iv.val].end,
                                    min_bp,
                                )
                            })
                            .count()
                    }
                }
                None => 0,
            })
            .collect()
    }

    fn any_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> Vec<bool> {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        let min_bp = min_overlap.unwrap_or(0);
        self.regions
            .iter()
            .map(|region| match index.get_chr_overlapper(&region.chr) {
                Some(lapper) => {
                    if min_bp <= 0 {
                        lapper.find_iter(region.start, region.end).next().is_some()
                    } else {
                        lapper
                            .find_iter(region.start, region.end)
                            .any(|iv| {
                                meets_min_overlap(
                                    region.start,
                                    region.end,
                                    other.regions[iv.val].start,
                                    other.regions[iv.val].end,
                                    min_bp,
                                )
                            })
                    }
                }
                None => false,
            })
            .collect()
    }

    fn find_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> Vec<Vec<usize>> {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        let min_bp = min_overlap.unwrap_or(0);
        self.regions
            .iter()
            .map(|region| match index.get_chr_overlapper(&region.chr) {
                Some(lapper) => {
                    if min_bp <= 0 {
                        lapper
                            .find_iter(region.start, region.end)
                            .map(|iv| iv.val)
                            .collect()
                    } else {
                        lapper
                            .find_iter(region.start, region.end)
                            .filter(|iv| {
                                meets_min_overlap(
                                    region.start,
                                    region.end,
                                    other.regions[iv.val].start,
                                    other.regions[iv.val].end,
                                    min_bp,
                                )
                            })
                            .map(|iv| iv.val)
                            .collect()
                    }
                }
                None => Vec::new(),
            })
            .collect()
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
        let result = a.subset_by_overlaps(&b, None);
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
        let counts = a.count_overlaps(&b, None);
        assert_eq!(counts, vec![2]);
    }

    #[test]
    fn test_any_overlaps() {
        let a = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
        ]);
        let b = RegionSet::from(vec![make_region("chr1", 150, 250)]);
        let any = a.any_overlaps(&b, None);
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
        let indices = a.find_overlaps(&b, None);
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
        let any = a.any_overlaps(&b, None);
        assert_eq!(any, vec![true, false]);
    }

    #[test]
    fn test_empty_other() {
        let a = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let b = RegionSet::from(vec![]);
        assert_eq!(a.count_overlaps(&b, None), vec![0]);
        assert_eq!(a.any_overlaps(&b, None), vec![false]);
        assert_eq!(a.find_overlaps(&b, None), vec![vec![] as Vec<usize>]);
        assert_eq!(a.subset_by_overlaps(&b, None).regions.len(), 0);
    }

    #[test]
    fn test_min_overlap_filters() {
        // a: [100, 200) = 100bp region
        // b: [190, 300) = overlaps a by 10bp
        let a = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let b = RegionSet::from(vec![make_region("chr1", 190, 300)]);

        // With no min_overlap, should find the overlap
        assert_eq!(a.count_overlaps(&b, None), vec![1]);
        assert_eq!(a.any_overlaps(&b, None), vec![true]);
        assert_eq!(a.find_overlaps(&b, None), vec![vec![0]]);
        assert_eq!(a.subset_by_overlaps(&b, None).regions.len(), 1);

        // With min_overlap=10, 10bp overlap passes
        assert_eq!(a.count_overlaps(&b, Some(10)), vec![1]);
        assert_eq!(a.any_overlaps(&b, Some(10)), vec![true]);

        // With min_overlap=11, 10bp overlap fails
        assert_eq!(a.count_overlaps(&b, Some(11)), vec![0]);
        assert_eq!(a.any_overlaps(&b, Some(11)), vec![false]);
        assert_eq!(a.find_overlaps(&b, Some(11)), vec![vec![] as Vec<usize>]);
        assert_eq!(a.subset_by_overlaps(&b, Some(11)).regions.len(), 0);
    }

    #[test]
    fn test_min_overlap_zero_same_as_none() {
        let a = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let b = RegionSet::from(vec![make_region("chr1", 199, 300)]);

        // 1bp overlap — should pass for both None and Some(0)
        assert_eq!(a.count_overlaps(&b, None), vec![1]);
        assert_eq!(a.count_overlaps(&b, Some(0)), vec![1]);
    }
}
