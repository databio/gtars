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
/// number of overlapping base pairs required. `None` or `Some(0)` or `Some(1)`
/// means any overlap at all (matching R's `minoverlap = 1L` default).
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

/// Compute actual overlap in base pairs between two regions.
#[inline]
fn overlap_bp(a_start: u32, a_end: u32, b_start: u32, b_end: u32) -> i64 {
    a_end.min(b_end) as i64 - a_start.max(b_start) as i64
}

impl RegionSetOverlaps for RegionSet {
    fn subset_by_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> RegionSet {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        let min_bp = min_overlap.unwrap_or(0);
        if min_bp <= 1 {
            // AIList hits guarantee >= 1bp overlap, so just check existence
            let flags = index.any_query_overlaps(self);
            let kept: Vec<_> = self
                .regions
                .iter()
                .zip(flags)
                .filter_map(|(r, hit)| if hit { Some(r.clone()) } else { None })
                .collect();
            RegionSet::from(kept)
        } else {
            let all_hits = index.find_query_overlaps(self);
            let kept: Vec<_> = self
                .regions
                .iter()
                .zip(all_hits.iter())
                .filter_map(|(r, hits)| {
                    let has_hit = hits.iter().any(|&idx| {
                        let b = &other.regions[idx];
                        overlap_bp(r.start, r.end, b.start, b.end) >= min_bp as i64
                    });
                    if has_hit { Some(r.clone()) } else { None }
                })
                .collect();
            RegionSet::from(kept)
        }
    }

    fn count_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> Vec<usize> {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        let min_bp = min_overlap.unwrap_or(0);
        if min_bp <= 1 {
            return index.count_query_overlaps(self);
        }
        // Iterate in-place to avoid materializing all hit indices
        self.regions
            .iter()
            .map(|region| match index.get_chr_overlapper(&region.chr) {
                Some(lapper) => lapper
                    .find_iter(region.start, region.end)
                    .filter(|iv| {
                        let b = &other.regions[iv.val];
                        overlap_bp(region.start, region.end, b.start, b.end) >= min_bp as i64
                    })
                    .count(),
                None => 0,
            })
            .collect()
    }

    fn any_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> Vec<bool> {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        let min_bp = min_overlap.unwrap_or(0);
        if min_bp <= 1 {
            return index.any_query_overlaps(self);
        }
        // Iterate in-place — short-circuits on first qualifying hit
        self.regions
            .iter()
            .map(|region| match index.get_chr_overlapper(&region.chr) {
                Some(lapper) => lapper
                    .find_iter(region.start, region.end)
                    .any(|iv| {
                        let b = &other.regions[iv.val];
                        overlap_bp(region.start, region.end, b.start, b.end) >= min_bp as i64
                    }),
                None => false,
            })
            .collect()
    }

    fn find_overlaps(&self, other: &RegionSet, min_overlap: Option<i32>) -> Vec<Vec<usize>> {
        let index = build_indexed_overlapper(other, OverlapperType::AIList);
        let min_bp = min_overlap.unwrap_or(0);
        if min_bp <= 1 {
            return index.find_query_overlaps(self);
        }
        let all_hits = index.find_query_overlaps(self);
        self.regions
            .iter()
            .zip(all_hits.into_iter())
            .map(|(r, hits)| {
                hits.into_iter()
                    .filter(|&idx| {
                        let b = &other.regions[idx];
                        overlap_bp(r.start, r.end, b.start, b.end) >= min_bp as i64
                    })
                    .collect()
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

    #[test]
    fn test_min_overlap_one_takes_fast_path() {
        // Some(1) should behave identically to None — both take the fast path
        // since AIList hits already guarantee >= 1bp overlap
        let a = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
        ]);
        let b = RegionSet::from(vec![
            make_region("chr1", 199, 300), // 1bp overlap with first
        ]);
        assert_eq!(a.count_overlaps(&b, Some(1)), vec![1, 0]);
        assert_eq!(a.any_overlaps(&b, Some(1)), vec![true, false]);
        assert_eq!(a.find_overlaps(&b, Some(1)), vec![vec![0], vec![]]);
        assert_eq!(a.subset_by_overlaps(&b, Some(1)).regions.len(), 1);
    }

    #[test]
    fn test_min_overlap_all_four_methods() {
        // Ensure min_overlap filtering works consistently across all methods
        // a: [100, 200) = 100bp; b: [180, 300) = overlaps by 20bp
        let a = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let b = RegionSet::from(vec![make_region("chr1", 180, 300)]);

        // min_overlap=20: exact boundary — should pass
        assert_eq!(a.count_overlaps(&b, Some(20)), vec![1]);
        assert_eq!(a.any_overlaps(&b, Some(20)), vec![true]);
        assert_eq!(a.find_overlaps(&b, Some(20)), vec![vec![0]]);
        assert_eq!(a.subset_by_overlaps(&b, Some(20)).regions.len(), 1);

        // min_overlap=21: one past boundary — should fail
        assert_eq!(a.count_overlaps(&b, Some(21)), vec![0]);
        assert_eq!(a.any_overlaps(&b, Some(21)), vec![false]);
        assert_eq!(a.find_overlaps(&b, Some(21)), vec![vec![] as Vec<usize>]);
        assert_eq!(a.subset_by_overlaps(&b, Some(21)).regions.len(), 0);
    }

    #[test]
    fn test_min_overlap_multi_chrom() {
        // Filtering should work across chromosomes
        let a = RegionSet::from(vec![
            make_region("chr1", 100, 200), // 100bp region
            make_region("chr2", 100, 200), // 100bp region
        ]);
        let b = RegionSet::from(vec![
            make_region("chr1", 150, 300), // 50bp overlap with chr1 region
            make_region("chr2", 190, 300), // 10bp overlap with chr2 region
        ]);

        // min_overlap=30: chr1 passes (50bp), chr2 fails (10bp)
        assert_eq!(a.count_overlaps(&b, Some(30)), vec![1, 0]);
        assert_eq!(a.any_overlaps(&b, Some(30)), vec![true, false]);
        assert_eq!(a.subset_by_overlaps(&b, Some(30)).regions.len(), 1);
        assert_eq!(a.subset_by_overlaps(&b, Some(30)).regions[0].chr, "chr1");
    }

    #[test]
    fn test_min_overlap_exceeds_all() {
        // min_overlap larger than any possible overlap — everything filtered
        let a = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
        ]);
        let b = RegionSet::from(vec![
            make_region("chr1", 150, 250),
            make_region("chr1", 350, 450),
        ]);
        // Each overlap is 50bp; require 100bp
        assert_eq!(a.count_overlaps(&b, Some(100)), vec![0, 0]);
        assert_eq!(a.any_overlaps(&b, Some(100)), vec![false, false]);
        assert_eq!(a.find_overlaps(&b, Some(100)), vec![vec![] as Vec<usize>, vec![] as Vec<usize>]);
        assert_eq!(a.subset_by_overlaps(&b, Some(100)).regions.len(), 0);
    }

    #[test]
    fn test_min_overlap_negative_same_as_none() {
        // Negative min_overlap should be treated as "any overlap"
        let a = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let b = RegionSet::from(vec![make_region("chr1", 199, 300)]);

        assert_eq!(a.count_overlaps(&b, Some(-5)), vec![1]);
        assert_eq!(a.any_overlaps(&b, Some(-5)), vec![true]);
    }

    #[test]
    fn test_count_overlaps_min_overlap_partial_filter() {
        // Multiple hits on one region, min_overlap filters some but not all
        let a = RegionSet::from(vec![make_region("chr1", 100, 300)]); // 200bp region
        let b = RegionSet::from(vec![
            make_region("chr1", 90, 110),  // 10bp overlap
            make_region("chr1", 150, 250), // 100bp overlap (fully contained)
            make_region("chr1", 290, 400), // 10bp overlap
        ]);

        // min_overlap=50: only the middle hit (100bp) passes
        assert_eq!(a.count_overlaps(&b, Some(50)), vec![1]);
        // min_overlap=5: all three pass
        assert_eq!(a.count_overlaps(&b, Some(5)), vec![3]);
    }
}
