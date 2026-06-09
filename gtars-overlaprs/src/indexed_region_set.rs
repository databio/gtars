//! A RegionSet paired with a pre-built overlap index for efficient queries.
//!
//! This module provides [`IndexedRegionSet`], a wrapper type that pairs a [`RegionSet`]
//! with a pre-built [`MultiChromOverlapper`] index. This enables the "build once, query many"
//! pattern: build the index once, then run multiple overlap queries without rebuilding.
//!
//! # Examples
//!
//! ```rust
//! use gtars_overlaprs::IndexedRegionSet;
//! use gtars_core::models::{Region, RegionSet};
//!
//! // Create a reference set of genomic intervals
//! let reference = RegionSet::from(vec![
//!     Region { chr: "chr1".to_string(), start: 100, end: 200, rest: None },
//!     Region { chr: "chr1".to_string(), start: 300, end: 400, rest: None },
//!     Region { chr: "chr2".to_string(), start: 500, end: 600, rest: None },
//! ]);
//!
//! // Build an indexed version for efficient queries
//! let indexed = IndexedRegionSet::new(reference);
//!
//! // Query multiple times without rebuilding the index
//! let query1 = RegionSet::from(vec![
//!     Region { chr: "chr1".to_string(), start: 150, end: 250, rest: None },
//! ]);
//! let query2 = RegionSet::from(vec![
//!     Region { chr: "chr2".to_string(), start: 550, end: 650, rest: None },
//! ]);
//!
//! // All queries reuse the same index
//! let counts1 = indexed.count_overlaps(&query1, None);
//! let counts2 = indexed.count_overlaps(&query2, None);
//! ```

use std::ops::Deref;

use gtars_core::models::region_set::IntervalSetOps;
use gtars_core::models::{Region, RegionSet};

use crate::multi_chrom_overlapper::{build_indexed_overlapper, MultiChromOverlapper};
use crate::OverlapperType;

/// A RegionSet paired with a pre-built overlap index for efficient queries.
///
/// `IndexedRegionSet` solves the index-reuse problem: build the index once,
/// then run multiple overlap queries without rebuilding.
///
/// # When to use
///
/// Use `IndexedRegionSet` when you need to:
/// - Query the same reference set multiple times
/// - Perform multiple different operations (count, find, subset) on the same reference
/// - Optimize performance by avoiding repeated index construction
///
/// For one-off operations, the inherent set-algebra methods on [`RegionSet`]
/// (and the [`IntervalSetOps`] trait) are simpler since they build no index.
///
/// # Examples
///
/// ```rust
/// use gtars_overlaprs::IndexedRegionSet;
/// use gtars_core::models::{Region, RegionSet};
///
/// let reference = RegionSet::from(vec![
///     Region { chr: "chr1".to_string(), start: 100, end: 200, rest: None },
///     Region { chr: "chr1".to_string(), start: 300, end: 400, rest: None },
/// ]);
///
/// let indexed = IndexedRegionSet::new(reference);
///
/// // Build once, query many times
/// let q1 = RegionSet::from(vec![
///     Region { chr: "chr1".to_string(), start: 150, end: 350, rest: None },
/// ]);
///
/// let counts = indexed.count_overlaps(&q1, None);
/// let any = indexed.any_overlaps(&q1, None);
/// let indices = indexed.find_overlaps(&q1, None);
/// ```
pub struct IndexedRegionSet {
    /// The source RegionSet that was indexed. This is the *single* place the
    /// source regions live; the inner [`MultiChromOverlapper`] is source-free.
    source: RegionSet,
    /// The pre-built source-free overlap index.
    index: MultiChromOverlapper<u32, ()>,
}

impl IndexedRegionSet {
    /// Create a new IndexedRegionSet from a RegionSet.
    ///
    /// Builds an overlap index using the default AIList algorithm. Use
    /// [`with_overlapper_type`](Self::with_overlapper_type) to select a different algorithm.
    /// The index stores indices back into the original RegionSet, so all Region data is preserved.
    ///
    /// # Arguments
    ///
    /// * `regions` - The RegionSet to index
    ///
    /// # Examples
    ///
    /// ```rust
    /// use gtars_overlaprs::IndexedRegionSet;
    /// use gtars_core::models::{Region, RegionSet};
    ///
    /// let rs = RegionSet::from(vec![
    ///     Region { chr: "chr1".to_string(), start: 100, end: 200, rest: None },
    /// ]);
    /// let indexed = IndexedRegionSet::new(rs);
    /// ```
    pub fn new(regions: RegionSet) -> Self {
        Self::with_overlapper_type(regions, OverlapperType::AIList)
    }

    /// Create with a specific overlapper type (AIList or Bits).
    ///
    /// # Arguments
    ///
    /// * `regions` - The RegionSet to index
    /// * `overlapper_type` - The type of overlap data structure to use
    ///
    /// # Examples
    ///
    /// ```rust
    /// use gtars_overlaprs::{IndexedRegionSet, OverlapperType};
    /// use gtars_core::models::{Region, RegionSet};
    ///
    /// let rs = RegionSet::from(vec![
    ///     Region { chr: "chr1".to_string(), start: 100, end: 200, rest: None },
    /// ]);
    /// let indexed = IndexedRegionSet::with_overlapper_type(rs, OverlapperType::Bits);
    /// ```
    pub fn with_overlapper_type(regions: RegionSet, overlapper_type: OverlapperType) -> Self {
        let index = build_indexed_overlapper(&regions, overlapper_type);
        Self {
            source: regions,
            index,
        }
    }

    /// Build a lookup from a region's `(chr, start, end)` coordinates to the
    /// list of source indices with those coordinates. Used to translate the
    /// source-free index's coordinate hits back into source indices for the
    /// index-returning query methods.
    fn coord_lookup(&self) -> std::collections::HashMap<(&str, u32, u32), Vec<usize>> {
        let mut map: std::collections::HashMap<(&str, u32, u32), Vec<usize>> =
            std::collections::HashMap::new();
        for (i, r) in self.source.regions.iter().enumerate() {
            map.entry((r.chr.as_str(), r.start, r.end)).or_default().push(i);
        }
        map
    }

    /// Access the underlying RegionSet.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use gtars_overlaprs::IndexedRegionSet;
    /// use gtars_core::models::{Region, RegionSet};
    ///
    /// let rs = RegionSet::from(vec![
    ///     Region { chr: "chr1".to_string(), start: 100, end: 200, rest: None },
    /// ]);
    /// let indexed = IndexedRegionSet::new(rs);
    /// assert_eq!(indexed.regions().len(), 1);
    /// ```
    pub fn regions(&self) -> &RegionSet {
        &self.source
    }

    /// Consume and return the underlying RegionSet.
    ///
    /// This is useful when you're done with the index and want to reclaim
    /// the original RegionSet without cloning.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use gtars_overlaprs::IndexedRegionSet;
    /// use gtars_core::models::{Region, RegionSet};
    ///
    /// let rs = RegionSet::from(vec![
    ///     Region { chr: "chr1".to_string(), start: 100, end: 200, rest: None },
    /// ]);
    /// let indexed = IndexedRegionSet::new(rs);
    /// let rs_back = indexed.into_regions();
    /// ```
    pub fn into_regions(self) -> RegionSet {
        self.source
    }

    // ---------------------------------------------------------------
    // Query methods (unified min_overlap-carrying family)
    // ---------------------------------------------------------------

    /// Returns regions from self that overlap any region in `query`.
    ///
    /// Returns the *source* regions (with their original `rest` metadata),
    /// deduplicated and returned in source index order.
    pub fn intersect_all(&self, query: &RegionSet) -> RegionSet {
        let mut hits: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        let hit_lists = self.find_overlaps(query, None);
        for list in hit_lists {
            for idx in list {
                hits.insert(idx);
            }
        }
        let kept: Vec<Region> = hits
            .into_iter()
            .map(|idx| self.source.regions[idx].clone())
            .collect();
        RegionSet::from(kept)
    }

    /// Return only the source regions that overlap at least one region in
    /// `query`, honoring an optional minimum overlap in base pairs.
    pub fn subset_by_overlaps(&self, query: &RegionSet, min_overlap: Option<i32>) -> RegionSet {
        let mut hits: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        for list in self.find_overlaps(query, min_overlap) {
            for idx in list {
                hits.insert(idx);
            }
        }
        let kept: Vec<Region> = hits
            .into_iter()
            .map(|idx| self.source.regions[idx].clone())
            .collect();
        RegionSet::from(kept)
    }

    /// For each region in `query`, count how many indexed regions overlap it,
    /// honoring an optional minimum overlap in base pairs.
    pub fn count_overlaps(&self, query: &RegionSet, min_overlap: Option<i32>) -> Vec<usize> {
        self.index.count_overlaps(query, min_overlap)
    }

    /// For each region in `query`, whether any indexed region overlaps it,
    /// honoring an optional minimum overlap in base pairs.
    pub fn any_overlaps(&self, query: &RegionSet, min_overlap: Option<i32>) -> Vec<bool> {
        self.index.any_overlaps(query, min_overlap)
    }

    /// For each region in `query`, return the indices of overlapping regions in
    /// the source RegionSet, honoring an optional minimum overlap in base pairs.
    pub fn find_overlaps(&self, query: &RegionSet, min_overlap: Option<i32>) -> Vec<Vec<usize>> {
        let lookup = self.coord_lookup();
        let region_hits = self.index.find_overlaps_regions(query, min_overlap);
        region_hits
            .into_iter()
            .map(|regions| {
                let mut idxs: Vec<usize> = Vec::new();
                for r in regions {
                    if let Some(src_idxs) = lookup.get(&(r.chr.as_str(), r.start, r.end)) {
                        idxs.extend_from_slice(src_idxs);
                    }
                }
                idxs.sort_unstable();
                idxs.dedup();
                idxs
            })
            .collect()
    }

    /// Find the nearest source region for each region in `other`.
    ///
    /// Returns `(self_index, other_index, distance)` tuples where `self_index`
    /// references the source RegionSet order.
    pub fn closest(&self, other: &RegionSet) -> Vec<(usize, usize, i64)> {
        self.source.closest(other)
    }

    /// Union of the source regions with `other`.
    pub fn union(&self, other: &RegionSet) -> RegionSet {
        self.source.union(other)
    }

    /// Cluster the source regions within `max_gap`.
    pub fn cluster(&self, max_gap: u32) -> Vec<u32> {
        self.source.cluster(max_gap)
    }
}

// ---------------------------------------------------------------
// IntervalSetOps facade: index-native set algebra via the inner MCO
// ---------------------------------------------------------------

impl IntervalSetOps for IndexedRegionSet {
    fn setdiff(&self, other: &RegionSet) -> RegionSet {
        self.index.setdiff(other)
    }

    fn intersect(&self, other: &RegionSet) -> RegionSet {
        self.index.intersect(other)
    }

    fn jaccard(&self, other: &RegionSet) -> f64 {
        self.index.jaccard(other)
    }

    fn coverage(&self, other: &RegionSet) -> f64 {
        self.index.coverage(other)
    }

    fn overlap_coefficient(&self, other: &RegionSet) -> f64 {
        self.index.overlap_coefficient(other)
    }

    fn closest(&self, other: &RegionSet) -> Vec<(usize, usize, i64)> {
        self.source.closest(other)
    }
}

// ---------------------------------------------------------------
// Deref to RegionSet for convenience
// ---------------------------------------------------------------

impl Deref for IndexedRegionSet {
    type Target = RegionSet;

    /// Dereference to the underlying RegionSet.
    ///
    /// This allows direct access to RegionSet methods like `iter()`, `iter_chroms()`, etc.
    fn deref(&self) -> &Self::Target {
        &self.source
    }
}

// ---------------------------------------------------------------
// From/Into conversions
// ---------------------------------------------------------------

impl From<RegionSet> for IndexedRegionSet {
    fn from(regions: RegionSet) -> Self {
        Self::new(regions)
    }
}

impl From<IndexedRegionSet> for RegionSet {
    fn from(indexed: IndexedRegionSet) -> Self {
        indexed.into_regions()
    }
}

// ---------------------------------------------------------------
// Tests
// ---------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_core::models::Region;

    fn make_region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    #[test]
    fn test_new_and_regions() {
        let rs = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
        ]);
        let indexed = IndexedRegionSet::new(rs.clone());
        assert_eq!(indexed.regions().len(), rs.len());
    }

    #[test]
    fn test_with_overlapper_type() {
        let rs = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let indexed_ailist = IndexedRegionSet::new(rs.clone());
        let indexed_bits = IndexedRegionSet::with_overlapper_type(rs, OverlapperType::Bits);

        // Both should work identically for queries
        let query = RegionSet::from(vec![make_region("chr1", 150, 250)]);
        assert_eq!(
            indexed_ailist.count_overlaps(&query, None),
            indexed_bits.count_overlaps(&query, None)
        );
    }

    #[test]
    fn test_into_regions() {
        let rs = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let indexed = IndexedRegionSet::new(rs.clone());
        let rs_back = indexed.into_regions();
        assert_eq!(rs_back.len(), rs.len());
    }

    #[test]
    fn test_intersect_all() {
        let reference = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
            make_region("chr2", 500, 600),
        ]);
        let indexed = IndexedRegionSet::new(reference);

        let query = RegionSet::from(vec![
            make_region("chr1", 150, 250), // overlaps first
            make_region("chr2", 550, 650), // overlaps third
        ]);

        let result = indexed.intersect_all(&query);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0].start, 100);
        assert_eq!(result.regions[0].end, 200);
        assert_eq!(result.regions[1].start, 500);
        assert_eq!(result.regions[1].end, 600);
    }

    #[test]
    fn test_count_overlaps() {
        let reference = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 150, 250),
            make_region("chr1", 300, 400),
        ]);
        let indexed = IndexedRegionSet::new(reference);

        let query = RegionSet::from(vec![make_region("chr1", 180, 220)]);
        let counts = indexed.count_overlaps(&query, None);
        assert_eq!(counts, vec![2]); // First two reference regions overlap
    }

    #[test]
    fn test_any_overlaps() {
        let reference = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let indexed = IndexedRegionSet::new(reference);

        let query = RegionSet::from(vec![
            make_region("chr1", 150, 250), // overlaps
            make_region("chr1", 300, 400), // no overlap
        ]);

        let any = indexed.any_overlaps(&query, None);
        assert_eq!(any, vec![true, false]);
    }

    #[test]
    fn test_find_overlaps() {
        let reference = RegionSet::from(vec![
            make_region("chr1", 100, 200), // index 0
            make_region("chr1", 300, 400), // index 1
        ]);
        let indexed = IndexedRegionSet::new(reference);

        let query = RegionSet::from(vec![make_region("chr1", 150, 350)]);

        let indices = indexed.find_overlaps(&query, None);
        assert_eq!(indices.len(), 1);
        let mut sorted = indices[0].clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1]);
    }

    #[test]
    fn test_deref_to_regionset() {
        let rs = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr2", 300, 400),
        ]);
        let indexed = IndexedRegionSet::new(rs);

        // Deref allows calling RegionSet methods directly
        assert_eq!(indexed.len(), 2);
        assert!(!indexed.is_empty());

        // iter_chroms() is a RegionSet method
        let chroms: Vec<_> = indexed.iter_chroms().collect();
        assert_eq!(chroms.len(), 2);
    }

    #[test]
    fn test_from_into_conversions() {
        let rs = RegionSet::from(vec![make_region("chr1", 100, 200)]);

        // From<RegionSet>
        let indexed: IndexedRegionSet = rs.clone().into();
        assert_eq!(indexed.len(), 1);

        // From<IndexedRegionSet>
        let rs_back: RegionSet = indexed.into();
        assert_eq!(rs_back.len(), 1);
    }

    #[test]
    fn test_empty_indexed_regionset() {
        let rs = RegionSet::from(vec![]);
        let indexed = IndexedRegionSet::new(rs);

        assert!(indexed.is_empty());
        assert_eq!(indexed.len(), 0);

        let query = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        assert_eq!(indexed.count_overlaps(&query, None), vec![0]);
        assert_eq!(indexed.any_overlaps(&query, None), vec![false]);
        assert_eq!(indexed.find_overlaps(&query, None), vec![vec![] as Vec<usize>]);
        assert_eq!(indexed.intersect_all(&query).regions.len(), 0);
    }

    #[test]
    fn test_empty_query() {
        let reference = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let indexed = IndexedRegionSet::new(reference);

        let query = RegionSet::from(vec![]);
        assert_eq!(indexed.count_overlaps(&query, None), Vec::<usize>::new());
        assert_eq!(indexed.any_overlaps(&query, None), Vec::<bool>::new());
        assert_eq!(indexed.find_overlaps(&query, None), Vec::<Vec<usize>>::new());
        assert_eq!(indexed.intersect_all(&query).regions.len(), 0);
    }

    #[test]
    fn test_multi_chrom() {
        let reference = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr2", 100, 200),
            make_region("chr3", 100, 200),
        ]);
        let indexed = IndexedRegionSet::new(reference);

        let query = RegionSet::from(vec![
            make_region("chr1", 150, 250),
            make_region("chr2", 150, 250),
            make_region("chr4", 150, 250), // nonexistent chrom
        ]);

        let counts = indexed.count_overlaps(&query, None);
        assert_eq!(counts, vec![1, 1, 0]);

        let any = indexed.any_overlaps(&query, None);
        assert_eq!(any, vec![true, true, false]);
    }

}
