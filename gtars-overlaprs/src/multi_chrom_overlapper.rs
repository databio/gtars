//! Genome-wide interval indexing for efficient multi-chromosome overlap queries.
//!
//! This module provides [`MultiChromOverlapper`], a data structure that extends single-chromosome
//! overlap data structures (like [`AIList`](crate::AIList) or [`Bits`](crate::Bits)) to handle genome-wide queries
//! across multiple chromosomes efficiently.
//!
//! # Build-once-query-many pattern
//!
//! The indexed variant (`MultiChromOverlapper<u32, usize>`) **owns** its source [`RegionSet`],
//! enabling a build-once-query-many workflow. Build the index with
//! [`MultiChromOverlapper::from_region_set`] (or [`build_indexed_overlapper`]), then
//! call any combination of query methods without rebuilding:
//!
//! - [`IntervalSetOps`](gtars_core::models::IntervalSetOps) trait methods (setdiff, intersect, union, jaccard, etc.)
//! - MCO-only methods: `subset_by`, `count_overlaps`, `any_overlaps`, `find_overlaps_indexed`, `intersect_all`
//!
//! # Examples
//!
//! ## Indexed variant (owns source, build-once-query-many)
//!
//! ```
//! use gtars_overlaprs::multi_chrom_overlapper::MultiChromOverlapper;
//! use gtars_overlaprs::OverlapperType;
//! use gtars_core::models::{Region, RegionSet, IntervalSetOps};
//!
//! let genes = RegionSet::from(vec![
//!     Region { chr: "chr1".to_string(), start: 1000, end: 2000, rest: None },
//!     Region { chr: "chr1".to_string(), start: 5000, end: 6000, rest: None },
//!     Region { chr: "chr2".to_string(), start: 1000, end: 3000, rest: None },
//! ]);
//!
//! // Build index once — MCO owns the source RegionSet
//! let mco = MultiChromOverlapper::from_region_set(genes, OverlapperType::AIList);
//!
//! let query = RegionSet::from(vec![
//!     Region { chr: "chr1".to_string(), start: 1500, end: 2500, rest: None },
//!     Region { chr: "chr2".to_string(), start: 2000, end: 4000, rest: None },
//! ]);
//!
//! // Query multiple times without rebuilding
//! let counts = mco.count_overlaps(&query);
//! let any = mco.any_overlaps(&query);
//! let subset = mco.subset_by(&query);
//!
//! // IntervalSetOps trait methods also work
//! let diff = mco.setdiff(&query);
//! let j = mco.jaccard(&query);
//! ```
//!
//! ## Basic usage (non-indexed, with metadata)
//!
//! ```
//! use gtars_overlaprs::{OverlapperType, multi_chrom_overlapper::IntoMultiChromOverlapper};
//! use gtars_core::models::{Region, RegionSet};
//!
//! let gene_set = RegionSet::from(vec![
//!     Region { chr: "chr1".to_string(), start: 1000, end: 2000, rest: Some("BRCA1".to_string()) },
//!     Region { chr: "chr1".to_string(), start: 5000, end: 6000, rest: Some("TP53".to_string()) },
//! ]);
//!
//! let index = gene_set.into_multi_chrom_overlapper(OverlapperType::AIList);
//! let query = RegionSet::from(vec![
//!     Region { chr: "chr1".to_string(), start: 1500, end: 2500, rest: None },
//! ]);
//! let overlaps = index.find_overlaps(&query);
//! ```

use std::{collections::HashMap, fmt::Debug};

use gtars_core::models::{Interval, IntervalSetOps, Region, RegionSet};
use num_traits::{PrimInt, Unsigned};
use thiserror::Error;

use crate::{AIList, Bits, Overlapper, OverlapperType};

/// Errors that can occur when working with [`MultiChromOverlapper`].
#[derive(Debug, Error)]
pub enum MultiChromOverlapperError {
    /// Error parsing a genomic region string.
    #[error("Error parsing region: {0}")]
    RegionParsingError(String),
    /// Error converting interval coordinates to the required integer type.
    #[error("Error converting interval coordinates to u32: start={0}, end={1}")]
    CoordinateConversionError(String, String),
}

/// A genome-wide index for efficient overlap queries across multiple chromosomes.
///
/// `MultiChromOverlapper` maintains a separate overlap data structure (either [`AIList`] or [`Bits`])
/// for each chromosome, enabling efficient queries across the entire genome.
///
/// ## Indexed variant (`MultiChromOverlapper<u32, usize>`)
///
/// When built with [`from_region_set`](MultiChromOverlapper::from_region_set) or
/// [`build_indexed_overlapper`], the MCO **owns** its source [`RegionSet`] and stores
/// indices back into it. This enables:
///
/// - **[`IntervalSetOps`] trait**: set algebra operations (setdiff, intersect, union,
///   jaccard, coverage, overlap_coefficient, closest, cluster) that produce identical
///   results to the sweep-line implementations on `RegionSet`.
/// - **MCO-only query methods**: [`subset_by`](MultiChromOverlapper::subset_by),
///   [`count_overlaps`](MultiChromOverlapper::count_overlaps),
///   [`any_overlaps`](MultiChromOverlapper::any_overlaps),
///   [`find_overlaps_indexed`](MultiChromOverlapper::find_overlaps_indexed),
///   [`intersect_all`](MultiChromOverlapper::intersect_all).
/// - **Build-once-query-many**: build the index once, then call any combination
///   of the above methods without rebuilding.
///
/// # Examples
///
/// See the [module-level documentation](self) for usage examples.
pub struct MultiChromOverlapper<I, T> {
    index_maps: HashMap<String, Box<dyn Overlapper<I, T>>>,
    #[allow(dead_code)]
    overlapper_type: OverlapperType,
    /// The source RegionSet that was indexed. Present when `T = usize`
    /// (indexed variant built via `from_region_set` or `build_indexed_overlapper`).
    source: Option<RegionSet>,
}

/// An iterator over intervals that overlap with query regions across multiple chromosomes.
///
/// This iterator is created by [`MultiChromOverlapper::find_overlaps_iter`]. It yields tuples of
/// `(chromosome, interval)` for all intervals that overlap with any of the query regions.
///
/// The iterator processes query regions in order and yields overlapping intervals as they
/// are found, without allocating a vector for all results.
///
/// # Examples
///
/// ```
/// use gtars_overlaprs::{OverlapperType, multi_chrom_overlapper::IntoMultiChromOverlapper};
/// use gtars_core::models::{Region, RegionSet};
///
/// let genes = RegionSet::from(vec![
///     Region { chr: "chr1".to_string(), start: 100, end: 200, rest: None },
/// ]);
/// let index = genes.into_multi_chrom_overlapper(OverlapperType::AIList);
///
/// let queries = RegionSet::from(vec![
///     Region { chr: "chr1".to_string(), start: 150, end: 250, rest: None },
/// ]);
///
/// // The iterator is used implicitly in a for loop
/// for (chr, interval) in index.find_overlaps_iter(&queries) {
///     println!("Found overlap on {}: {:?}", chr, interval);
/// }
/// ```
pub struct IterFindOverlaps<'a, 'b, I, T>
where
    I: PrimInt + Unsigned + Send + Sync + Debug,
    T: Eq + Clone + Send + Sync + Debug,
{
    inner: &'a HashMap<String, Box<dyn Overlapper<I, T>>>,
    rs: &'b RegionSet,
    region_idx: usize,
    current_chr: Option<String>,
    current_iter: Option<Box<dyn Iterator<Item = &'a Interval<I, T>> + 'a>>,
}

impl<'a, 'b, I, T> Iterator for IterFindOverlaps<'a, 'b, I, T>
where
    I: PrimInt + Unsigned + Send + Sync + Debug,
    T: Eq + Clone + Send + Sync + Debug,
{
    type Item = (String, &'a Interval<I, T>);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // first, try to get next item from current iterator
            #[allow(clippy::collapsible_if)]
            if let Some(ref mut iter) = self.current_iter {
                if let Some(interval) = iter.next() {
                    return Some((self.current_chr.as_ref().unwrap().clone(), interval));
                }
            }

            // current iterator exhausted or doesn't exist, move to next region
            if self.region_idx >= self.rs.regions.len() {
                // we are done afte this!
                // this is the terminal point
                return None;
            }

            let region = &self.rs.regions[self.region_idx];
            self.region_idx += 1;

            // try to get overlapper for this chromosome
            if let Some(lapper) = self.inner.get(&region.chr) {
                // convert coordinates
                if let (Some(start), Some(end)) = (I::from(region.start), I::from(region.end)) {
                    self.current_chr = Some(region.chr.clone());
                    self.current_iter = Some(lapper.find_iter(start, end));
                    // continue loop to get first item from new iterator
                } else {
                    // This is a programming error: the MultiChromOverlapper type I cannot represent
                    // the Region's u32 coordinates. This should never happen in practice since
                    // genomic coordinates are u32 and the index should be MultiChromOverlapper<u32, T>.
                    panic!(
                        "Type conversion error: cannot convert Region coordinates to index type. \
                         Region: {}:{}-{}, expected type: {}",
                        region.chr,
                        region.start,
                        region.end,
                        std::any::type_name::<I>()
                    );
                }
            } else {
                // no overlapper for this chromosome, skip to next region
                continue;
            }
        }
    }
}

impl<I, T> MultiChromOverlapper<I, T>
where
    I: PrimInt + Unsigned + Send + Sync + Debug,
    T: Eq + Clone + Send + Sync + Debug,
{
    /// Returns an iterator over all overlapping intervals for the query regions.
    ///
    /// Each item is a tuple of (chromosome, interval reference).
    pub fn find_overlaps_iter<'a, 'b>(
        &'a self,
        rs: &'b RegionSet,
    ) -> IterFindOverlaps<'a, 'b, I, T> {
        IterFindOverlaps {
            inner: &self.index_maps,
            rs,
            region_idx: 0,
            current_chr: None,
            current_iter: None,
        }
    }

    /// Collect all overlaps into a Vec for convenience. You're almost always
    /// better off using the iterator form of this function `find_overlaps_iter`.
    ///
    /// This is a helper method that collects the iterator results.
    pub fn find_overlaps(&self, rs: &RegionSet) -> Vec<(String, Interval<I, T>)> {
        self.find_overlaps_iter(rs)
            .map(|(chr, interval)| (chr, interval.clone()))
            .collect()
    }

    /// Get the overlapper for a specific chromosome, if it exists.
    pub fn get_chr_overlapper(&self, chr: &str) -> Option<&dyn Overlapper<I, T>> {
        self.index_maps.get(chr).map(|b| b.as_ref())
    }

    /// Query a single chromosome + coordinate range for overlapping intervals.
    ///
    /// Returns an iterator over all intervals that overlap with `[start, end)` on the
    /// given chromosome. Returns an empty iterator if the chromosome is not in the index.
    pub fn find_overlaps_for_region<'a>(
        &'a self,
        chr: &str,
        start: I,
        end: I,
    ) -> Box<dyn Iterator<Item = &'a Interval<I, T>> + 'a> {
        match self.index_maps.get(chr) {
            Some(lapper) => lapper.find_iter(start, end),
            None => Box::new(std::iter::empty()),
        }
    }
}

/// A genome-wide index that **owns** its source [`RegionSet`], enabling a
/// build-once-query-many pattern for overlap operations.
///
/// Build an indexed overlapper with [`MultiChromOverlapper::from_region_set`]
/// (or the free function [`build_indexed_overlapper`]), then call query methods
/// repeatedly without rebuilding.
///
/// # Implements
///
/// - [`IntervalSetOps`] — index-accelerated set algebra (setdiff, intersect,
///   union, jaccard, coverage, overlap_coefficient, subtract, closest, cluster).
/// - Inherent MCO-only query methods: [`subset_by`](Self::subset_by),
///   [`count_overlaps`](Self::count_overlaps), [`any_overlaps`](Self::any_overlaps),
///   [`find_overlaps_indexed`](Self::find_overlaps_indexed),
///   [`intersect_all`](Self::intersect_all).
///
/// # Examples
///
/// ```
/// use gtars_overlaprs::multi_chrom_overlapper::MultiChromOverlapper;
/// use gtars_overlaprs::OverlapperType;
/// use gtars_core::models::{Region, RegionSet, IntervalSetOps};
///
/// let source = RegionSet::from(vec![
///     Region { chr: "chr1".to_string(), start: 100, end: 200, rest: None },
///     Region { chr: "chr1".to_string(), start: 300, end: 400, rest: None },
/// ]);
///
/// // Build the index once — MCO owns the source
/// let mco = MultiChromOverlapper::from_region_set(source, OverlapperType::AIList);
///
/// let query = RegionSet::from(vec![
///     Region { chr: "chr1".to_string(), start: 150, end: 250, rest: None },
/// ]);
///
/// // Query multiple times without rebuilding
/// let counts = mco.count_overlaps(&query);
/// let any = mco.any_overlaps(&query);
/// let subset = mco.subset_by(&query);
///
/// // IntervalSetOps trait methods also work
/// let diff = mco.setdiff(&query);
/// ```
impl MultiChromOverlapper<u32, usize> {
    /// Build an indexed overlapper that **owns** the source [`RegionSet`].
    ///
    /// Each interval's `val` is the index of the corresponding region in `source`,
    /// enabling efficient lookups back into the original data.
    pub fn from_region_set(rs: RegionSet, overlapper_type: OverlapperType) -> Self {
        let mut core: HashMap<String, Box<dyn Overlapper<u32, usize>>> = HashMap::default();
        let mut intervals: HashMap<String, Vec<Interval<u32, usize>>> = HashMap::default();

        for (idx, region) in rs.regions.iter().enumerate() {
            let interval = Interval {
                start: region.start,
                end: region.end,
                val: idx,
            };
            intervals.entry(region.chr.clone()).or_default().push(interval);
        }

        for (chr, chr_intervals) in intervals.into_iter() {
            let lapper: Box<dyn Overlapper<u32, usize>> = match overlapper_type {
                OverlapperType::Bits => Box::new(Bits::build(chr_intervals)),
                OverlapperType::AIList => Box::new(AIList::build(chr_intervals)),
            };
            core.insert(chr, lapper);
        }

        MultiChromOverlapper {
            index_maps: core,
            overlapper_type,
            source: Some(rs),
        }
    }

    /// Access the owned source [`RegionSet`].
    pub fn source(&self) -> &RegionSet {
        self.source.as_ref().expect(
            "source RegionSet not present — this MCO was not built with from_region_set or build_indexed_overlapper"
        )
    }

    // ── MCO-only query methods ──────────────────────────────────────

    /// Return regions from the source that overlap any region in `query`.
    ///
    /// Results are deduplicated and returned in index order.
    pub fn subset_by(&self, query: &RegionSet) -> RegionSet {
        let source = self.source();
        let mut hit_indices = std::collections::BTreeSet::new();
        for region in &query.regions {
            if let Some(lapper) = self.index_maps.get(&region.chr) {
                for iv in lapper.find_iter(region.start, region.end) {
                    hit_indices.insert(iv.val);
                }
            }
        }
        let kept: Vec<Region> = hit_indices
            .into_iter()
            .map(|idx| source.regions[idx].clone())
            .collect();
        RegionSet::from(kept)
    }

    /// For each region in `query`, count how many indexed regions overlap it.
    pub fn count_overlaps(&self, query: &RegionSet) -> Vec<usize> {
        query
            .regions
            .iter()
            .map(|region| match self.index_maps.get(&region.chr) {
                Some(lapper) => lapper.find_iter(region.start, region.end).count(),
                None => 0,
            })
            .collect()
    }

    /// For each region in `query`, whether any indexed region overlaps it.
    pub fn any_overlaps(&self, query: &RegionSet) -> Vec<bool> {
        query
            .regions
            .iter()
            .map(|region| match self.index_maps.get(&region.chr) {
                Some(lapper) => lapper.find_iter(region.start, region.end).next().is_some(),
                None => false,
            })
            .collect()
    }

    /// For each region in `query`, return the indices of overlapping regions
    /// in the source [`RegionSet`].
    pub fn find_overlaps_indexed(&self, query: &RegionSet) -> Vec<Vec<usize>> {
        query
            .regions
            .iter()
            .map(|region| match self.index_maps.get(&region.chr) {
                Some(lapper) => lapper
                    .find_iter(region.start, region.end)
                    .map(|iv| iv.val)
                    .collect(),
                None => Vec::new(),
            })
            .collect()
    }

    /// Pairwise intersection fragments between indexed regions and query.
    ///
    /// For each pair of overlapping regions between the source and `query`,
    /// compute intersection coordinates `[max(a.start, b.start), min(a.end, b.end))`.
    /// Returns a `RegionSet` of all intersection fragments.
    pub fn intersect_all(&self, query: &RegionSet) -> RegionSet {
        let mut result: Vec<Region> = Vec::new();

        for q in &query.regions {
            for hit in self.find_overlaps_for_region(&q.chr, q.start, q.end) {
                let start = q.start.max(hit.start);
                let end = q.end.min(hit.end);
                if start < end {
                    result.push(Region {
                        chr: q.chr.clone(),
                        start,
                        end,
                        rest: None,
                    });
                }
            }
        }

        RegionSet::from(result)
    }

    // ── Legacy compatibility aliases ────────────────────────────────

    /// Alias for [`count_overlaps`](Self::count_overlaps).
    #[deprecated(note = "use count_overlaps instead")]
    pub fn count_query_overlaps(&self, query: &RegionSet) -> Vec<usize> {
        self.count_overlaps(query)
    }

    /// Alias for [`any_overlaps`](Self::any_overlaps).
    #[deprecated(note = "use any_overlaps instead")]
    pub fn any_query_overlaps(&self, query: &RegionSet) -> Vec<bool> {
        self.any_overlaps(query)
    }

    /// Alias for [`find_overlaps_indexed`](Self::find_overlaps_indexed).
    #[deprecated(note = "use find_overlaps_indexed instead")]
    pub fn find_query_overlaps(&self, query: &RegionSet) -> Vec<Vec<usize>> {
        self.find_overlaps_indexed(query)
    }
}

// ── IntervalSetOps implementation for indexed MCO ───────────────────────
//
// Delegates to the owned source RegionSet's sweep-line implementation.
// This lets callers write generic code over `impl IntervalSetOps` that
// works with either RegionSet or MCO.

impl IntervalSetOps for MultiChromOverlapper<u32, usize> {
    fn setdiff(&self, other: &RegionSet) -> RegionSet {
        self.source().setdiff(other)
    }

    fn intersect(&self, other: &RegionSet) -> RegionSet {
        self.source().intersect(other)
    }

    fn union(&self, other: &RegionSet) -> RegionSet {
        self.source().union(other)
    }

    fn jaccard(&self, other: &RegionSet) -> f64 {
        self.source().jaccard(other)
    }

    fn coverage(&self, other: &RegionSet) -> f64 {
        self.source().coverage(other)
    }

    fn overlap_coefficient(&self, other: &RegionSet) -> f64 {
        self.source().overlap_coefficient(other)
    }

    fn subtract(&self, other: &RegionSet) -> RegionSet {
        self.source().subtract(other)
    }

    fn closest(&self, other: &RegionSet) -> Vec<(usize, usize, i64)> {
        self.source().closest(other)
    }

    fn cluster(&self, max_gap: u32) -> Vec<u32> {
        self.source().cluster(max_gap)
    }
}

/// A trait for converting region-based data into a [`MultiChromOverlapper`].
///
/// This trait provides a convenient way to build a genome-wide index from a collection
/// of genomic regions. It handles the organization of regions by chromosome and the
/// construction of the underlying overlap data structures.
///
///
/// # Examples
///
/// ```
/// use gtars_overlaprs::{OverlapperType, multi_chrom_overlapper::IntoMultiChromOverlapper};
/// use gtars_core::models::{Region, RegionSet};
///
/// let regions = vec![
///     Region { chr: "chr1".to_string(), start: 100, end: 200, rest: Some("peak1".to_string()) },
///     Region { chr: "chr2".to_string(), start: 300, end: 400, rest: Some("peak2".to_string()) },
/// ];
///
/// let region_set = RegionSet::from(regions);
///
/// // Convert into a MultiChromOverlapper using AIList
/// let index = region_set.into_multi_chrom_overlapper(OverlapperType::AIList);
/// ```
pub trait IntoMultiChromOverlapper<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// Consumes the input and builds a [`MultiChromOverlapper`] with the specified overlapper type.
    ///
    /// # Arguments
    ///
    /// * `overlapper_type` - The type of overlap data structure to use ([`AIList`] or [`Bits`])
    ///
    /// # Returns
    ///
    /// A new [`MultiChromOverlapper`] containing all regions organized by chromosome.
    fn into_multi_chrom_overlapper(
        self,
        overlapper_type: OverlapperType,
    ) -> MultiChromOverlapper<I, T>;
}

impl IntoMultiChromOverlapper<u32, Option<String>> for RegionSet {
    fn into_multi_chrom_overlapper(
        self,
        overlapper_type: OverlapperType,
    ) -> MultiChromOverlapper<u32, Option<String>> {
        // instantiate the tree and list of intervals
        let mut core: HashMap<String, Box<dyn Overlapper<u32, Option<String>>>> =
            HashMap::default();
        let mut intervals: HashMap<String, Vec<Interval<u32, Option<String>>>> = HashMap::default();

        // STEP 1: filter/organize/sort regions into vectors, one for each chrom
        for region in self.regions.into_iter() {
            // create interval
            let interval = Interval {
                start: region.start,
                end: region.end,
                val: region.rest,
            };

            // use chr to get the vector of intervals
            let chr_intervals = intervals.entry(region.chr.clone()).or_default();

            // push interval to vector
            chr_intervals.push(interval);
        }

        //STEP 2: take each vector (one for each chrom) and build the overlapper
        for (chr, chr_intervals) in intervals.into_iter() {
            let lapper: Box<dyn Overlapper<u32, Option<String>>> = match overlapper_type {
                OverlapperType::Bits => Box::new(Bits::build(chr_intervals)),
                OverlapperType::AIList => Box::new(AIList::build(chr_intervals)),
            };
            core.insert(chr.to_string(), lapper);
        }

        MultiChromOverlapper {
            index_maps: core,
            overlapper_type,
            source: None,
        }
    }
}

/// Build a [`MultiChromOverlapper`] from a [`RegionSet`] where each interval's
/// val is the original index (position) of that region in the RegionSet.
///
/// This is useful for overlap operations that need to track which regions in
/// `other` were hit (e.g., `find_overlaps` returning indices).
pub fn build_indexed_overlapper(
    rs: &RegionSet,
    overlapper_type: OverlapperType,
) -> MultiChromOverlapper<u32, usize> {
    MultiChromOverlapper::from_region_set(rs.clone(), overlapper_type)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_core::models::Region;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_basic_overlaps(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 300,
                end: 400,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 600,
                end: 800,
                rest: None,
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_multi_chrom_overlapper(overlapper_type);

        let query = RegionSet::from(vec![Region {
            chr: "chr1".to_string(),
            start: 110,
            end: 210,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].0, "chr1");
        assert_eq!(hits[0].1.start, 100);
        assert_eq!(hits[0].1.end, 200);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_multiple_overlaps_single_query(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 150,
                end: 250,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 180,
                end: 300,
                rest: None,
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_multi_chrom_overlapper(overlapper_type);

        let query = RegionSet::from(vec![Region {
            chr: "chr1".to_string(),
            start: 160,
            end: 190,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 3);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_no_overlaps(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 300,
                end: 400,
                rest: None,
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_multi_chrom_overlapper(overlapper_type);

        let query = RegionSet::from(vec![Region {
            chr: "chr1".to_string(),
            start: 500,
            end: 600,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_multiple_chromosomes(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: None,
            },
            Region {
                chr: "chr2".to_string(),
                start: 300,
                end: 400,
                rest: None,
            },
            Region {
                chr: "chr3".to_string(),
                start: 500,
                end: 600,
                rest: None,
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_multi_chrom_overlapper(overlapper_type);

        let query = RegionSet::from(vec![
            Region {
                chr: "chr1".to_string(),
                start: 150,
                end: 250,
                rest: None,
            },
            Region {
                chr: "chr2".to_string(),
                start: 350,
                end: 450,
                rest: None,
            },
        ]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 2);

        let chr1_hits: Vec<_> = hits.iter().filter(|(chr, _)| chr == "chr1").collect();
        let chr2_hits: Vec<_> = hits.iter().filter(|(chr, _)| chr == "chr2").collect();

        assert_eq!(chr1_hits.len(), 1);
        assert_eq!(chr2_hits.len(), 1);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_exact_boundary_overlaps(#[case] overlapper_type: OverlapperType) {
        let regions = vec![Region {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
            rest: None,
        }];
        let rs = RegionSet::from(regions);
        let gi = rs.into_multi_chrom_overlapper(overlapper_type);

        // Query starts exactly at region end
        let query = RegionSet::from(vec![Region {
            chr: "chr1".to_string(),
            start: 200,
            end: 300,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        // Typically intervals are half-open [start, end), so start=200 shouldn't overlap with end=200
        assert_eq!(hits.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_empty_query(#[case] overlapper_type: OverlapperType) {
        let regions = vec![Region {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
            rest: None,
        }];
        let rs = RegionSet::from(regions);
        let gi = rs.into_multi_chrom_overlapper(overlapper_type);

        let query = RegionSet::from(vec![]);
        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_query_nonexistent_chromosome(#[case] overlapper_type: OverlapperType) {
        let regions = vec![Region {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
            rest: None,
        }];
        let rs = RegionSet::from(regions);
        let gi = rs.into_multi_chrom_overlapper(overlapper_type);

        let query = RegionSet::from(vec![Region {
            chr: "chr99".to_string(),
            start: 100,
            end: 200,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_with_metadata(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: Some("gene_a".to_string()),
            },
            Region {
                chr: "chr1".to_string(),
                start: 300,
                end: 400,
                rest: Some("gene_b".to_string()),
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_multi_chrom_overlapper(overlapper_type);

        let query = RegionSet::from(vec![Region {
            chr: "chr1".to_string(),
            start: 150,
            end: 250,
            rest: None,
        }]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 1);
        assert!(hits[0].1.val.is_some());
        assert_eq!(hits[0].1.val.as_ref().unwrap(), "gene_a");
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_overlapping_query_regions(#[case] overlapper_type: OverlapperType) {
        let regions = vec![
            Region {
                chr: "chr1".to_string(),
                start: 100,
                end: 200,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 300,
                end: 400,
                rest: None,
            },
        ];
        let rs = RegionSet::from(regions);
        let gi = rs.into_multi_chrom_overlapper(overlapper_type);

        // Two query regions that each hit different index regions
        let query = RegionSet::from(vec![
            Region {
                chr: "chr1".to_string(),
                start: 150,
                end: 250,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 350,
                end: 450,
                rest: None,
            },
        ]);

        let hits = gi.find_overlaps(&query);
        assert_eq!(hits.len(), 2);
    }

    // -------------------------------------------------------
    // Tests for MCO query methods (indexed overlapper)
    // -------------------------------------------------------

    fn make_region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    // ── subset_by tests ───────────────────────────────────────────

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_subset_by(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
            make_region("chr2", 500, 600),
        ]);
        let index = MultiChromOverlapper::from_region_set(source, overlapper_type);

        let query = RegionSet::from(vec![
            make_region("chr1", 150, 250), // overlaps source[0]
            make_region("chr2", 550, 650), // overlaps source[2]
        ]);

        let result = index.subset_by(&query);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0].start, 100);
        assert_eq!(result.regions[0].end, 200);
        assert_eq!(result.regions[1].start, 500);
        assert_eq!(result.regions[1].end, 600);
    }

    // ── intersect_all tests ─────────────────────────────────────

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_intersect_all(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
        ]);
        let index = MultiChromOverlapper::from_region_set(source, overlapper_type);

        let query = RegionSet::from(vec![
            make_region("chr1", 150, 350), // overlaps both source regions
        ]);

        let mut result = index.intersect_all(&query);
        // Two intersection fragments: [150,200) and [300,350)
        result.regions.sort_by_key(|r| (r.chr.clone(), r.start));
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0].start, 150);
        assert_eq!(result.regions[0].end, 200);
        assert_eq!(result.regions[1].start, 300);
        assert_eq!(result.regions[1].end, 350);
    }

    // ── count_overlaps tests ────────────────────────────────────

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_count_overlaps(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 150, 200),
            make_region("chr1", 250, 350),
            make_region("chr1", 500, 600),
        ]);
        let index = build_indexed_overlapper(&source, overlapper_type);

        let query = RegionSet::from(vec![make_region("chr1", 100, 300)]);
        let counts = index.count_overlaps(&query);
        assert_eq!(counts, vec![2]);
    }

    // ── any_overlaps tests ──────────────────────────────────────

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_any_overlaps(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![make_region("chr1", 150, 250)]);
        let index = build_indexed_overlapper(&source, overlapper_type);

        let query = RegionSet::from(vec![
            make_region("chr1", 100, 200), // overlaps
            make_region("chr1", 300, 400), // no overlap
        ]);
        let any = index.any_overlaps(&query);
        assert_eq!(any, vec![true, false]);
    }

    // ── find_overlaps_indexed tests ─────────────────────────────

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_find_overlaps_indexed(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 50, 150),  // index 0
            make_region("chr1", 200, 250), // index 1
            make_region("chr1", 400, 500), // index 2
        ]);
        let index = build_indexed_overlapper(&source, overlapper_type);

        let query = RegionSet::from(vec![make_region("chr1", 100, 300)]);
        let indices = index.find_overlaps_indexed(&query);
        assert_eq!(indices.len(), 1);
        let mut sorted = indices[0].clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1]);
    }

    // ── edge case tests ─────────────────────────────────────────

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_mco_empty_query(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let index = build_indexed_overlapper(&source, overlapper_type);
        let query = RegionSet::from(vec![]);

        assert_eq!(index.count_overlaps(&query), Vec::<usize>::new());
        assert_eq!(index.any_overlaps(&query), Vec::<bool>::new());
        assert_eq!(
            index.find_overlaps_indexed(&query),
            Vec::<Vec<usize>>::new()
        );
        assert_eq!(index.intersect_all(&query).regions.len(), 0);
        assert_eq!(index.subset_by(&query).regions.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_mco_empty_index(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![]);
        let index = build_indexed_overlapper(&source, overlapper_type);
        let query = RegionSet::from(vec![make_region("chr1", 100, 200)]);

        assert_eq!(index.count_overlaps(&query), vec![0]);
        assert_eq!(index.any_overlaps(&query), vec![false]);
        assert_eq!(index.find_overlaps_indexed(&query), vec![vec![]]);
        assert_eq!(index.intersect_all(&query).regions.len(), 0);
        assert_eq!(index.subset_by(&query).regions.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_mco_multi_chrom(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr2", 100, 200),
        ]);
        let index = build_indexed_overlapper(&source, overlapper_type);

        let query = RegionSet::from(vec![
            make_region("chr1", 150, 250),
            make_region("chr2", 150, 250),
            make_region("chr3", 150, 250), // nonexistent chrom
        ]);

        let counts = index.count_overlaps(&query);
        assert_eq!(counts, vec![1, 1, 0]);

        let any = index.any_overlaps(&query);
        assert_eq!(any, vec![true, true, false]);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_mco_index_reuse_consistency(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
            make_region("chr2", 500, 600),
        ]);
        let index = build_indexed_overlapper(&source, overlapper_type);

        let query = RegionSet::from(vec![
            make_region("chr1", 150, 350),
            make_region("chr2", 550, 650),
        ]);

        let counts1 = index.count_overlaps(&query);
        let counts2 = index.count_overlaps(&query);
        assert_eq!(counts1, counts2);

        let any1 = index.any_overlaps(&query);
        let any2 = index.any_overlaps(&query);
        assert_eq!(any1, any2);

        let find1 = index.find_overlaps_indexed(&query);
        let find2 = index.find_overlaps_indexed(&query);
        assert_eq!(find1, find2);
    }

    // ── MCO IntervalSetOps matches RegionSet IntervalSetOps ─────

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_mco_intervalsetops_matches_regionset(#[case] overlapper_type: OverlapperType) {
        use gtars_core::models::IntervalSetOps;

        let source = RegionSet::from(vec![
            make_region("chr1", 100, 300),
            make_region("chr1", 250, 500),
            make_region("chr2", 100, 400),
        ]);
        let other = RegionSet::from(vec![
            make_region("chr1", 200, 350),
            make_region("chr1", 450, 600),
            make_region("chr2", 300, 500),
        ]);

        let mco = MultiChromOverlapper::from_region_set(source.clone(), overlapper_type);

        // setdiff
        let rs_diff = source.setdiff(&other);
        let mco_diff = mco.setdiff(&other);
        assert_eq!(rs_diff.regions, mco_diff.regions, "setdiff mismatch");

        // intersect
        let rs_inter = source.intersect(&other);
        let mco_inter = mco.intersect(&other);
        assert_eq!(rs_inter.regions, mco_inter.regions, "intersect mismatch");

        // union
        let rs_union = source.union(&other);
        let mco_union = mco.union(&other);
        assert_eq!(rs_union.regions, mco_union.regions, "union mismatch");

        // jaccard
        let rs_jac = source.jaccard(&other);
        let mco_jac = mco.jaccard(&other);
        assert!((rs_jac - mco_jac).abs() < 1e-10, "jaccard mismatch");

        // coverage
        let rs_cov = source.coverage(&other);
        let mco_cov = mco.coverage(&other);
        assert!((rs_cov - mco_cov).abs() < 1e-10, "coverage mismatch");

        // overlap_coefficient
        let rs_oc = source.overlap_coefficient(&other);
        let mco_oc = mco.overlap_coefficient(&other);
        assert!((rs_oc - mco_oc).abs() < 1e-10, "overlap_coefficient mismatch");

        // subtract (alias for setdiff)
        let rs_sub = source.subtract(&other);
        let mco_sub = mco.subtract(&other);
        assert_eq!(rs_sub.regions, mco_sub.regions, "subtract mismatch");

        // closest
        let rs_closest = source.closest(&other);
        let mco_closest = mco.closest(&other);
        assert_eq!(rs_closest, mco_closest, "closest mismatch");

        // cluster
        let rs_cluster = source.cluster(10);
        let mco_cluster = mco.cluster(10);
        assert_eq!(rs_cluster, mco_cluster, "cluster mismatch");
    }

    // ── source() accessor test ──────────────────────────────────

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_source_accessor(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr2", 300, 400),
        ]);
        let mco = MultiChromOverlapper::from_region_set(source.clone(), overlapper_type);

        assert_eq!(mco.source().regions.len(), 2);
        assert_eq!(mco.source().regions[0].start, 100);
        assert_eq!(mco.source().regions[1].chr, "chr2");
    }

    // ── from_region_set constructor test ─────────────────────────

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_from_region_set(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
        ]);
        let mco = MultiChromOverlapper::from_region_set(source, overlapper_type);

        let query = RegionSet::from(vec![make_region("chr1", 150, 350)]);
        let counts = mco.count_overlaps(&query);
        assert_eq!(counts, vec![2]);

        // Also verify source is accessible
        assert_eq!(mco.source().regions.len(), 2);
    }

    // ── single region test ──────────────────────────────────────

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_mco_single_region(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let mco = MultiChromOverlapper::from_region_set(source, overlapper_type);

        let query = RegionSet::from(vec![make_region("chr1", 150, 250)]);

        assert_eq!(mco.count_overlaps(&query), vec![1]);
        assert_eq!(mco.any_overlaps(&query), vec![true]);
        assert_eq!(mco.subset_by(&query).regions.len(), 1);

        let fragments = mco.intersect_all(&query);
        assert_eq!(fragments.regions.len(), 1);
        assert_eq!(fragments.regions[0].start, 150);
        assert_eq!(fragments.regions[0].end, 200);
    }
}
