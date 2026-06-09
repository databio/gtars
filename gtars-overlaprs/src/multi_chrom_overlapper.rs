//! Genome-wide interval indexing for efficient multi-chromosome overlap queries.
//!
//! This module provides [`MultiChromOverlapper`](crate::multi_chrom_overlapper::MultiChromOverlapper), a data structure that extends single-chromosome
//! overlap data structures (like [`AIList`](crate::AIList) or [`Bits`](crate::Bits)) to handle genome-wide queries
//! across multiple chromosomes efficiently.
//!
//!
//! # Examples
//!
//! ## Basic Usage
//!
//! ```
//! use gtars_overlaprs::{OverlapperType, multi_chrom_overlapper::IntoMultiChromOverlapper};
//! use gtars_core::models::{Region, RegionSet};
//!
//! // Create a genome-wide set of intervals (e.g., gene annotations)
//! let genes = vec![
//!     Region { chr: "chr1".to_string(), start: 1000, end: 2000, rest: Some("BRCA1".to_string()) },
//!     Region { chr: "chr1".to_string(), start: 5000, end: 6000, rest: Some("TP53".to_string()) },
//!     Region { chr: "chr2".to_string(), start: 1000, end: 3000, rest: Some("EGFR".to_string()) },
//! ];
//!
//! let gene_set = RegionSet::from(genes);
//! let multi_chrom_overlapper = gene_set.into_multi_chrom_overlapper(OverlapperType::AIList);
//!
//! // Query multiple regions across different chromosomes
//! let query_regions = RegionSet::from(vec![
//!     Region { chr: "chr1".to_string(), start: 1500, end: 2500, rest: None },
//!     Region { chr: "chr2".to_string(), start: 2000, end: 4000, rest: None },
//! ]);
//!
//! let overlaps = multi_chrom_overlapper.find_overlaps(&query_regions);
//! println!("Found {} overlapping features", overlaps.len());
//! ```

use std::{collections::HashMap, fmt::Debug};

use gtars_core::models::{Interval, Region, RegionSet};
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
/// ## Source-free indexed variant (`MultiChromOverlapper<u32, ()>`)
///
/// `MultiChromOverlapper` is **internal index machinery**, not a public
/// "indexed region set" — that role belongs to [`IndexedRegionSet`](crate::IndexedRegionSet).
/// When built from a [`RegionSet`] (via [`build_indexed_overlapper`] or
/// [`from_region_set`](MultiChromOverlapper::from_region_set)) it is
/// `MultiChromOverlapper<u32, ()>`: it stores **no** source RegionSet and **no**
/// back-reference payload. Intervals are enumerated straight from the index via
/// [`Overlapper::iter`], and per-chromosome reduced regions are reconstructed on
/// demand by `reduced_by_chr`.
///
/// This enables:
///
/// - **[`IntervalSetOps`](gtars_core::models::IntervalSetOps) trait**: index-native
///   set algebra computed from the index (not delegated to a stored source).
/// - **MCO query methods**: [`subset_by`](MultiChromOverlapper::subset_by),
///   [`count_overlaps`](MultiChromOverlapper::count_overlaps),
///   [`any_overlaps`](MultiChromOverlapper::any_overlaps),
///   [`find_overlaps_regions`](MultiChromOverlapper::find_overlaps_regions),
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
        }
    }
}

/// Build a source-free [`MultiChromOverlapper`] index from a [`RegionSet`].
///
/// The resulting `MultiChromOverlapper<u32, ()>` stores no source RegionSet and
/// no back-reference payload — intervals carry `()`. Regions are reconstructed
/// from the index on demand via [`Overlapper::iter`].
pub fn build_indexed_overlapper(
    rs: &RegionSet,
    overlapper_type: OverlapperType,
) -> MultiChromOverlapper<u32, ()> {
    let mut core: HashMap<String, Box<dyn Overlapper<u32, ()>>> = HashMap::default();
    let mut intervals: HashMap<String, Vec<Interval<u32, ()>>> = HashMap::default();

    for region in rs.regions.iter() {
        let interval = Interval {
            start: region.start,
            end: region.end,
            val: (),
        };
        let chr_intervals = intervals.entry(region.chr.clone()).or_default();
        chr_intervals.push(interval);
    }

    for (chr, chr_intervals) in intervals.into_iter() {
        let lapper: Box<dyn Overlapper<u32, ()>> = match overlapper_type {
            OverlapperType::Bits => Box::new(Bits::build(chr_intervals)),
            OverlapperType::AIList => Box::new(AIList::build(chr_intervals)),
        };
        core.insert(chr, lapper);
    }

    MultiChromOverlapper {
        index_maps: core,
        overlapper_type,
    }
}

impl MultiChromOverlapper<u32, ()> {
    /// Build a new source-free MCO index from a [`RegionSet`].
    ///
    /// This enables the build-once-query-many pattern: build the index once,
    /// then call any of the query methods or `IntervalSetOps` methods repeatedly
    /// without rebuilding. The source RegionSet is not retained.
    pub fn from_region_set(rs: RegionSet, overlapper_type: OverlapperType) -> Self {
        build_indexed_overlapper(&rs, overlapper_type)
    }

    /// Reconstruct a [`RegionSet`] by enumerating every interval in the index.
    ///
    /// Regions are reconstructed directly from the index (via [`Overlapper::iter`]),
    /// not from a stored source. The resulting order is the per-chromosome stored
    /// order of the underlying overlappers and is not guaranteed to match any
    /// original input order. Reconstructed regions carry no `rest` metadata.
    pub fn to_region_set(&self) -> RegionSet {
        let mut regions: Vec<Region> = Vec::new();
        let mut chrs: Vec<&String> = self.index_maps.keys().collect();
        chrs.sort();
        for chr in chrs {
            let lapper = &self.index_maps[chr];
            for iv in lapper.iter() {
                regions.push(Region {
                    chr: chr.clone(),
                    start: iv.start,
                    end: iv.end,
                    rest: None,
                });
            }
        }
        RegionSet::from(regions)
    }

    /// Union of the index-reconstructed regions with `other`.
    pub fn union(&self, other: &RegionSet) -> RegionSet {
        self.to_region_set().union(other)
    }

    /// Cluster the index-reconstructed regions within `max_gap`.
    pub fn cluster(&self, max_gap: u32) -> Vec<u32> {
        self.to_region_set().cluster(max_gap)
    }

    /// Produce sorted, reduced (merged) regions per chromosome directly from the
    /// index, without reconstructing a full `RegionSet`.
    fn reduced_by_chr(&self) -> HashMap<String, Vec<Region>> {
        let mut result: HashMap<String, Vec<Region>> = HashMap::new();

        for (chr, overlapper) in &self.index_maps {
            // Collect intervals from the overlapper and sort by start.
            let mut intervals: Vec<(u32, u32)> =
                overlapper.iter().map(|iv| (iv.start, iv.end)).collect();
            intervals.sort_unstable();

            if intervals.is_empty() {
                continue;
            }

            // Merge overlapping/adjacent intervals inline.
            let mut merged: Vec<Region> = Vec::new();
            let (mut cur_start, mut cur_end) = intervals[0];

            for &(s, e) in &intervals[1..] {
                if s <= cur_end {
                    cur_end = cur_end.max(e);
                } else {
                    merged.push(Region {
                        chr: chr.clone(),
                        start: cur_start,
                        end: cur_end,
                        rest: None,
                    });
                    cur_start = s;
                    cur_end = e;
                }
            }
            merged.push(Region {
                chr: chr.clone(),
                start: cur_start,
                end: cur_end,
                rest: None,
            });

            result.insert(chr.clone(), merged);
        }

        result
    }

    /// Return reconstructed regions from the index that overlap any region in
    /// `query`, honoring an optional minimum overlap in base pairs.
    ///
    /// Because the index is source-free, the returned regions are reconstructed
    /// from the index and carry no `rest` metadata. Results are deduplicated and
    /// returned sorted by chromosome and start.
    pub fn subset_by(&self, query: &RegionSet) -> RegionSet {
        self.subset_by_overlaps(query, None)
    }

    /// Return reconstructed overlapping regions, honoring `min_overlap` (bp).
    pub fn subset_by_overlaps(&self, query: &RegionSet, min_overlap: Option<i32>) -> RegionSet {
        let min_bp = min_overlap.unwrap_or(0);
        let mut hits: std::collections::BTreeSet<(String, u32, u32)> =
            std::collections::BTreeSet::new();
        for region in &query.regions {
            for iv in self.find_overlaps_for_region(&region.chr, region.start, region.end) {
                if min_bp > 1
                    && overlap_bp(region.start, region.end, iv.start, iv.end) < min_bp as i64
                {
                    continue;
                }
                hits.insert((region.chr.clone(), iv.start, iv.end));
            }
        }
        let kept: Vec<Region> = hits
            .into_iter()
            .map(|(chr, start, end)| Region {
                chr,
                start,
                end,
                rest: None,
            })
            .collect();
        RegionSet::from(kept)
    }

    /// For each region in `query`, count how many indexed regions overlap it.
    ///
    /// Optionally filter by minimum overlap in base pairs.
    pub fn count_overlaps(&self, query: &RegionSet, min_overlap: Option<i32>) -> Vec<usize> {
        let min_bp = min_overlap.unwrap_or(0);
        query
            .regions
            .iter()
            .map(|region| {
                self.find_overlaps_for_region(&region.chr, region.start, region.end)
                    .filter(|iv| {
                        min_bp <= 1
                            || overlap_bp(region.start, region.end, iv.start, iv.end)
                                >= min_bp as i64
                    })
                    .count()
            })
            .collect()
    }

    /// For each region in `query`, whether any indexed region overlaps it.
    ///
    /// Optionally filter by minimum overlap in base pairs.
    pub fn any_overlaps(&self, query: &RegionSet, min_overlap: Option<i32>) -> Vec<bool> {
        let min_bp = min_overlap.unwrap_or(0);
        query
            .regions
            .iter()
            .map(|region| {
                self.find_overlaps_for_region(&region.chr, region.start, region.end)
                    .any(|iv| {
                        min_bp <= 1
                            || overlap_bp(region.start, region.end, iv.start, iv.end)
                                >= min_bp as i64
                    })
            })
            .collect()
    }

    /// For each region in `query`, return the reconstructed indexed regions that
    /// overlap it, honoring an optional minimum overlap in base pairs.
    ///
    /// Because the index is source-free, overlapping regions are reconstructed
    /// from the index (carrying no `rest` metadata) rather than returned as
    /// indices into a stored source.
    pub fn find_overlaps_regions(
        &self,
        query: &RegionSet,
        min_overlap: Option<i32>,
    ) -> Vec<Vec<Region>> {
        let min_bp = min_overlap.unwrap_or(0);
        query
            .regions
            .iter()
            .map(|region| {
                self.find_overlaps_for_region(&region.chr, region.start, region.end)
                    .filter(|iv| {
                        min_bp <= 1
                            || overlap_bp(region.start, region.end, iv.start, iv.end)
                                >= min_bp as i64
                    })
                    .map(|iv| Region {
                        chr: region.chr.clone(),
                        start: iv.start,
                        end: iv.end,
                        rest: None,
                    })
                    .collect()
            })
            .collect()
    }

    /// Return all indexed regions (reconstructed) that overlap any region in
    /// `query`, deduplicated and sorted. Equivalent to [`subset_by`](Self::subset_by).
    pub fn intersect_all(&self, query: &RegionSet) -> RegionSet {
        self.subset_by(query)
    }
}

/// Compute actual overlap in base pairs between two regions.
#[inline]
fn overlap_bp(a_start: u32, a_end: u32, b_start: u32, b_end: u32) -> i64 {
    a_end.min(b_end) as i64 - a_start.max(b_start) as i64
}

// ── Index-native IntervalSetOps implementation for MCO ──────────────────
//
// Uses per-chromosome sweep-line algorithms operating directly on the index
// (via `reduced_by_chr`), avoiding the cost of reconstructing a full
// RegionSet and never consulting a stored source (there is none).

use gtars_core::models::region_set::{sweep_intersect_chr, sweep_setdiff_chr};
use gtars_core::models::IntervalSetOps;

impl IntervalSetOps for MultiChromOverlapper<u32, ()> {
    fn setdiff(&self, other: &RegionSet) -> RegionSet {
        let self_by_chr = self.reduced_by_chr();

        let b = other.reduce();
        let mut b_by_chr: HashMap<String, Vec<Region>> = HashMap::new();
        for r in b.regions {
            b_by_chr.entry(r.chr.clone()).or_default().push(r);
        }

        // Iterate chromosomes in sorted order for deterministic output.
        let mut chrs: Vec<&String> = self_by_chr.keys().collect();
        chrs.sort();

        let mut result: Vec<Region> = Vec::new();
        for chr in chrs {
            let a_regions = &self_by_chr[chr];
            let empty = vec![];
            let b_regions = b_by_chr.get(chr.as_str()).unwrap_or(&empty);
            result.extend(sweep_setdiff_chr(chr, a_regions, b_regions));
        }

        RegionSet::from(result)
    }

    fn intersect(&self, other: &RegionSet) -> RegionSet {
        let self_by_chr = self.reduced_by_chr();

        let b = other.reduce();
        let mut b_by_chr: HashMap<String, Vec<Region>> = HashMap::new();
        for r in b.regions {
            b_by_chr.entry(r.chr.clone()).or_default().push(r);
        }

        let mut chrs: Vec<&String> = self_by_chr.keys().collect();
        chrs.sort();

        let mut result: Vec<Region> = Vec::new();
        for chr in chrs {
            let a_regions = &self_by_chr[chr];
            if let Some(b_regions) = b_by_chr.get(chr.as_str()) {
                result.extend(sweep_intersect_chr(chr, a_regions, b_regions));
            }
        }

        RegionSet::from(result)
    }

    fn jaccard(&self, other: &RegionSet) -> f64 {
        let self_by_chr = self.reduced_by_chr();
        let self_bp: u32 = self_by_chr
            .values()
            .flat_map(|rs| rs.iter())
            .map(|r| r.end - r.start)
            .sum();
        let b = other.reduce();
        let other_bp = b.nucleotides_length();

        let mut b_by_chr: HashMap<String, Vec<Region>> = HashMap::new();
        for r in b.regions {
            b_by_chr.entry(r.chr.clone()).or_default().push(r);
        }

        let intersection_bp = intersect_bp(&self_by_chr, &b_by_chr);
        let union_bp = self_bp + other_bp - intersection_bp;
        if union_bp == 0 {
            return 0.0;
        }
        intersection_bp as f64 / union_bp as f64
    }

    fn coverage(&self, other: &RegionSet) -> f64 {
        let self_by_chr = self.reduced_by_chr();
        let self_bp: u32 = self_by_chr
            .values()
            .flat_map(|rs| rs.iter())
            .map(|r| r.end - r.start)
            .sum();
        if self_bp == 0 {
            return 0.0;
        }
        let diff_bp = self.setdiff(other).nucleotides_length();
        1.0 - (diff_bp as f64 / self_bp as f64)
    }

    fn overlap_coefficient(&self, other: &RegionSet) -> f64 {
        let self_by_chr = self.reduced_by_chr();
        let self_bp: u32 = self_by_chr
            .values()
            .flat_map(|rs| rs.iter())
            .map(|r| r.end - r.start)
            .sum();
        let b = other.reduce();
        let other_bp = b.nucleotides_length();
        let min_bp = self_bp.min(other_bp);
        if min_bp == 0 {
            return 0.0;
        }

        let mut b_by_chr: HashMap<String, Vec<Region>> = HashMap::new();
        for r in b.regions {
            b_by_chr.entry(r.chr.clone()).or_default().push(r);
        }

        let intersection_bp = intersect_bp(&self_by_chr, &b_by_chr);
        intersection_bp as f64 / min_bp as f64
    }

    fn closest(&self, other: &RegionSet) -> Vec<(usize, usize, i64)> {
        // `closest` indices reference the order produced by `to_region_set()`
        // (the reconstructed ordering), not any original source order.
        self.to_region_set().closest(other)
    }
}

/// Sum of intersected base pairs across chromosomes, computed per-chromosome
/// from the sweep-line intersection.
fn intersect_bp(
    self_by_chr: &HashMap<String, Vec<Region>>,
    b_by_chr: &HashMap<String, Vec<Region>>,
) -> u32 {
    let mut chrs: Vec<&String> = self_by_chr.keys().collect();
    chrs.sort();
    let mut total = 0u32;
    for chr in chrs {
        if let Some(b_regions) = b_by_chr.get(chr.as_str()) {
            for r in sweep_intersect_chr(chr, &self_by_chr[chr], b_regions) {
                total += r.end - r.start;
            }
        }
    }
    total
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
    // Tests for source-free MCO query methods + IntervalSetOps
    // -------------------------------------------------------

    fn make_region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    fn coords(rs: &RegionSet) -> Vec<(String, u32, u32)> {
        let mut v: Vec<(String, u32, u32)> = rs
            .regions
            .iter()
            .map(|r| (r.chr.clone(), r.start, r.end))
            .collect();
        v.sort();
        v
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_intersect_all(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
            make_region("chr2", 500, 600),
        ]);
        let index = build_indexed_overlapper(&source, overlapper_type);

        let query = RegionSet::from(vec![
            make_region("chr1", 150, 250), // overlaps source[0]
            make_region("chr2", 550, 650), // overlaps source[2]
        ]);

        let result = index.intersect_all(&query);
        assert_eq!(
            coords(&result),
            vec![
                ("chr1".to_string(), 100, 200),
                ("chr2".to_string(), 500, 600),
            ]
        );
    }

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
        assert_eq!(index.count_overlaps(&query, None), vec![2]);
    }

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
        assert_eq!(index.any_overlaps(&query, None), vec![true, false]);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_find_overlaps_regions(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 50, 150),
            make_region("chr1", 200, 250),
            make_region("chr1", 400, 500),
        ]);
        let index = build_indexed_overlapper(&source, overlapper_type);

        let query = RegionSet::from(vec![make_region("chr1", 100, 300)]);
        let regions = index.find_overlaps_regions(&query, None);
        assert_eq!(regions.len(), 1);
        let mut got: Vec<(u32, u32)> = regions[0].iter().map(|r| (r.start, r.end)).collect();
        got.sort();
        assert_eq!(got, vec![(50, 150), (200, 250)]);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_min_overlap_filtering(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![make_region("chr1", 100, 110)]);
        let index = build_indexed_overlapper(&source, overlapper_type);

        // 5bp overlap with query [105, 200)
        let query = RegionSet::from(vec![make_region("chr1", 105, 200)]);
        assert_eq!(index.count_overlaps(&query, Some(5)), vec![1]);
        assert_eq!(index.count_overlaps(&query, Some(6)), vec![0]);
        assert_eq!(index.any_overlaps(&query, Some(6)), vec![false]);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_mco_empty_query(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let index = build_indexed_overlapper(&source, overlapper_type);
        let query = RegionSet::from(vec![]);

        assert_eq!(index.count_overlaps(&query, None), Vec::<usize>::new());
        assert_eq!(index.any_overlaps(&query, None), Vec::<bool>::new());
        assert_eq!(index.find_overlaps_regions(&query, None).len(), 0);
        assert_eq!(index.intersect_all(&query).regions.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_mco_empty_index(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![]);
        let index = build_indexed_overlapper(&source, overlapper_type);
        let query = RegionSet::from(vec![make_region("chr1", 100, 200)]);

        assert_eq!(index.count_overlaps(&query, None), vec![0]);
        assert_eq!(index.any_overlaps(&query, None), vec![false]);
        assert_eq!(index.intersect_all(&query).regions.len(), 0);
    }

    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_to_region_set_roundtrip(#[case] overlapper_type: OverlapperType) {
        let source = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr2", 300, 400),
            make_region("chr1", 500, 600),
        ]);
        let index = build_indexed_overlapper(&source, overlapper_type);
        let reconstructed = index.to_region_set();
        assert_eq!(coords(&reconstructed), coords(&source));
    }

    /// Cross-check that the index-native IntervalSetOps on the source-free MCO
    /// agrees with RegionSet's direct implementation for all 6 methods.
    #[rstest]
    #[case(OverlapperType::AIList)]
    #[case(OverlapperType::Bits)]
    fn test_mco_intervalsetops_matches_regionset(#[case] overlapper_type: OverlapperType) {
        // Includes overlapping self-intervals on chr1 to exercise the merge
        // inside reduced_by_chr().
        let self_rs = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 180, 260), // overlaps previous
            make_region("chr1", 400, 500),
            make_region("chr2", 100, 200),
        ]);
        let other = RegionSet::from(vec![
            make_region("chr1", 150, 220),
            make_region("chr1", 450, 470),
            make_region("chr2", 50, 150),
            make_region("chr2", 300, 400),
            make_region("chr3", 0, 100),
        ]);

        let index = build_indexed_overlapper(&self_rs, overlapper_type);

        assert_eq!(
            coords(&index.setdiff(&other)),
            coords(&self_rs.setdiff(&other)),
            "setdiff mismatch"
        );
        assert_eq!(
            coords(&index.intersect(&other)),
            coords(&self_rs.intersect(&other)),
            "intersect mismatch"
        );
        assert!(
            (index.jaccard(&other) - self_rs.jaccard(&other)).abs() < 1e-9,
            "jaccard mismatch"
        );
        assert!(
            (index.coverage(&other) - self_rs.coverage(&other)).abs() < 1e-9,
            "coverage mismatch"
        );
        assert!(
            (index.overlap_coefficient(&other) - self_rs.overlap_coefficient(&other)).abs() < 1e-9,
            "overlap_coefficient mismatch"
        );

        // closest: indices reference reconstructed ordering, so compare the
        // distances (which are order-independent) as a set.
        let mut mco_dists: Vec<i64> = index.closest(&other).iter().map(|t| t.2).collect();
        let mut rs_dists: Vec<i64> = self_rs.closest(&other).iter().map(|t| t.2).collect();
        mco_dists.sort();
        rs_dists.sort();
        assert_eq!(mco_dists, rs_dists, "closest distances mismatch");
    }
}
