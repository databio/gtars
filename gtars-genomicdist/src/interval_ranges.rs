//! Interval set algebra operations for genomic region sets.
//!
//! Provides GenomicRanges/IRanges-style operations: trim, promoters, reduce,
//! setdiff, and pintersect. All operations use 0-based half-open coordinates
//! (BED convention) and are strand-unaware.

use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};
use gtars_overlaprs::OverlapperType;
use gtars_overlaprs::multi_chrom_overlapper::IntoMultiChromOverlapper;

use crate::models::{SortedRegionSet, Strand, StrandedRegionSet};

/// Interval set algebra operations on genomic region sets.
///
/// Modeled after R's GenomicRanges/IRanges package. All functions return new
/// `RegionSet` instances (immutable pattern) with 0-based half-open coordinates.
///
/// **Note:** All operations produce regions with `rest: None`. Metadata from
/// the `rest` field of input regions is not preserved — operations like `reduce()`
/// merge multiple regions, so there is no unambiguous `rest` to carry forward.
pub trait IntervalRanges {
    /// Clamp regions to chromosome boundaries.
    ///
    /// Regions extending past chromosome ends are trimmed to `[0, chrom_size)`.
    /// Regions on chromosomes not present in `chrom_sizes` are dropped.
    /// Empty regions (start > end after clamping) are dropped; zero-width
    /// regions where start == end are kept.
    fn trim(&self, chrom_sizes: &HashMap<String, u32>) -> RegionSet;

    /// Generate promoter regions relative to each region's start position.
    ///
    /// For each region, produces `[start - upstream, start + downstream)`.
    /// Uses saturating subtraction at coordinate 0.
    fn promoters(&self, upstream: u32, downstream: u32) -> RegionSet;

    /// Merge overlapping and adjacent intervals per chromosome.
    ///
    /// Sorts by (chr, start), then sweeps to merge intervals where
    /// `next.start <= current.end`. Returns a minimal set of non-overlapping regions.
    fn reduce(&self) -> RegionSet;

    /// Subtract one region set from another (set difference).
    ///
    /// Removes portions of `self` that overlap with `other`. Both inputs are
    /// reduced internally before subtraction. Operates per-chromosome with a
    /// sweep-line algorithm.
    ///
    /// # Example
    /// ```text
    /// A: chr1 100–200
    /// B: chr1 120–140, chr1 160–180
    /// setdiff(A, B): chr1 100–120, chr1 140–160, chr1 180–200
    /// ```
    fn setdiff(&self, other: &RegionSet) -> RegionSet;

    /// Pairwise intersection of two region sets by index position.
    ///
    /// For each pair at the same index, computes `[max(a.start, b.start), min(a.end, b.end))`.
    /// Pairs with no overlap or mismatched chromosomes are dropped.
    /// If the region sets differ in length, intersections are computed only up
    /// to the length of the shorter set; excess regions are ignored.
    ///
    /// Note: this intersects by index position (1st with 1st, 2nd with 2nd, etc.),
    /// **not** by genomic overlap across all regions.
    ///
    /// # Example
    /// ```text
    /// A: chr1 0–10, chr1 20–30
    /// B: chr1 5–15, chr1 25–35
    /// pintersect(A, B): chr1 5–10, chr1 25–30
    /// ```
    fn pintersect(&self, other: &RegionSet) -> RegionSet;

    /// Combine two region sets without merging overlapping intervals.
    ///
    /// Clones regions from both sets into a single `RegionSet`. No sorting,
    /// deduplication, or merging is performed.
    fn concat(&self, other: &RegionSet) -> RegionSet;

    /// Merge two region sets into a minimal non-overlapping set.
    ///
    /// Equivalent to `self.concat(other).reduce()`.
    fn union(&self, other: &RegionSet) -> RegionSet;

    /// Nucleotide-level Jaccard similarity between two region sets.
    ///
    /// Computes `|intersection| / |union|` in base pairs, where both sets are
    /// reduced first. Returns 0.0 if the union has zero base pairs.
    fn jaccard(&self, other: &RegionSet) -> f64;

    /// Fraction of self's base pairs covered by other.
    ///
    /// Both sets are reduced first. Returns
    /// `1.0 - (setdiff_bp / self_bp)`, or 0.0 if self has zero base pairs.
    /// Result is in [0.0, 1.0].
    fn coverage(&self, other: &RegionSet) -> f64;

    /// Overlap coefficient between two region sets.
    ///
    /// Computes `|intersection_bp| / min(|self_bp|, |other_bp|)` after reducing
    /// both sets. Returns 0.0 if either set has zero base pairs.
    /// Result is in [0.0, 1.0].
    fn overlap_coefficient(&self, other: &RegionSet) -> f64;

    /// Shift all regions by a fixed offset.
    ///
    /// Adds `offset` to both start and end of every region. Negative offsets
    /// use saturating subtraction at coordinate 0.
    fn shift(&self, offset: i64) -> RegionSet;

    /// Generate flanking regions.
    ///
    /// - `use_start = true`: flank upstream of start → `[start - width, start)`
    /// - `use_start = false`: flank downstream of end → `[end, end + width)`
    /// - `both = true`: flank on both sides → `[anchor - width, anchor + width)`
    ///
    /// Uses saturating arithmetic at coordinate 0.
    fn flank(&self, width: u32, use_start: bool, both: bool) -> RegionSet;

    /// Resize regions to a fixed width, anchored at start, end, or center.
    ///
    /// - `"start"`: keep start, set end = start + width
    /// - `"end"`: keep end, set start = end - width
    /// - `"center"`: keep midpoint, expand symmetrically
    fn resize(&self, width: u32, fix: &str) -> RegionSet;

    /// Narrow each region by specifying a relative sub-range within it.
    ///
    /// Parameters are 1-based relative positions within each region (matching
    /// GenomicRanges convention). Exactly two of the three must be provided.
    fn narrow(&self, start: Option<u32>, end: Option<u32>, width: Option<u32>) -> RegionSet;

    /// Break all regions into non-overlapping disjoint pieces.
    ///
    /// Every boundary in the input becomes a boundary in the output. The result
    /// tiles the covered positions exactly, with no overlaps.
    fn disjoin(&self) -> RegionSet;

    /// Return the gaps between regions per chromosome, bounded by chromosome sizes.
    ///
    /// Reduces the input first, then emits intervals that tile the peak-free
    /// regions of each chromosome listed in `chrom_sizes`:
    ///
    /// - a **leading gap** from position 0 to the first region's start
    ///   (omitted if the first region starts at 0),
    /// - an **inter-region gap** between each consecutive pair of reduced
    ///   regions,
    /// - a **trailing gap** from the last region's end to the chromosome
    ///   size (omitted if the last region already reaches the chromosome
    ///   end, or extends past it due to assembly mismatch),
    /// - a **full-chromosome gap** `0..chrom_size` for any chromosome in
    ///   `chrom_sizes` that has no regions at all.
    ///
    /// Regions on chromosomes not present in `chrom_sizes` are skipped.
    /// Regions that extend past the stated chromosome size are clipped to
    /// `chrom_size` when computing the trailing gap, matching the
    /// clipping behavior of `trim()`.
    fn gaps(&self, chrom_sizes: &HashMap<String, u32>) -> RegionSet;

    /// Range-level intersection of two region sets.
    ///
    /// Returns the set of positions covered by *both* sets. Both inputs are
    /// reduced before computing the intersection.
    fn intersect(&self, other: &RegionSet) -> RegionSet;

    /// All-vs-all genomic intersection.
    ///
    /// For each pair of overlapping regions between `self` and `other`,
    /// compute intersection coordinates `[max(a.start, b.start), min(a.end, b.end))`.
    /// Returns a `RegionSet` of all intersection fragments.
    ///
    /// Unlike `pintersect` (which pairs by index position), this finds ALL
    /// overlapping pairs across the two sets using an AIList index.
    fn intersect_all(&self, other: &RegionSet) -> RegionSet;

    /// Genomic subtraction (alias for `setdiff`).
    ///
    /// Remove portions of `self` that overlap with `other`.
    /// Both inputs are reduced before subtraction.
    fn subtract(&self, other: &RegionSet) -> RegionSet;

    /// Find the nearest region in `other` for each region in `self`.
    ///
    /// Returns a list of `(self_index, other_index, distance)` tuples.
    /// Distance is 0 for overlapping regions. Regions on chromosomes
    /// absent in `other` are omitted from results.
    fn closest(&self, other: &RegionSet) -> Vec<(usize, usize, i64)>;

    /// Cluster nearby regions.
    ///
    /// Assign a cluster ID to each region. Regions within `max_gap`
    /// distance on the same chromosome are assigned the same cluster.
    /// Returns cluster IDs in original region order.
    fn cluster(&self, max_gap: u32) -> Vec<u32>;
}

impl IntervalRanges for RegionSet {
    fn trim(&self, chrom_sizes: &HashMap<String, u32>) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .filter_map(|r| {
                let chrom_size = chrom_sizes.get(&r.chr)?;
                let start = r.start.min(*chrom_size);
                let end = r.end.min(*chrom_size);
                if start > end {
                    None
                } else {
                    Some(Region {
                        chr: r.chr.clone(),
                        start,
                        end,
                        rest: None,
                    })
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    fn promoters(&self, upstream: u32, downstream: u32) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| Region {
                chr: r.chr.clone(),
                start: r.start.saturating_sub(upstream),
                end: r.start.saturating_add(downstream),
                rest: None,
            })
            .collect();
        RegionSet::from(regions)
    }

    fn reduce(&self) -> RegionSet {
        if self.regions.is_empty() {
            return RegionSet::from(Vec::<Region>::new());
        }

        // Clone and sort via SortedRegionSet (moves, sorts in place — no double clone)
        let sorted = SortedRegionSet::new(RegionSet::from(self.regions.clone()));
        let regions = &sorted.0.regions;

        let mut merged: Vec<Region> = Vec::new();
        let mut current = regions[0].clone();

        for r in &regions[1..] {
            if r.chr == current.chr && r.start <= current.end {
                // Overlapping or adjacent -- extend
                current.end = current.end.max(r.end);
            } else {
                merged.push(Region {
                    chr: current.chr.clone(),
                    start: current.start,
                    end: current.end,
                    rest: None,
                });
                current = r.clone();
            }
        }
        merged.push(Region {
            chr: current.chr,
            start: current.start,
            end: current.end,
            rest: None,
        });

        RegionSet::from(merged)
    }

    fn setdiff(&self, other: &RegionSet) -> RegionSet {
        let a = self.reduce();
        let b = other.reduce();

        // Group b regions by chromosome (already sorted from reduce)
        let mut b_by_chr: HashMap<String, Vec<&Region>> = HashMap::new();
        for r in &b.regions {
            b_by_chr.entry(r.chr.clone()).or_default().push(r);
        }

        let mut result: Vec<Region> = Vec::new();

        // Group a by chromosome and process with a cursor over b
        let mut a_chr_start = 0;
        while a_chr_start < a.regions.len() {
            let chr = &a.regions[a_chr_start].chr;
            let mut a_chr_end = a_chr_start;
            while a_chr_end < a.regions.len() && a.regions[a_chr_end].chr == *chr {
                a_chr_end += 1;
            }

            let empty_vec = vec![];
            let b_chr = b_by_chr.get(chr.as_str()).unwrap_or(&empty_vec);
            let mut b_idx = 0;

            for a_region in &a.regions[a_chr_start..a_chr_end] {
                // Advance b cursor past intervals that end before this region starts
                while b_idx < b_chr.len() && b_chr[b_idx].end <= a_region.start {
                    b_idx += 1;
                }

                let mut pos = a_region.start;
                let mut j = b_idx;

                while j < b_chr.len() && b_chr[j].start < a_region.end && pos < a_region.end {
                    if b_chr[j].start > pos {
                        // Gap before this subtraction interval
                        result.push(Region {
                            chr: chr.clone(),
                            start: pos,
                            end: b_chr[j].start,
                            rest: None,
                        });
                    }
                    pos = pos.max(b_chr[j].end);
                    j += 1;
                }

                // Remaining tail after all subtraction intervals
                if pos < a_region.end {
                    result.push(Region {
                        chr: chr.clone(),
                        start: pos,
                        end: a_region.end,
                        rest: None,
                    });
                }
            }

            a_chr_start = a_chr_end;
        }

        RegionSet::from(result)
    }

    fn pintersect(&self, other: &RegionSet) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .zip(other.regions.iter())
            .map(|(a, b)| {
                if a.chr != b.chr {
                    // Different chromosomes — zero-width interval at a's start
                    return Region {
                        chr: a.chr.clone(),
                        start: a.start,
                        end: a.start,
                        rest: None,
                    };
                }
                let start = a.start.max(b.start);
                let end = a.end.min(b.end);
                if start >= end {
                    // No overlap — zero-width interval
                    Region {
                        chr: a.chr.clone(),
                        start,
                        end: start,
                        rest: None,
                    }
                } else {
                    Region {
                        chr: a.chr.clone(),
                        start,
                        end,
                        rest: None,
                    }
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    fn concat(&self, other: &RegionSet) -> RegionSet {
        let mut regions = self.regions.clone();
        regions.extend(other.regions.iter().cloned());
        RegionSet::from(regions)
    }

    fn union(&self, other: &RegionSet) -> RegionSet {
        self.concat(other).reduce()
    }

    fn jaccard(&self, other: &RegionSet) -> f64 {
        let a_reduced = self.reduce();
        let b_reduced = other.reduce();
        let a_bp = a_reduced.nucleotides_length();
        let b_bp = b_reduced.nucleotides_length();
        let union_bp = a_reduced.union(&b_reduced).nucleotides_length();
        if union_bp == 0 {
            return 0.0;
        }
        let intersection_bp = a_bp + b_bp - union_bp;
        intersection_bp as f64 / union_bp as f64
    }

    fn coverage(&self, other: &RegionSet) -> f64 {
        let self_reduced = self.reduce();
        let self_bp = self_reduced.nucleotides_length();
        if self_bp == 0 {
            return 0.0;
        }
        let diff = self_reduced.setdiff(other);
        let diff_bp = diff.nucleotides_length();
        1.0 - (diff_bp as f64 / self_bp as f64)
    }

    fn overlap_coefficient(&self, other: &RegionSet) -> f64 {
        let a_bp = self.reduce().nucleotides_length();
        let b_bp = other.reduce().nucleotides_length();
        let min_bp = a_bp.min(b_bp);
        if min_bp == 0 {
            return 0.0;
        }
        let union_bp = self.union(other).nucleotides_length();
        let intersection_bp = a_bp + b_bp - union_bp;
        intersection_bp as f64 / min_bp as f64
    }

    fn shift(&self, offset: i64) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                let start = (r.start as i64 + offset).max(0) as u32;
                let end = (r.end as i64 + offset).max(start as i64) as u32;
                Region {
                    chr: r.chr.clone(),
                    start,
                    end,
                    rest: None,
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    fn flank(&self, width: u32, use_start: bool, both: bool) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                if both {
                    let anchor = if use_start { r.start } else { r.end };
                    Region {
                        chr: r.chr.clone(),
                        start: anchor.saturating_sub(width),
                        end: anchor.saturating_add(width),
                        rest: None,
                    }
                } else if use_start {
                    Region {
                        chr: r.chr.clone(),
                        start: r.start.saturating_sub(width),
                        end: r.start,
                        rest: None,
                    }
                } else {
                    Region {
                        chr: r.chr.clone(),
                        start: r.end,
                        end: r.end.saturating_add(width),
                        rest: None,
                    }
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    fn resize(&self, width: u32, fix: &str) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| match fix {
                "end" => Region {
                    chr: r.chr.clone(),
                    start: r.end.saturating_sub(width),
                    end: r.end,
                    rest: None,
                },
                "center" => {
                    let mid = r.start + (r.end - r.start) / 2;
                    let half = width / 2;
                    Region {
                        chr: r.chr.clone(),
                        start: mid.saturating_sub(half),
                        end: mid.saturating_sub(half).saturating_add(width),
                        rest: None,
                    }
                }
                _ => Region {
                    // "start" or default
                    chr: r.chr.clone(),
                    start: r.start,
                    end: r.start.saturating_add(width),
                    rest: None,
                },
            })
            .collect();
        RegionSet::from(regions)
    }

    fn narrow(&self, start: Option<u32>, end: Option<u32>, width: Option<u32>) -> RegionSet {
        let regions: Vec<Region> = self
            .regions
            .iter()
            .map(|r| {
                let region_width = r.end - r.start;
                // Parameters are 1-based relative positions within the region
                let (rel_start, rel_end) = match (start, end, width) {
                    (Some(s), Some(e), None) => (s.saturating_sub(1), e),
                    (Some(s), None, Some(w)) => (s.saturating_sub(1), s.saturating_sub(1) + w),
                    (None, Some(e), Some(w)) => (e.saturating_sub(w), e),
                    (None, None, Some(w)) => (0, w),
                    _ => (0, region_width),
                };
                let abs_start = r.start + rel_start.min(region_width);
                let abs_end = r.start + rel_end.min(region_width);
                Region {
                    chr: r.chr.clone(),
                    start: abs_start,
                    end: abs_end.max(abs_start),
                    rest: None,
                }
            })
            .collect();
        RegionSet::from(regions)
    }

    fn disjoin(&self) -> RegionSet {
        if self.regions.is_empty() {
            return RegionSet::from(Vec::<Region>::new());
        }

        let sorted = SortedRegionSet::new(RegionSet::from(self.regions.clone()));
        let regions = &sorted.0.regions;

        let mut result: Vec<Region> = Vec::new();

        // Process per chromosome
        let mut i = 0;
        while i < regions.len() {
            let chr = &regions[i].chr;
            let mut chr_end = i;
            while chr_end < regions.len() && regions[chr_end].chr == *chr {
                chr_end += 1;
            }

            // Collect all boundary events for this chromosome
            let mut events: Vec<u32> = Vec::with_capacity((chr_end - i) * 2);
            for r in &regions[i..chr_end] {
                events.push(r.start);
                events.push(r.end);
            }
            events.sort_unstable();
            events.dedup();

            // Walk the sorted, deduplicated boundary events directly.
            // Between consecutive events, emit a disjoint piece if
            // any input region covers that span.
            //
            // We use a sweep-line with a coverage counter: +1 at each
            // region start, -1 at each region end. Since events[] is
            // already sorted and deduped, we just need to know the
            // coverage change at each event position.
            //
            // Build a map: position → net coverage delta
            let mut delta_map: std::collections::HashMap<u32, i32> =
                std::collections::HashMap::with_capacity((chr_end - i) * 2);
            for r in &regions[i..chr_end] {
                *delta_map.entry(r.start).or_insert(0) += 1;
                *delta_map.entry(r.end).or_insert(0) -= 1;
            }

            // Walk events in order; emit a piece between consecutive
            // events whenever coverage > 0.
            let mut coverage: i32 = 0;
            for ei in 0..events.len() {
                if let Some(&d) = delta_map.get(&events[ei]) {
                    coverage += d;
                }
                // Emit interval [events[ei], events[ei+1]) if covered
                if coverage > 0 && ei + 1 < events.len() {
                    result.push(Region {
                        chr: chr.clone(),
                        start: events[ei],
                        end: events[ei + 1],
                        rest: None,
                    });
                }
            }

            i = chr_end;
        }

        RegionSet::from(result)
    }

    fn gaps(&self, chrom_sizes: &HashMap<String, u32>) -> RegionSet {
        let reduced = self.reduce();

        // Group reduced regions by chromosome so we can emit per-chrom gaps
        // and also detect chromosomes with zero regions (full-chrom gaps).
        let mut by_chr: HashMap<&str, Vec<&Region>> = HashMap::new();
        for r in &reduced.regions {
            // Skip chromosomes we don't have a size for — can't bound trailing gaps,
            // and including leading gaps for unknown-size chromosomes is misleading.
            if chrom_sizes.contains_key(&r.chr) {
                by_chr.entry(r.chr.as_str()).or_default().push(r);
            }
        }

        let mut result: Vec<Region> = Vec::new();

        // Emit gaps for every chromosome named in chrom_sizes, not just those
        // present in the input — this way a chromosome with zero regions
        // contributes a full-chromosome gap, matching bedtools complement.
        for (chr_name, &chrom_size) in chrom_sizes.iter() {
            if chrom_size == 0 {
                continue;
            }

            match by_chr.get(chr_name.as_str()) {
                None => {
                    // No regions on this chromosome — whole chromosome is a gap.
                    result.push(Region {
                        chr: chr_name.clone(),
                        start: 0,
                        end: chrom_size,
                        rest: None,
                    });
                }
                Some(regions) => {
                    // Leading gap from 0 to the first region's start.
                    if regions[0].start > 0 {
                        // Leading gap is clipped to chrom_size as a safety net:
                        // if the first region starts past chrom_size (assembly
                        // mismatch) we still produce a valid [0, chrom_size) gap.
                        let lead_end = regions[0].start.min(chrom_size);
                        result.push(Region {
                            chr: chr_name.clone(),
                            start: 0,
                            end: lead_end,
                            rest: None,
                        });
                    }

                    // Inter-region gaps.
                    for pair in regions.windows(2) {
                        let gap_start = pair[0].end;
                        let gap_end = pair[1].start;
                        if gap_start < gap_end {
                            // Clip both bounds to chrom_size so the whole emitted
                            // gap lies within the chromosome.
                            let cs = chrom_size;
                            let clipped_start = gap_start.min(cs);
                            let clipped_end = gap_end.min(cs);
                            if clipped_start < clipped_end {
                                result.push(Region {
                                    chr: chr_name.clone(),
                                    start: clipped_start,
                                    end: clipped_end,
                                    rest: None,
                                });
                            }
                        }
                    }

                    // Trailing gap from last region's end to chrom_size.
                    let last_end = regions[regions.len() - 1].end;
                    if last_end < chrom_size {
                        result.push(Region {
                            chr: chr_name.clone(),
                            start: last_end,
                            end: chrom_size,
                            rest: None,
                        });
                    }
                }
            }
        }

        // Karyotypic chromosome ordering so output is stable across runs.
        result.sort_by(|a, b| {
            crate::utils::chrom_karyotype_key(&a.chr)
                .cmp(&crate::utils::chrom_karyotype_key(&b.chr))
                .then(a.start.cmp(&b.start))
        });

        RegionSet::from(result)
    }

    fn intersect(&self, other: &RegionSet) -> RegionSet {
        let a = self.reduce();
        let b = other.reduce();

        // Group b regions by chromosome
        let mut b_by_chr: HashMap<String, Vec<&Region>> = HashMap::new();
        for r in &b.regions {
            b_by_chr.entry(r.chr.clone()).or_default().push(r);
        }

        let mut result: Vec<Region> = Vec::new();

        // Process a by chromosome groups
        let mut a_i = 0;
        while a_i < a.regions.len() {
            let chr = &a.regions[a_i].chr;
            let mut a_end = a_i;
            while a_end < a.regions.len() && a.regions[a_end].chr == *chr {
                a_end += 1;
            }

            if let Some(b_chr) = b_by_chr.get(chr.as_str()) {
                let mut b_idx = 0;
                for a_region in &a.regions[a_i..a_end] {
                    // Advance b past regions that end before a starts
                    while b_idx < b_chr.len() && b_chr[b_idx].end <= a_region.start {
                        b_idx += 1;
                    }
                    // Collect overlaps
                    let mut j = b_idx;
                    while j < b_chr.len() && b_chr[j].start < a_region.end {
                        let start = a_region.start.max(b_chr[j].start);
                        let end = a_region.end.min(b_chr[j].end);
                        if start < end {
                            result.push(Region {
                                chr: chr.clone(),
                                start,
                                end,
                                rest: None,
                            });
                        }
                        j += 1;
                    }
                }
            }

            a_i = a_end;
        }

        RegionSet::from(result)
    }

    fn intersect_all(&self, other: &RegionSet) -> RegionSet {
        if self.regions.is_empty() || other.regions.is_empty() {
            return RegionSet::from(Vec::<Region>::new());
        }

        let index = other.clone().into_multi_chrom_overlapper(OverlapperType::AIList);
        let mut result: Vec<Region> = Vec::new();

        for a in &self.regions {
            for hit in index.find_overlaps_for_region(&a.chr, a.start, a.end) {
                let start = a.start.max(hit.start);
                let end = a.end.min(hit.end);
                if start < end {
                    result.push(Region {
                        chr: a.chr.clone(),
                        start,
                        end,
                        rest: None,
                    });
                }
            }
        }

        RegionSet::from(result)
    }

    fn subtract(&self, other: &RegionSet) -> RegionSet {
        self.setdiff(other)
    }

    fn closest(&self, other: &RegionSet) -> Vec<(usize, usize, i64)> {
        if other.regions.is_empty() {
            return Vec::new();
        }

        // Group other by chromosome, keeping original indices, sorted by start
        let mut other_by_chr: HashMap<String, Vec<(usize, &Region)>> = HashMap::new();
        for (idx, r) in other.regions.iter().enumerate() {
            other_by_chr
                .entry(r.chr.clone())
                .or_default()
                .push((idx, r));
        }
        for v in other_by_chr.values_mut() {
            v.sort_by_key(|(_, r)| r.start);
        }

        // Precompute max region width per chromosome for left-side early termination
        let mut max_width_by_chr: HashMap<String, u32> = HashMap::new();
        for (chr, candidates) in &other_by_chr {
            let max_w = candidates.iter().map(|(_, r)| r.end - r.start).max().unwrap_or(0);
            max_width_by_chr.insert(chr.clone(), max_w);
        }

        let mut result: Vec<(usize, usize, i64)> = Vec::new();

        for (self_idx, a) in self.regions.iter().enumerate() {
            let Some(candidates) = other_by_chr.get(&a.chr) else {
                continue; // skip regions on chromosomes absent in other
            };

            // Binary search for insertion point based on a.start
            let ins = candidates
                .binary_search_by_key(&a.start, |(_, r)| r.start)
                .unwrap_or_else(|x| x);

            let gap_dist = |a_reg: &Region, b_reg: &Region| -> i64 {
                if a_reg.start < b_reg.end && b_reg.start < a_reg.end {
                    0i64
                } else if b_reg.end <= a_reg.start {
                    (a_reg.start as i64) - (b_reg.end as i64)
                } else {
                    (b_reg.start as i64) - (a_reg.end as i64)
                }
            };

            let max_width = *max_width_by_chr.get(&a.chr).unwrap_or(&0) as i64;

            let mut best_other_idx = 0usize;
            let mut best_dist = i64::MAX;

            let mut left_done = ins == 0;
            let mut right_done = ins >= candidates.len();
            let mut li = if ins > 0 { ins - 1 } else { 0 };
            let mut ri = ins;

            while !left_done || !right_done {
                if !right_done {
                    let (other_idx, b) = candidates[ri];
                    let dist = gap_dist(a, b);
                    if dist.abs() < best_dist.abs() {
                        best_dist = dist;
                        best_other_idx = other_idx;
                    }
                    if best_dist == 0 { break; }
                    ri += 1;
                    if ri >= candidates.len() || (b.start as i64 - a.end as i64 > best_dist.abs()) {
                        right_done = true;
                    }
                }

                if !left_done {
                    let (other_idx, b) = candidates[li];
                    let dist = gap_dist(a, b);
                    if dist.abs() < best_dist.abs() {
                        best_dist = dist;
                        best_other_idx = other_idx;
                    }
                    if best_dist == 0 { break; }
                    if li == 0 || (a.start as i64 - b.start as i64 > best_dist.abs() + max_width) {
                        left_done = true;
                    } else {
                        li -= 1;
                    }
                }
            }

            result.push((self_idx, best_other_idx, best_dist));
        }

        result
    }

    fn cluster(&self, max_gap: u32) -> Vec<u32> {
        if self.regions.is_empty() {
            return Vec::new();
        }

        let n = self.regions.len();
        let mut result = vec![0u32; n];

        // Create sorted indices
        let mut sorted_indices: Vec<usize> = (0..n).collect();
        sorted_indices.sort_by(|&i, &j| {
            self.regions[i]
                .chr
                .cmp(&self.regions[j].chr)
                .then(self.regions[i].start.cmp(&self.regions[j].start))
                .then(self.regions[i].end.cmp(&self.regions[j].end))
        });

        let mut cluster_id: u32 = 0;
        let mut cluster_end = self.regions[sorted_indices[0]].end;
        let mut current_chr = &self.regions[sorted_indices[0]].chr;
        result[sorted_indices[0]] = cluster_id;

        for &idx in &sorted_indices[1..] {
            let r = &self.regions[idx];
            if r.chr != *current_chr
                || r.start > cluster_end.saturating_add(max_gap)
            {
                // New cluster
                cluster_id += 1;
                cluster_end = r.end;
                current_chr = &r.chr;
            } else {
                // Same cluster - extend end
                cluster_end = cluster_end.max(r.end);
            }
            result[idx] = cluster_id;
        }

        result
    }
}

// ── Strand-aware operations on StrandedRegionSet ─────────────────────────

impl StrandedRegionSet {
    /// Clip regions to chromosome boundaries, preserving strand information.
    pub fn trim(&self, chrom_sizes: &HashMap<String, u32>) -> StrandedRegionSet {
        let mut regions = Vec::new();
        let mut strands = Vec::new();
        for (r, s) in self.inner.regions.iter().zip(self.strands.iter()) {
            if let Some(&chrom_size) = chrom_sizes.get(&r.chr) {
                let start = r.start.min(chrom_size);
                let end = r.end.min(chrom_size);
                if start < end {
                    regions.push(Region {
                        chr: r.chr.clone(),
                        start,
                        end,
                        rest: None,
                    });
                    strands.push(*s);
                }
                // Drop zero-width regions (start == end after clipping)
            } else {
                // Chromosome not in chrom_sizes — keep region as-is (no trimming)
                regions.push(r.clone());
                strands.push(*s);
            }
        }
        StrandedRegionSet {
            inner: RegionSet::from(regions),
            strands,
        }
    }

    /// Strand-aware promoter computation. Returns an unstranded `RegionSet`.
    ///
    /// - Plus / Unstranded: `[start - upstream, start + downstream)`
    /// - Minus: `[end - downstream, end + upstream)`
    pub fn promoters(&self, upstream: u32, downstream: u32) -> RegionSet {
        self.promoters_stranded(upstream, downstream).inner
    }

    /// Like `promoters()` but preserves strand information, so a subsequent
    /// strand-aware `reduce()` merges only same-strand promoters.
    pub fn promoters_stranded(&self, upstream: u32, downstream: u32) -> StrandedRegionSet {
        let regions: Vec<Region> = self
            .inner
            .regions
            .iter()
            .zip(self.strands.iter())
            .map(|(r, strand)| match strand {
                Strand::Minus => Region {
                    chr: r.chr.clone(),
                    start: r.end.saturating_sub(downstream),
                    end: r.end.saturating_add(upstream),
                    rest: None,
                },
                _ => Region {
                    chr: r.chr.clone(),
                    start: r.start.saturating_sub(upstream),
                    end: r.start.saturating_add(downstream),
                    rest: None,
                },
            })
            .collect();
        StrandedRegionSet {
            inner: RegionSet::from(regions),
            strands: self.strands.clone(),
        }
    }

    /// Strand-aware reduce: merge overlapping/adjacent intervals only within
    /// the same (chr, strand) group.
    pub fn reduce(&self) -> StrandedRegionSet {
        if self.inner.regions.is_empty() {
            return StrandedRegionSet::new(
                RegionSet::from(Vec::<Region>::new()),
                Vec::new(),
            );
        }

        // Build (region, strand) pairs and sort by (chr, strand, start)
        let mut pairs: Vec<(&Region, &Strand)> = self
            .inner
            .regions
            .iter()
            .zip(self.strands.iter())
            .collect();
        pairs.sort_by(|(a, sa), (b, sb)| {
            a.chr
                .cmp(&b.chr)
                .then(strand_ord(**sa).cmp(&strand_ord(**sb)))
                .then(a.start.cmp(&b.start))
        });

        let mut merged_regions: Vec<Region> = Vec::new();
        let mut merged_strands: Vec<Strand> = Vec::new();

        let (mut cur_r, mut cur_s) = (pairs[0].0.clone(), *pairs[0].1);
        for &(r, s) in &pairs[1..] {
            if r.chr == cur_r.chr && *s == cur_s && r.start <= cur_r.end {
                cur_r.end = cur_r.end.max(r.end);
            } else {
                merged_regions.push(Region {
                    chr: cur_r.chr.clone(),
                    start: cur_r.start,
                    end: cur_r.end,
                    rest: None,
                });
                merged_strands.push(cur_s);
                cur_r = r.clone();
                cur_s = *s;
            }
        }
        merged_regions.push(Region {
            chr: cur_r.chr,
            start: cur_r.start,
            end: cur_r.end,
            rest: None,
        });
        merged_strands.push(cur_s);

        StrandedRegionSet::new(RegionSet::from(merged_regions), merged_strands)
    }

    /// Strand-aware setdiff: subtract `other` from `self`, matching only
    /// within the same (chr, strand) group.
    pub fn setdiff(&self, other: &StrandedRegionSet) -> StrandedRegionSet {
        let a = self.reduce();
        let b = other.reduce();

        // Group b by (chr, strand)
        let mut b_map: HashMap<(String, Strand), Vec<&Region>> = HashMap::new();
        for (r, s) in b.inner.regions.iter().zip(b.strands.iter()) {
            b_map
                .entry((r.chr.clone(), *s))
                .or_default()
                .push(r);
        }

        let mut result_regions: Vec<Region> = Vec::new();
        let mut result_strands: Vec<Strand> = Vec::new();

        // Process a by (chr, strand) groups
        let mut i = 0;
        while i < a.inner.regions.len() {
            let chr = &a.inner.regions[i].chr;
            let strand = a.strands[i];

            // Find the extent of this (chr, strand) group
            let mut j = i;
            while j < a.inner.regions.len()
                && a.inner.regions[j].chr == *chr
                && a.strands[j] == strand
            {
                j += 1;
            }

            let empty_vec = vec![];
            let b_chr_strand = b_map
                .get(&(chr.clone(), strand))
                .unwrap_or(&empty_vec);
            let mut b_idx = 0;

            for a_region in &a.inner.regions[i..j] {
                while b_idx < b_chr_strand.len()
                    && b_chr_strand[b_idx].end <= a_region.start
                {
                    b_idx += 1;
                }

                let mut pos = a_region.start;
                let mut k = b_idx;

                while k < b_chr_strand.len()
                    && b_chr_strand[k].start < a_region.end
                    && pos < a_region.end
                {
                    if b_chr_strand[k].start > pos {
                        result_regions.push(Region {
                            chr: chr.clone(),
                            start: pos,
                            end: b_chr_strand[k].start,
                            rest: None,
                        });
                        result_strands.push(strand);
                    }
                    pos = pos.max(b_chr_strand[k].end);
                    k += 1;
                }

                if pos < a_region.end {
                    result_regions.push(Region {
                        chr: chr.clone(),
                        start: pos,
                        end: a_region.end,
                        rest: None,
                    });
                    result_strands.push(strand);
                }
            }

            i = j;
        }

        StrandedRegionSet::new(RegionSet::from(result_regions), result_strands)
    }
}

/// Ordering key for Strand so sorting groups by (chr, strand, start).
fn strand_ord(s: Strand) -> u8 {
    match s {
        Strand::Plus => 0,
        Strand::Minus => 1,
        Strand::Unstranded => 2,
    }
}

// =========================================================================
// Batch pairwise Jaccard
// =========================================================================

/// Compute intersection bp between two *already-reduced, sorted* RegionSets
/// using a linear merge-walk. O(n+m) with no allocations.
fn merge_intersection_bp(a: &RegionSet, b: &RegionSet) -> u64 {
    let a_regions = &a.regions;
    let b_regions = &b.regions;
    let mut ai = 0usize;
    let mut bi = 0usize;
    let mut intersection: u64 = 0;

    while ai < a_regions.len() && bi < b_regions.len() {
        let ra = &a_regions[ai];
        let rb = &b_regions[bi];

        // Skip ahead if on different chromosomes
        match ra.chr.cmp(&rb.chr) {
            std::cmp::Ordering::Less => {
                ai += 1;
                continue;
            }
            std::cmp::Ordering::Greater => {
                bi += 1;
                continue;
            }
            std::cmp::Ordering::Equal => {}
        }

        // Same chromosome — check overlap
        if ra.start < rb.end && rb.start < ra.end {
            let ovl_start = ra.start.max(rb.start);
            let ovl_end = ra.end.min(rb.end);
            intersection += (ovl_end - ovl_start) as u64;
        }

        // Advance whichever interval ends first
        if ra.end <= rb.end {
            ai += 1;
        } else {
            bi += 1;
        }
    }

    intersection
}

/// Pairwise Jaccard similarity using merge-walk intersection (Algo B).
///
/// Pre-reduces each set once, then computes all N*(N-1)/2 pairs using a linear
/// merge-walk. Returns a flat Vec<f64> of length N*N in row-major order
/// (symmetric matrix with 1.0 on the diagonal).
pub fn pairwise_jaccard(sets: &[RegionSet]) -> Vec<f64> {
    let n = sets.len();
    let mut matrix = vec![0.0f64; n * n];

    // Pre-reduce all sets and cache their bp counts
    let reduced: Vec<RegionSet> = sets.iter().map(|s| s.reduce()).collect();
    let bp: Vec<u64> = reduced.iter().map(|s| s.nucleotides_length() as u64).collect();

    // Diagonal
    for i in 0..n {
        matrix[i * n + i] = 1.0;
    }

    // Upper triangle
    for i in 0..n {
        for j in (i + 1)..n {
            let intersection = merge_intersection_bp(&reduced[i], &reduced[j]);
            let union = bp[i] + bp[j] - intersection;
            let jac = if union == 0 { 0.0 } else { intersection as f64 / union as f64 };
            matrix[i * n + j] = jac;
            matrix[j * n + i] = jac;
        }
    }

    matrix
}

// --- Indexed operations on RegionSetList ---
//
// These let callers operate on pairs by index without cloning full RegionSets
// across an FFI boundary (wasm, Python, R).

use gtars_core::models::RegionSetList;

/// Indexed pair operations on a RegionSetList.
pub trait RegionSetListOps {
    fn pintersect_at(&self, i: usize, j: usize) -> Option<RegionSet>;
    fn pintersect_count(&self, i: usize, j: usize) -> Option<u32>;
    fn jaccard_at(&self, i: usize, j: usize) -> Option<f64>;
    fn union_at(&self, i: usize, j: usize) -> Option<RegionSet>;
    fn setdiff_at(&self, i: usize, j: usize) -> Option<RegionSet>;
    fn region_count(&self, i: usize) -> Option<u32>;
    fn union_except(&self, skip: usize) -> Option<RegionSet>;
    /// Compute union-of-all and all N union-except results in O(n) unions
    /// using prefix/suffix arrays. Returns (full_union, vec_of_union_except).
    fn bulk_union_except(&self) -> Option<(RegionSet, Vec<RegionSet>)>;
    /// Fold all sets into a single union.
    fn union_all(&self) -> Option<RegionSet>;
    /// Fold all sets into a single intersection.
    fn intersect_all(&self) -> Option<RegionSet>;
}

impl RegionSetListOps for RegionSetList {
    fn pintersect_at(&self, i: usize, j: usize) -> Option<RegionSet> {
        let a = self.get(i)?;
        let b = self.get(j)?;
        Some(a.pintersect(b))
    }

    fn pintersect_count(&self, i: usize, j: usize) -> Option<u32> {
        self.pintersect_at(i, j).map(|rs| rs.len() as u32)
    }

    fn jaccard_at(&self, i: usize, j: usize) -> Option<f64> {
        let a = self.get(i)?;
        let b = self.get(j)?;
        Some(a.jaccard(b))
    }

    fn union_at(&self, i: usize, j: usize) -> Option<RegionSet> {
        let a = self.get(i)?;
        let b = self.get(j)?;
        Some(a.union(b))
    }

    fn setdiff_at(&self, i: usize, j: usize) -> Option<RegionSet> {
        let a = self.get(i)?;
        let b = self.get(j)?;
        Some(a.setdiff(b))
    }

    fn region_count(&self, i: usize) -> Option<u32> {
        self.get(i).map(|rs| rs.len() as u32)
    }

    fn union_except(&self, skip: usize) -> Option<RegionSet> {
        let n = self.len();
        if n < 2 || skip >= n { return None; }
        let first = if skip == 0 { 1 } else { 0 };
        let mut acc = self.get(first)?.clone();
        for k in (first + 1)..n {
            if k == skip { continue; }
            if let Some(other) = self.get(k) {
                acc = acc.union(other);
            }
        }
        Some(acc)
    }

    fn bulk_union_except(&self) -> Option<(RegionSet, Vec<RegionSet>)> {
        let n = self.len();
        if n < 2 { return None; }

        // prefix[i] = union(set[0]..=set[i])
        let mut prefix = Vec::with_capacity(n);
        prefix.push(self.get(0)?.clone());
        for i in 1..n {
            let prev = &prefix[i - 1];
            prefix.push(prev.union(self.get(i)?));
        }

        // suffix[i] = union(set[i]..=set[n-1]), built incrementally from right
        let mut suffix = vec![None; n];
        suffix[n - 1] = Some(self.get(n - 1)?.clone());
        for i in (0..n - 1).rev() {
            suffix[i] = Some(self.get(i)?.union(suffix[i + 1].as_ref().unwrap()));
        }
        let suffix: Vec<RegionSet> = suffix.into_iter().map(|s| s.unwrap()).collect();

        let full_union = prefix[n - 1].clone();

        // union_except[i] = union(prefix[i-1], suffix[i+1])
        let mut results = Vec::with_capacity(n);
        for i in 0..n {
            let except = match (i > 0, i < n - 1) {
                (false, true) => suffix[1].clone(),
                (true, false) => prefix[i - 1].clone(),
                (true, true) => prefix[i - 1].union(&suffix[i + 1]),
                (false, false) => unreachable!(), // n >= 2
            };
            results.push(except);
        }

        Some((full_union, results))
    }

    fn union_all(&self) -> Option<RegionSet> {
        let n = self.len();
        if n == 0 { return None; }
        let mut acc = self.get(0)?.clone();
        for i in 1..n {
            if let Some(other) = self.get(i) {
                acc = acc.union(other);
            }
        }
        Some(acc)
    }

    fn intersect_all(&self) -> Option<RegionSet> {
        let n = self.len();
        if n == 0 { return None; }
        let mut acc = self.get(0)?.clone();
        for i in 1..n {
            if let Some(other) = self.get(i) {
                acc = acc.intersect(other);
            }
        }
        Some(acc)
    }
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

    fn make_region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    fn make_regionset(regions: Vec<(&str, u32, u32)>) -> RegionSet {
        let regions: Vec<Region> = regions
            .into_iter()
            .map(|(chr, start, end)| make_region(chr, start, end))
            .collect();
        RegionSet::from(regions)
    }

    // ── trim tests ──────────────────────────────────────────────────────

    #[rstest]
    fn test_trim_clamps_past_chrom_end() {
        let rs = make_regionset(vec![("chr1", 90, 150)]);
        let chrom_sizes: HashMap<String, u32> =
            [("chr1".to_string(), 100)].into_iter().collect();
        let trimmed = rs.trim(&chrom_sizes);
        assert_eq!(trimmed.regions.len(), 1);
        assert_eq!(trimmed.regions[0].start, 90);
        assert_eq!(trimmed.regions[0].end, 100);
    }

    #[rstest]
    fn test_trim_drops_unknown_chrom() {
        let rs = make_regionset(vec![("chrX", 0, 50), ("chr1", 10, 20)]);
        let chrom_sizes: HashMap<String, u32> =
            [("chr1".to_string(), 100)].into_iter().collect();
        let trimmed = rs.trim(&chrom_sizes);
        assert_eq!(trimmed.regions.len(), 1);
        assert_eq!(trimmed.regions[0].chr, "chr1");
    }

    #[rstest]
    fn test_trim_keeps_zero_width_after_clamp() {
        // chr1:100-200 clamped to chrom_size=100 -> chr1:100-100 (zero-width, kept)
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let chrom_sizes: HashMap<String, u32> =
            [("chr1".to_string(), 100)].into_iter().collect();
        let trimmed = rs.trim(&chrom_sizes);
        assert_eq!(trimmed.regions.len(), 1);
        assert_eq!(trimmed.regions[0].start, 100);
        assert_eq!(trimmed.regions[0].end, 100);
    }

    #[rstest]
    fn test_trim_drops_inverted_after_clamp() {
        // chr1:150-100 is inverted (start > end), should be dropped
        let rs = make_regionset(vec![("chr1", 150, 100)]);
        let chrom_sizes: HashMap<String, u32> =
            [("chr1".to_string(), 200)].into_iter().collect();
        let trimmed = rs.trim(&chrom_sizes);
        assert_eq!(trimmed.regions.len(), 0);
    }

    #[rstest]
    fn test_trim_no_change_when_within_bounds() {
        let rs = make_regionset(vec![("chr1", 10, 50)]);
        let chrom_sizes: HashMap<String, u32> =
            [("chr1".to_string(), 100)].into_iter().collect();
        let trimmed = rs.trim(&chrom_sizes);
        assert_eq!(trimmed.regions.len(), 1);
        assert_eq!(trimmed.regions[0].start, 10);
        assert_eq!(trimmed.regions[0].end, 50);
    }

    // ── promoters tests ─────────────────────────────────────────────────

    #[rstest]
    fn test_promoters_standard() {
        let rs = make_regionset(vec![("chr1", 1000, 2000)]);
        let result = rs.promoters(500, 200);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0].start, 500);
        assert_eq!(result.regions[0].end, 1200);
    }

    #[rstest]
    fn test_promoters_saturating_sub_at_zero() {
        let rs = make_regionset(vec![("chr1", 100, 500)]);
        let result = rs.promoters(200, 50);
        assert_eq!(result.regions[0].start, 0); // 100 - 200 saturates to 0
        assert_eq!(result.regions[0].end, 150);
    }

    #[rstest]
    fn test_promoters_at_origin() {
        let rs = make_regionset(vec![("chr1", 0, 100)]);
        let result = rs.promoters(500, 200);
        assert_eq!(result.regions[0].start, 0);
        assert_eq!(result.regions[0].end, 200);
    }

    // ── reduce tests ────────────────────────────────────────────────────

    #[rstest]
    fn test_reduce_dummy_bed() {
        // dummy.bed has: chr1:2-6, chr1:4-7, chr1:5-9, chr1:7-12
        // All overlap -> single merged region chr1:2-12
        let path = get_test_path("dummy.bed");
        let rs = RegionSet::try_from(path.to_str().unwrap()).unwrap();
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 1);
        assert_eq!(reduced.regions[0].chr, "chr1");
        assert_eq!(reduced.regions[0].start, 2);
        assert_eq!(reduced.regions[0].end, 12);
    }

    #[rstest]
    fn test_reduce_non_overlapping() {
        let rs = make_regionset(vec![("chr1", 0, 5), ("chr1", 10, 15), ("chr1", 20, 25)]);
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 3);
    }

    #[rstest]
    fn test_reduce_adjacent_merged() {
        // Adjacent intervals (end == start of next) should merge
        let rs = make_regionset(vec![("chr1", 0, 10), ("chr1", 10, 20)]);
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 1);
        assert_eq!(reduced.regions[0].start, 0);
        assert_eq!(reduced.regions[0].end, 20);
    }

    #[rstest]
    fn test_reduce_multi_chrom() {
        let rs = make_regionset(vec![
            ("chr1", 0, 10),
            ("chr1", 5, 15),
            ("chr2", 0, 10),
            ("chr2", 20, 30),
        ]);
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 3);
        // chr1:0-15, chr2:0-10, chr2:20-30
        assert_eq!(reduced.regions[0], make_region("chr1", 0, 15));
        assert_eq!(reduced.regions[1], make_region("chr2", 0, 10));
        assert_eq!(reduced.regions[2], make_region("chr2", 20, 30));
    }

    #[rstest]
    fn test_reduce_lexicographic_chrom_order() {
        // BED convention uses lexicographic chromosome order
        let rs = make_regionset(vec![
            ("chr10", 0, 10),
            ("chr2", 0, 10),
            ("chr1", 0, 10),
            ("chrX", 0, 10),
            ("chrM", 0, 10),
            ("chrY", 0, 10),
        ]);
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 6);
        assert_eq!(reduced.regions[0].chr, "chr1");
        assert_eq!(reduced.regions[1].chr, "chr10");
        assert_eq!(reduced.regions[2].chr, "chr2");
        assert_eq!(reduced.regions[3].chr, "chrM");
        assert_eq!(reduced.regions[4].chr, "chrX");
        assert_eq!(reduced.regions[5].chr, "chrY");
    }

    #[rstest]
    fn test_reduce_lexicographic_order_merges_correctly() {
        // Overlapping regions on chr2 and chr10 — lexicographic sort
        // groups them correctly before merging
        let rs = make_regionset(vec![
            ("chr10", 5, 15),
            ("chr2", 0, 10),
            ("chr10", 0, 8),
            ("chr2", 5, 20),
        ]);
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 2);
        assert_eq!(reduced.regions[0], make_region("chr10", 0, 15));
        assert_eq!(reduced.regions[1], make_region("chr2", 0, 20));
    }

    #[rstest]
    fn test_reduce_empty() {
        let rs = RegionSet::from(Vec::<Region>::new());
        let reduced = rs.reduce();
        assert_eq!(reduced.regions.len(), 0);
    }

    // ── setdiff tests ───────────────────────────────────────────────────

    #[rstest]
    fn test_setdiff_middle_subtraction() {
        // [0,10) - [3,7) = [0,3) + [7,10)
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 3, 7)]);
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 0, 3));
        assert_eq!(result.regions[1], make_region("chr1", 7, 10));
    }

    #[rstest]
    fn test_setdiff_complete_subtraction() {
        let a = make_regionset(vec![("chr1", 3, 7)]);
        let b = make_regionset(vec![("chr1", 0, 10)]);
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 0);
    }

    #[rstest]
    fn test_setdiff_no_overlap() {
        let a = make_regionset(vec![("chr1", 0, 5)]);
        let b = make_regionset(vec![("chr1", 10, 20)]);
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 0, 5));
    }

    #[rstest]
    fn test_setdiff_multi_chrom() {
        let a = make_regionset(vec![("chr1", 0, 10), ("chr2", 0, 10)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 0, 5));
        assert_eq!(result.regions[1], make_region("chr2", 0, 10));
    }

    #[rstest]
    fn test_setdiff_from_bed_files() {
        // dummy.bed: chr1:2-6, chr1:4-7, chr1:5-9, chr1:7-12 -> reduces to chr1:2-12
        // dummy_b.bed: chr1:3-5, chr1:8-10
        // setdiff: chr1:2-12 - chr1:3-5 - chr1:8-10 = chr1:2-3 + chr1:5-8 + chr1:10-12
        let path_a = get_test_path("dummy.bed");
        let path_b = get_test_path("dummy_b.bed");
        let a = RegionSet::try_from(path_a.to_str().unwrap()).unwrap();
        let b = RegionSet::try_from(path_b.to_str().unwrap()).unwrap();
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 3);
        assert_eq!(result.regions[0], make_region("chr1", 2, 3));
        assert_eq!(result.regions[1], make_region("chr1", 5, 8));
        assert_eq!(result.regions[2], make_region("chr1", 10, 12));
    }

    #[rstest]
    fn test_setdiff_multiple_subtractions() {
        // [0,20) - [2,5) - [8,12) - [15,18) = [0,2) + [5,8) + [12,15) + [18,20)
        let a = make_regionset(vec![("chr1", 0, 20)]);
        let b = make_regionset(vec![("chr1", 2, 5), ("chr1", 8, 12), ("chr1", 15, 18)]);
        let result = a.setdiff(&b);
        assert_eq!(result.regions.len(), 4);
        assert_eq!(result.regions[0], make_region("chr1", 0, 2));
        assert_eq!(result.regions[1], make_region("chr1", 5, 8));
        assert_eq!(result.regions[2], make_region("chr1", 12, 15));
        assert_eq!(result.regions[3], make_region("chr1", 18, 20));
    }

    // ── pintersect tests ────────────────────────────────────────────────

    #[rstest]
    fn test_pintersect_overlapping_pair() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let result = a.pintersect(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 5, 10));
    }

    #[rstest]
    fn test_pintersect_no_overlap_zero_width() {
        let a = make_regionset(vec![("chr1", 0, 5)]);
        let b = make_regionset(vec![("chr1", 10, 20)]);
        let result = a.pintersect(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0].start, 10);
        assert_eq!(result.regions[0].end, 10);
    }

    #[rstest]
    fn test_pintersect_chrom_mismatch_zero_width() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr2", 0, 10)]);
        let result = a.pintersect(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0].chr, "chr1");
        assert_eq!(result.regions[0].start, 0);
        assert_eq!(result.regions[0].end, 0);
    }

    #[rstest]
    fn test_pintersect_multiple_pairs() {
        let a = make_regionset(vec![
            ("chr1", 0, 10),
            ("chr1", 20, 30),
            ("chr2", 0, 100),
        ]);
        let b = make_regionset(vec![
            ("chr1", 5, 15),
            ("chr1", 25, 35),
            ("chr2", 50, 60),
        ]);
        let result = a.pintersect(&b);
        assert_eq!(result.regions.len(), 3);
        assert_eq!(result.regions[0], make_region("chr1", 5, 10));
        assert_eq!(result.regions[1], make_region("chr1", 25, 30));
        assert_eq!(result.regions[2], make_region("chr2", 50, 60));
    }

    #[rstest]
    fn test_pintersect_contained() {
        // b fully inside a
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 30, 70)]);
        let result = a.pintersect(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 30, 70));
    }

    // ── StrandedRegionSet promoters tests ────────────────────────────

    fn make_stranded(regions: Vec<(&str, u32, u32)>, strands: Vec<Strand>) -> StrandedRegionSet {
        StrandedRegionSet::new(make_regionset(regions), strands)
    }

    #[rstest]
    fn test_stranded_promoters_plus() {
        // Plus-strand gene at [1000, 5000): promoter 100bp upstream of start
        let srs = make_stranded(vec![("chr1", 1000, 5000)], vec![Strand::Plus]);
        let result = srs.promoters(100, 0);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 900, 1000));
    }

    #[rstest]
    fn test_stranded_promoters_minus() {
        // Minus-strand gene at [3000, 8000): promoter 100bp upstream of end
        let srs = make_stranded(vec![("chr2", 3000, 8000)], vec![Strand::Minus]);
        let result = srs.promoters(100, 0);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr2", 8000, 8100));
    }

    #[rstest]
    fn test_stranded_promoters_unstranded_matches_plus() {
        // Unstranded should behave like plus-strand (same as original behavior)
        let srs = make_stranded(vec![("chr1", 1000, 5000)], vec![Strand::Unstranded]);
        let result = srs.promoters(100, 0);
        assert_eq!(result.regions[0], make_region("chr1", 900, 1000));
    }

    #[rstest]
    fn test_stranded_promoters_with_downstream() {
        // Minus-strand with both upstream and downstream
        let srs = make_stranded(vec![("chr1", 3000, 8000)], vec![Strand::Minus]);
        let result = srs.promoters(200, 50);
        // [end - downstream, end + upstream) = [7950, 8200)
        assert_eq!(result.regions[0], make_region("chr1", 7950, 8200));
    }

    #[rstest]
    fn test_stranded_promoters_saturating_at_zero() {
        // Plus-strand gene near origin
        let srs = make_stranded(vec![("chr1", 50, 500)], vec![Strand::Plus]);
        let result = srs.promoters(200, 0);
        assert_eq!(result.regions[0].start, 0); // saturates
        assert_eq!(result.regions[0].end, 50);
    }

    // ── StrandedRegionSet reduce tests ───────────────────────────────

    #[rstest]
    fn test_stranded_reduce_keeps_opposite_strand_separate() {
        // Two overlapping regions on opposite strands should NOT merge
        let srs = make_stranded(
            vec![("chr1", 100, 300), ("chr1", 200, 400)],
            vec![Strand::Plus, Strand::Minus],
        );
        let reduced = srs.reduce();
        assert_eq!(reduced.len(), 2);
    }

    #[rstest]
    fn test_stranded_reduce_merges_same_strand() {
        // Two overlapping regions on the same strand should merge
        let srs = make_stranded(
            vec![("chr1", 100, 300), ("chr1", 200, 400)],
            vec![Strand::Plus, Strand::Plus],
        );
        let reduced = srs.reduce();
        assert_eq!(reduced.len(), 1);
        assert_eq!(reduced.inner.regions[0], make_region("chr1", 100, 400));
    }

    #[rstest]
    fn test_stranded_reduce_empty() {
        let srs = StrandedRegionSet::new(
            RegionSet::from(Vec::<Region>::new()),
            Vec::new(),
        );
        let reduced = srs.reduce();
        assert_eq!(reduced.len(), 0);
    }

    // ── StrandedRegionSet setdiff tests ──────────────────────────────

    #[rstest]
    fn test_stranded_setdiff_same_strand_subtracts() {
        let a = make_stranded(vec![("chr1", 0, 100)], vec![Strand::Plus]);
        let b = make_stranded(vec![("chr1", 30, 70)], vec![Strand::Plus]);
        let result = a.setdiff(&b);
        assert_eq!(result.len(), 2);
        assert_eq!(result.inner.regions[0], make_region("chr1", 0, 30));
        assert_eq!(result.inner.regions[1], make_region("chr1", 70, 100));
    }

    #[rstest]
    fn test_stranded_setdiff_different_strand_no_subtraction() {
        // Opposite strand: should NOT subtract
        let a = make_stranded(vec![("chr1", 0, 100)], vec![Strand::Plus]);
        let b = make_stranded(vec![("chr1", 30, 70)], vec![Strand::Minus]);
        let result = a.setdiff(&b);
        assert_eq!(result.len(), 1);
        assert_eq!(result.inner.regions[0], make_region("chr1", 0, 100));
    }

    // ── concat tests ────────────────────────────────────────────────────

    #[rstest]
    fn test_concat() {
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let result = a.concat(&b);
        assert_eq!(result.regions.len(), 3);
    }

    #[rstest]
    fn test_concat_preserves_regions_and_order() {
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let result = a.concat(&b);
        assert_eq!(result.regions[0], make_region("chr1", 0, 10));
        assert_eq!(result.regions[1], make_region("chr1", 20, 30));
        assert_eq!(result.regions[2], make_region("chr1", 5, 15));
    }

    #[rstest]
    fn test_concat_both_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = RegionSet::from(Vec::<Region>::new());
        let result = a.concat(&b);
        assert_eq!(result.regions.len(), 0);
    }

    #[rstest]
    fn test_concat_one_empty() {
        let nonempty = make_regionset(vec![("chr1", 0, 10)]);
        let empty = RegionSet::from(Vec::<Region>::new());
        // empty + nonempty
        let r1 = empty.concat(&nonempty);
        assert_eq!(r1.regions.len(), 1);
        assert_eq!(r1.regions[0], make_region("chr1", 0, 10));
        // nonempty + empty
        let r2 = nonempty.concat(&empty);
        assert_eq!(r2.regions.len(), 1);
        assert_eq!(r2.regions[0], make_region("chr1", 0, 10));
    }

    #[rstest]
    fn test_concat_multi_chromosome() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr2", 50, 60)]);
        let result = a.concat(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 0, 10));
        assert_eq!(result.regions[1], make_region("chr2", 50, 60));
    }

    #[rstest]
    fn test_concat_does_not_merge_overlapping() {
        let a = make_regionset(vec![("chr1", 0, 20)]);
        let b = make_regionset(vec![("chr1", 10, 30)]);
        let result = a.concat(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 0, 20));
        assert_eq!(result.regions[1], make_region("chr1", 10, 30));
    }

    #[rstest]
    fn test_concat_preserves_duplicates() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 0, 10)]);
        let result = a.concat(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], result.regions[1]);
    }

    // ── union tests ─────────────────────────────────────────────────────

    #[rstest]
    fn test_union_from_bed_files() {
        // dummy.bed reduces to chr1:2-12, dummy_b.bed has chr1:3-5, chr1:8-10
        // All fall within chr1:2-12, so union = chr1:2-12
        let path_a = get_test_path("dummy.bed");
        let path_b = get_test_path("dummy_b.bed");
        let a = RegionSet::try_from(path_a.to_str().unwrap()).unwrap();
        let b = RegionSet::try_from(path_b.to_str().unwrap()).unwrap();
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 2, 12));
    }

    #[rstest]
    fn test_union_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 20, 30)]);
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 0, 10));
        assert_eq!(result.regions[1], make_region("chr1", 20, 30));
    }

    #[rstest]
    fn test_union_overlapping() {
        let a = make_regionset(vec![("chr1", 0, 15)]);
        let b = make_regionset(vec![("chr1", 10, 25)]);
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 0, 25));
    }

    #[rstest]
    fn test_union_adjacent() {
        // [0,10) and [10,20) are adjacent; reduce merges them (10 <= 10)
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 10, 20)]);
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 0, 20));
    }

    #[rstest]
    fn test_union_one_gap() {
        // [0,10) and [11,20) have a 1bp gap at position 10 -- stay separate
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 11, 20)]);
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 0, 10));
        assert_eq!(result.regions[1], make_region("chr1", 11, 20));
    }

    #[rstest]
    fn test_union_both_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = RegionSet::from(Vec::<Region>::new());
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 0);
    }

    #[rstest]
    fn test_union_one_empty() {
        let nonempty = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]);
        let empty = RegionSet::from(Vec::<Region>::new());
        let result = nonempty.union(&empty);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 0, 10));
        assert_eq!(result.regions[1], make_region("chr1", 20, 30));
    }

    #[rstest]
    fn test_union_multi_chromosome() {
        // chr1 and chr2 intervals never merge; output sorted lexicographically
        let a = make_regionset(vec![("chr2", 0, 10)]);
        let b = make_regionset(vec![("chr1", 0, 10)]);
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 0, 10));
        assert_eq!(result.regions[1], make_region("chr2", 0, 10));
    }

    #[rstest]
    fn test_union_contained() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 20, 50)]);
        let result = a.union(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 0, 100));
    }

    #[rstest]
    fn test_union_identical() {
        // union(A, A) = A (idempotent)
        let a = make_regionset(vec![("chr1", 10, 20), ("chr1", 30, 40)]);
        let result = a.union(&a);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 10, 20));
        assert_eq!(result.regions[1], make_region("chr1", 30, 40));
    }

    // ── jaccard tests ───────────────────────────────────────────────────

    #[rstest]
    fn test_jaccard() {
        // dummy.bed reduces to chr1:2-12 (10bp)
        // dummy_b.bed reduces to chr1:3-5 + chr1:8-10 (4bp)
        // Union = chr1:2-12 (10bp). Intersection = 10 + 4 - 10 = 4bp.
        // Jaccard = 4/10 = 0.4
        let path_a = get_test_path("dummy.bed");
        let path_b = get_test_path("dummy_b.bed");
        let a = RegionSet::try_from(path_a.to_str().unwrap()).unwrap();
        let b = RegionSet::try_from(path_b.to_str().unwrap()).unwrap();
        let j = a.jaccard(&b);
        assert!((j - 0.4).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_identical() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        assert!((a.jaccard(&a) - 1.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 20, 30)]);
        assert!((a.jaccard(&b)).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = RegionSet::from(Vec::<Region>::new());
        assert!((a.jaccard(&b)).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_partial_overlap_inline() {
        // A=[0,10) 10bp, B=[5,15) 10bp, union=[0,15) 15bp, intersect=5bp
        // J = 5/15 = 1/3
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let j = a.jaccard(&b);
        assert!((j - 1.0 / 3.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_contained() {
        // A=[0,100) 100bp, B=[20,50) 30bp, union=[0,100) 100bp, intersect=30bp
        // J = 30/100 = 0.3
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 20, 50)]);
        let j = a.jaccard(&b);
        assert!((j - 0.3).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_multi_chromosome() {
        // A: chr1:[0,20) + chr2:[0,10) = 30bp
        // B: chr1:[10,30) + chr2:[5,15) = 30bp
        // union: chr1:[0,30) + chr2:[0,15) = 45bp
        // intersect: 30+30-45 = 15bp  (chr1:[10,20)=10bp + chr2:[5,10)=5bp)
        // J = 15/45 = 1/3
        let a = make_regionset(vec![("chr1", 0, 20), ("chr2", 0, 10)]);
        let b = make_regionset(vec![("chr1", 10, 30), ("chr2", 5, 15)]);
        let j = a.jaccard(&b);
        assert!((j - 1.0 / 3.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_different_chromosomes() {
        // Same coordinates but different chroms -> no shared positions
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr2", 0, 100)]);
        assert!((a.jaccard(&b)).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_symmetry() {
        let a = make_regionset(vec![("chr1", 0, 15), ("chr1", 30, 50)]);
        let b = make_regionset(vec![("chr1", 10, 40)]);
        assert!((a.jaccard(&b) - b.jaccard(&a)).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_overlapping_input() {
        // A has self-overlapping regions; reduces to [0,25) = 25bp
        // B = [5,20) = 15bp
        // union = [0,25) = 25bp, intersect = 25+15-25 = 15bp
        // J = 15/25 = 0.6
        let a = make_regionset(vec![("chr1", 0, 15), ("chr1", 10, 25)]);
        let b = make_regionset(vec![("chr1", 5, 20)]);
        let j = a.jaccard(&b);
        assert!((j - 0.6).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_one_empty() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = RegionSet::from(Vec::<Region>::new());
        assert!((a.jaccard(&b)).abs() < 1e-10);
    }

    #[rstest]
    fn test_jaccard_adjacent() {
        // [0,10) covers positions 0-9, [10,20) covers 10-19: no shared positions
        // union via reduce: [0,20) = 20bp, intersect = 10+10-20 = 0bp
        // J = 0/20 = 0.0
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 10, 20)]);
        assert!((a.jaccard(&b)).abs() < 1e-10);
    }

    // ── coverage tests ──────────────────────────────────────────────────

    #[rstest]
    fn test_coverage_identical() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        assert!((a.coverage(&a) - 1.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_coverage_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 20, 30)]);
        assert!(a.coverage(&b).abs() < 1e-10);
    }

    #[rstest]
    fn test_coverage_partial_overlap() {
        // A=[0,10) 10bp, B=[5,15) 10bp
        // setdiff(A, B) = [0,5) = 5bp, coverage = 1 - 5/10 = 0.5
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        assert!((a.coverage(&b) - 0.5).abs() < 1e-10);
    }

    #[rstest]
    fn test_coverage_subset() {
        // A is a subset of B => coverage(A, B) = 1.0
        let a = make_regionset(vec![("chr1", 20, 50)]);
        let b = make_regionset(vec![("chr1", 0, 100)]);
        assert!((a.coverage(&b) - 1.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_coverage_superset() {
        // A=[0,100) 100bp, B=[20,50) 30bp
        // setdiff(A, B) = [0,20) + [50,100) = 70bp, coverage = 1 - 70/100 = 0.3
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 20, 50)]);
        assert!((a.coverage(&b) - 0.3).abs() < 1e-10);
    }

    #[rstest]
    fn test_coverage_empty_self() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = make_regionset(vec![("chr1", 0, 100)]);
        assert!(a.coverage(&b).abs() < 1e-10);
    }

    #[rstest]
    fn test_coverage_asymmetry() {
        // A=[0,10), B=[5,15)
        // coverage(A,B) = 0.5, coverage(B,A) = 0.5 (symmetric in this case)
        // Use asymmetric example: A=[0,100), B=[20,50)
        // coverage(A,B) = 0.3, coverage(B,A) = 1.0
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 20, 50)]);
        let cov_ab = a.coverage(&b);
        let cov_ba = b.coverage(&a);
        assert!((cov_ab - 0.3).abs() < 1e-10);
        assert!((cov_ba - 1.0).abs() < 1e-10);
        assert!((cov_ab - cov_ba).abs() > 0.1); // they differ
    }

    // ── overlap_coefficient tests ──────────────────────────────────────

    #[rstest]
    fn test_overlap_coefficient_identical() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        assert!((a.overlap_coefficient(&a) - 1.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_overlap_coefficient_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 20, 30)]);
        assert!(a.overlap_coefficient(&b).abs() < 1e-10);
    }

    #[rstest]
    fn test_overlap_coefficient_partial_overlap() {
        // A=[0,10) 10bp, B=[5,15) 10bp
        // intersection = 5bp, min = 10bp, oc = 5/10 = 0.5
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        assert!((a.overlap_coefficient(&b) - 0.5).abs() < 1e-10);
    }

    #[rstest]
    fn test_overlap_coefficient_subset() {
        // A=[20,50) 30bp inside B=[0,100) 100bp
        // intersection = 30bp, min = 30bp, oc = 1.0
        let a = make_regionset(vec![("chr1", 20, 50)]);
        let b = make_regionset(vec![("chr1", 0, 100)]);
        assert!((a.overlap_coefficient(&b) - 1.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_overlap_coefficient_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = make_regionset(vec![("chr1", 0, 100)]);
        assert!(a.overlap_coefficient(&b).abs() < 1e-10);
    }

    #[rstest]
    fn test_overlap_coefficient_symmetry() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 20, 50)]);
        assert!((a.overlap_coefficient(&b) - b.overlap_coefficient(&a)).abs() < 1e-10);
    }

    #[rstest]
    fn test_overlap_coefficient_both_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = RegionSet::from(Vec::<Region>::new());
        assert!(a.overlap_coefficient(&b).abs() < 1e-10);
    }

    // ── intersect tests ────────────────────────────────────────────────

    #[rstest]
    fn test_intersect_overlapping() {
        let a = make_regionset(vec![("chr1", 100, 200)]);
        let b = make_regionset(vec![("chr1", 150, 250)]);
        let result = a.intersect(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 150, 200));
    }

    #[rstest]
    fn test_intersect_no_overlap() {
        let a = make_regionset(vec![("chr1", 100, 200)]);
        let b = make_regionset(vec![("chr1", 300, 400)]);
        let result = a.intersect(&b);
        assert_eq!(result.regions.len(), 0);
    }

    #[rstest]
    fn test_intersect_one_vs_multiple() {
        // One region in self overlaps multiple in other
        let a = make_regionset(vec![("chr1", 100, 300)]);
        let b = make_regionset(vec![
            ("chr1", 120, 150),
            ("chr1", 200, 250),
            ("chr1", 400, 500),
        ]);
        let result = a.intersect(&b);
        assert_eq!(result.regions.len(), 2);
        let mut sorted: Vec<_> = result.regions.clone();
        sorted.sort_by_key(|r| (r.start, r.end));
        assert_eq!(sorted[0], make_region("chr1", 120, 150));
        assert_eq!(sorted[1], make_region("chr1", 200, 250));
    }

    #[rstest]
    fn test_intersect_multi_chrom() {
        let a = make_regionset(vec![("chr1", 100, 200), ("chr2", 100, 200)]);
        let b = make_regionset(vec![("chr1", 150, 250), ("chr2", 150, 250)]);
        let result = a.intersect(&b);
        assert_eq!(result.regions.len(), 2);
        assert_eq!(result.regions[0], make_region("chr1", 150, 200));
        assert_eq!(result.regions[1], make_region("chr2", 150, 200));
    }

    #[rstest]
    fn test_intersect_empty_inputs() {
        let empty = RegionSet::from(Vec::<Region>::new());
        let a = make_regionset(vec![("chr1", 100, 200)]);
        assert_eq!(a.intersect(&empty).regions.len(), 0);
        assert_eq!(empty.intersect(&a).regions.len(), 0);
    }

    #[rstest]
    fn test_intersect_contained() {
        // b fully inside a
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 30, 70)]);
        let result = a.intersect(&b);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 30, 70));
    }

    // ── subtract tests (alias) ─────────────────────────────────────────

    #[rstest]
    fn test_subtract_same_as_setdiff() {
        let a = make_regionset(vec![("chr1", 0, 20)]);
        let b = make_regionset(vec![("chr1", 5, 15)]);
        let setdiff_result = a.setdiff(&b);
        let subtract_result = a.subtract(&b);
        assert_eq!(setdiff_result.regions, subtract_result.regions);
    }

    // ── closest tests ──────────────────────────────────────────────────

    #[rstest]
    fn test_closest_overlapping() {
        let a = make_regionset(vec![("chr1", 100, 200)]);
        let b = make_regionset(vec![("chr1", 150, 250)]);
        let result = a.closest(&b);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], (0, 0, 0));
    }

    #[rstest]
    fn test_closest_upstream() {
        let a = make_regionset(vec![("chr1", 200, 300)]);
        let b = make_regionset(vec![("chr1", 100, 150)]);
        let result = a.closest(&b);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], (0, 0, 50)); // distance = 200 - 150 = 50
    }

    #[rstest]
    fn test_closest_downstream() {
        let a = make_regionset(vec![("chr1", 100, 150)]);
        let b = make_regionset(vec![("chr1", 200, 300)]);
        let result = a.closest(&b);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], (0, 0, 50)); // distance = 200 - 150 = 50
    }

    #[rstest]
    fn test_closest_choose_nearer() {
        let a = make_regionset(vec![("chr1", 100, 110)]);
        let b = make_regionset(vec![("chr1", 50, 60), ("chr1", 115, 120)]);
        let result = a.closest(&b);
        assert_eq!(result.len(), 1);
        // distance to (50,60) = 100-60 = 40, distance to (115,120) = 115-110 = 5
        assert_eq!(result[0].1, 1); // closer is index 1
        assert_eq!(result[0].2, 5);
    }

    #[rstest]
    fn test_closest_multi_chrom_no_cross() {
        let a = make_regionset(vec![("chr1", 100, 200), ("chr2", 100, 200)]);
        let b = make_regionset(vec![("chr1", 300, 400)]);
        let result = a.closest(&b);
        // chr2 region has no candidates in b, should be omitted
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, 0); // only the chr1 region
    }

    #[rstest]
    fn test_closest_absent_chrom_omitted() {
        let a = make_regionset(vec![("chr2", 100, 200)]);
        let b = make_regionset(vec![("chr1", 100, 200)]);
        let result = a.closest(&b);
        assert_eq!(result.len(), 0);
    }

    #[rstest]
    fn test_closest_nearest_far_left_of_insertion() {
        let a = make_regionset(vec![("chr1", 1000, 1010)]);
        let b = make_regionset(vec![
            ("chr1", 0, 999),     // dist = 1 (nearest)
            ("chr1", 1, 2),       // dist = 998
            ("chr1", 3, 4),       // dist = 996
            ("chr1", 5, 6),       // dist = 994
            ("chr1", 1020, 2000), // dist = 10
        ]);
        let result = a.closest(&b);
        assert_eq!(result[0].2, 1);
    }

    #[rstest]
    fn test_closest_nearest_far_right_of_insertion() {
        let a = make_regionset(vec![("chr1", 500, 510)]);
        let b = make_regionset(vec![
            ("chr1", 0, 490),     // dist = 10
            ("chr1", 520, 521),   // dist = 10
            ("chr1", 522, 523),   // dist = 12
            ("chr1", 524, 525),   // dist = 14
            ("chr1", 526, 527),   // dist = 16
            ("chr1", 528, 529),   // dist = 18
            ("chr1", 530, 531),   // dist = 20
            ("chr1", 508, 509),   // dist = 0 (overlap!)
        ]);
        // After sorting by start, the overlapping region should be found
        let result = a.closest(&b);
        assert_eq!(result[0].2, 0);
    }

    #[rstest]
    fn test_closest_variable_density() {
        let a = make_regionset(vec![("chr1", 5000, 5010)]);
        let mut regions: Vec<(&str, u32, u32)> = Vec::new();
        // 20 tiny regions clustered at the start
        for i in 0..20 {
            regions.push(("chr1", i * 2, i * 2 + 1));
        }
        // One region that's actually close
        regions.push(("chr1", 4990, 4999)); // dist = 1
        // Some regions after
        regions.push(("chr1", 6000, 6001));
        let b = make_regionset(regions);
        let result = a.closest(&b);
        assert_eq!(result[0].2, 1);
    }

    #[rstest]
    fn test_closest_single_candidate() {
        let a = make_regionset(vec![("chr1", 100, 200)]);
        let b = make_regionset(vec![("chr1", 500, 600)]);
        let result = a.closest(&b);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].2, 300);
    }

    #[rstest]
    fn test_closest_wide_region_far_left() {
        let a = make_regionset(vec![("chr1", 1000, 1010)]);
        let b = make_regionset(vec![
            ("chr1", 10, 998),    // dist = 2 (wide region, close end)
            ("chr1", 500, 501),   // dist = 499
            ("chr1", 600, 601),   // dist = 399
            ("chr1", 700, 701),   // dist = 299
            ("chr1", 1050, 1060), // dist = 40
        ]);
        let result = a.closest(&b);
        assert_eq!(result[0].2, 2);
    }

    // ── cluster tests ──────────────────────────────────────────────────

    #[rstest]
    fn test_cluster_non_overlapping_separate() {
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30), ("chr1", 40, 50)]);
        let ids = a.cluster(0);
        assert_eq!(ids.len(), 3);
        assert_ne!(ids[0], ids[1]);
        assert_ne!(ids[1], ids[2]);
    }

    #[rstest]
    fn test_cluster_overlapping_same() {
        let a = make_regionset(vec![("chr1", 0, 15), ("chr1", 10, 25)]);
        let ids = a.cluster(0);
        assert_eq!(ids[0], ids[1]);
    }

    #[rstest]
    fn test_cluster_adjacent_same() {
        // end == start of next => gap is 0 => same cluster with max_gap=0
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 10, 20)]);
        let ids = a.cluster(0);
        assert_eq!(ids[0], ids[1]);
    }

    #[rstest]
    fn test_cluster_within_max_gap() {
        // gap of 5bp, with max_gap=5 they should cluster
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 15, 25)]);
        let ids = a.cluster(5);
        assert_eq!(ids[0], ids[1]);
    }

    #[rstest]
    fn test_cluster_multi_chrom() {
        let a = make_regionset(vec![("chr1", 0, 10), ("chr2", 0, 10)]);
        let ids = a.cluster(0);
        assert_ne!(ids[0], ids[1]);
    }

    #[rstest]
    fn test_cluster_single_region() {
        let a = make_regionset(vec![("chr1", 100, 200)]);
        let ids = a.cluster(0);
        assert_eq!(ids, vec![0]);
    }

    #[rstest]
    fn test_cluster_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let ids = a.cluster(0);
        assert!(ids.is_empty());
    }

    #[rstest]
    fn test_cluster_original_order() {
        // Regions given out of sorted order; cluster IDs should map back to original indices
        let a = make_regionset(vec![
            ("chr1", 100, 110), // idx 0 -> after sort: second
            ("chr1", 0, 10),    // idx 1 -> after sort: first
        ]);
        let ids = a.cluster(0);
        assert_eq!(ids.len(), 2);
        // They're on separate clusters. idx 1 (chr1:0-10) should get cluster 0,
        // idx 0 (chr1:100-110) should get cluster 1
        assert_eq!(ids[0], 1); // original idx 0 maps to cluster 1
        assert_eq!(ids[1], 0); // original idx 1 maps to cluster 0
    }

    // ── gaps tests ──────────────────────────────────────────────────────

    fn chrom_sizes(entries: &[(&str, u32)]) -> HashMap<String, u32> {
        entries
            .iter()
            .map(|(c, s)| (c.to_string(), *s))
            .collect()
    }

    #[test]
    fn test_gaps_basic() {
        // Three peaks on chr1 with gaps between them; leading + trailing
        // gaps also present.
        let rs = make_regionset(vec![
            ("chr1", 10, 20),
            ("chr1", 30, 40),
            ("chr1", 50, 60),
        ]);
        let cs = chrom_sizes(&[("chr1", 100)]);
        let result = rs.gaps(&cs);
        let gaps: Vec<(&str, u32, u32)> = result
            .regions
            .iter()
            .map(|r| (r.chr.as_str(), r.start, r.end))
            .collect();
        assert_eq!(
            gaps,
            vec![
                ("chr1", 0, 10),   // leading
                ("chr1", 20, 30),  // between peak 1 and 2
                ("chr1", 40, 50),  // between peak 2 and 3
                ("chr1", 60, 100), // trailing
            ]
        );
    }

    #[test]
    fn test_gaps_peak_at_origin_no_leading() {
        // First peak starts at 0 — no leading gap.
        let rs = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]);
        let cs = chrom_sizes(&[("chr1", 100)]);
        let result = rs.gaps(&cs);
        let gaps: Vec<(u32, u32)> = result
            .regions
            .iter()
            .map(|r| (r.start, r.end))
            .collect();
        assert_eq!(gaps, vec![(10, 20), (30, 100)]);
    }

    #[test]
    fn test_gaps_peak_at_chrom_end_no_trailing() {
        // Last peak ends at chrom_size — no trailing gap.
        let rs = make_regionset(vec![("chr1", 10, 20), ("chr1", 80, 100)]);
        let cs = chrom_sizes(&[("chr1", 100)]);
        let result = rs.gaps(&cs);
        let gaps: Vec<(u32, u32)> = result
            .regions
            .iter()
            .map(|r| (r.start, r.end))
            .collect();
        assert_eq!(gaps, vec![(0, 10), (20, 80)]);
    }

    #[test]
    fn test_gaps_peak_past_chrom_end_clipped() {
        // Last peak extends past chrom_size — should be clipped, no trailing.
        let rs = make_regionset(vec![("chr1", 10, 20), ("chr1", 80, 150)]);
        let cs = chrom_sizes(&[("chr1", 100)]);
        let result = rs.gaps(&cs);
        let gaps: Vec<(u32, u32)> = result
            .regions
            .iter()
            .map(|r| (r.start, r.end))
            .collect();
        assert_eq!(gaps, vec![(0, 10), (20, 80)]);
    }

    #[test]
    fn test_gaps_empty_regionset_populated_chrom_sizes() {
        // No regions, but chrom_sizes has entries — emit whole-chrom gaps.
        let rs = RegionSet::from(Vec::<Region>::new());
        let cs = chrom_sizes(&[("chr1", 100), ("chr2", 50)]);
        let result = rs.gaps(&cs);
        let mut gaps: Vec<(String, u32, u32)> = result
            .regions
            .iter()
            .map(|r| (r.chr.clone(), r.start, r.end))
            .collect();
        gaps.sort();
        assert_eq!(
            gaps,
            vec![
                ("chr1".to_string(), 0, 100),
                ("chr2".to_string(), 0, 50),
            ]
        );
    }

    #[test]
    fn test_gaps_empty_regionset_empty_chrom_sizes() {
        let rs = RegionSet::from(Vec::<Region>::new());
        let cs: HashMap<String, u32> = HashMap::new();
        let result = rs.gaps(&cs);
        assert!(result.regions.is_empty());
    }

    #[test]
    fn test_gaps_chromosome_not_in_chrom_sizes_skipped() {
        // Peak on chr2 with no chr2 entry in chrom_sizes — should be ignored.
        let rs = make_regionset(vec![("chr1", 10, 20), ("chr2", 5, 15)]);
        let cs = chrom_sizes(&[("chr1", 100)]);
        let result = rs.gaps(&cs);
        // Only chr1 gaps emitted.
        for r in &result.regions {
            assert_eq!(r.chr, "chr1");
        }
    }

    #[test]
    fn test_gaps_full_chrom_gap_for_unrepresented_chrom() {
        // chrom_sizes has chr2 but input has no chr2 peaks — emit whole chr2.
        let rs = make_regionset(vec![("chr1", 10, 20)]);
        let cs = chrom_sizes(&[("chr1", 100), ("chr2", 200)]);
        let result = rs.gaps(&cs);
        let chr2_gaps: Vec<(u32, u32)> = result
            .regions
            .iter()
            .filter(|r| r.chr == "chr2")
            .map(|r| (r.start, r.end))
            .collect();
        assert_eq!(chr2_gaps, vec![(0, 200)]);
    }

    #[test]
    fn test_gaps_overlapping_peaks_reduced() {
        // Overlapping peaks get merged by reduce() before gap computation.
        let rs = make_regionset(vec![
            ("chr1", 10, 30),
            ("chr1", 25, 40), // overlaps with previous
            ("chr1", 50, 60),
        ]);
        let cs = chrom_sizes(&[("chr1", 100)]);
        let result = rs.gaps(&cs);
        let gaps: Vec<(u32, u32)> = result
            .regions
            .iter()
            .map(|r| (r.start, r.end))
            .collect();
        // After reduce: [10,40], [50,60] → gaps: [0,10], [40,50], [60,100]
        assert_eq!(gaps, vec![(0, 10), (40, 50), (60, 100)]);
    }

    #[test]
    fn test_gaps_karyotypic_ordering() {
        // Output should be karyotypically ordered regardless of chrom_sizes insertion order.
        let rs = make_regionset(vec![
            ("chr2", 10, 20),
            ("chr1", 10, 20),
            ("chr10", 10, 20),
        ]);
        let cs = chrom_sizes(&[("chr10", 100), ("chr1", 100), ("chr2", 100)]);
        let result = rs.gaps(&cs);
        // First gap on each chr; collect chr names in order of appearance.
        let order: Vec<&str> = result
            .regions
            .iter()
            .map(|r| r.chr.as_str())
            .scan("", |prev, chr| {
                if chr != *prev {
                    *prev = chr;
                    Some(chr)
                } else {
                    Some("") // repeat, skip
                }
            })
            .filter(|s| !s.is_empty())
            .collect();
        assert_eq!(order, vec!["chr1", "chr2", "chr10"]);
    }

    #[test]
    fn test_resize_center_large_coords() {
        // Regression: coordinates near u32::MAX/2 should not overflow in midpoint calc
        let rs = make_regionset(vec![("chr1", 2_000_000_000, 2_100_000_000)]);
        let resized = rs.resize(1000, "center");
        assert_eq!(resized.regions.len(), 1);
        let r = &resized.regions[0];
        let expected_mid = 2_000_000_000u32 + (2_100_000_000u32 - 2_000_000_000u32) / 2;
        assert_eq!(r.start, expected_mid - 500);
        assert_eq!(r.end, expected_mid - 500 + 1000);
    }

    #[test]
    fn test_narrow_start_zero() {
        // Regression: narrow with start=0 should not panic from underflow
        let rs = make_regionset(vec![("chr1", 100, 200)]);
        let narrowed = rs.narrow(Some(0), Some(50), None);
        assert_eq!(narrowed.regions.len(), 1);
        // start=0, saturating_sub(1) = 0, so abs_start = 100 + 0 = 100
    }

    // ── pairwise_jaccard tests ──────────────────────────────────────────

    #[rstest]
    fn test_pairwise_jaccard_matches_single() {
        // Pairwise matrix entries should match individual jaccard() calls
        let a = make_regionset(vec![("chr1", 0, 20), ("chr2", 0, 10)]);
        let b = make_regionset(vec![("chr1", 10, 30), ("chr2", 5, 15)]);
        let c = make_regionset(vec![("chr1", 50, 60)]);

        let matrix = super::pairwise_jaccard(&[a.clone(), b.clone(), c.clone()]);
        assert_eq!(matrix.len(), 9); // 3x3

        // Compare each pair against the single jaccard
        let eps = 1e-10;
        assert!((matrix[0 * 3 + 1] - a.jaccard(&b)).abs() < eps);
        assert!((matrix[0 * 3 + 2] - a.jaccard(&c)).abs() < eps);
        assert!((matrix[1 * 3 + 2] - b.jaccard(&c)).abs() < eps);
    }

    #[rstest]
    fn test_pairwise_jaccard_diagonal() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 50, 150)]);
        let matrix = super::pairwise_jaccard(&[a, b]);
        assert!((matrix[0] - 1.0).abs() < 1e-10); // [0,0]
        assert!((matrix[3] - 1.0).abs() < 1e-10); // [1,1]
    }

    #[rstest]
    fn test_pairwise_jaccard_symmetric() {
        let a = make_regionset(vec![("chr1", 0, 15), ("chr1", 30, 50)]);
        let b = make_regionset(vec![("chr1", 10, 40)]);
        let c = make_regionset(vec![("chr2", 0, 100)]);
        let matrix = super::pairwise_jaccard(&[a, b, c]);
        let n = 3;
        for i in 0..n {
            for j in 0..n {
                assert!((matrix[i * n + j] - matrix[j * n + i]).abs() < 1e-10);
            }
        }
    }

    #[rstest]
    fn test_pairwise_jaccard_empty_sets() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = RegionSet::from(Vec::<Region>::new());
        let matrix = super::pairwise_jaccard(&[a, b]);
        // Empty sets: intersection=0, union=0 -> jaccard=0
        // Diagonal: 1.0 (by convention)
        assert!((matrix[0] - 1.0).abs() < 1e-10);
        assert!((matrix[3] - 1.0).abs() < 1e-10);
        assert!(matrix[1].abs() < 1e-10);
        assert!(matrix[2].abs() < 1e-10);
    }

    #[rstest]
    fn test_pairwise_jaccard_single_set() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let matrix = super::pairwise_jaccard(&[a]);
        assert_eq!(matrix.len(), 1);
        assert!((matrix[0] - 1.0).abs() < 1e-10);
    }

    #[rstest]
    fn test_pairwise_jaccard_no_overlap_different_chroms() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr2", 0, 100)]);
        let matrix = super::pairwise_jaccard(&[a, b]);
        assert!(matrix[1].abs() < 1e-10); // [0,1]
        assert!(matrix[2].abs() < 1e-10); // [1,0]
    }

    #[rstest]
    fn test_pairwise_jaccard_identical_sets() {
        let a = make_regionset(vec![("chr1", 0, 50), ("chr2", 10, 30)]);
        let b = make_regionset(vec![("chr1", 0, 50), ("chr2", 10, 30)]);
        let matrix = super::pairwise_jaccard(&[a, b]);
        assert!((matrix[1] - 1.0).abs() < 1e-10); // identical -> jaccard=1
    }

    #[rstest]
    fn test_pairwise_jaccard_overlapping_input_regions() {
        // Input regions that need reducing (self-overlapping)
        let a = make_regionset(vec![("chr1", 0, 15), ("chr1", 10, 25)]); // reduces to [0,25)
        let b = make_regionset(vec![("chr1", 5, 20)]);
        let matrix = super::pairwise_jaccard(&[a.clone(), b.clone()]);
        let expected = a.jaccard(&b);
        assert!((matrix[1] - expected).abs() < 1e-10);
    }

    #[rstest]
    fn test_merge_intersection_bp_basic() {
        // Already-reduced sets: a=[0,10) b=[5,15) -> intersection=[5,10)=5bp
        let a = make_regionset(vec![("chr1", 0, 10)]).reduce();
        let b = make_regionset(vec![("chr1", 5, 15)]).reduce();
        assert_eq!(super::merge_intersection_bp(&a, &b), 5);
    }

    #[rstest]
    fn test_merge_intersection_bp_multi_chrom() {
        let a = make_regionset(vec![("chr1", 0, 20), ("chr2", 0, 10)]).reduce();
        let b = make_regionset(vec![("chr1", 10, 30), ("chr2", 5, 15)]).reduce();
        // chr1: [10,20) = 10bp, chr2: [5,10) = 5bp -> total 15bp
        assert_eq!(super::merge_intersection_bp(&a, &b), 15);
    }

    #[rstest]
    fn test_merge_intersection_bp_no_overlap() {
        let a = make_regionset(vec![("chr1", 0, 10)]).reduce();
        let b = make_regionset(vec![("chr1", 20, 30)]).reduce();
        assert_eq!(super::merge_intersection_bp(&a, &b), 0);
    }

    #[rstest]
    fn test_merge_intersection_bp_adjacent() {
        // Half-open [0,10) and [10,20) do NOT overlap
        let a = make_regionset(vec![("chr1", 0, 10)]).reduce();
        let b = make_regionset(vec![("chr1", 10, 20)]).reduce();
        assert_eq!(super::merge_intersection_bp(&a, &b), 0);
    }

    #[rstest]
    fn test_merge_intersection_bp_contained() {
        // b is fully inside a
        let a = make_regionset(vec![("chr1", 0, 100)]).reduce();
        let b = make_regionset(vec![("chr1", 20, 50)]).reduce();
        assert_eq!(super::merge_intersection_bp(&a, &b), 30);
    }

    #[rstest]
    fn test_merge_intersection_bp_multiple_overlaps() {
        // a has two intervals, b overlaps both
        let a = make_regionset(vec![("chr1", 0, 10), ("chr1", 20, 30)]).reduce();
        let b = make_regionset(vec![("chr1", 5, 25)]).reduce();
        // overlap with [0,10): [5,10)=5bp, overlap with [20,30): [20,25)=5bp
        assert_eq!(super::merge_intersection_bp(&a, &b), 10);
    }

    // ── RegionSetListOps tests ─────────────────────────────────────────

    fn make_rsl(sets: Vec<RegionSet>) -> RegionSetList {
        RegionSetList::from(sets)
    }

    #[rstest]
    fn test_rsl_pintersect_count() {
        let a = make_regionset(vec![("chr1", 0, 100), ("chr1", 200, 300)]);
        let b = make_regionset(vec![("chr1", 50, 150), ("chr1", 250, 350)]);
        let rsl = make_rsl(vec![a, b]);
        // a has 2 regions, b has 2 regions, both overlap each other
        let count = rsl.pintersect_count(0, 1).unwrap();
        assert_eq!(count, 2); // both pairs overlap
    }

    #[rstest]
    fn test_rsl_pintersect_count_no_overlap() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr1", 100, 200)]);
        let rsl = make_rsl(vec![a, b]);
        // paired by index: chr1:0-10 vs chr1:100-200 → no genomic overlap,
        // but pintersect produces a zero-width region (start=end) per pair
        let count = rsl.pintersect_count(0, 1).unwrap();
        assert_eq!(count, 1); // zero-width region still counted
    }

    #[rstest]
    fn test_rsl_pintersect_count_out_of_bounds() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let rsl = make_rsl(vec![a]);
        assert!(rsl.pintersect_count(0, 5).is_none());
    }

    #[rstest]
    fn test_rsl_jaccard_at() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 0, 100)]);
        let rsl = make_rsl(vec![a, b]);
        let j = rsl.jaccard_at(0, 1).unwrap();
        assert!((j - 1.0).abs() < 1e-9, "identical sets should have jaccard=1.0");
    }

    #[rstest]
    fn test_rsl_jaccard_at_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 200, 300)]);
        let rsl = make_rsl(vec![a, b]);
        let j = rsl.jaccard_at(0, 1).unwrap();
        assert!((j - 0.0).abs() < 1e-9, "disjoint sets should have jaccard=0.0");
    }

    #[rstest]
    fn test_rsl_union_at() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 50, 150)]);
        let rsl = make_rsl(vec![a, b]);
        let u = rsl.union_at(0, 1).unwrap();
        assert_eq!(u.regions.len(), 1);
        assert_eq!(u.regions[0].start, 0);
        assert_eq!(u.regions[0].end, 150);
    }

    #[rstest]
    fn test_rsl_setdiff_at() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 50, 150)]);
        let rsl = make_rsl(vec![a, b]);
        let diff = rsl.setdiff_at(0, 1).unwrap();
        // a minus b: chr1:0-50
        assert_eq!(diff.regions.len(), 1);
        assert_eq!(diff.regions[0].start, 0);
        assert_eq!(diff.regions[0].end, 50);
    }

    #[rstest]
    fn test_rsl_region_count() {
        let a = make_regionset(vec![("chr1", 0, 100), ("chr1", 200, 300)]);
        let b = make_regionset(vec![("chr1", 50, 150)]);
        let rsl = make_rsl(vec![a, b]);
        assert_eq!(rsl.region_count(0).unwrap(), 2);
        assert_eq!(rsl.region_count(1).unwrap(), 1);
        assert!(rsl.region_count(5).is_none());
    }

    #[rstest]
    fn test_rsl_union_all() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 50, 200)]);
        let c = make_regionset(vec![("chr1", 150, 300)]);
        let rsl = make_rsl(vec![a, b, c]);
        let u = rsl.union_all().unwrap();
        assert_eq!(u.regions.len(), 1);
        assert_eq!(u.regions[0].start, 0);
        assert_eq!(u.regions[0].end, 300);
    }

    #[rstest]
    fn test_rsl_union_all_empty() {
        let rsl = make_rsl(vec![]);
        assert!(rsl.union_all().is_none());
    }

    #[rstest]
    fn test_rsl_union_all_single() {
        let a = make_regionset(vec![("chr1", 10, 50)]);
        let rsl = make_rsl(vec![a]);
        let u = rsl.union_all().unwrap();
        assert_eq!(u.regions.len(), 1);
        assert_eq!(u.regions[0].start, 10);
        assert_eq!(u.regions[0].end, 50);
    }

    #[rstest]
    fn test_rsl_intersect_all() {
        // Three overlapping sets — intersection is the region shared by all three
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 30, 200)]);
        let c = make_regionset(vec![("chr1", 60, 150)]);
        let rsl = make_rsl(vec![a, b, c]);
        let inter = rsl.intersect_all().unwrap();
        assert_eq!(inter.regions.len(), 1);
        assert_eq!(inter.regions[0].start, 60);
        assert_eq!(inter.regions[0].end, 100);
    }

    #[rstest]
    fn test_rsl_intersect_all_disjoint() {
        let a = make_regionset(vec![("chr1", 0, 50)]);
        let b = make_regionset(vec![("chr1", 100, 200)]);
        let rsl = make_rsl(vec![a, b]);
        let inter = rsl.intersect_all().unwrap();
        assert_eq!(inter.regions.len(), 0);
    }

    #[rstest]
    fn test_rsl_intersect_all_empty() {
        let rsl = make_rsl(vec![]);
        assert!(rsl.intersect_all().is_none());
    }

    #[rstest]
    fn test_rsl_intersect_all_different_sizes() {
        // This is the case where pintersect would give wrong results:
        // sets have different numbers of regions, but share genomic coverage
        let a = make_regionset(vec![("chr1", 0, 100), ("chr1", 200, 300)]);
        let b = make_regionset(vec![("chr1", 50, 250)]);
        let rsl = make_rsl(vec![a, b]);
        let inter = rsl.intersect_all().unwrap();
        // Shared coverage: [50,100) and [200,250)
        assert_eq!(inter.regions.len(), 2);
        assert_eq!(inter.regions[0].start, 50);
        assert_eq!(inter.regions[0].end, 100);
        assert_eq!(inter.regions[1].start, 200);
        assert_eq!(inter.regions[1].end, 250);
    }

    #[rstest]
    fn test_rsl_union_except() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 200, 300)]);
        let c = make_regionset(vec![("chr1", 400, 500)]);
        let rsl = make_rsl(vec![a, b, c]);
        // union_except(1) = union of sets 0 and 2 (skip set 1)
        let ue = rsl.union_except(1).unwrap();
        assert_eq!(ue.regions.len(), 2);
        assert_eq!(ue.regions[0].start, 0);
        assert_eq!(ue.regions[0].end, 100);
        assert_eq!(ue.regions[1].start, 400);
        assert_eq!(ue.regions[1].end, 500);
    }

    #[rstest]
    fn test_rsl_union_except_too_small() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let rsl = make_rsl(vec![a]);
        assert!(rsl.union_except(0).is_none());
    }

    #[rstest]
    fn test_rsl_bulk_union_except_n2() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 200, 300)]);
        let rsl = make_rsl(vec![a, b]);
        let (full_union, excepts) = rsl.bulk_union_except().unwrap();

        // Full union covers both regions
        assert_eq!(full_union.regions.len(), 2);

        // except[0] = union of everything except set 0 = set 1
        assert_eq!(excepts.len(), 2);
        assert_eq!(excepts[0].regions.len(), 1);
        assert_eq!(excepts[0].regions[0].start, 200);
        assert_eq!(excepts[0].regions[0].end, 300);

        // except[1] = union of everything except set 1 = set 0
        assert_eq!(excepts[1].regions.len(), 1);
        assert_eq!(excepts[1].regions[0].start, 0);
        assert_eq!(excepts[1].regions[0].end, 100);
    }

    #[rstest]
    fn test_rsl_bulk_union_except_n3() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let b = make_regionset(vec![("chr1", 200, 300)]);
        let c = make_regionset(vec![("chr1", 400, 500)]);
        let rsl = make_rsl(vec![a, b, c]);
        let (full_union, excepts) = rsl.bulk_union_except().unwrap();

        assert_eq!(full_union.regions.len(), 3);
        assert_eq!(excepts.len(), 3);

        // except[0] = union(b, c)
        assert_eq!(excepts[0].regions.len(), 2);
        assert_eq!(excepts[0].regions[0].start, 200);
        assert_eq!(excepts[0].regions[1].start, 400);

        // except[1] = union(a, c)
        assert_eq!(excepts[1].regions.len(), 2);
        assert_eq!(excepts[1].regions[0].start, 0);
        assert_eq!(excepts[1].regions[1].start, 400);

        // except[2] = union(a, b)
        assert_eq!(excepts[2].regions.len(), 2);
        assert_eq!(excepts[2].regions[0].start, 0);
        assert_eq!(excepts[2].regions[1].start, 200);
    }

    #[rstest]
    fn test_rsl_bulk_union_except_too_small() {
        let a = make_regionset(vec![("chr1", 0, 100)]);
        let rsl = make_rsl(vec![a]);
        assert!(rsl.bulk_union_except().is_none());

        let rsl_empty = make_rsl(vec![]);
        assert!(rsl_empty.bulk_union_except().is_none());
    }

    #[rstest]
    fn test_rsl_bulk_union_except_matches_union_except() {
        // Verify bulk algorithm produces same results as individual union_except calls
        let a = make_regionset(vec![("chr1", 0, 100), ("chr2", 50, 200)]);
        let b = make_regionset(vec![("chr1", 80, 180), ("chr2", 100, 300)]);
        let c = make_regionset(vec![("chr1", 150, 250)]);
        let d = make_regionset(vec![("chr2", 0, 150)]);
        let rsl = make_rsl(vec![a, b, c, d]);

        let (_, bulk_excepts) = rsl.bulk_union_except().unwrap();

        for i in 0..4 {
            let individual = rsl.union_except(i).unwrap();
            assert_eq!(
                bulk_excepts[i].regions.len(),
                individual.regions.len(),
                "region count mismatch at index {}",
                i
            );
            for (j, (bulk_r, indiv_r)) in bulk_excepts[i]
                .regions
                .iter()
                .zip(individual.regions.iter())
                .enumerate()
            {
                assert_eq!(
                    bulk_r.chr, indiv_r.chr,
                    "chr mismatch at except[{}][{}]", i, j
                );
                assert_eq!(
                    bulk_r.start, indiv_r.start,
                    "start mismatch at except[{}][{}]", i, j
                );
                assert_eq!(
                    bulk_r.end, indiv_r.end,
                    "end mismatch at except[{}][{}]", i, j
                );
            }
        }
    }
}
