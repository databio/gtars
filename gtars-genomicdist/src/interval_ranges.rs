//! Interval set algebra operations for genomic region sets.
//!
//! Provides GenomicRanges/IRanges-style operations: trim, promoters, reduce,
//! setdiff, and pintersect. All operations use 0-based half-open coordinates
//! (BED convention) and are strand-unaware.

use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};

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
            .filter_map(|(a, b)| {
                if a.chr != b.chr {
                    return None;
                }
                let start = a.start.max(b.start);
                let end = a.end.min(b.end);
                if start > end {
                    None
                } else {
                    Some(Region {
                        chr: a.chr.clone(),
                        start,
                        end,
                        rest: None,
                    })
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
        let a_bp = self.reduce().nucleotides_length();
        let b_bp = other.reduce().nucleotides_length();
        let union_bp = self.union(other).nucleotides_length();
        if union_bp == 0 {
            return 0.0;
        }
        let intersection_bp = a_bp + b_bp - union_bp;
        intersection_bp as f64 / union_bp as f64
    }
}

// ── Strand-aware operations on StrandedRegionSet ─────────────────────────

impl StrandedRegionSet {
    /// Strand-aware promoter computation. Returns an unstranded `RegionSet`.
    ///
    /// - Plus / Unstranded: `[start - upstream, start + downstream)`
    /// - Minus: `[end - downstream, end + upstream)`
    pub fn promoters(&self, upstream: u32, downstream: u32) -> RegionSet {
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
        RegionSet::from(regions)
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
    fn test_pintersect_no_overlap_dropped() {
        let a = make_regionset(vec![("chr1", 0, 5)]);
        let b = make_regionset(vec![("chr1", 10, 20)]);
        let result = a.pintersect(&b);
        assert_eq!(result.regions.len(), 0);
    }

    #[rstest]
    fn test_pintersect_chrom_mismatch_dropped() {
        let a = make_regionset(vec![("chr1", 0, 10)]);
        let b = make_regionset(vec![("chr2", 0, 10)]);
        let result = a.pintersect(&b);
        assert_eq!(result.regions.len(), 0);
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
}
