//! Index-dependent interval operations and strand-aware extensions.
//!
//! The `IntervalRanges` trait provides operations that require an overlap index
//! (e.g., AIList). Structural operations (trim, reduce, shift, etc.) and
//! sweep-line set operations (setdiff, intersect, union, jaccard, etc.) have
//! moved to `gtars_core::models::RegionSet` (inherent methods) and
//! `gtars_core::models::IntervalSetOps` (trait), respectively.

use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};
use gtars_overlaprs::OverlapperType;
use gtars_overlaprs::multi_chrom_overlapper::IntoMultiChromOverlapper;

use crate::models::{Strand, StrandedRegionSet};

/// Index-dependent interval operations on genomic region sets.
///
/// Contains only operations that require an overlap index (AIList).
/// For structural operations, use inherent methods on `RegionSet`.
/// For sweep-line set operations, use `IntervalSetOps` from `gtars_core`.
pub trait IntervalRanges {
    /// All-vs-all genomic intersection using an overlap index.
    ///
    /// For each pair of overlapping regions between `self` and `other`,
    /// compute intersection coordinates `[max(a.start, b.start), min(a.end, b.end))`.
    /// Returns a `RegionSet` of all intersection fragments.
    fn intersect_all(&self, other: &RegionSet) -> RegionSet;
}

impl IntervalRanges for RegionSet {
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

    // ── StrandedRegionSet promoters tests ────────────────────────────

    fn make_stranded(regions: Vec<(&str, u32, u32)>, strands: Vec<Strand>) -> StrandedRegionSet {
        StrandedRegionSet::new(make_regionset(regions), strands)
    }

    #[rstest]
    fn test_stranded_promoters_plus() {
        let srs = make_stranded(vec![("chr1", 1000, 5000)], vec![Strand::Plus]);
        let result = srs.promoters(100, 0);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr1", 900, 1000));
    }

    #[rstest]
    fn test_stranded_promoters_minus() {
        let srs = make_stranded(vec![("chr2", 3000, 8000)], vec![Strand::Minus]);
        let result = srs.promoters(100, 0);
        assert_eq!(result.regions.len(), 1);
        assert_eq!(result.regions[0], make_region("chr2", 8000, 8100));
    }

    #[rstest]
    fn test_stranded_promoters_unstranded_matches_plus() {
        let srs = make_stranded(vec![("chr1", 1000, 5000)], vec![Strand::Unstranded]);
        let result = srs.promoters(100, 0);
        assert_eq!(result.regions[0], make_region("chr1", 900, 1000));
    }

    #[rstest]
    fn test_stranded_promoters_with_downstream() {
        let srs = make_stranded(vec![("chr1", 3000, 8000)], vec![Strand::Minus]);
        let result = srs.promoters(200, 50);
        assert_eq!(result.regions[0], make_region("chr1", 7950, 8200));
    }

    #[rstest]
    fn test_stranded_promoters_saturating_at_zero() {
        let srs = make_stranded(vec![("chr1", 50, 500)], vec![Strand::Plus]);
        let result = srs.promoters(200, 0);
        assert_eq!(result.regions[0].start, 0);
        assert_eq!(result.regions[0].end, 50);
    }

    // ── StrandedRegionSet reduce tests ───────────────────────────────

    #[rstest]
    fn test_stranded_reduce_keeps_opposite_strand_separate() {
        let srs = make_stranded(
            vec![("chr1", 100, 300), ("chr1", 200, 400)],
            vec![Strand::Plus, Strand::Minus],
        );
        let reduced = srs.reduce();
        assert_eq!(reduced.len(), 2);
    }

    #[rstest]
    fn test_stranded_reduce_merges_same_strand() {
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
        let a = make_stranded(vec![("chr1", 0, 100)], vec![Strand::Plus]);
        let b = make_stranded(vec![("chr1", 30, 70)], vec![Strand::Minus]);
        let result = a.setdiff(&b);
        assert_eq!(result.len(), 1);
        assert_eq!(result.inner.regions[0], make_region("chr1", 0, 100));
    }

    // ── intersect_all tests ──────────────────────────────────────────

    #[rstest]
    fn test_intersect_all_basic() {
        let a = make_regionset(vec![("chr1", 0, 20)]);
        let b = make_regionset(vec![("chr1", 10, 30), ("chr1", 15, 25)]);
        let result = a.intersect_all(&b);
        // Should find two intersection fragments: [10,20) and [15,20)
        assert_eq!(result.regions.len(), 2);
    }

    #[rstest]
    fn test_intersect_all_empty() {
        let a = RegionSet::from(Vec::<Region>::new());
        let b = make_regionset(vec![("chr1", 0, 10)]);
        assert_eq!(a.intersect_all(&b).regions.len(), 0);
        assert_eq!(b.intersect_all(&a).regions.len(), 0);
    }
}
