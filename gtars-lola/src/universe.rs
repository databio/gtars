//! Universe validation, user set redefinition, and restricted universe construction.

use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};

/// Diagnostic report from universe appropriateness check.
#[derive(Debug, Clone)]
pub struct UniverseReport {
    /// Per-user-set diagnostics.
    pub user_set_reports: Vec<UserSetReport>,
}

/// Per-user-set diagnostics from universe check.
#[derive(Debug, Clone)]
pub struct UserSetReport {
    /// Index of the user set.
    pub user_set_index: usize,
    /// Total regions in the user set.
    pub total_regions: usize,
    /// Number of user regions that overlap at least one universe region.
    pub regions_in_universe: usize,
    /// Fraction of user regions overlapping the universe (0.0–1.0).
    pub coverage: f64,
    /// Number of user regions with many-to-many mappings
    /// (overlapping more than one universe region).
    pub many_to_many_count: usize,
    /// Warning messages.
    pub warnings: Vec<String>,
}

/// Check whether the universe is appropriate for the given user sets.
///
/// For each user set, reports:
/// - What fraction of user set regions overlap at least one universe region
/// - Whether there are many-to-many mappings (user region overlaps >1 universe region)
/// - Warnings if coverage is low
pub fn check_universe_appropriateness(
    user_sets: &[RegionSet],
    universe: &RegionSet,
) -> UniverseReport {
    // Build a per-chromosome sorted index of universe regions for efficient overlap queries
    let uni_index = build_chrom_index(universe);

    let mut user_set_reports = Vec::with_capacity(user_sets.len());

    for (us_idx, user_set) in user_sets.iter().enumerate() {
        let total = user_set.regions.len();
        let mut in_universe = 0usize;
        let mut many_to_many = 0usize;

        for region in &user_set.regions {
            let overlap_count = count_overlapping_regions(
                &region.chr,
                region.start,
                region.end,
                &uni_index,
            );
            if overlap_count > 0 {
                in_universe += 1;
            }
            if overlap_count > 1 {
                many_to_many += 1;
            }
        }

        let coverage = if total > 0 {
            in_universe as f64 / total as f64
        } else {
            0.0
        };

        let mut warnings = Vec::new();
        if coverage < 0.5 {
            warnings.push(format!(
                "User set {}: only {:.1}% of regions overlap the universe. \
                 Consider using a more appropriate universe.",
                us_idx,
                coverage * 100.0
            ));
        } else if coverage < 0.9 {
            warnings.push(format!(
                "User set {}: {:.1}% of regions overlap the universe. \
                 Some regions may not be represented.",
                us_idx,
                coverage * 100.0
            ));
        }
        if many_to_many > 0 {
            warnings.push(format!(
                "User set {}: {} regions overlap multiple universe regions (many-to-many). \
                 Consider using redefine_user_sets() to eliminate artifacts.",
                us_idx, many_to_many
            ));
        }

        user_set_reports.push(UserSetReport {
            user_set_index: us_idx,
            total_regions: total,
            regions_in_universe: in_universe,
            coverage,
            many_to_many_count: many_to_many,
            warnings,
        });
    }

    UniverseReport { user_set_reports }
}

/// Redefine user sets in terms of universe regions.
///
/// For each user set, finds all universe regions that overlap any user region
/// and returns those universe regions as the new user set. This eliminates
/// many-to-many mapping artifacts.
///
/// This is the Rust equivalent of R LOLA's `redefineUserSets()`.
pub fn redefine_user_sets(user_sets: &[RegionSet], universe: &RegionSet) -> Vec<RegionSet> {
    let uni_index = build_chrom_index(universe);

    user_sets
        .iter()
        .map(|user_set| {
            let mut new_regions: Vec<Region> = Vec::new();
            let mut seen = std::collections::HashSet::new();

            for user_region in &user_set.regions {
                let overlapping = find_overlapping_regions(
                    &user_region.chr,
                    user_region.start,
                    user_region.end,
                    &uni_index,
                );

                for uni_region in overlapping {
                    // Deduplicate by (chr, start, end)
                    let key = (uni_region.chr.clone(), uni_region.start, uni_region.end);
                    if seen.insert(key) {
                        new_regions.push(uni_region);
                    }
                }
            }

            RegionSet::from(new_regions)
        })
        .collect()
}

/// Build a restricted universe from user sets.
///
/// Merges all regions from all user sets, then merges overlapping intervals
/// to produce a disjoint set. Used for differential enrichment analysis.
pub fn build_restricted_universe(user_sets: &[RegionSet]) -> RegionSet {
    // Collect all regions from all user sets
    let mut all_regions: Vec<Region> = Vec::new();
    for us in user_sets {
        all_regions.extend(us.regions.iter().cloned());
    }

    // Merge overlapping intervals per chromosome
    let merged = merge_regions(&all_regions);
    RegionSet::from(merged)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

type ChromIndex<'a> = HashMap<String, Vec<(u32, u32, usize)>>; // (start, end, index into universe.regions)

/// Build a per-chromosome sorted index of regions for overlap queries.
fn build_chrom_index(region_set: &RegionSet) -> ChromIndex<'_> {
    let mut index: ChromIndex = HashMap::new();

    for (i, region) in region_set.regions.iter().enumerate() {
        index
            .entry(region.chr.clone())
            .or_default()
            .push((region.start, region.end, i));
    }

    // Sort each chromosome's entries by start position
    for entries in index.values_mut() {
        entries.sort_by_key(|&(start, _, _)| start);
    }

    index
}

/// Count how many regions in the index overlap the given interval.
fn count_overlapping_regions(
    chrom: &str,
    start: u32,
    end: u32,
    index: &ChromIndex,
) -> usize {
    let entries = match index.get(chrom) {
        Some(e) => e,
        None => return 0,
    };

    let mut count = 0;
    for &(s, e, _) in entries {
        if s >= end {
            break; // sorted by start, no more overlaps possible
        }
        if e > start {
            count += 1;
        }
    }
    count
}

/// Find all regions in the index that overlap the given interval.
fn find_overlapping_regions(
    chrom: &str,
    start: u32,
    end: u32,
    index: &ChromIndex,
) -> Vec<Region> {
    let entries = match index.get(chrom) {
        Some(e) => e,
        None => return Vec::new(),
    };

    let mut result = Vec::new();
    for &(s, e, _) in entries {
        if s >= end {
            break;
        }
        if e > start {
            result.push(Region {
                chr: chrom.to_string(),
                start: s,
                end: e,
                rest: None,
            });
        }
    }
    result
}

/// Merge overlapping intervals across all chromosomes.
/// Returns a sorted, disjoint set of regions.
fn merge_regions(regions: &[Region]) -> Vec<Region> {
    // Group by chromosome
    let mut by_chr: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    for r in regions {
        by_chr
            .entry(r.chr.clone())
            .or_default()
            .push((r.start, r.end));
    }

    let mut merged: Vec<Region> = Vec::new();

    let mut chroms: Vec<String> = by_chr.keys().cloned().collect();
    chroms.sort();

    for chr in &chroms {
        let intervals = by_chr.get_mut(chr.as_str()).unwrap();
        intervals.sort_by_key(|&(s, _)| s);

        let mut it = intervals.iter();
        if let Some(&(mut cur_start, mut cur_end)) = it.next() {
            for &(s, e) in it {
                if s <= cur_end {
                    // Overlapping or adjacent — extend
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
        }
    }

    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    fn rs(regions: Vec<Region>) -> RegionSet {
        RegionSet::from(regions)
    }

    // -------------------------------------------------------------------
    // check_universe_appropriateness
    // -------------------------------------------------------------------

    #[test]
    fn test_check_universe_full_coverage() {
        let universe = rs(vec![
            r("chr1", 0, 1000),
            r("chr1", 2000, 3000),
        ]);
        let user = rs(vec![
            r("chr1", 100, 200),   // inside [0,1000)
            r("chr1", 2500, 2600), // inside [2000,3000)
        ]);

        let report = check_universe_appropriateness(&[user], &universe);
        assert_eq!(report.user_set_reports.len(), 1);

        let ur = &report.user_set_reports[0];
        assert_eq!(ur.total_regions, 2);
        assert_eq!(ur.regions_in_universe, 2);
        assert!((ur.coverage - 1.0).abs() < 1e-10);
        assert_eq!(ur.many_to_many_count, 0);
        assert!(ur.warnings.is_empty());
    }

    #[test]
    fn test_check_universe_low_coverage() {
        let universe = rs(vec![r("chr1", 0, 100)]);
        let user = rs(vec![
            r("chr1", 50, 80),    // overlaps
            r("chr1", 500, 600),  // does NOT overlap
            r("chr1", 700, 800),  // does NOT overlap
        ]);

        let report = check_universe_appropriateness(&[user], &universe);
        let ur = &report.user_set_reports[0];

        assert_eq!(ur.regions_in_universe, 1);
        assert!((ur.coverage - 1.0 / 3.0).abs() < 0.01);
        assert!(!ur.warnings.is_empty()); // should warn about low coverage
    }

    #[test]
    fn test_check_universe_many_to_many() {
        // User region overlaps two universe regions
        let universe = rs(vec![
            r("chr1", 100, 200),
            r("chr1", 150, 250),
        ]);
        let user = rs(vec![r("chr1", 120, 220)]); // overlaps both

        let report = check_universe_appropriateness(&[user], &universe);
        let ur = &report.user_set_reports[0];

        assert_eq!(ur.many_to_many_count, 1);
        assert!(ur.warnings.iter().any(|w| w.contains("many-to-many")));
    }

    // -------------------------------------------------------------------
    // redefine_user_sets
    // -------------------------------------------------------------------

    #[test]
    fn test_redefine_user_sets_basic() {
        let universe = rs(vec![
            r("chr1", 100, 200),
            r("chr1", 300, 400),
            r("chr1", 500, 600),
        ]);
        // User region overlaps first two universe regions
        let user = rs(vec![r("chr1", 150, 350)]);

        let redefined = redefine_user_sets(&[user], &universe);
        assert_eq!(redefined.len(), 1);

        let new_set = &redefined[0];
        assert_eq!(new_set.regions.len(), 2);
        assert_eq!(new_set.regions[0].start, 100);
        assert_eq!(new_set.regions[0].end, 200);
        assert_eq!(new_set.regions[1].start, 300);
        assert_eq!(new_set.regions[1].end, 400);
    }

    #[test]
    fn test_redefine_user_sets_dedup() {
        // Two user regions overlap the same universe region
        let universe = rs(vec![r("chr1", 100, 300)]);
        let user = rs(vec![
            r("chr1", 120, 150),
            r("chr1", 200, 250),
        ]);

        let redefined = redefine_user_sets(&[user], &universe);
        assert_eq!(redefined[0].regions.len(), 1); // deduplicated
    }

    #[test]
    fn test_redefine_user_sets_no_overlap() {
        let universe = rs(vec![r("chr1", 100, 200)]);
        let user = rs(vec![r("chr1", 500, 600)]);

        let redefined = redefine_user_sets(&[user], &universe);
        assert_eq!(redefined[0].regions.len(), 0);
    }

    // -------------------------------------------------------------------
    // build_restricted_universe
    // -------------------------------------------------------------------

    #[test]
    fn test_build_restricted_universe_basic() {
        let user0 = rs(vec![
            r("chr1", 100, 200),
            r("chr1", 300, 400),
        ]);
        let user1 = rs(vec![
            r("chr1", 150, 250), // overlaps user0[0]
            r("chr2", 100, 200),
        ]);

        let restricted = build_restricted_universe(&[user0, user1]);

        // chr1: [100,200) + [150,250) merge to [100,250); [300,400) stays
        // chr2: [100,200)
        assert_eq!(restricted.regions.len(), 3);

        let chr1: Vec<&Region> = restricted
            .regions
            .iter()
            .filter(|r| r.chr == "chr1")
            .collect();
        assert_eq!(chr1.len(), 2);
        assert_eq!(chr1[0].start, 100);
        assert_eq!(chr1[0].end, 250); // merged
        assert_eq!(chr1[1].start, 300);
        assert_eq!(chr1[1].end, 400);

        let chr2: Vec<&Region> = restricted
            .regions
            .iter()
            .filter(|r| r.chr == "chr2")
            .collect();
        assert_eq!(chr2.len(), 1);
    }

    #[test]
    fn test_build_restricted_universe_empty() {
        let restricted = build_restricted_universe(&[]);
        assert!(restricted.regions.is_empty());
    }

    #[test]
    fn test_build_restricted_universe_disjoint() {
        // Non-overlapping regions should stay separate
        let user = rs(vec![
            r("chr1", 100, 200),
            r("chr1", 300, 400),
            r("chr1", 500, 600),
        ]);

        let restricted = build_restricted_universe(&[user]);
        assert_eq!(restricted.regions.len(), 3);
    }

    // -------------------------------------------------------------------
    // merge_regions
    // -------------------------------------------------------------------

    #[test]
    fn test_merge_regions_overlapping() {
        let regions = vec![
            r("chr1", 100, 300),
            r("chr1", 200, 400),
            r("chr1", 350, 500),
        ];
        let merged = merge_regions(&regions);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].end, 500);
    }

    #[test]
    fn test_merge_regions_adjacent() {
        // Adjacent intervals (end == start of next) should merge
        let regions = vec![
            r("chr1", 100, 200),
            r("chr1", 200, 300),
        ];
        let merged = merge_regions(&regions);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].end, 300);
    }

    #[test]
    fn test_merge_regions_multi_chrom() {
        let regions = vec![
            r("chr2", 100, 200),
            r("chr1", 100, 200),
            r("chr1", 150, 250),
        ];
        let merged = merge_regions(&regions);
        assert_eq!(merged.len(), 2); // chr1: [100,250), chr2: [100,200)
        assert_eq!(merged[0].chr, "chr1");
        assert_eq!(merged[1].chr, "chr2");
    }
}
