//! Core enrichment engine: contingency tables, Fisher's exact test, ranking.

use statrs::distribution::{DiscreteCDF, Hypergeometric};

use gtars_core::models::RegionSet;
use gtars_igd::igd::Igd;

use crate::errors::LolaError;
use crate::models::{ContingencyTable, Direction, LolaConfig, LolaResult};

impl ContingencyTable {
    /// Compute the one-sided p-value using Fisher's exact test.
    ///
    /// Uses the hypergeometric distribution where:
    /// - Population size N = a + b + c + d
    /// - Number of success states K = a + b
    /// - Number of draws n = a + c
    /// - Observed successes k = a
    pub fn fisher_pvalue(&self, direction: Direction) -> f64 {
        let n_pop = self.a + self.b + self.c + self.d;
        let k_success = self.a + self.b;
        let n_draws = self.a + self.c;

        // Edge cases
        if n_pop == 0 || k_success == 0 || n_draws == 0 {
            return 1.0;
        }
        if k_success > n_pop || n_draws > n_pop {
            return 1.0;
        }

        let hyper = match Hypergeometric::new(n_pop, k_success, n_draws) {
            Ok(h) => h,
            Err(_) => return 1.0,
        };

        match direction {
            Direction::Enrichment => {
                // P(X >= a): 1 - P(X <= a-1)
                if self.a == 0 {
                    1.0
                } else {
                    1.0 - hyper.cdf(self.a - 1)
                }
            }
            Direction::Depletion => {
                // P(X <= a)
                hyper.cdf(self.a)
            }
        }
    }

    /// Compute the odds ratio: (a*d) / (b*c).
    /// Returns f64::INFINITY if denominator is zero and numerator > 0.
    /// Returns 0.0 if numerator is zero.
    pub fn odds_ratio(&self) -> f64 {
        let num = self.a as f64 * self.d as f64;
        let den = self.b as f64 * self.c as f64;

        if den == 0.0 {
            if num > 0.0 {
                f64::INFINITY
            } else {
                0.0 // 0/0 case
            }
        } else {
            num / den
        }
    }

    /// Compute -log10(p-value).
    pub fn p_value_log(&self, direction: Direction) -> f64 {
        let p = self.fisher_pvalue(direction);
        if p <= 0.0 {
            f64::INFINITY
        } else {
            -p.log10()
        }
    }
}

/// Run LOLA enrichment analysis.
///
/// For each (user_set, db_set) pair, builds a contingency table and runs
/// Fisher's exact test. Results are ranked by p-value, odds ratio, and support.
///
/// # Arguments
/// * `igd` - The IGD index containing all database region sets
/// * `user_sets` - User region sets to test
/// * `universe` - The universe (background) region set
/// * `config` - LOLA configuration (min_overlap, direction)
pub fn run_lola(
    igd: &Igd,
    user_sets: &[RegionSet],
    universe: &RegionSet,
    config: &LolaConfig,
) -> Result<Vec<LolaResult>, LolaError> {
    let n_db = igd.num_files();
    if n_db == 0 {
        return Err(LolaError::EmptyDatabase);
    }

    let universe_size = universe.regions.len() as u64;
    if universe_size == 0 {
        return Err(LolaError::EmptyUniverse);
    }

    // Step 1: Query universe against IGD (done once for all user sets)
    let universe_hits = igd.count_set_overlaps(universe, config.min_overlap);

    // Step 2: For each user set, query against IGD and build contingency tables
    let mut all_results: Vec<LolaResult> = Vec::new();

    for (us_idx, user_set) in user_sets.iter().enumerate() {
        let user_set_size = user_set.regions.len() as u64;

        let user_hits = igd.count_set_overlaps(user_set, config.min_overlap);

        let mut user_results: Vec<LolaResult> = Vec::with_capacity(n_db);

        for db_idx in 0..n_db {
            let a = user_hits[db_idx];
            let b = universe_hits[db_idx].saturating_sub(a);
            let c = user_set_size.saturating_sub(a);
            let d = universe_size.saturating_sub(a + b + c);

            let ct = ContingencyTable { a, b, c, d };
            let pv_log = ct.p_value_log(config.direction);
            let or = ct.odds_ratio();

            let filename = if db_idx < igd.file_info.len() {
                igd.file_info[db_idx].filename.clone()
            } else {
                String::new()
            };

            user_results.push(LolaResult {
                user_set: us_idx,
                db_set: db_idx,
                p_value_log: pv_log,
                odds_ratio: or,
                support: a,
                rnk_pv: 0,
                rnk_or: 0,
                rnk_sup: 0,
                max_rnk: 0,
                mean_rnk: 0.0,
                b,
                c,
                d,
                q_value: None,
                filename,
            });
        }

        // Step 3: Rank results within this user set
        rank_results(&mut user_results);

        all_results.extend(user_results);
    }

    Ok(all_results)
}

/// Assign ranks to results by p_value_log (descending), odds_ratio (descending),
/// and support (descending). Then compute max_rnk and mean_rnk.
fn rank_results(results: &mut [LolaResult]) {
    let n = results.len();
    if n == 0 {
        return;
    }

    // Rank by p_value_log (descending — highest -log10(p) = most significant gets rank 1)
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        results[b]
            .p_value_log
            .partial_cmp(&results[a].p_value_log)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for (rank, &idx) in indices.iter().enumerate() {
        results[idx].rnk_pv = rank + 1;
    }

    // Rank by odds_ratio (descending)
    indices.sort_by(|&a, &b| {
        results[b]
            .odds_ratio
            .partial_cmp(&results[a].odds_ratio)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for (rank, &idx) in indices.iter().enumerate() {
        results[idx].rnk_or = rank + 1;
    }

    // Rank by support (descending)
    indices.sort_by(|&a, &b| results[b].support.cmp(&results[a].support));
    for (rank, &idx) in indices.iter().enumerate() {
        results[idx].rnk_sup = rank + 1;
    }

    // Compute combined ranks
    for r in results.iter_mut() {
        r.max_rnk = r.rnk_pv.max(r.rnk_or).max(r.rnk_sup);
        r.mean_rnk = (r.rnk_pv + r.rnk_or + r.rnk_sup) as f64 / 3.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_core::models::Region;
    use gtars_igd::igd::Igd;

    // -------------------------------------------------------------------
    // ContingencyTable unit tests
    // -------------------------------------------------------------------

    #[test]
    fn test_contingency_odds_ratio() {
        let ct = ContingencyTable {
            a: 10,
            b: 20,
            c: 30,
            d: 40,
        };
        // (10*40) / (20*30) = 400/600 = 0.6667
        let or = ct.odds_ratio();
        assert!((or - 0.6667).abs() < 0.001);
    }

    #[test]
    fn test_contingency_odds_ratio_zero_denom() {
        // b=0: denominator is 0, numerator > 0 → infinity
        let ct = ContingencyTable {
            a: 10,
            b: 0,
            c: 5,
            d: 100,
        };
        assert_eq!(ct.odds_ratio(), f64::INFINITY);

        // a=0, d=0: numerator is 0 → 0
        let ct2 = ContingencyTable {
            a: 0,
            b: 5,
            c: 0,
            d: 0,
        };
        assert_eq!(ct2.odds_ratio(), 0.0);
    }

    #[test]
    fn test_fisher_enrichment_significant() {
        // Strong enrichment: a is much larger than expected
        let ct = ContingencyTable {
            a: 50,
            b: 10,
            c: 5,
            d: 1000,
        };
        let p = ct.fisher_pvalue(Direction::Enrichment);
        assert!(p < 0.001, "Expected significant enrichment, got p={}", p);
    }

    #[test]
    fn test_fisher_enrichment_not_significant() {
        // No enrichment: a is small relative to expected
        let ct = ContingencyTable {
            a: 1,
            b: 100,
            c: 100,
            d: 1000,
        };
        let p = ct.fisher_pvalue(Direction::Enrichment);
        assert!(p > 0.05, "Expected non-significant, got p={}", p);
    }

    #[test]
    fn test_fisher_depletion() {
        // Strong depletion: a is much smaller than expected
        let ct = ContingencyTable {
            a: 1,
            b: 100,
            c: 100,
            d: 10,
        };
        let p = ct.fisher_pvalue(Direction::Depletion);
        assert!(p < 0.05, "Expected significant depletion, got p={}", p);
    }

    #[test]
    fn test_fisher_edge_cases() {
        // Empty population
        let ct = ContingencyTable {
            a: 0,
            b: 0,
            c: 0,
            d: 0,
        };
        assert_eq!(ct.fisher_pvalue(Direction::Enrichment), 1.0);

        // a=0, enrichment test → p should be high (no overlap at all)
        let ct2 = ContingencyTable {
            a: 0,
            b: 50,
            c: 50,
            d: 100,
        };
        assert_eq!(ct2.fisher_pvalue(Direction::Enrichment), 1.0);
    }

    #[test]
    fn test_p_value_log() {
        // Use moderate values so p-value doesn't underflow to 0
        let ct = ContingencyTable {
            a: 5,
            b: 15,
            c: 10,
            d: 100,
        };
        let pvl = ct.p_value_log(Direction::Enrichment);
        assert!(pvl > 0.0, "Expected positive -log10(p), got {}", pvl);
        // Cross-check: -log10(p) should match manual calculation
        let p = ct.fisher_pvalue(Direction::Enrichment);
        let expected = -p.log10();
        assert!(
            (pvl - expected).abs() < 1e-10,
            "pvl={} expected={}",
            pvl,
            expected
        );
    }

    #[test]
    fn test_p_value_log_extreme() {
        // Very significant result → p-value underflows to 0 → -log10(0) = inf
        let ct = ContingencyTable {
            a: 50,
            b: 10,
            c: 5,
            d: 1000,
        };
        let pvl = ct.p_value_log(Direction::Enrichment);
        assert_eq!(pvl, f64::INFINITY);
    }

    // -------------------------------------------------------------------
    // Ranking tests
    // -------------------------------------------------------------------

    #[test]
    fn test_ranking() {
        let mut results = vec![
            LolaResult {
                user_set: 0,
                db_set: 0,
                p_value_log: 5.0,
                odds_ratio: 2.0,
                support: 100,
                rnk_pv: 0,
                rnk_or: 0,
                rnk_sup: 0,
                max_rnk: 0,
                mean_rnk: 0.0,
                b: 0,
                c: 0,
                d: 0,
                q_value: None,
                filename: "a.bed".into(),
            },
            LolaResult {
                user_set: 0,
                db_set: 1,
                p_value_log: 10.0,
                odds_ratio: 1.0,
                support: 200,
                rnk_pv: 0,
                rnk_or: 0,
                rnk_sup: 0,
                max_rnk: 0,
                mean_rnk: 0.0,
                b: 0,
                c: 0,
                d: 0,
                q_value: None,
                filename: "b.bed".into(),
            },
            LolaResult {
                user_set: 0,
                db_set: 2,
                p_value_log: 3.0,
                odds_ratio: 5.0,
                support: 50,
                rnk_pv: 0,
                rnk_or: 0,
                rnk_sup: 0,
                max_rnk: 0,
                mean_rnk: 0.0,
                b: 0,
                c: 0,
                d: 0,
                q_value: None,
                filename: "c.bed".into(),
            },
        ];

        rank_results(&mut results);

        // p_value_log descending: 10.0 > 5.0 > 3.0 → ranks 1, 2, 3
        assert_eq!(results[0].rnk_pv, 2); // db_set 0: pvl=5.0
        assert_eq!(results[1].rnk_pv, 1); // db_set 1: pvl=10.0
        assert_eq!(results[2].rnk_pv, 3); // db_set 2: pvl=3.0

        // odds_ratio descending: 5.0 > 2.0 > 1.0 → ranks 1, 2, 3
        assert_eq!(results[0].rnk_or, 2); // or=2.0
        assert_eq!(results[1].rnk_or, 3); // or=1.0
        assert_eq!(results[2].rnk_or, 1); // or=5.0

        // support descending: 200 > 100 > 50 → ranks 1, 2, 3
        assert_eq!(results[0].rnk_sup, 2);
        assert_eq!(results[1].rnk_sup, 1);
        assert_eq!(results[2].rnk_sup, 3);

        // max_rnk: max of all three ranks
        assert_eq!(results[0].max_rnk, 2); // max(2,2,2) = 2
        assert_eq!(results[1].max_rnk, 3); // max(1,3,1) = 3
        assert_eq!(results[2].max_rnk, 3); // max(3,1,3) = 3

        // mean_rnk
        assert!((results[0].mean_rnk - 2.0).abs() < 1e-10);
        assert!((results[1].mean_rnk - 5.0 / 3.0).abs() < 1e-10);
        assert!((results[2].mean_rnk - 7.0 / 3.0).abs() < 1e-10);
    }

    // -------------------------------------------------------------------
    // Integration: run_lola
    // -------------------------------------------------------------------

    fn make_region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    fn make_region_set(regions: Vec<Region>) -> RegionSet {
        RegionSet::from(regions)
    }

    #[test]
    fn test_run_lola_basic() {
        // Set up: 2 DB sets
        // db0: regions at chr1:[100,200), chr1:[500,600)
        // db1: regions at chr1:[150,250)
        let sets = vec![
            (
                "db0.bed".to_string(),
                vec![
                    ("chr1".to_string(), 100, 200),
                    ("chr1".to_string(), 500, 600),
                ],
            ),
            (
                "db1.bed".to_string(),
                vec![("chr1".to_string(), 150, 250)],
            ),
        ];
        let igd = Igd::from_region_sets(sets);

        // User set: one region overlapping db0[0] and db1[0]
        let user = make_region_set(vec![make_region("chr1", 120, 180)]);

        // Universe: covers a broader range
        let universe = make_region_set(vec![
            make_region("chr1", 50, 250),
            make_region("chr1", 450, 650),
            make_region("chr1", 700, 800),
            make_region("chr1", 900, 1000),
            make_region("chr1", 1100, 1200),
        ]);

        let config = LolaConfig::default();
        let results = run_lola(&igd, &[user], &universe, &config).unwrap();

        assert_eq!(results.len(), 2); // 1 user set × 2 db sets

        // Both db sets should have support > 0
        let r0 = results.iter().find(|r| r.db_set == 0).unwrap();
        let r1 = results.iter().find(|r| r.db_set == 1).unwrap();

        assert_eq!(r0.support, 1); // user overlaps db0[0] (one pairwise overlap)
        assert_eq!(r1.support, 1); // user overlaps db1[0]

        // Ranks should be assigned
        assert!(r0.rnk_pv >= 1);
        assert!(r1.rnk_pv >= 1);
    }

    #[test]
    fn test_run_lola_no_overlap() {
        let sets = vec![(
            "db0.bed".to_string(),
            vec![("chr1".to_string(), 100, 200)],
        )];
        let igd = Igd::from_region_sets(sets);

        // User set doesn't overlap db at all
        let user = make_region_set(vec![make_region("chr1", 500, 600)]);
        let universe = make_region_set(vec![
            make_region("chr1", 50, 250),
            make_region("chr1", 450, 650),
        ]);

        let results = run_lola(&igd, &[user], &universe, &LolaConfig::default()).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].support, 0);
        assert_eq!(results[0].p_value_log, 0.0); // -log10(1.0) = 0
    }

    #[test]
    fn test_run_lola_multiple_user_sets() {
        let sets = vec![(
            "db0.bed".to_string(),
            vec![
                ("chr1".to_string(), 100, 200),
                ("chr1".to_string(), 300, 400),
            ],
        )];
        let igd = Igd::from_region_sets(sets);

        let user0 = make_region_set(vec![make_region("chr1", 150, 180)]);
        let user1 = make_region_set(vec![
            make_region("chr1", 150, 180),
            make_region("chr1", 350, 380),
        ]);

        let universe = make_region_set(vec![
            make_region("chr1", 50, 250),
            make_region("chr1", 250, 450),
            make_region("chr1", 500, 600),
        ]);

        let results =
            run_lola(&igd, &[user0, user1], &universe, &LolaConfig::default()).unwrap();

        assert_eq!(results.len(), 2); // 2 user sets × 1 db set

        let r0 = results.iter().find(|r| r.user_set == 0).unwrap();
        let r1 = results.iter().find(|r| r.user_set == 1).unwrap();

        assert_eq!(r0.support, 1); // user0 has 1 overlapping region
        assert_eq!(r1.support, 2); // user1 has 2 overlapping regions
    }

    #[test]
    fn test_run_lola_pairwise_counting() {
        // One query region overlaps 3 DB regions in same file → support = 3
        let sets = vec![(
            "db0.bed".to_string(),
            vec![
                ("chr1".to_string(), 100, 200),
                ("chr1".to_string(), 120, 220),
                ("chr1".to_string(), 140, 240),
            ],
        )];
        let igd = Igd::from_region_sets(sets);

        let user = make_region_set(vec![make_region("chr1", 150, 190)]);
        let universe = make_region_set(vec![
            make_region("chr1", 50, 300),
            make_region("chr1", 400, 500),
        ]);

        let results = run_lola(&igd, &[user], &universe, &LolaConfig::default()).unwrap();
        assert_eq!(results[0].support, 3);
    }

    #[test]
    fn test_run_lola_empty_db() {
        let igd = Igd::from_region_sets(Vec::<(String, Vec<(String, i32, i32)>)>::new());
        let user = make_region_set(vec![make_region("chr1", 100, 200)]);
        let universe = make_region_set(vec![make_region("chr1", 50, 300)]);

        let result = run_lola(&igd, &[user], &universe, &LolaConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn test_run_lola_empty_universe() {
        let sets = vec![(
            "db0.bed".to_string(),
            vec![("chr1".to_string(), 100, 200)],
        )];
        let igd = Igd::from_region_sets(sets);
        let user = make_region_set(vec![make_region("chr1", 100, 200)]);
        let universe = make_region_set(vec![]);

        let result = run_lola(&igd, &[user], &universe, &LolaConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn test_run_lola_contingency_table_arithmetic() {
        // Verify the contingency table values are computed correctly
        //
        // Setup:
        // - Universe: 10 regions, 2 overlap db0
        // - User: 3 regions, 1 overlaps db0
        // Expected: a=1, b=2-1=1, c=3-1=2, d=10-1-1-2=6
        let sets = vec![(
            "db0.bed".to_string(),
            vec![
                ("chr1".to_string(), 100, 200),
                ("chr1".to_string(), 300, 400),
            ],
        )];
        let igd = Igd::from_region_sets(sets);

        let user = make_region_set(vec![
            make_region("chr1", 150, 180),  // overlaps db0[0]
            make_region("chr1", 500, 600),  // no overlap
            make_region("chr1", 700, 800),  // no overlap
        ]);

        let universe = make_region_set(vec![
            make_region("chr1", 50, 250),   // overlaps db0[0]
            make_region("chr1", 250, 450),  // overlaps db0[1]
            make_region("chr1", 450, 550),
            make_region("chr1", 550, 650),
            make_region("chr1", 650, 750),
            make_region("chr1", 750, 850),
            make_region("chr1", 850, 950),
            make_region("chr1", 950, 1050),
            make_region("chr1", 1050, 1150),
            make_region("chr1", 1150, 1250),
        ]);

        let results = run_lola(&igd, &[user], &universe, &LolaConfig::default()).unwrap();
        let r = &results[0];

        assert_eq!(r.support, 1, "a (support)");
        assert_eq!(r.b, 1, "b = universe_hits - a = 2 - 1");
        assert_eq!(r.c, 2, "c = user_size - a = 3 - 1");
        assert_eq!(r.d, 6, "d = universe_size - a - b - c = 10 - 1 - 1 - 2");
    }

    #[test]
    fn test_run_lola_with_min_overlap() {
        let sets = vec![(
            "db0.bed".to_string(),
            vec![("chr1".to_string(), 100, 200)],
        )];
        let igd = Igd::from_region_sets(sets);

        // User region [190, 250) — only 10bp overlap with db0 [100, 200)
        let user = make_region_set(vec![make_region("chr1", 190, 250)]);
        let universe = make_region_set(vec![
            make_region("chr1", 50, 300),
            make_region("chr1", 400, 500),
        ]);

        // min_overlap=1: should count
        let config1 = LolaConfig {
            min_overlap: 1,
            ..Default::default()
        };
        let r1 = run_lola(&igd, &[user.clone()], &universe, &config1).unwrap();
        assert_eq!(r1[0].support, 1);

        // min_overlap=50: should NOT count (only 10bp overlap)
        let config50 = LolaConfig {
            min_overlap: 50,
            ..Default::default()
        };
        let r50 = run_lola(&igd, &[user], &universe, &config50).unwrap();
        assert_eq!(r50[0].support, 0);
    }
}
