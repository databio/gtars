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
                // P(X >= a) = P(X > a-1) = sf(a-1)
                // Using sf() (survival function) avoids catastrophic cancellation
                // that occurs with 1.0 - cdf(a-1) when the CDF is near 1.0.
                if self.a == 0 {
                    1.0
                } else {
                    hyper.sf(self.a - 1)
                }
            }
            Direction::Depletion => {
                // P(X <= a)
                hyper.cdf(self.a)
            }
        }
    }

    /// Compute the conditional maximum likelihood estimate (CMLE) of the odds
    /// ratio, matching R's `fisher.test()$estimate`.
    ///
    /// This finds the noncentrality parameter ω of the Fisher noncentral
    /// hypergeometric distribution such that E[X | ω] equals the observed
    /// cell count `a`.  When the observed `a` is at the boundary of its
    /// support the MLE is 0 or ∞.
    pub fn odds_ratio(&self) -> f64 {
        let m = self.a + self.c; // column 1 total
        let n = self.b + self.d; // column 2 total
        let k = self.a + self.b; // row 1 total
        let x = self.a;

        let lo = if k > n { k - n } else { 0 };
        let hi = k.min(m);

        // Boundary / degenerate cases
        if lo == hi {
            return f64::NAN; // only one possible table
        }
        if x == lo {
            return 0.0;
        }
        if x == hi {
            return f64::INFINITY;
        }

        // Precompute log-densities of the central hypergeometric distribution
        // using pure recurrence (no lgamma calls) for maximum precision.
        // d(y+1)/d(y) = (m-y)(k-y) / ((y+1)(n-k+y+1))
        // We set logdc[0] = 0 since only relative values matter (they cancel
        // in the normalised mean computation).
        let support_size = (hi - lo + 1) as usize;
        let mut logdc = Vec::with_capacity(support_size);
        logdc.push(0.0); // logdc[0] = log(1) = 0 (arbitrary reference)
        for i in 1..support_size {
            let y = lo + i as u64 - 1; // previous support value
            let log_ratio = ((m - y) as f64).ln() + ((k - y) as f64).ln()
                - ((y + 1) as f64).ln()
                - ((n - k + y + 1) as f64).ln();
            logdc.push(logdc[i - 1] + log_ratio);
        }

        // Mean of the noncentral hypergeometric for a given omega.
        // Uses compensated (Kahan) summation for precision.
        let mean_nhyper = |omega: f64| -> f64 {
            if omega == 0.0 {
                return lo as f64;
            }
            if omega.is_infinite() {
                return hi as f64;
            }
            let log_omega = omega.ln();
            let log_vals: Vec<f64> = logdc
                .iter()
                .enumerate()
                .map(|(i, &ld)| ld + (lo as f64 + i as f64) * log_omega)
                .collect();
            let max_log = log_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

            // Kahan summation for both sum and weighted sum
            let mut sum = 0.0f64;
            let mut sum_c = 0.0f64;
            let mut wsum = 0.0f64;
            let mut wsum_c = 0.0f64;
            for (i, &lv) in log_vals.iter().enumerate() {
                let w = (lv - max_log).exp();
                let y = lo as f64 + i as f64;

                let yw = y * w - wsum_c;
                let wt = wsum + yw;
                wsum_c = (wt - wsum) - yw;
                wsum = wt;

                let sw = w - sum_c;
                let st = sum + sw;
                sum_c = (st - sum) - sw;
                sum = st;
            }
            wsum / sum
        };

        let xf = x as f64;
        let mu1 = mean_nhyper(1.0);

        if (mu1 - xf).abs() < 1e-12 {
            return 1.0;
        }

        // Brent's method (bisection with ITP-like speedup)
        // E[X | omega] is monotonically increasing in omega, so a unique root exists.
        if mu1 > xf {
            // omega < 1: search [0, 1]
            brent(|t| mean_nhyper(t) - xf, 0.0, 1.0, 1e-8, 100)
        } else {
            // omega > 1: reparameterise as t = 1/omega, search t in (eps, 1)
            let t = brent(
                |t| mean_nhyper(1.0 / t) - xf,
                f64::EPSILON,
                1.0,
                1e-8,
                100,
            );
            1.0 / t
        }
    }

    /// Compute -log10(p-value).
    ///
    /// Adds 1e-322 floor to avoid infinity, matching R LOLA's
    /// `-log10(pValueLog + 10^-322)` behavior. Caps at ~322.
    pub fn p_value_log(&self, direction: Direction) -> f64 {
        let p = self.fisher_pvalue(direction);
        -(p + 1e-322_f64).log10()
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

    // Step 1: Query universe against IGD (binary per region)
    // count_region_hits gives the number of universe regions overlapping each DB set.
    let universe_hits = igd.count_region_hits(universe, config.min_overlap);

    // Step 2: For each user set, query against IGD and build contingency tables
    let mut all_results: Vec<LolaResult> = Vec::new();

    for (us_idx, user_set) in user_sets.iter().enumerate() {
        let user_set_size = user_set.regions.len() as u64;

        // Binary per-region overlap counting (matches R LOLA's countOverlaps semantics)
        let user_hits = igd.count_region_hits(user_set, config.min_overlap);

        let mut user_results: Vec<LolaResult> = Vec::with_capacity(n_db);

        for db_idx in 0..n_db {
            let a = user_hits[db_idx];

            // Compute b, c, d as signed values, matching R LOLA behavior.
            // R LOLA passes negative values through and warns.
            let b = universe_hits[db_idx] as i64 - a as i64;
            let c = user_set_size as i64 - a as i64;
            let d = universe_size as i64 - a as i64 - b - c;

            let has_negative = b < 0 || c < 0 || d < 0;
            if has_negative {
                eprintln!(
                    "Warning: negative contingency value for db_set {} (user_set {}). \
                     This means your user sets contain regions outside the universe.",
                    db_idx, us_idx
                );
            }

            // Only compute Fisher test when all values are non-negative.
            // For negative rows, return p_value_log=0 and odds_ratio=NaN.
            let (pv_log, or) = if has_negative {
                (0.0, f64::NAN)
            } else {
                let ct = ContingencyTable {
                    a,
                    b: b as u64,
                    c: c as u64,
                    d: d as u64,
                };
                (ct.p_value_log(config.direction), ct.odds_ratio())
            };

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
                collection: None,
                description: None,
                cell_type: None,
                tissue: None,
                antibody: None,
                treatment: None,
                data_source: None,
                db_set_size: 0,
            });
        }

        // Step 3: Rank results within this user set
        rank_results(&mut user_results);

        all_results.extend(user_results);
    }

    // Sort by descending pValueLog, then ascending meanRnk (matches R LOLA output order)
    all_results.sort_by(|a, b| {
        b.p_value_log
            .partial_cmp(&a.p_value_log)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                a.mean_rnk
                    .partial_cmp(&b.mean_rnk)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    Ok(all_results)
}

/// Compare two f64 values for tie-breaking equality.
/// NaN == NaN, and normal values use bitwise equality (exact match).
fn f64_tied(a: f64, b: f64) -> bool {
    if a.is_nan() && b.is_nan() {
        return true;
    }
    a.to_bits() == b.to_bits()
}

/// Assign min-rank to a pre-sorted index array (ties.method="min").
/// `get_val` extracts the ranking key; `set_rank` writes the rank back.
fn assign_min_ranks_f64(
    indices: &[usize],
    results: &mut [LolaResult],
    get_val: impl Fn(&LolaResult) -> f64,
    mut set_rank: impl FnMut(&mut LolaResult, usize),
) {
    if indices.is_empty() {
        return;
    }
    let mut rank = 1usize; // rank of the current tie group
    set_rank(&mut results[indices[0]], 1);
    for i in 1..indices.len() {
        let prev = get_val(&results[indices[i - 1]]);
        let curr = get_val(&results[indices[i]]);
        if !f64_tied(prev, curr) {
            rank = i + 1; // new group starts at 1-based position
        }
        set_rank(&mut results[indices[i]], rank);
    }
}

/// Assign min-rank for u64 values.
fn assign_min_ranks_u64(
    indices: &[usize],
    results: &mut [LolaResult],
    get_val: impl Fn(&LolaResult) -> u64,
    mut set_rank: impl FnMut(&mut LolaResult, usize),
) {
    if indices.is_empty() {
        return;
    }
    let mut rank = 1usize;
    set_rank(&mut results[indices[0]], 1);
    for i in 1..indices.len() {
        if get_val(&results[indices[i - 1]]) != get_val(&results[indices[i]]) {
            rank = i + 1;
        }
        set_rank(&mut results[indices[i]], rank);
    }
}

/// Assign ranks to results by p_value_log (descending), odds_ratio (descending),
/// and support (descending). Then compute max_rnk and mean_rnk.
fn rank_results(results: &mut [LolaResult]) {
    let n = results.len();
    if n == 0 {
        return;
    }

    // Helper: assign min-rank (ties.method="min") after sorting indices descending.
    // All tied values get the rank of the first occurrence in sorted order.
    let mut indices: Vec<usize> = (0..n).collect();

    // Rank by p_value_log (descending — highest -log10(p) = most significant gets rank 1)
    indices.sort_by(|&a, &b| {
        results[b]
            .p_value_log
            .partial_cmp(&results[a].p_value_log)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    assign_min_ranks_f64(&indices, results, |r| r.p_value_log, |r, v| r.rnk_pv = v);

    // Rank by odds_ratio (descending, NaN gets worst rank)
    indices.sort_by(|&a, &b| {
        let ra = results[a].odds_ratio;
        let rb = results[b].odds_ratio;
        match (ra.is_nan(), rb.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater,
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => rb.partial_cmp(&ra).unwrap_or(std::cmp::Ordering::Equal),
        }
    });
    assign_min_ranks_f64(&indices, results, |r| r.odds_ratio, |r, v| r.rnk_or = v);

    // Rank by support (descending)
    indices.sort_by(|&a, &b| results[b].support.cmp(&results[a].support));
    assign_min_ranks_u64(&indices, results, |r| r.support, |r, v| r.rnk_sup = v);

    // Compute combined ranks
    for r in results.iter_mut() {
        r.max_rnk = r.rnk_pv.max(r.rnk_or).max(r.rnk_sup);
        r.mean_rnk = (r.rnk_pv + r.rnk_or + r.rnk_sup) as f64 / 3.0;
    }
}

/// Brent's root-finding method on the interval [a, b].
///
/// `f` must be continuous and `f(a)` and `f(b)` must have opposite signs
/// (or one must be zero).  Returns the root to within `tol`.
fn brent<F: Fn(f64) -> f64>(f: F, mut a: f64, mut b: f64, tol: f64, max_iter: usize) -> f64 {
    let mut fa = f(a);
    let mut fb = f(b);

    // If one endpoint is already a root, return it.
    if fa.abs() < tol {
        return a;
    }
    if fb.abs() < tol {
        return b;
    }

    // Ensure fa and fb have opposite signs; if not, fall back to midpoint.
    if fa * fb > 0.0 {
        return (a + b) / 2.0;
    }

    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;

    for _ in 0..max_iter {
        if fb * fc > 0.0 {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        let tol1 = 2.0 * f64::EPSILON * b.abs() + 0.5 * tol;
        let m = 0.5 * (c - b);

        if m.abs() <= tol1 || fb == 0.0 {
            return b;
        }

        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            // Attempt inverse quadratic interpolation
            let s = fb / fa;
            let (p, q) = if (a - c).abs() < f64::EPSILON {
                let p = 2.0 * m * s;
                let q = 1.0 - s;
                (p, q)
            } else {
                let q_val = fa / fc;
                let r = fb / fc;
                let p = s * (2.0 * m * q_val * (q_val - r) - (b - a) * (r - 1.0));
                let q = (q_val - 1.0) * (r - 1.0) * (s - 1.0);
                (p, q)
            };

            let (p, q) = if p > 0.0 { (p, -q) } else { (-p, q) };

            if 2.0 * p < (3.0 * m * q - (tol1 * q).abs()).min(e * q) {
                e = d;
                d = p / q;
            } else {
                d = m;
                e = m;
            }
        } else {
            d = m;
            e = m;
        }

        a = b;
        fa = fb;

        if d.abs() > tol1 {
            b += d;
        } else {
            b += if m > 0.0 { tol1 } else { -tol1 };
        }
        fb = f(b);
    }

    b
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
    fn test_contingency_odds_ratio_cmle() {
        // Compare CMLE against R's fisher.test()$estimate
        // R: fisher.test(matrix(c(10,30,20,40), nrow=2))$estimate = 0.6693434
        let ct = ContingencyTable {
            a: 10,
            b: 20,
            c: 30,
            d: 40,
        };
        let or = ct.odds_ratio();
        assert!(
            (or - 0.6693434).abs() < 0.001,
            "CMLE odds ratio should be ~0.6693, got {}",
            or
        );

    }

    #[test]
    fn test_contingency_odds_ratio_boundary() {
        // a at upper boundary (a = hi = min(k, m)) → Infinity
        let ct = ContingencyTable {
            a: 10,
            b: 0,
            c: 5,
            d: 100,
        };
        assert_eq!(ct.odds_ratio(), f64::INFINITY);

        // a at lower boundary (a = lo = max(0, k-n)) → 0
        let ct2 = ContingencyTable {
            a: 0,
            b: 5,
            c: 10,
            d: 100,
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
        // Very significant result — sf() avoids precision loss
        let ct = ContingencyTable {
            a: 50,
            b: 10,
            c: 5,
            d: 1000,
        };
        let pvl = ct.p_value_log(Direction::Enrichment);
        // Should be a large finite value, not infinity
        assert!(pvl > 30.0, "expected large -log10(p), got {}", pvl);
        assert!(pvl.is_finite(), "expected finite, got {}", pvl);
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
                ..Default::default()
            },
            LolaResult {
                user_set: 0,
                db_set: 1,
                p_value_log: 10.0,
                odds_ratio: 1.0,
                support: 200,
                filename: "b.bed".into(),
                ..Default::default()
            },
            LolaResult {
                user_set: 0,
                db_set: 2,
                p_value_log: 3.0,
                odds_ratio: 5.0,
                support: 50,
                filename: "c.bed".into(),
                ..Default::default()
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
    fn test_run_lola_binary_counting() {
        // One query region overlaps 3 DB regions in same file.
        // Support should be 1 (binary per query region), matching R LOLA semantics.
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
        // 1 query region overlaps db0 → support = 1 (not 3)
        assert_eq!(results[0].support, 1);
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
    fn test_empty_user_set() {
        // User set with 0 regions should return results with 0 support
        let sets = vec![(
            "db0.bed".to_string(),
            vec![("chr1".to_string(), 100, 200)],
        )];
        let igd = Igd::from_region_sets(sets);

        let user = make_region_set(vec![]); // empty user set
        let universe = make_region_set(vec![
            make_region("chr1", 50, 250),
            make_region("chr1", 400, 500),
        ]);

        let results = run_lola(&igd, &[user], &universe, &LolaConfig::default()).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].support, 0);
        assert_eq!(results[0].p_value_log, 0.0); // no overlap → p=1.0, -log10(1)=0
        // c should be 0 (user_set_size - a = 0 - 0)
        assert_eq!(results[0].c, 0);
    }

    #[test]
    fn test_single_region_user_set() {
        // User set with exactly 1 region
        let sets = vec![(
            "db0.bed".to_string(),
            vec![
                ("chr1".to_string(), 100, 200),
                ("chr1".to_string(), 300, 400),
            ],
        )];
        let igd = Igd::from_region_sets(sets);

        let user = make_region_set(vec![make_region("chr1", 150, 180)]);
        let universe = make_region_set(vec![
            make_region("chr1", 50, 250),
            make_region("chr1", 250, 450),
            make_region("chr1", 500, 600),
            make_region("chr1", 700, 800),
            make_region("chr1", 900, 1000),
        ]);

        let results = run_lola(&igd, &[user], &universe, &LolaConfig::default()).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].support, 1); // single region overlaps db0[0]
        // a=1, user_set_size=1, so c = 1-1 = 0
        assert_eq!(results[0].c, 0);
    }

    #[test]
    fn test_universe_equals_user_set() {
        // Edge case where universe == user set (same regions)
        let sets = vec![(
            "db0.bed".to_string(),
            vec![
                ("chr1".to_string(), 100, 200),
                ("chr1".to_string(), 500, 600),
            ],
        )];
        let igd = Igd::from_region_sets(sets);

        let regions = vec![
            make_region("chr1", 50, 250),
            make_region("chr1", 450, 650),
            make_region("chr1", 700, 800),
        ];
        let user = make_region_set(regions.clone());
        let universe = make_region_set(regions);

        let config = LolaConfig::default();
        let results = run_lola(&igd, &[user], &universe, &config).unwrap();
        assert_eq!(results.len(), 1);
        // When user == universe, a should equal universe_hits (all universe hits are user hits)
        // b = universe_hits - a = 0
        assert_eq!(results[0].b, 0);
        // c = user_size - a = universe_size - a
        // d should be 0 since universe_size - a - b - c = universe_size - a - 0 - (universe_size - a) = 0
        assert_eq!(results[0].d, 0);
    }

    #[test]
    fn test_large_universe_tiny_user_set() {
        // Imbalanced case: large universe, small user set
        let sets = vec![(
            "db0.bed".to_string(),
            vec![("chr1".to_string(), 100, 200)],
        )];
        let igd = Igd::from_region_sets(sets);

        // 5 user regions, 1 overlaps db
        let user = make_region_set(vec![
            make_region("chr1", 150, 180),
            make_region("chr2", 100, 200),
            make_region("chr3", 100, 200),
            make_region("chr4", 100, 200),
            make_region("chr5", 100, 200),
        ]);

        // 100 universe regions (only first overlaps db)
        let mut uni_regions = vec![make_region("chr1", 50, 250)];
        for i in 1..100 {
            uni_regions.push(make_region("chr1", (1000 + i * 200) as u32, (1000 + i * 200 + 100) as u32));
        }
        let universe = make_region_set(uni_regions);

        let results = run_lola(&igd, &[user], &universe, &LolaConfig::default()).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].support, 1); // only 1 user region overlaps
        // universe_size=100, user_size=5
        // a=1, b=universe_hits-1, c=5-1=4
        assert_eq!(results[0].c, 4);
        // d = 100 - a - b - c; should be large since universe is much bigger
        assert!(results[0].d > 90, "d should be large, got {}", results[0].d);
    }

    #[test]
    fn test_min_overlap_boundary() {
        // Test regions that overlap by exactly min_overlap bp vs min_overlap-1 bp
        let sets = vec![(
            "db0.bed".to_string(),
            vec![("chr1".to_string(), 100, 200)], // db region: [100, 200)
        )];
        let igd = Igd::from_region_sets(sets);

        let universe = make_region_set(vec![
            make_region("chr1", 0, 500),
            make_region("chr1", 600, 700),
            make_region("chr1", 800, 900),
        ]);

        // User region [190, 210) overlaps db [100,200) by exactly 10bp
        let user_10bp = make_region_set(vec![make_region("chr1", 190, 210)]);

        // min_overlap=10: exactly 10bp overlap → should count
        let config10 = LolaConfig {
            min_overlap: 10,
            ..Default::default()
        };
        let r10 = run_lola(&igd, &[user_10bp.clone()], &universe, &config10).unwrap();
        assert_eq!(r10[0].support, 1, "10bp overlap with min_overlap=10 should count");

        // min_overlap=11: only 10bp overlap → should NOT count
        let config11 = LolaConfig {
            min_overlap: 11,
            ..Default::default()
        };
        let r11 = run_lola(&igd, &[user_10bp], &universe, &config11).unwrap();
        assert_eq!(r11[0].support, 0, "10bp overlap with min_overlap=11 should not count");
    }

    #[test]
    fn test_run_lola_user_not_subset_of_universe_passes_negative() {
        // User set has regions NOT in the universe — should warn and pass through negative b
        let sets = vec![(
            "db0.bed".to_string(),
            vec![("chr1".to_string(), 100, 200)],
        )];
        let igd = Igd::from_region_sets(sets);

        // User has 3 regions that overlap db0, but universe only has 1 that overlaps db0.
        // This means user_hits (3) > universe_hits (1) — b = 1 - 3 = -2.
        let user = make_region_set(vec![
            make_region("chr1", 110, 120),
            make_region("chr1", 130, 140),
            make_region("chr1", 150, 160),
        ]);
        // Universe has only 1 region overlapping db0 and 1 not overlapping
        let universe = make_region_set(vec![
            make_region("chr1", 110, 120),
            make_region("chr1", 500, 600),
        ]);

        let result = run_lola(&igd, &[user], &universe, &LolaConfig::default());
        assert!(result.is_ok(), "Should succeed with negative values, not error");
        let results = result.unwrap();
        assert!(results[0].b < 0, "b should be negative, got {}", results[0].b);
        assert_eq!(results[0].p_value_log, 0.0, "p_value_log should be 0 for negative row");
        assert!(results[0].odds_ratio.is_nan(), "odds_ratio should be NaN for negative row");
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

    #[test]
    fn test_run_lola_depletion_direction() {
        // Run with Direction::Depletion — p-values should differ from enrichment
        let sets = vec![(
            "db0.bed".to_string(),
            vec![
                ("chr1".to_string(), 100, 200),
                ("chr1".to_string(), 300, 400),
            ],
        )];
        let igd = Igd::from_region_sets(sets);

        let user = make_region_set(vec![make_region("chr1", 150, 180)]);
        let universe = make_region_set(vec![
            make_region("chr1", 50, 250),
            make_region("chr1", 250, 450),
            make_region("chr1", 500, 600),
            make_region("chr1", 700, 800),
            make_region("chr1", 900, 1000),
        ]);

        let config_enrich = LolaConfig {
            direction: Direction::Enrichment,
            ..Default::default()
        };
        let config_deplete = LolaConfig {
            direction: Direction::Depletion,
            ..Default::default()
        };

        let r_enrich = run_lola(&igd, &[user.clone()], &universe, &config_enrich).unwrap();
        let r_deplete = run_lola(&igd, &[user], &universe, &config_deplete).unwrap();

        // The p-value logs should differ between enrichment and depletion
        assert_ne!(
            r_enrich[0].p_value_log, r_deplete[0].p_value_log,
            "Enrichment and depletion p-values should differ"
        );
    }

    #[test]
    fn test_odds_ratio_nan_ranking() {
        // When odds_ratio is NaN, it should get worst (highest) rank
        let mut results = vec![
            LolaResult {
                user_set: 0,
                db_set: 0,
                p_value_log: 5.0,
                odds_ratio: f64::NAN,
                support: 10,
                filename: "nan.bed".into(),
                ..Default::default()
            },
            LolaResult {
                user_set: 0,
                db_set: 1,
                p_value_log: 3.0,
                odds_ratio: 2.0,
                support: 20,
                filename: "real.bed".into(),
                ..Default::default()
            },
        ];

        rank_results(&mut results);

        // NaN odds ratio should get worst rank (2 out of 2)
        assert_eq!(results[0].rnk_or, 2, "NaN odds ratio should get worst rank");
        assert_eq!(results[1].rnk_or, 1, "Real odds ratio should get rank 1");
    }
}
