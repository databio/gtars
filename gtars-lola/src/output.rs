//! Output formatting and FDR correction.

use std::io::Write;
use std::path::Path;

use crate::models::LolaResult;

/// Apply Benjamini-Hochberg FDR correction to LOLA results.
///
/// Computes q-values (adjusted p-values) for each user set independently.
/// Modifies results in place, setting the `q_value` field.
pub fn apply_fdr_correction(results: &mut [LolaResult]) {
    if results.is_empty() {
        return;
    }

    // Group by user_set and apply BH correction within each group
    let mut max_user_set = 0;
    for r in results.iter() {
        if r.user_set > max_user_set {
            max_user_set = r.user_set;
        }
    }

    for us in 0..=max_user_set {
        // Collect indices for this user set
        let mut indices: Vec<usize> = results
            .iter()
            .enumerate()
            .filter(|(_, r)| r.user_set == us)
            .map(|(i, _)| i)
            .collect();

        if indices.is_empty() {
            continue;
        }

        let n = indices.len();

        // Sort by p-value (ascending = by -log10(p) descending)
        indices.sort_by(|&a, &b| {
            results[b]
                .p_value_log
                .partial_cmp(&results[a].p_value_log)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Convert -log10(p) back to p-values
        let p_values: Vec<f64> = indices
            .iter()
            .map(|&i| {
                let pvl = results[i].p_value_log;
                if pvl == f64::INFINITY {
                    0.0 // p-value underflowed
                } else {
                    10.0_f64.powf(-pvl)
                }
            })
            .collect();

        // Benjamini-Hochberg procedure:
        // Sort p-values ascending, then q[i] = min(p[i] * n / rank, q[i+1])
        // working from the largest p-value down.

        // indices is already sorted by ascending p-value (descending -log10(p))
        // So p_values[0] is smallest, p_values[n-1] is largest

        let mut q_values = vec![0.0f64; n];
        q_values[n - 1] = (p_values[n - 1] * n as f64 / n as f64).min(1.0);

        for i in (0..n - 1).rev() {
            let rank = i + 1; // 1-based rank
            let q = p_values[i] * n as f64 / rank as f64;
            q_values[i] = q.min(q_values[i + 1]).min(1.0);
        }

        // Assign q-values back
        for (j, &idx) in indices.iter().enumerate() {
            results[idx].q_value = Some(q_values[j]);
        }
    }
}

/// Write LOLA results as TSV matching R LOLA's `writeCombinedEnrichment` format.
pub fn write_results_tsv<W: Write>(
    writer: &mut W,
    results: &[LolaResult],
) -> std::io::Result<()> {
    // Header
    writeln!(
        writer,
        "userSet\tdbSet\tpValueLog\toddsRatio\tsupport\trnkPV\trnkOR\trnkSup\tmaxRnk\tmeanRnk\tb\tc\td\tqValue\tfilename"
    )?;

    for r in results {
        let qv = r
            .q_value
            .map(|q| format!("{:.6e}", q))
            .unwrap_or_else(|| "NA".to_string());

        writeln!(
            writer,
            "{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}",
            r.user_set + 1, // 1-based for R compatibility
            r.db_set + 1,
            r.p_value_log,
            r.odds_ratio,
            r.support,
            r.rnk_pv,
            r.rnk_or,
            r.rnk_sup,
            r.max_rnk,
            r.mean_rnk,
            r.b,
            r.c,
            r.d,
            qv,
            r.filename
        )?;
    }

    Ok(())
}

/// Write results to a TSV file on disk.
pub fn write_results_to_file(
    path: &Path,
    results: &[LolaResult],
) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    write_results_tsv(&mut file, results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::LolaResult;

    fn make_result(user_set: usize, db_set: usize, p_value_log: f64) -> LolaResult {
        LolaResult {
            user_set,
            db_set,
            p_value_log,
            odds_ratio: 1.0,
            support: 10,
            rnk_pv: 1,
            rnk_or: 1,
            rnk_sup: 1,
            max_rnk: 1,
            mean_rnk: 1.0,
            b: 5,
            c: 5,
            d: 100,
            q_value: None,
            filename: format!("file{}.bed", db_set),
        }
    }

    #[test]
    fn test_fdr_basic() {
        // 3 results for same user set, with known p-values
        let mut results = vec![
            make_result(0, 0, 3.0),  // p = 0.001
            make_result(0, 1, 2.0),  // p = 0.01
            make_result(0, 2, 1.0),  // p = 0.1
        ];

        apply_fdr_correction(&mut results);

        // All should have q_values set
        for r in &results {
            assert!(r.q_value.is_some());
        }

        // q-values should be >= p-values
        for r in &results {
            let p = 10.0_f64.powf(-r.p_value_log);
            assert!(
                r.q_value.unwrap() >= p - 1e-10,
                "q={} should be >= p={}",
                r.q_value.unwrap(),
                p
            );
        }

        // q-values should be monotonically non-decreasing when sorted by p-value ascending
        // (i.e., by -log10(p) descending)
        let q0 = results[0].q_value.unwrap(); // most significant
        let q1 = results[1].q_value.unwrap();
        let q2 = results[2].q_value.unwrap(); // least significant
        assert!(q0 <= q1 + 1e-10);
        assert!(q1 <= q2 + 1e-10);
    }

    #[test]
    fn test_fdr_multiple_user_sets() {
        let mut results = vec![
            make_result(0, 0, 5.0),
            make_result(0, 1, 2.0),
            make_result(1, 0, 3.0),
            make_result(1, 1, 1.0),
        ];

        apply_fdr_correction(&mut results);

        // All should have q-values
        for r in &results {
            assert!(r.q_value.is_some());
        }
    }

    #[test]
    fn test_fdr_single_result() {
        let mut results = vec![make_result(0, 0, 5.0)];
        apply_fdr_correction(&mut results);
        assert!(results[0].q_value.is_some());
        // Single result: q = p (no correction needed)
        let p = 10.0_f64.powf(-5.0);
        assert!((results[0].q_value.unwrap() - p).abs() < 1e-10);
    }

    #[test]
    fn test_fdr_empty() {
        let mut results: Vec<LolaResult> = vec![];
        apply_fdr_correction(&mut results); // should not panic
    }

    #[test]
    fn test_fdr_identical_pvalues() {
        // All results have the same p-value (tie handling)
        let mut results = vec![
            make_result(0, 0, 3.0), // p = 0.001
            make_result(0, 1, 3.0), // p = 0.001
            make_result(0, 2, 3.0), // p = 0.001
            make_result(0, 3, 3.0), // p = 0.001
        ];

        apply_fdr_correction(&mut results);

        // All should have q-values set
        for r in &results {
            assert!(r.q_value.is_some(), "q_value should be set");
        }

        // With identical p-values, BH correction should produce q-values
        // that are >= p-values and <= 1.0
        let p = 10.0_f64.powf(-3.0);
        for r in &results {
            let q = r.q_value.unwrap();
            assert!(q >= p - 1e-10, "q={} should be >= p={}", q, p);
            assert!(q <= 1.0, "q={} should be <= 1.0", q);
        }
    }

    #[test]
    fn test_fdr_preserves_order() {
        // After FDR correction, q-values should be monotonically non-decreasing
        // when results are sorted by p-value ascending (i.e., -log10(p) descending)
        let mut results = vec![
            make_result(0, 0, 10.0), // most significant
            make_result(0, 1, 7.0),
            make_result(0, 2, 5.0),
            make_result(0, 3, 3.0),
            make_result(0, 4, 2.0),
            make_result(0, 5, 1.0),
            make_result(0, 6, 0.5),
            make_result(0, 7, 0.1), // least significant
        ];

        apply_fdr_correction(&mut results);

        // Sort by p_value_log descending (= p-value ascending)
        let mut sorted = results.clone();
        sorted.sort_by(|a, b| {
            b.p_value_log
                .partial_cmp(&a.p_value_log)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // q-values should be monotonically non-decreasing
        for i in 1..sorted.len() {
            let q_prev = sorted[i - 1].q_value.unwrap();
            let q_curr = sorted[i].q_value.unwrap();
            assert!(
                q_curr >= q_prev - 1e-10,
                "q-values not monotonic: q[{}]={} < q[{}]={}",
                i,
                q_curr,
                i - 1,
                q_prev
            );
        }
    }

    #[test]
    fn test_fdr_single_very_significant() {
        // One very significant result among many non-significant ones
        let mut results = vec![
            make_result(0, 0, 20.0), // very significant (p ~ 1e-20)
            make_result(0, 1, 0.1),  // not significant
            make_result(0, 2, 0.05), // not significant
            make_result(0, 3, 0.01), // not significant
            make_result(0, 4, 0.0),  // p = 1.0 (no signal)
        ];

        apply_fdr_correction(&mut results);

        // The very significant result should still be significant after correction
        let best = results.iter().find(|r| r.db_set == 0).unwrap();
        assert!(
            best.q_value.unwrap() < 0.05,
            "Very significant result should remain significant after FDR, q={}",
            best.q_value.unwrap()
        );

        // The non-significant results should have q-values >= their p-values
        for r in results.iter().filter(|r| r.db_set != 0) {
            let p = 10.0_f64.powf(-r.p_value_log);
            assert!(
                r.q_value.unwrap() >= p - 1e-10,
                "q={} should be >= p={} for db_set={}",
                r.q_value.unwrap(),
                p,
                r.db_set
            );
        }
    }

    #[test]
    fn test_write_tsv() {
        let mut results = vec![
            make_result(0, 0, 5.0),
            make_result(0, 1, 2.0),
        ];
        apply_fdr_correction(&mut results);

        let mut buf = Vec::new();
        write_results_tsv(&mut buf, &results).unwrap();
        let output = String::from_utf8(buf).unwrap();

        // Check header
        assert!(output.starts_with("userSet\tdbSet\tpValueLog\t"));
        // Check data lines
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 3); // header + 2 data rows
        // Check 1-based indexing
        assert!(lines[1].starts_with("1\t1\t"));
        assert!(lines[2].starts_with("1\t2\t"));
    }

    #[test]
    fn test_write_tsv_no_qvalue() {
        let results = vec![make_result(0, 0, 5.0)]; // q_value is None
        let mut buf = Vec::new();
        write_results_tsv(&mut buf, &results).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("NA")); // q_value should be NA
    }
}
