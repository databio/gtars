//! Output formatting, FDR correction, and annotation.

use std::io::Write;

use crate::database::RegionDB;
use crate::models::LolaResult;

/// Attach DB metadata (collection, description, cellType, etc.) to results.
///
/// Uses each result's `db_set` index to look up annotations from the RegionDB.
/// Also fills in `db_set_size` from the original region sets.
pub fn annotate_results(results: &mut [LolaResult], db: &RegionDB) {
    for r in results.iter_mut() {
        if r.db_set < db.region_anno.len() {
            let anno = &db.region_anno[r.db_set];
            r.collection = anno.collection.clone();
            // Truncate description to 80 chars (matches R LOLA behavior)
            r.description = anno.description.as_ref().map(|d| d.chars().take(80).collect());
            r.cell_type = anno.cell_type.clone();
            r.tissue = anno.tissue.clone();
            r.antibody = anno.antibody.clone();
            r.treatment = anno.treatment.clone();
            r.data_source = anno.data_source.clone();
        }
        if r.db_set < db.region_sets.len() {
            r.db_set_size = db.region_sets[r.db_set].regions.len() as u64;
        }
    }
}

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

/// Column-oriented representation of LOLA results.
///
/// Each field is a parallel Vec — row `i` across all fields describes one result.
/// Bindings should convert this to their native columnar type (JS object, PyDict,
/// R data.frame) rather than reimplementing the row→column pivot.
#[derive(Debug, Clone)]
pub struct LolaColumnar {
    pub user_set: Vec<usize>,
    pub db_set: Vec<usize>,
    pub p_value_log: Vec<f64>,
    pub odds_ratio: Vec<f64>,
    pub support: Vec<u64>,
    pub rnk_pv: Vec<usize>,
    pub rnk_or: Vec<usize>,
    pub rnk_sup: Vec<usize>,
    pub max_rnk: Vec<usize>,
    pub mean_rnk: Vec<f64>,
    pub b: Vec<i64>,
    pub c: Vec<i64>,
    pub d: Vec<i64>,
    pub q_value: Vec<Option<f64>>,
    pub filename: Vec<String>,
    pub collection: Vec<Option<String>>,
    pub description: Vec<Option<String>>,
    pub cell_type: Vec<Option<String>>,
    pub tissue: Vec<Option<String>>,
    pub antibody: Vec<Option<String>>,
    pub treatment: Vec<Option<String>>,
    pub data_source: Vec<Option<String>>,
    pub db_set_size: Vec<u64>,
}

/// Convert a slice of LolaResults into column-oriented vectors.
pub fn results_to_columns(results: &[LolaResult]) -> LolaColumnar {
    let n = results.len();
    let mut c = LolaColumnar {
        user_set: Vec::with_capacity(n),
        db_set: Vec::with_capacity(n),
        p_value_log: Vec::with_capacity(n),
        odds_ratio: Vec::with_capacity(n),
        support: Vec::with_capacity(n),
        rnk_pv: Vec::with_capacity(n),
        rnk_or: Vec::with_capacity(n),
        rnk_sup: Vec::with_capacity(n),
        max_rnk: Vec::with_capacity(n),
        mean_rnk: Vec::with_capacity(n),
        b: Vec::with_capacity(n),
        c: Vec::with_capacity(n),
        d: Vec::with_capacity(n),
        q_value: Vec::with_capacity(n),
        filename: Vec::with_capacity(n),
        collection: Vec::with_capacity(n),
        description: Vec::with_capacity(n),
        cell_type: Vec::with_capacity(n),
        tissue: Vec::with_capacity(n),
        antibody: Vec::with_capacity(n),
        treatment: Vec::with_capacity(n),
        data_source: Vec::with_capacity(n),
        db_set_size: Vec::with_capacity(n),
    };
    for r in results {
        c.user_set.push(r.user_set);
        c.db_set.push(r.db_set);
        c.p_value_log.push(r.p_value_log);
        c.odds_ratio.push(r.odds_ratio);
        c.support.push(r.support);
        c.rnk_pv.push(r.rnk_pv);
        c.rnk_or.push(r.rnk_or);
        c.rnk_sup.push(r.rnk_sup);
        c.max_rnk.push(r.max_rnk);
        c.mean_rnk.push(r.mean_rnk);
        c.b.push(r.b);
        c.c.push(r.c);
        c.d.push(r.d);
        c.q_value.push(r.q_value);
        c.filename.push(r.filename.clone());
        c.collection.push(r.collection.clone());
        c.description.push(r.description.clone());
        c.cell_type.push(r.cell_type.clone());
        c.tissue.push(r.tissue.clone());
        c.antibody.push(r.antibody.clone());
        c.treatment.push(r.treatment.clone());
        c.data_source.push(r.data_source.clone());
        c.db_set_size.push(r.db_set_size);
    }
    c
}

/// Write LOLA results as TSV matching R LOLA's `writeCombinedEnrichment` format.
pub fn write_results_tsv<W: Write>(
    writer: &mut W,
    results: &[LolaResult],
) -> std::io::Result<()> {
    // Header
    writeln!(
        writer,
        "userSet\tdbSet\tcollection\tpValueLog\toddsRatio\tsupport\t\
         rnkPV\trnkOR\trnkSup\tmaxRnk\tmeanRnk\tb\tc\td\t\
         description\tcellType\ttissue\tantibody\ttreatment\tdataSource\t\
         filename\tqValue\tsize"
    )?;

    for r in results {
        let qv = r
            .q_value
            .map(|q| format!("{:.6e}", q))
            .unwrap_or_else(|| "NA".to_string());
        writeln!(
            writer,
            "{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t\
             {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.user_set + 1, // 1-based for R compatibility
            r.db_set + 1,
            r.collection.as_deref().unwrap_or(""),
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
            r.description.as_deref().unwrap_or(""),
            r.cell_type.as_deref().unwrap_or(""),
            r.tissue.as_deref().unwrap_or(""),
            r.antibody.as_deref().unwrap_or(""),
            r.treatment.as_deref().unwrap_or(""),
            r.data_source.as_deref().unwrap_or(""),
            r.filename,
            qv,
            r.db_set_size,
        )?;
    }

    Ok(())
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
            collection: None,
            description: None,
            cell_type: None,
            tissue: None,
            antibody: None,
            treatment: None,
            data_source: None,
            db_set_size: 0,
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
        assert!(output.starts_with("userSet\tdbSet\tcollection\tpValueLog\t"));
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

    #[test]
    fn test_results_to_columns_basic() {
        let results = vec![
            make_result(0, 0, 3.0),
            make_result(1, 2, 5.0),
        ];
        let c = results_to_columns(&results);

        assert_eq!(c.user_set.len(), 2);
        assert_eq!(c.user_set, vec![0, 1]);
        assert_eq!(c.db_set, vec![0, 2]);
        assert_eq!(c.p_value_log, vec![3.0, 5.0]);
        assert_eq!(c.odds_ratio, vec![1.0, 1.0]);
        assert_eq!(c.support, vec![10, 10]);
        assert_eq!(c.b, vec![5, 5]);
        assert_eq!(c.c, vec![5, 5]);
        assert_eq!(c.d, vec![100, 100]);
        assert_eq!(c.filename, vec!["file0.bed", "file2.bed"]);
        assert_eq!(c.db_set_size, vec![0, 0]);
        // empty strings → None
        assert_eq!(c.collection, vec![None, None]);
        assert_eq!(c.description, vec![None, None]);
        assert_eq!(c.cell_type, vec![None, None]);
        assert_eq!(c.tissue, vec![None, None]);
        assert_eq!(c.antibody, vec![None, None]);
        assert_eq!(c.treatment, vec![None, None]);
        assert_eq!(c.data_source, vec![None, None]);
    }

    #[test]
    fn test_results_to_columns_empty() {
        let c = results_to_columns(&[]);
        assert!(c.user_set.is_empty());
        assert!(c.filename.is_empty());
    }

    #[test]
    fn test_results_to_columns_with_metadata() {
        let mut r = make_result(0, 0, 1.0);
        r.collection = Some("ENCODE".to_string());
        r.cell_type = Some("K562".to_string());
        r.tissue = None; // stays None
        r.q_value = Some(0.05);

        let c = results_to_columns(&[r]);
        assert_eq!(c.collection, vec![Some("ENCODE".to_string())]);
        assert_eq!(c.cell_type, vec![Some("K562".to_string())]);
        assert_eq!(c.tissue, vec![None]);
        assert_eq!(c.q_value, vec![Some(0.05)]);
    }
}
