//! Signal matrix overlap and summary statistics.
//!
//! Ports R's `calcSummarySignal` from GenomicDistributions: overlaps query
//! regions with a signal matrix (TSV of region × condition signal values),
//! aggregates by MAX per query region, and computes Tukey boxplot statistics.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use gtars_core::models::{Region, RegionSet};
use gtars_overlaprs::traits::{Interval, Overlapper};
use gtars_overlaprs::AIList;

use crate::errors::GtarsGenomicDistError;

/// A matrix of signal values across genomic regions and conditions.
///
/// Rows are genomic regions, columns are conditions (e.g. cell types).
/// Loaded from a TSV file where the first column is `chr_start_end`.
pub struct SignalMatrix {
    /// The genomic regions corresponding to each row
    pub regions: RegionSet,
    /// Condition/cell-type names (column headers, excluding first column)
    pub condition_names: Vec<String>,
    /// Signal values: outer vec = rows, inner vec = conditions (same order as condition_names)
    pub values: Vec<Vec<f64>>,
}

/// Result of [`calc_summary_signal`].
pub struct SignalSummaryResult {
    /// Per-query-region signal summaries.
    /// Each entry: (query region label as `chr_start_end`, max signal per condition).
    pub signal_matrix: Vec<(String, Vec<f64>)>,
    /// Boxplot statistics per condition, in same order as `condition_names`.
    pub matrix_stats: Vec<ConditionStats>,
    /// Condition names, for labeling.
    pub condition_names: Vec<String>,
}

/// Tukey boxplot statistics for a single condition (matches R's `boxplot.stats`).
#[derive(Debug, Clone)]
pub struct ConditionStats {
    pub condition: String,
    pub lower_whisker: f64,
    pub lower_hinge: f64,
    pub median: f64,
    pub upper_hinge: f64,
    pub upper_whisker: f64,
}

impl SignalMatrix {
    /// Load a signal matrix from a TSV file.
    ///
    /// Format: first column is region ID as `chr_start_end`, remaining columns
    /// are condition names with float signal values. Rows where the region ID
    /// splits into more than 3 parts on `_` are skipped.
    pub fn from_tsv<P: AsRef<Path>>(path: P) -> Result<Self, GtarsGenomicDistError> {
        let file = File::open(path.as_ref())?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Header line
        let header = lines
            .next()
            .ok_or_else(|| {
                GtarsGenomicDistError::SignalMatrixError("Empty signal matrix file".into())
            })?
            .map_err(|e| {
                GtarsGenomicDistError::SignalMatrixError(format!("Reading header: {}", e))
            })?;

        let header_fields: Vec<&str> = header.split('\t').collect();
        if header_fields.len() < 2 {
            return Err(GtarsGenomicDistError::SignalMatrixError(
                "Signal matrix must have at least 2 columns".into(),
            ));
        }
        let condition_names: Vec<String> = header_fields[1..].iter().map(|s| s.to_string()).collect();
        let n_conditions = condition_names.len();

        let mut regions = Vec::new();
        let mut values = Vec::new();

        for (line_num, line_result) in lines.enumerate() {
            let line = line_result.map_err(|e| {
                GtarsGenomicDistError::SignalMatrixError(format!("Line {}: {}", line_num + 2, e))
            })?;

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.is_empty() {
                continue;
            }

            // Parse region ID: split on '_', take last two as start/end, rest as chr
            let parts: Vec<&str> = fields[0].split('_').collect();
            if parts.len() < 3 || parts.len() > 3 {
                // Skip malformed rows (matching R's behavior of dropping >3 splits)
                continue;
            }

            let chr = parts[0].to_string();
            let start: u32 = match parts[1].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            let end: u32 = match parts[2].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };

            // Parse signal values
            if fields.len() < 1 + n_conditions {
                continue;
            }
            let row_values: Result<Vec<f64>, _> =
                fields[1..=n_conditions].iter().map(|s| s.parse::<f64>()).collect();
            let row_values = match row_values {
                Ok(v) => v,
                Err(_) => continue,
            };

            regions.push(Region {
                chr,
                start,
                end,
                rest: None,
            });
            values.push(row_values);
        }

        if regions.is_empty() {
            return Err(GtarsGenomicDistError::SignalMatrixError(
                "No valid rows in signal matrix".into(),
            ));
        }

        Ok(SignalMatrix {
            regions: RegionSet::from(regions),
            condition_names,
            values,
        })
    }
}

/// Overlap query regions with a signal matrix and compute per-condition statistics.
///
/// For each query region that overlaps at least one signal matrix region, takes
/// the MAX signal value per condition across all overlapping rows. Then computes
/// Tukey boxplot statistics per condition on the aggregated values.
///
/// Uses per-chromosome AIList indexes for O(Q × log S) overlap queries, where
/// Q = number of query regions and S = signal regions per chromosome.
pub fn calc_summary_signal(
    query: &RegionSet,
    signal_matrix: &SignalMatrix,
) -> Result<SignalSummaryResult, GtarsGenomicDistError> {
    let n_conditions = signal_matrix.condition_names.len();

    // Step 1: Build per-chromosome AIList with row indices in val
    let mut chr_intervals: HashMap<String, Vec<Interval<u32, usize>>> = HashMap::new();
    for (idx, region) in signal_matrix.regions.regions.iter().enumerate() {
        chr_intervals
            .entry(region.chr.clone())
            .or_default()
            .push(Interval {
                start: region.start,
                end: region.end,
                val: idx,
            });
    }
    let chr_ailist: HashMap<String, AIList<u32, usize>> = chr_intervals
        .into_iter()
        .map(|(chr, intervals)| (chr, AIList::build(intervals)))
        .collect();

    // Step 2: For each query region, find overlapping signal rows and take MAX per condition
    let mut signal_rows: Vec<(String, Vec<f64>)> = Vec::new();

    for query_region in &query.regions {
        if let Some(ailist) = chr_ailist.get(&query_region.chr) {
            let mut max_vals: Option<Vec<f64>> = None;

            for hit in ailist.find_iter(query_region.start, query_region.end) {
                let signal_values = &signal_matrix.values[hit.val];
                match &mut max_vals {
                    Some(existing) => {
                        for (ci, val) in signal_values.iter().enumerate() {
                            if *val > existing[ci] {
                                existing[ci] = *val;
                            }
                        }
                    }
                    None => {
                        max_vals = Some(signal_values.clone());
                    }
                }
            }

            if let Some(vals) = max_vals {
                let label = format!(
                    "{}_{}_{}",
                    query_region.chr, query_region.start, query_region.end
                );
                signal_rows.push((label, vals));
            }
        }
    }

    // Step 3: Compute boxplot stats per condition
    let matrix_stats = if signal_rows.is_empty() {
        Vec::new()
    } else {
        let mut condition_columns: Vec<Vec<f64>> = vec![Vec::new(); n_conditions];
        for (_label, vals) in &signal_rows {
            for (ci, v) in vals.iter().enumerate() {
                condition_columns[ci].push(*v);
            }
        }

        condition_columns
            .iter_mut()
            .enumerate()
            .map(|(ci, col)| {
                let mut stats = boxplot_stats(col);
                stats.condition = signal_matrix.condition_names[ci].clone();
                stats
            })
            .collect()
    };

    Ok(SignalSummaryResult {
        signal_matrix: signal_rows,
        matrix_stats,
        condition_names: signal_matrix.condition_names.clone(),
    })
}

/// Compute Tukey boxplot statistics matching R's `boxplot.stats` / `fivenum`.
///
/// Sorts the data in place, computes hinges using R's `fivenum` algorithm
/// (median of each half, including the median element for odd n), then
/// computes whiskers as the most extreme data points within 1.5 × IQR fences.
fn boxplot_stats(data: &mut [f64]) -> ConditionStats {
    data.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = data.len();

    if n == 0 {
        return ConditionStats {
            condition: String::new(),
            lower_whisker: 0.0,
            lower_hinge: 0.0,
            median: 0.0,
            upper_hinge: 0.0,
            upper_whisker: 0.0,
        };
    }

    let median = fivenum_median(data);

    // R's fivenum: lower half includes median for odd n, upper half always starts at mid
    let mid = n / 2;
    let lower_half = if n % 2 == 0 {
        &data[..mid]
    } else {
        &data[..=mid]
    };
    let upper_half = &data[mid..];
    let lower_hinge = fivenum_median(lower_half);
    let upper_hinge = fivenum_median(upper_half);

    let iqr = upper_hinge - lower_hinge;
    let lower_fence = lower_hinge - 1.5 * iqr;
    let upper_fence = upper_hinge + 1.5 * iqr;

    let lower_whisker = data
        .iter()
        .copied()
        .find(|&x| x >= lower_fence)
        .unwrap_or(lower_hinge);
    let upper_whisker = data
        .iter()
        .rev()
        .copied()
        .find(|&x| x <= upper_fence)
        .unwrap_or(upper_hinge);

    ConditionStats {
        condition: String::new(),
        lower_whisker,
        lower_hinge,
        median,
        upper_hinge,
        upper_whisker,
    }
}

/// Median of a sorted slice, matching R's fivenum internal logic.
fn fivenum_median(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_test_signal_matrix() -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "V1\tcond_A\tcond_B\tcond_C").unwrap();
        writeln!(f, "chr1_100_200\t0.5\t0.3\t0.1").unwrap();
        writeln!(f, "chr1_150_250\t0.2\t0.8\t0.4").unwrap();
        writeln!(f, "chr1_300_400\t0.9\t0.1\t0.7").unwrap();
        writeln!(f, "chr2_100_200\t0.3\t0.6\t0.2").unwrap();
        f.flush().unwrap();
        f
    }

    #[rstest]
    fn test_signal_matrix_load() {
        let f = write_test_signal_matrix();
        let sm = SignalMatrix::from_tsv(f.path()).unwrap();

        assert_eq!(sm.condition_names, vec!["cond_A", "cond_B", "cond_C"]);
        assert_eq!(sm.values.len(), 4);
        assert_eq!(sm.regions.len(), 4);
        assert!((sm.values[0][0] - 0.5).abs() < 1e-9);
    }

    #[rstest]
    fn test_signal_matrix_skips_malformed_rows() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "V1\tcond_A").unwrap();
        writeln!(f, "chr1_100_200\t0.5").unwrap();
        writeln!(f, "chr1_random_100_200\t0.3").unwrap(); // >3 parts, skipped
        writeln!(f, "bad_row\t0.1").unwrap(); // <3 parts, skipped
        f.flush().unwrap();

        let sm = SignalMatrix::from_tsv(f.path()).unwrap();
        assert_eq!(sm.values.len(), 1);
    }

    #[rstest]
    fn test_boxplot_stats_known_values() {
        // R: boxplot.stats(c(1,2,3,4,5))
        // $stats = [1, 2, 3, 4, 5] (lower_whisker, Q1, median, Q3, upper_whisker)
        let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let stats = boxplot_stats(&mut data);

        assert!((stats.median - 3.0).abs() < 1e-9);
        assert!((stats.lower_hinge - 2.0).abs() < 1e-9);
        assert!((stats.upper_hinge - 4.0).abs() < 1e-9);
        assert!((stats.lower_whisker - 1.0).abs() < 1e-9);
        assert!((stats.upper_whisker - 5.0).abs() < 1e-9);
    }

    #[rstest]
    fn test_boxplot_stats_even_count() {
        // R: boxplot.stats(c(1,2,3,4,5,6))
        // $stats = [1, 2, 3.5, 5, 6]
        let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let stats = boxplot_stats(&mut data);

        assert!((stats.median - 3.5).abs() < 1e-9);
        assert!((stats.lower_hinge - 2.0).abs() < 1e-9);
        assert!((stats.upper_hinge - 5.0).abs() < 1e-9);
        assert!((stats.lower_whisker - 1.0).abs() < 1e-9);
        assert!((stats.upper_whisker - 6.0).abs() < 1e-9);
    }

    #[rstest]
    fn test_boxplot_stats_with_outliers() {
        // R: boxplot.stats(c(1,2,3,4,5,100))
        // IQR = 5-2 = 3, fence = 5 + 4.5 = 9.5
        // $stats = [1, 2, 3.5, 5, 5]  (100 is outlier)
        let mut data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 100.0];
        let stats = boxplot_stats(&mut data);

        assert!((stats.upper_whisker - 5.0).abs() < 1e-9);
    }

    #[rstest]
    fn test_calc_summary_signal_end_to_end() {
        let f = write_test_signal_matrix();
        let sm = SignalMatrix::from_tsv(f.path()).unwrap();

        // Query regions:
        // chr1:120-180 overlaps chr1_100_200 and chr1_150_250
        // chr1:350-380 overlaps chr1_300_400
        // chr2:500-600 has no overlaps
        let query = RegionSet::from(vec![
            Region {
                chr: "chr1".to_string(),
                start: 120,
                end: 180,
                rest: None,
            },
            Region {
                chr: "chr1".to_string(),
                start: 350,
                end: 380,
                rest: None,
            },
            Region {
                chr: "chr2".to_string(),
                start: 500,
                end: 600,
                rest: None,
            },
        ]);

        let result = calc_summary_signal(&query, &sm).unwrap();

        // 2 query regions had overlaps
        assert_eq!(result.signal_matrix.len(), 2);

        // chr1:120-180 → max of rows chr1_100_200 and chr1_150_250
        // cond_A: max(0.5, 0.2) = 0.5
        // cond_B: max(0.3, 0.8) = 0.8
        // cond_C: max(0.1, 0.4) = 0.4
        let (label0, vals0) = &result.signal_matrix[0];
        assert_eq!(label0, "chr1_120_180");
        assert!((vals0[0] - 0.5).abs() < 1e-9);
        assert!((vals0[1] - 0.8).abs() < 1e-9);
        assert!((vals0[2] - 0.4).abs() < 1e-9);

        // chr1:350-380 → chr1_300_400 only
        // cond_A: 0.9, cond_B: 0.1, cond_C: 0.7
        let (label1, vals1) = &result.signal_matrix[1];
        assert_eq!(label1, "chr1_350_380");
        assert!((vals1[0] - 0.9).abs() < 1e-9);
        assert!((vals1[1] - 0.1).abs() < 1e-9);
        assert!((vals1[2] - 0.7).abs() < 1e-9);

        // Boxplot stats computed on 2-element columns
        assert_eq!(result.matrix_stats.len(), 3);
        assert_eq!(result.matrix_stats[0].condition, "cond_A");
    }

    #[rstest]
    fn test_calc_summary_signal_no_overlaps() {
        let f = write_test_signal_matrix();
        let sm = SignalMatrix::from_tsv(f.path()).unwrap();

        let query = RegionSet::from(vec![Region {
            chr: "chr3".to_string(),
            start: 100,
            end: 200,
            rest: None,
        }]);

        let result = calc_summary_signal(&query, &sm).unwrap();
        assert!(result.signal_matrix.is_empty());
        assert!(result.matrix_stats.is_empty());
    }

    #[rstest]
    fn test_fivenum_median() {
        assert!((fivenum_median(&[1.0, 2.0, 3.0]) - 2.0).abs() < 1e-9);
        assert!((fivenum_median(&[1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-9);
        assert!((fivenum_median(&[5.0]) - 5.0).abs() < 1e-9);
        assert!((fivenum_median(&[]) - 0.0).abs() < 1e-9);
    }
}
