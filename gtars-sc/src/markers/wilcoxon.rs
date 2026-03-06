//! Wilcoxon rank-sum test for differential expression / differential accessibility.
//!
//! Implements Seurat-style `FindAllMarkers`:
//! 1. For each cluster vs rest, for each gene/feature:
//!    - Filter by min_pct and min_log2fc
//!    - Wilcoxon rank-sum test (normal approximation)
//! 2. BH (Benjamini-Hochberg) correction across all tests
//!
//! Operates on the sparse log-normalized matrix directly via CSR conversion.

use anyhow::Result;

use crate::types::{FeatureMatrix, MarkerConfig, MarkerResult};

/// Find marker genes/features for all clusters using Wilcoxon rank-sum test.
///
/// Tests each cluster against all other cells. Filters by `min_pct` (minimum
/// fraction of cells expressing the feature) and `min_log2fc`. Applies BH
/// correction for multiple testing.
pub fn find_all_markers(
    matrix: &FeatureMatrix,
    clusters: &[u32],
    config: &MarkerConfig,
) -> Result<Vec<MarkerResult>> {
    let n_cells = matrix.n_cells();
    let n_features = matrix.n_features();

    if n_cells == 0 || n_features == 0 {
        return Ok(Vec::new());
    }

    let n_clusters = *clusters.iter().max().unwrap_or(&0) as usize + 1;

    // Build cluster membership lists
    let mut cluster_cells: Vec<Vec<usize>> = vec![Vec::new(); n_clusters];
    for (i, &c) in clusters.iter().enumerate() {
        cluster_cells[c as usize].push(i);
    }

    // Convert to CSR for efficient per-gene (row) iteration
    let csr = matrix.matrix.to_csr();

    let mut all_results: Vec<MarkerResult> = Vec::new();

    for gene_idx in 0..n_features {
        // Extract expression values for this gene across all cells
        let row = csr.outer_view(gene_idx);
        let gene_values: Vec<(usize, f64)> = match row {
            Some(rv) => rv.iter().map(|(col, &val)| (col, val)).collect(),
            None => Vec::new(),
        };

        // Build a dense-ish lookup: cell_idx → value (0.0 if absent)
        let mut values = vec![0.0f64; n_cells];
        for &(col, val) in &gene_values {
            values[col] = val;
        }

        for cluster_id in 0..n_clusters {
            let in_cells = &cluster_cells[cluster_id];
            if in_cells.is_empty() {
                continue;
            }

            // Collect values for in-group and out-group
            let in_values: Vec<f64> = in_cells.iter().map(|&c| values[c]).collect();
            let out_values: Vec<f64> = (0..n_cells)
                .filter(|c| clusters[*c] as usize != cluster_id)
                .map(|c| values[c])
                .collect();

            if out_values.is_empty() {
                continue;
            }

            let n1 = in_values.len();
            let n2 = out_values.len();

            // pct expressing (value > 0)
            let pct_in = in_values.iter().filter(|&&v| v > 0.0).count() as f64 / n1 as f64;
            let pct_out = out_values.iter().filter(|&&v| v > 0.0).count() as f64 / n2 as f64;

            // Filter: at least min_pct expressing in one group
            if pct_in.max(pct_out) < config.min_pct {
                continue;
            }

            // avg_log2fc: Seurat-compatible log fold change.
            // Reverse log1p to get linear-scale values, compute means, then log back.
            // This avoids Jensen's inequality bias from averaging on the log scale.
            //   fc = log2(mean(expm1(x_in)) + 1) - log2(mean(expm1(x_out)) + 1)
            let mean_linear_in =
                in_values.iter().map(|&v| v.exp_m1()).sum::<f64>() / n1 as f64;
            let mean_linear_out =
                out_values.iter().map(|&v| v.exp_m1()).sum::<f64>() / n2 as f64;
            let avg_log2fc = ((mean_linear_in + 1.0).ln() - (mean_linear_out + 1.0).ln())
                / std::f64::consts::LN_2;

            // Filter by fold change
            if config.only_positive && avg_log2fc < config.min_log2fc {
                continue;
            }
            if !config.only_positive && avg_log2fc.abs() < config.min_log2fc {
                continue;
            }

            // Wilcoxon rank-sum test
            let pval = wilcoxon_rank_sum(&in_values, &out_values);

            all_results.push(MarkerResult {
                gene: matrix.feature_names[gene_idx].clone(),
                cluster: cluster_id as u32,
                avg_log2fc,
                pval,
                pval_adj: pval, // will be corrected below
                pct_in,
                pct_out,
            });
        }
    }

    // BH correction
    benjamini_hochberg(&mut all_results);

    // Sort by cluster, then by p-value
    all_results.sort_by(|a, b| {
        a.cluster
            .cmp(&b.cluster)
            .then(a.pval.partial_cmp(&b.pval).unwrap_or(std::cmp::Ordering::Equal))
    });

    Ok(all_results)
}

/// Wilcoxon rank-sum test (Mann-Whitney U) with normal approximation.
///
/// Returns a two-sided p-value. Uses normal approximation with continuity
/// correction for n1, n2 > 10.
fn wilcoxon_rank_sum(group1: &[f64], group2: &[f64]) -> f64 {
    let n1 = group1.len();
    let n2 = group2.len();

    if n1 == 0 || n2 == 0 {
        return 1.0;
    }

    // Combine and rank
    let mut combined: Vec<(f64, bool)> = Vec::with_capacity(n1 + n2);
    for &v in group1 {
        combined.push((v, true)); // in-group
    }
    for &v in group2 {
        combined.push((v, false)); // out-group
    }
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    // Assign ranks (handle ties by averaging)
    let n = combined.len();
    let mut ranks = vec![0.0f64; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j < n && combined[j].0 == combined[i].0 {
            j += 1;
        }
        // Positions i..j have the same value → average rank
        let avg_rank = (i + 1 + j) as f64 / 2.0; // 1-based
        for k in i..j {
            ranks[k] = avg_rank;
        }
        i = j;
    }

    // Sum of ranks for group 1
    let r1: f64 = ranks
        .iter()
        .zip(combined.iter())
        .filter(|(_, c)| c.1)
        .map(|(r, _)| r)
        .sum();

    // U statistic
    let u = r1 - (n1 as f64 * (n1 as f64 + 1.0)) / 2.0;

    // Expected value and variance under null
    let expected = n1 as f64 * n2 as f64 / 2.0;

    // Variance with tie correction
    let n_f = n as f64;
    let mut tie_correction = 0.0f64;
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j < n && combined[j].0 == combined[i].0 {
            j += 1;
        }
        let t = (j - i) as f64;
        if t > 1.0 {
            tie_correction += t * t * t - t;
        }
        i = j;
    }

    let variance = (n1 as f64 * n2 as f64 / 12.0)
        * ((n_f + 1.0) - tie_correction / (n_f * (n_f - 1.0)));

    if variance <= 0.0 {
        return 1.0;
    }

    // Normal approximation with continuity correction
    let z = (u - expected).abs() - 0.5; // continuity correction
    let z = z / variance.sqrt();

    // Two-sided p-value
    2.0 * norm_cdf(-z.abs())
}

/// Normal CDF via error function approximation.
/// Abramowitz and Stegun formula 7.1.26, max error ~1.5e-7.
fn norm_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2))
}

fn erf(x: f64) -> f64 {
    let a = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * a);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    let result = 1.0 - poly * (-a * a).exp();
    if x >= 0.0 {
        result
    } else {
        -result
    }
}

/// Benjamini-Hochberg FDR correction (in-place on MarkerResult.pval_adj).
fn benjamini_hochberg(results: &mut [MarkerResult]) {
    let n = results.len();
    if n == 0 {
        return;
    }

    // Sort by p-value descending
    let mut indexed: Vec<(usize, f64)> = results.iter().enumerate().map(|(i, r)| (i, r.pval)).collect();
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    // Compute adjusted p-values with cumulative minimum
    let mut cummin = f64::INFINITY;
    for (rank_from_top, &(orig_idx, pval)) in indexed.iter().enumerate() {
        let rank = n - rank_from_top; // 1-based rank from bottom
        let adjusted = (pval * n as f64 / rank as f64).min(1.0);
        cummin = cummin.min(adjusted);
        results[orig_idx].pval_adj = cummin;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::{CsMat, TriMat};

    /// Build a simple FeatureMatrix from dense data.
    fn make_feature_matrix(
        data: &[Vec<f64>],
        feature_names: Vec<String>,
    ) -> FeatureMatrix {
        use crate::types::FeatureType;
        let n_features = data.len();
        let n_cells = data[0].len();
        let mut tri = TriMat::new((n_features, n_cells));
        for (r, row) in data.iter().enumerate() {
            for (c, &val) in row.iter().enumerate() {
                if val != 0.0 {
                    tri.add_triplet(r, c, val);
                }
            }
        }
        let mat: CsMat<f64> = tri.to_csc();
        FeatureMatrix {
            matrix: mat,
            feature_names: feature_names.clone(),
            feature_ids: feature_names,
            cell_ids: (0..n_cells).map(|i| format!("cell_{i}")).collect(),
            feature_type: FeatureType::Gene,
        }
    }

    #[test]
    fn test_planted_markers() {
        // 4 genes, 6 cells, 2 clusters
        // Gene A: high in cluster 0, low in cluster 1
        // Gene B: high in cluster 1, low in cluster 0
        // Gene C: expressed equally → not a marker
        // Gene D: not expressed → not a marker
        let fm = make_feature_matrix(
            &[
                vec![5.0, 5.0, 5.0, 0.1, 0.1, 0.1], // Gene A — marker for cluster 0
                vec![0.1, 0.1, 0.1, 5.0, 5.0, 5.0], // Gene B — marker for cluster 1
                vec![2.0, 2.0, 2.0, 2.0, 2.0, 2.0], // Gene C — not DE
                vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0], // Gene D — not expressed
            ],
            vec!["A".into(), "B".into(), "C".into(), "D".into()],
        );
        let clusters = vec![0, 0, 0, 1, 1, 1];
        let config = MarkerConfig {
            min_pct: 0.0,
            min_log2fc: 0.0,
            only_positive: true,
        };

        let markers = find_all_markers(&fm, &clusters, &config).unwrap();

        // Should find A as marker for cluster 0, B as marker for cluster 1
        let cluster0_markers: Vec<&str> = markers
            .iter()
            .filter(|m| m.cluster == 0)
            .map(|m| m.gene.as_str())
            .collect();
        let cluster1_markers: Vec<&str> = markers
            .iter()
            .filter(|m| m.cluster == 1)
            .map(|m| m.gene.as_str())
            .collect();

        assert!(cluster0_markers.contains(&"A"), "Gene A should be marker for cluster 0");
        assert!(cluster1_markers.contains(&"B"), "Gene B should be marker for cluster 1");
    }

    #[test]
    fn test_bh_correction() {
        let mut results = vec![
            MarkerResult {
                gene: "a".into(),
                cluster: 0,
                avg_log2fc: 1.0,
                pval: 0.01,
                pval_adj: 0.0,
                pct_in: 1.0,
                pct_out: 0.5,
            },
            MarkerResult {
                gene: "b".into(),
                cluster: 0,
                avg_log2fc: 1.0,
                pval: 0.04,
                pval_adj: 0.0,
                pct_in: 1.0,
                pct_out: 0.5,
            },
            MarkerResult {
                gene: "c".into(),
                cluster: 0,
                avg_log2fc: 1.0,
                pval: 0.03,
                pval_adj: 0.0,
                pct_in: 1.0,
                pct_out: 0.5,
            },
        ];

        benjamini_hochberg(&mut results);

        // Sorted p-values: 0.01, 0.03, 0.04
        // BH adjusted: 0.01*3/1=0.03, 0.03*3/2=0.045, 0.04*3/3=0.04
        // With cummin (from largest): 0.04, 0.04 (min of 0.045, 0.04), 0.03
        // So: a→0.03, b→0.04, c→0.04
        assert!((results[0].pval_adj - 0.03).abs() < 1e-10, "a adj: {}", results[0].pval_adj);
        assert!((results[1].pval_adj - 0.04).abs() < 1e-10, "b adj: {}", results[1].pval_adj);
        assert!((results[2].pval_adj - 0.04).abs() < 1e-10, "c adj: {}", results[2].pval_adj);
    }

    #[test]
    fn test_wilcoxon_identical_groups() {
        let g1 = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let g2 = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let p = wilcoxon_rank_sum(&g1, &g2);
        // Identical distributions → p should be close to 1.0
        assert!(p > 0.5, "identical groups should have high p-value, got {p}");
    }

    #[test]
    fn test_wilcoxon_separated_groups() {
        let g1 = vec![10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0];
        let g2 = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let p = wilcoxon_rank_sum(&g1, &g2);
        // Well-separated → p should be very small
        assert!(p < 0.01, "separated groups should have low p-value, got {p}");
    }

    #[test]
    fn test_norm_cdf() {
        // Standard normal CDF at 0 should be 0.5
        assert!((norm_cdf(0.0) - 0.5).abs() < 1e-6);
        // CDF at -inf → 0, +inf → 1
        assert!(norm_cdf(-10.0) < 1e-10);
        assert!((norm_cdf(10.0) - 1.0).abs() < 1e-10);
        // CDF at 1.96 ≈ 0.975
        assert!((norm_cdf(1.96) - 0.975).abs() < 0.001);
    }
}
