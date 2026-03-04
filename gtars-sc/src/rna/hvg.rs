use anyhow::{Result, bail};

use crate::types::FeatureMatrix;

/// Find highly variable features using Seurat-style VST.
///
/// Matches Seurat's `FindVariableFeatures(selection.method = "vst")`:
/// 1. Compute per-gene mean and variance
/// 2. Fit smoothed log10(variance) ~ log10(mean) via binned averaging (approximates loess)
/// 3. Standardize each value: clip((x - mean) / sqrt(expected_var), clip_max)
/// 4. Compute variance of standardized values per gene
/// 5. Select top n genes by standardized variance
pub fn find_variable_features(
    matrix: &FeatureMatrix,
    n_features: usize,
) -> Result<Vec<String>> {
    let n_genes = matrix.n_features();
    let n_cells = matrix.n_cells();

    if n_features >= n_genes {
        return Ok(matrix.feature_names.clone());
    }
    if n_genes == 0 || n_cells == 0 {
        bail!("empty matrix");
    }

    // --- Pass 1: per-gene mean and variance ---
    let mut sums = vec![0.0f64; n_genes];
    let mut sq_sums = vec![0.0f64; n_genes];
    let mut nnz_per_gene = vec![0u32; n_genes];

    for col_idx in 0..n_cells {
        if let Some(col) = matrix.matrix.outer_view(col_idx) {
            for (row, &val) in col.iter() {
                sums[row] += val;
                sq_sums[row] += val * val;
                nnz_per_gene[row] += 1;
            }
        }
    }

    let n_cells_f = n_cells as f64;
    let mut means = vec![0.0f64; n_genes];
    let mut variances = vec![0.0f64; n_genes];
    for i in 0..n_genes {
        means[i] = sums[i] / n_cells_f;
        // Sample variance: (sum_sq - n*mean^2) / (n-1)
        variances[i] = (sq_sums[i] - n_cells_f * means[i] * means[i]) / (n_cells_f - 1.0);
        if variances[i] < 0.0 {
            variances[i] = 0.0;
        }
    }

    // --- Fit smoothed log10(variance) ~ log10(mean) via local regression ---
    // Matches Seurat: loess(log10(variance) ~ log10(mean), span = 0.3)
    let not_const: Vec<bool> = variances.iter().map(|&v| v > 0.0).collect();
    let log10_means: Vec<f64> = means.iter().map(|&m| (m.max(1e-300)).log10()).collect();
    let log10_vars: Vec<f64> = variances.iter().map(|&v| (v.max(1e-300)).log10()).collect();

    // Collect non-constant genes for fitting
    let mut fit_x: Vec<f64> = Vec::new();
    let mut fit_y: Vec<f64> = Vec::new();
    let mut fit_idx: Vec<usize> = Vec::new();
    for i in 0..n_genes {
        if not_const[i] {
            fit_x.push(log10_means[i]);
            fit_y.push(log10_vars[i]);
            fit_idx.push(i);
        }
    }

    // Evaluate loess at grid points, then interpolate per-gene (O(n_grid * n_fit) instead of O(n^2))
    let n_grid = 200;
    let (x_min, x_max) = fit_x.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
        (lo.min(v), hi.max(v))
    });
    let grid_step = if (x_max - x_min).abs() < 1e-15 { 1.0 } else { (x_max - x_min) / (n_grid - 1) as f64 };
    let grid_x: Vec<f64> = (0..n_grid).map(|i| x_min + i as f64 * grid_step).collect();
    let grid_y = loess_fit_at(&fit_x, &fit_y, &grid_x, 0.3);

    let mut expected_var = vec![1e-10f64; n_genes];
    for &i in &fit_idx {
        let lm = log10_means[i];
        let expected_lv = interpolate(&grid_x, &grid_y, lm);
        expected_var[i] = 10.0_f64.powf(expected_lv);
    }

    // --- Pass 2: variance of standardized values (Seurat convention) ---
    // Seurat clips standardized values only from ABOVE: min(x, vmax), not clamped below.
    let clip_max = (n_cells_f).sqrt();
    let mut std_sq_sum = vec![0.0f64; n_genes];

    // Precompute standardized value for zero entries (only clip from above)
    let mut zero_std = vec![0.0f64; n_genes];
    for i in 0..n_genes {
        let ev_sqrt = expected_var[i].sqrt();
        if ev_sqrt > 1e-15 {
            zero_std[i] = ((-means[i]) / ev_sqrt).min(clip_max);
        }
    }

    // Accumulate squared standardized values from nonzero entries
    for col_idx in 0..n_cells {
        if let Some(col) = matrix.matrix.outer_view(col_idx) {
            for (row, &val) in col.iter() {
                let ev_sqrt = expected_var[row].sqrt();
                let s = if ev_sqrt > 1e-15 {
                    ((val - means[row]) / ev_sqrt).min(clip_max)
                } else {
                    0.0
                };
                std_sq_sum[row] += s * s;
            }
        }
    }

    // Add zero-entry contributions
    for i in 0..n_genes {
        let n_zeros = (n_cells as u32 - nnz_per_gene[i]) as f64;
        std_sq_sum[i] += n_zeros * zero_std[i] * zero_std[i];
    }

    // Seurat: variance = sum(std^2) / (n_cells - 1)
    let mut scores: Vec<(usize, f64)> = (0..n_genes)
        .map(|i| {
            let var_std = std_sq_sum[i] / (n_cells_f - 1.0);
            (i, var_std)
        })
        .collect();

    // Sort descending by standardized variance
    scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let hvgs: Vec<String> = scores
        .iter()
        .take(n_features)
        .map(|&(i, _)| matrix.feature_names[i].clone())
        .collect();

    Ok(hvgs)
}

/// Evaluate loess (local weighted linear regression, degree=1) at query points.
///
/// `data_x` and `data_y` are the training data; `query_x` are the points to predict.
/// Uses tricube kernel with the given span (fraction of data points in each window).
fn loess_fit_at(data_x: &[f64], data_y: &[f64], query_x: &[f64], span: f64) -> Vec<f64> {
    let n = data_x.len();
    if n == 0 {
        return vec![0.0; query_x.len()];
    }

    let window = ((span * n as f64).ceil() as usize).max(2).min(n);

    // Sort data by x
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| data_x[a].partial_cmp(&data_x[b]).unwrap());
    let sorted_x: Vec<f64> = order.iter().map(|&i| data_x[i]).collect();
    let sorted_y: Vec<f64> = order.iter().map(|&i| data_y[i]).collect();

    let mut fitted = Vec::with_capacity(query_x.len());

    for &xi in query_x {
        // Find window of nearest neighbors in sorted data
        let pos = sorted_x.partition_point(|&v| v < xi);
        let mut lo = pos.saturating_sub(window);
        let mut hi = (pos + window).min(n);

        while hi - lo > window {
            let d_lo = (sorted_x[lo] - xi).abs();
            let d_hi = (sorted_x[hi - 1] - xi).abs();
            if d_lo > d_hi {
                lo += 1;
            } else {
                hi -= 1;
            }
        }

        let max_dist = sorted_x[lo..hi]
            .iter()
            .map(|&v| (v - xi).abs())
            .fold(1e-15f64, f64::max);

        // Tricube-weighted linear regression
        let mut sw = 0.0f64;
        let mut swx = 0.0f64;
        let mut swy = 0.0f64;
        let mut swxx = 0.0f64;
        let mut swxy = 0.0f64;

        for j in lo..hi {
            let u = ((sorted_x[j] - xi).abs() / max_dist).min(1.0);
            let t = 1.0 - u * u * u;
            let w = t * t * t;
            sw += w;
            swx += w * sorted_x[j];
            swy += w * sorted_y[j];
            swxx += w * sorted_x[j] * sorted_x[j];
            swxy += w * sorted_x[j] * sorted_y[j];
        }

        let val = if sw < 1e-15 {
            0.0
        } else {
            let denom = sw * swxx - swx * swx;
            if denom.abs() < 1e-15 {
                swy / sw
            } else {
                let b = (sw * swxy - swx * swy) / denom;
                let a = (swy - b * swx) / sw;
                a + b * xi
            }
        };
        fitted.push(val);
    }

    fitted
}

/// Linear interpolation from sorted (x, y) control points.
fn interpolate(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    if xs.len() == 1 || x <= xs[0] {
        return ys[0];
    }
    if x >= xs[xs.len() - 1] {
        return ys[ys.len() - 1];
    }
    let pos = xs.partition_point(|&v| v < x);
    if pos == 0 {
        return ys[0];
    }
    let i = pos - 1;
    let frac = (x - xs[i]) / (xs[i + 1] - xs[i]);
    ys[i] + frac * (ys[i + 1] - ys[i])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::read_10x;
    use std::path::PathBuf;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/small_10x")
    }

    #[test]
    fn test_hvg_returns_requested_count() {
        let fm = read_10x(&test_data_dir()).unwrap();
        let hvgs = find_variable_features(&fm, 3).unwrap();
        assert_eq!(hvgs.len(), 3);
    }

    #[test]
    fn test_hvg_all_features_if_n_exceeds() {
        let fm = read_10x(&test_data_dir()).unwrap();
        let hvgs = find_variable_features(&fm, 100).unwrap();
        assert_eq!(hvgs.len(), 5); // all features returned
    }

    #[test]
    fn test_hvg_results_are_gene_names() {
        let fm = read_10x(&test_data_dir()).unwrap();
        let hvgs = find_variable_features(&fm, 2).unwrap();
        for name in &hvgs {
            assert!(fm.feature_names.contains(name), "unknown gene: {name}");
        }
    }
}
