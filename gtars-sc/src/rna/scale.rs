use std::collections::HashMap;

use anyhow::{Result, bail};
use ndarray::Array2;

use crate::types::FeatureMatrix;

/// Subset the sparse matrix to the given features, convert to dense,
/// center (mean=0), scale (std=1), and optionally clip values.
///
/// Returns a dense (n_features × n_cells) matrix.
///
/// If `regress` is provided (e.g. percent_mt per cell), simple OLS residuals
/// are computed before centering/scaling.
pub fn scale_data(
    matrix: &FeatureMatrix,
    features: &[String],
    clip: Option<f64>,
    regress: Option<&[f64]>,
) -> Result<Array2<f64>> {
    if features.is_empty() {
        bail!("no features to scale");
    }
    if matrix.n_cells() < 2 {
        bail!("scale_data requires at least 2 cells (got {})", matrix.n_cells());
    }

    // Map feature names to row indices
    let feature_indices: Vec<usize> = features
        .iter()
        .filter_map(|name| matrix.feature_names.iter().position(|n| n == name))
        .collect();

    if feature_indices.len() != features.len() {
        bail!(
            "some features not found in matrix ({} requested, {} found)",
            features.len(),
            feature_indices.len()
        );
    }

    let n_feats = feature_indices.len();
    let n_cells = matrix.n_cells();

    // Build O(1) lookup: old_row -> new_row
    let row_lookup: HashMap<usize, usize> = feature_indices
        .iter()
        .enumerate()
        .map(|(new, &old)| (old, new))
        .collect();

    // Convert sparse subset to dense
    let mut dense = Array2::<f64>::zeros((n_feats, n_cells));
    for col_idx in 0..n_cells {
        if let Some(col) = matrix.matrix.outer_view(col_idx) {
            for (old_row, &val) in col.iter() {
                if let Some(&new_row) = row_lookup.get(&old_row) {
                    dense[[new_row, col_idx]] = val;
                }
            }
        }
    }

    // Optional regression
    if let Some(covariate) = regress {
        if covariate.len() != n_cells {
            bail!("covariate length {} != n_cells {}", covariate.len(), n_cells);
        }
        regress_out(&mut dense, covariate);
    }

    // Center and scale each feature (row)
    for i in 0..n_feats {
        let mut row = dense.row_mut(i);
        let n = row.len() as f64;

        let mean = row.sum() / n;
        row -= mean;

        let var = row.iter().map(|&v| v * v).sum::<f64>() / (n - 1.0);
        let std = var.sqrt();
        if std > 1e-15 {
            row /= std;
        }

        // Clip
        if let Some(clip_val) = clip {
            for val in row.iter_mut() {
                *val = val.clamp(-clip_val, clip_val);
            }
        }
    }

    Ok(dense)
}

/// Simple OLS regression: for each feature, regress out a single covariate.
fn regress_out(data: &mut Array2<f64>, covariate: &[f64]) {
    let n = covariate.len() as f64;
    let cov_mean: f64 = covariate.iter().sum::<f64>() / n;
    let cov_var: f64 = covariate.iter().map(|&x| (x - cov_mean).powi(2)).sum::<f64>();

    if cov_var < 1e-15 {
        return; // constant covariate, nothing to regress
    }

    for i in 0..data.nrows() {
        let row = data.row(i);
        let row_mean = row.sum() / n;

        let cov_sum: f64 = row
            .iter()
            .zip(covariate.iter())
            .map(|(&y, &x)| (y - row_mean) * (x - cov_mean))
            .sum();

        let beta = cov_sum / cov_var;
        let intercept = row_mean - beta * cov_mean;

        let mut row_mut = data.row_mut(i);
        for (j, val) in row_mut.iter_mut().enumerate() {
            *val -= intercept + beta * covariate[j];
        }
    }
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
    fn test_scale_output_shape() {
        let fm = read_10x(&test_data_dir()).unwrap();
        let features = vec!["Gene1".to_string(), "Gene2".to_string()];
        let scaled = scale_data(&fm, &features, Some(10.0), None).unwrap();
        assert_eq!(scaled.shape(), &[2, 3]);
    }

    #[test]
    fn test_scale_centered() {
        let fm = read_10x(&test_data_dir()).unwrap();
        let features = vec!["Gene1".to_string(), "Gene4".to_string(), "Gene5".to_string()];
        let scaled = scale_data(&fm, &features, None, None).unwrap();

        // Each row should have mean ≈ 0
        for i in 0..scaled.nrows() {
            let row_mean: f64 = scaled.row(i).sum() / scaled.ncols() as f64;
            assert!(
                row_mean.abs() < 1e-10,
                "row {i} mean = {row_mean}, expected ~0"
            );
        }
    }

    #[test]
    fn test_scale_clipping() {
        let fm = read_10x(&test_data_dir()).unwrap();
        let features = vec!["Gene1".to_string()];
        let scaled = scale_data(&fm, &features, Some(0.5), None).unwrap();

        for &val in scaled.iter() {
            assert!(val.abs() <= 0.5 + 1e-10, "val {val} exceeds clip");
        }
    }

    #[test]
    fn test_scale_missing_feature_errors() {
        let fm = read_10x(&test_data_dir()).unwrap();
        let features = vec!["Nonexistent".to_string()];
        assert!(scale_data(&fm, &features, None, None).is_err());
    }
}
