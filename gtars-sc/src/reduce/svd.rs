use anyhow::{Result, bail};
use faer::Mat;
use faer::linalg::solvers::Svd;
use ndarray::{Array1, Array2};

use crate::types::SvdResult;

/// Compute truncated SVD of a dense matrix.
///
/// Input: (n_features × n_cells) matrix.
/// Returns: U (n_features × k), S (k), Vt (k × n_cells), variance_explained.
///
/// Uses faer's full thin SVD, then truncates to the top k components.
/// Variance explained is normalized to sum to 1.0 across the retained components
/// (matching Seurat/irlba convention).
pub fn truncated_svd(matrix: &Array2<f64>, n_components: usize) -> Result<SvdResult> {
    let (nrows, ncols) = matrix.dim();
    let k = n_components.min(nrows).min(ncols);

    if k == 0 {
        bail!("cannot compute SVD with 0 components");
    }

    // Convert ndarray → faer Mat<f64>
    let faer_mat = Mat::<f64>::from_fn(nrows, ncols, |i, j| matrix[[i, j]]);

    // Compute thin SVD: A = U * S * V^T
    let svd = Svd::new_thin(faer_mat.as_ref())
        .map_err(|e| anyhow::anyhow!("SVD failed: {e:?}"))?;

    let s_col = svd.S().column_vector();
    let full_k = s_col.nrows();
    let k = k.min(full_k);

    let s_vec: Vec<f64> = (0..k).map(|i| *s_col.get(i)).collect();

    // Variance explained: s_i^2 / sum(s_retained^2)
    let retained_var: f64 = s_vec.iter().map(|&sv| sv * sv).sum();
    let variance_explained: Vec<f64> = s_vec
        .iter()
        .map(|&sv| if retained_var > 0.0 { sv * sv / retained_var } else { 0.0 })
        .collect();

    // Extract U (nrows × k)
    let u_faer = svd.U();
    let mut u = Array2::<f64>::zeros((nrows, k));
    for i in 0..nrows {
        for j in 0..k {
            u[[i, j]] = *u_faer.get(i, j);
        }
    }

    // Extract Vt (k × ncols) from V
    let v_faer = svd.V();
    let mut vt = Array2::<f64>::zeros((k, ncols));
    for i in 0..k {
        for j in 0..ncols {
            vt[[i, j]] = *v_faer.get(j, i);
        }
    }

    let s = Array1::from_vec(s_vec);

    Ok(SvdResult {
        u,
        s,
        vt,
        variance_explained,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_truncated_svd_shapes() {
        let mat = array![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0],
        ];

        let result = truncated_svd(&mat, 2).unwrap();
        assert_eq!(result.u.dim(), (4, 2));
        assert_eq!(result.s.len(), 2);
        assert_eq!(result.vt.dim(), (2, 3));
        assert_eq!(result.variance_explained.len(), 2);
    }

    #[test]
    fn test_svd_variance_explained_sums_to_near_one() {
        let mat = array![
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ];

        let result = truncated_svd(&mat, 3).unwrap();
        let total: f64 = result.variance_explained.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-10,
            "variance explained total = {total}, expected ~1.0"
        );
    }

    #[test]
    fn test_svd_reconstruction() {
        let mat = array![
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
        ];

        let result = truncated_svd(&mat, 2).unwrap();

        // Reconstruct: A ≈ U * diag(S) * Vt
        let mut reconstructed = Array2::<f64>::zeros((3, 2));
        for i in 0..3 {
            for j in 0..2 {
                let mut val = 0.0;
                for k in 0..2 {
                    val += result.u[[i, k]] * result.s[k] * result.vt[[k, j]];
                }
                reconstructed[[i, j]] = val;
            }
        }

        for i in 0..3 {
            for j in 0..2 {
                assert!(
                    (reconstructed[[i, j]] - mat[[i, j]]).abs() < 1e-10,
                    "reconstruction error at ({i},{j}): got {}, expected {}",
                    reconstructed[[i, j]],
                    mat[[i, j]]
                );
            }
        }
    }
}
