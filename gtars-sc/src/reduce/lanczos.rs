//! Augmented Implicitly Restarted Lanczos Bidiagonalization (IRLBA)
//!
//! Implements the algorithm described in:
//!   Baglama, J. and Reichel, L. (2005).
//!   "Augmented implicitly restarted Lanczos bidiagonalization methods."
//!   SIAM Journal on Scientific Computing, 27(1), 19–42.
//!   doi:10.1137/04060593X

use anyhow::{Result, bail};
use faer::linalg::matmul::matmul;
use faer::linalg::solvers::Svd;
use faer::{Accum, Mat, Par};
use ndarray::{Array1, Array2};
use rand::Rng;

use crate::types::SvdResult;

/// Return the parallelism strategy for large matrix-vector products.
/// Uses rayon when the `parallel` feature is enabled, sequential otherwise.
fn matmul_par() -> Par {
    #[cfg(feature = "parallel")]
    {
        Par::rayon(0)
    }
    #[cfg(not(feature = "parallel"))]
    {
        Par::Seq
    }
}

/// Compute top-k truncated SVD using augmented implicitly restarted
/// Lanczos bidiagonalization (Baglama & Reichel 2005).
///
/// Input: (n_features × n_cells) matrix, number of components k.
/// Returns: U (n_features × k), S (k), Vt (k × n_cells), variance_explained.
///
/// Falls back to exact SVD when the matrix is too small for the iterative approach.
pub fn lanczos_svd(matrix: &Array2<f64>, n_components: usize) -> Result<SvdResult> {
    let (m, n) = matrix.dim();
    let k = n_components.min(m).min(n);

    if k == 0 {
        bail!("cannot compute SVD with 0 components");
    }

    let min_dim = m.min(n);

    // Work = subspace dimension. Must satisfy k < work < min_dim.
    // Match R's irlba default: work = k + 7. Using the same work dimension
    // as R keeps the restart path similar, improving correlation of later
    // singular vectors on near-degenerate spectra.
    let work = k + 7;

    if work >= min_dim || k >= min_dim {
        return crate::reduce::svd::truncated_svd(matrix, n_components);
    }

    let tol: f64 = 1e-5;
    let maxit: usize = 1000;
    let p = work;

    // Convert ndarray → faer
    let a = Mat::<f64>::from_fn(m, n, |i, j| matrix[[i, j]]);

    // Working matrices
    let mut w_basis = Mat::<f64>::zeros(m, p);
    let mut v_basis = Mat::<f64>::zeros(n, p + 1);
    let mut b_mat = Mat::<f64>::zeros(p, p); // bidiagonal (or augmented after restart)

    // Working column vectors as (rows × 1) Mat for matmul compatibility
    let mut w_col = Mat::<f64>::zeros(m, 1);
    let mut f_vec = Mat::<f64>::zeros(n, 1);

    // Random starting vector
    let mut rng = rand::rng();
    for i in 0..n {
        *v_basis.get_mut(i, 0) = rng.random::<f64>() - 0.5;
    }
    let norm = v_basis.as_ref().subcols(0, 1).norm_l2();
    if norm < f64::EPSILON {
        bail!("random start vector has near-zero norm");
    }
    for i in 0..n {
        *v_basis.get_mut(i, 0) /= norm;
    }

    let mut f_norm: f64 = 0.0;
    let mut k_restart: usize = 0;

    for _iter in 0..maxit {
        // --- Lanczos bidiagonalization from j = k_restart to p-1 ---
        for j in k_restart..p {
            // w = A * v_j
            {
                let v_j = v_basis.as_ref().subcols(j, 1);
                matmul(w_col.as_mut(), Accum::Replace, a.as_ref(), v_j, 1.0, matmul_par());
            }

            // Reorthogonalize w against W[:, 0..j]
            if j > 0 {
                reorthogonalize(&mut w_col, &w_basis, j);
            }

            // s = ||w||
            let s = w_col.as_ref().norm_l2();
            // Relative tolerance: scale with largest diagonal seen so far
            let max_diag = (0..=j).map(|idx| (*b_mat.get(idx, idx)).abs()).fold(0.0f64, f64::max);
            let breakdown_tol = 1e-12 * max_diag.max(1e-15);
            if s < breakdown_tol {
                // Lucky breakdown: w is in span(W). Set alpha=0, W[:,j]=0.
                *b_mat.get_mut(j, j) = 0.0;
                for i in 0..m {
                    *w_basis.get_mut(i, j) = 0.0;
                }
                // Generate random f for continuation
                for i in 0..n {
                    *f_vec.get_mut(i, 0) = rng.random::<f64>() - 0.5;
                }
                reorthogonalize(&mut f_vec, &v_basis, j + 1);
                f_norm = f_vec.as_ref().norm_l2();
                if j + 1 < p && f_norm > breakdown_tol {
                    *b_mat.get_mut(j, j + 1) = 0.0; // no coupling
                    for i in 0..n {
                        *v_basis.get_mut(i, j + 1) = *f_vec.get(i, 0) / f_norm;
                    }
                }
                continue;
            }
            *b_mat.get_mut(j, j) = s;

            // W[:, j] = w / s
            for i in 0..m {
                *w_basis.get_mut(i, j) = *w_col.get(i, 0) / s;
            }

            // f = A^T * W[:, j] - s * V[:, j]
            {
                let wj = w_basis.as_ref().subcols(j, 1);
                matmul(f_vec.as_mut(), Accum::Replace, a.transpose(), wj, 1.0, matmul_par());
            }
            for i in 0..n {
                *f_vec.get_mut(i, 0) -= s * *v_basis.get(i, j);
            }

            // Reorthogonalize f against V[:, 0..j+1]
            reorthogonalize(&mut f_vec, &v_basis, j + 1);

            // f_norm = ||f||
            f_norm = f_vec.as_ref().norm_l2();

            if j + 1 < p {
                *b_mat.get_mut(j, j + 1) = f_norm;

                // Use same relative tolerance for residual breakdown
                let res_tol = 1e-12 * (0..=j).map(|idx| (*b_mat.get(idx, idx)).abs()).fold(0.0f64, f64::max).max(1e-15);
                if f_norm < res_tol {
                    // Lucky breakdown: residual is zero. Generate random restart vector.
                    for i in 0..n {
                        *f_vec.get_mut(i, 0) = rng.random::<f64>() - 0.5;
                    }
                    reorthogonalize(&mut f_vec, &v_basis, j + 1);
                    f_norm = f_vec.as_ref().norm_l2();
                    if f_norm > res_tol {
                        for i in 0..n {
                            *v_basis.get_mut(i, j + 1) = *f_vec.get(i, 0) / f_norm;
                        }
                    }
                    // Reset f_norm to the actual (near-zero) residual for convergence check
                    f_norm = 0.0;
                } else {
                    // V[:, j+1] = f / f_norm
                    for i in 0..n {
                        *v_basis.get_mut(i, j + 1) = *f_vec.get(i, 0) / f_norm;
                    }
                }
            }
        }

        // --- SVD of B (p × p) ---
        let b_svd = Svd::new_thin(b_mat.as_ref())
            .map_err(|e| anyhow::anyhow!("SVD of B failed: {e:?}"))?;

        let s_b = b_svd.S().column_vector();
        let u_b = b_svd.U();
        let v_b = b_svd.V();

        // --- Convergence check ---
        let smax = *s_b.get(0);
        if smax == 0.0 {
            bail!("all singular values are zero");
        }
        let threshold = tol * smax;

        let mut n_converged = 0usize;
        for i in 0..k {
            let residual = f_norm * (*u_b.get(p - 1, i)).abs();
            if residual < threshold {
                n_converged += 1;
            }
        }

        if n_converged >= k {
            return extract_result(&w_basis, &v_basis, &b_svd, m, n, k, p);
        }

        // --- Augmented implicit restart ---
        // Match R's irlba: always restart with exactly k Ritz vectors.
        // This leaves (p - k) fresh Lanczos steps per restart (= 7 with
        // work = k + 7), matching R's convergence path.
        let k_r = k;

        // Rotate V basis: V_new[:, 0..k_r] = V[:, 0..p] * V_B[:, 0..k_r]
        let mut v_new = Mat::<f64>::zeros(n, k_r);
        {
            let v_sub = v_basis.as_ref().subcols(0, p);
            let vb_sub = v_b.subcols(0, k_r);
            matmul(v_new.as_mut(), Accum::Replace, v_sub, vb_sub, 1.0, Par::Seq);
        }
        for j in 0..k_r {
            for i in 0..n {
                *v_basis.get_mut(i, j) = *v_new.get(i, j);
            }
        }

        // V[:, k_r] = f / f_norm (augmentation vector)
        if f_norm > 0.0 {
            for i in 0..n {
                *v_basis.get_mut(i, k_r) = *f_vec.get(i, 0) / f_norm;
            }
        }

        // Rotate W basis: W_new[:, 0..k_r] = W[:, 0..p] * U_B[:, 0..k_r]
        let mut w_new = Mat::<f64>::zeros(m, k_r);
        {
            let w_sub = w_basis.as_ref().subcols(0, p);
            let ub_sub = u_b.subcols(0, k_r);
            matmul(w_new.as_mut(), Accum::Replace, w_sub, ub_sub, 1.0, Par::Seq);
        }
        for j in 0..k_r {
            for i in 0..m {
                *w_basis.get_mut(i, j) = *w_new.get(i, j);
            }
        }

        // Reconstruct B for next iteration:
        // B = diag(S_B[0..k_r]) with residuals in column k_r
        b_mat = Mat::<f64>::zeros(p, p);
        for l in 0..k_r {
            *b_mat.get_mut(l, l) = *s_b.get(l);
            *b_mat.get_mut(l, k_r) = f_norm * *u_b.get(p - 1, l);
        }

        k_restart = k_r;
    }

    bail!("Lanczos SVD did not converge in {maxit} iterations");
}

/// Reorthogonalize column vector `v` (stored as m×1 Mat) against the first
/// `ncols` columns of `basis`. Uses double Classical Gram-Schmidt (two passes)
/// to maintain numerical orthogonality — a single pass loses orthogonality
/// in finite precision, causing ghost eigenvalues and poor convergence of
/// later singular vectors.
fn reorthogonalize(v: &mut Mat<f64>, basis: &Mat<f64>, ncols: usize) {
    if ncols == 0 {
        return;
    }
    let q = basis.as_ref().subcols(0, ncols);
    let mut h = Mat::<f64>::zeros(ncols, 1);

    // First pass: h = Q^T * v; v = v - Q * h
    matmul(h.as_mut(), Accum::Replace, q.transpose(), v.as_ref(), 1.0, Par::Seq);
    matmul(v.as_mut(), Accum::Add, q, h.as_ref(), -1.0, Par::Seq);

    // Second pass: removes residual components lost to rounding
    matmul(h.as_mut(), Accum::Replace, q.transpose(), v.as_ref(), 1.0, Par::Seq);
    matmul(v.as_mut(), Accum::Add, q, h.as_ref(), -1.0, Par::Seq);
}

/// Extract final SVD result from converged Lanczos bases.
fn extract_result(
    w_basis: &Mat<f64>,
    v_basis: &Mat<f64>,
    b_svd: &Svd<f64>,
    m: usize,
    n: usize,
    k: usize,
    p: usize,
) -> Result<SvdResult> {
    let u_b = b_svd.U();
    let v_b = b_svd.V();
    let s_col = b_svd.S().column_vector();

    let s_vec: Vec<f64> = (0..k).map(|i| *s_col.get(i)).collect();

    let retained_var: f64 = s_vec.iter().map(|&sv| sv * sv).sum();
    let variance_explained: Vec<f64> = s_vec
        .iter()
        .map(|&sv| if retained_var > 0.0 { sv * sv / retained_var } else { 0.0 })
        .collect();

    // U = W * U_B[:, 0..k]  → (m × k)
    let mut u_result = Mat::<f64>::zeros(m, k);
    {
        let w_sub = w_basis.as_ref().subcols(0, p);
        let ub_sub = u_b.subcols(0, k);
        matmul(u_result.as_mut(), Accum::Replace, w_sub, ub_sub, 1.0, Par::Seq);
    }

    // V_out = V * V_B[:, 0..k]  → (n × k)
    let mut v_result = Mat::<f64>::zeros(n, k);
    {
        let v_sub = v_basis.as_ref().subcols(0, p);
        let vb_sub = v_b.subcols(0, k);
        matmul(v_result.as_mut(), Accum::Replace, v_sub, vb_sub, 1.0, Par::Seq);
    }

    // Convert to ndarray
    let mut u = Array2::<f64>::zeros((m, k));
    for i in 0..m {
        for j in 0..k {
            u[[i, j]] = *u_result.get(i, j);
        }
    }

    let mut vt = Array2::<f64>::zeros((k, n));
    for i in 0..k {
        for j in 0..n {
            vt[[i, j]] = *v_result.get(j, i);
        }
    }

    Ok(SvdResult {
        u,
        s: Array1::from_vec(s_vec),
        vt,
        variance_explained,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reduce::svd::truncated_svd;
    use ndarray::array;

    #[test]
    fn test_lanczos_shapes() {
        // Small matrix → falls back to exact SVD
        let mat = array![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0],
        ];
        let result = lanczos_svd(&mat, 2).unwrap();
        assert_eq!(result.u.dim(), (4, 2));
        assert_eq!(result.s.len(), 2);
        assert_eq!(result.vt.dim(), (2, 3));
        assert_eq!(result.variance_explained.len(), 2);
    }

    #[test]
    fn test_lanczos_variance_sums_to_one() {
        let mat = array![
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ];
        let result = lanczos_svd(&mat, 3).unwrap();
        let total: f64 = result.variance_explained.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-10,
            "variance explained total = {total}, expected ~1.0"
        );
    }

    #[test]
    fn test_lanczos_reconstruction_small() {
        // Small matrix (fallback path)
        let mat = array![
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
        ];
        let result = lanczos_svd(&mat, 2).unwrap();

        for i in 0..3 {
            for j in 0..2 {
                let mut val = 0.0;
                for kk in 0..2 {
                    val += result.u[[i, kk]] * result.s[kk] * result.vt[[kk, j]];
                }
                assert!(
                    (val - mat[[i, j]]).abs() < 1e-10,
                    "reconstruction error at ({i},{j}): got {val}, expected {}",
                    mat[[i, j]]
                );
            }
        }
    }

    #[test]
    fn test_lanczos_agrees_with_exact_svd() {
        {
            // Large enough matrix to exercise the Lanczos path
            let m = 100;
            let n = 80;
            let k = 10;

            let mut rng = rand::rng();
            let mat = Array2::from_shape_fn((m, n), |_| rng.random::<f64>());

            let exact = truncated_svd(&mat, k).unwrap();
            let lanczos = lanczos_svd(&mat, k).unwrap();

            // Singular values should agree within tolerance
            for i in 0..k {
                let rel_err = (exact.s[i] - lanczos.s[i]).abs() / exact.s[i].max(1e-15);
                assert!(
                    rel_err < 1e-3,
                    "singular value {i}: exact={}, lanczos={}, rel_err={rel_err}",
                    exact.s[i], lanczos.s[i]
                );
            }

            // Singular vectors should agree up to sign
            for j in 0..k {
                let dot: f64 = (0..m).map(|i| exact.u[[i, j]] * lanczos.u[[i, j]]).sum();
                assert!(
                    (dot.abs() - 1.0).abs() < 1e-2,
                    "left singular vector {j} disagreement: |dot| = {}",
                    dot.abs()
                );
            }

            // Variance explained should match
            for i in 0..k {
                let diff = (exact.variance_explained[i] - lanczos.variance_explained[i]).abs();
                assert!(
                    diff < 1e-3,
                    "variance_explained[{i}]: exact={}, lanczos={}",
                    exact.variance_explained[i], lanczos.variance_explained[i]
                );
            }
        }
    }

    #[test]
    fn test_lanczos_reconstruction_large() {
        {
            // True rank-2 matrix: A = u1*v1^T + u2*v2^T
            // Large enough for iterative path (50×40, k=2, work=9 < min(50,40)=40)
            let m = 50;
            let n = 40;
            let k = 2;

            let mut rng = rand::rng();
            let u1: Vec<f64> = (0..m).map(|_| rng.random::<f64>()).collect();
            let u2: Vec<f64> = (0..m).map(|_| rng.random::<f64>()).collect();
            let v1: Vec<f64> = (0..n).map(|_| rng.random::<f64>()).collect();
            let v2: Vec<f64> = (0..n).map(|_| rng.random::<f64>()).collect();

            let mut mat = Array2::<f64>::zeros((m, n));
            for i in 0..m {
                for j in 0..n {
                    mat[[i, j]] = 3.0 * u1[i] * v1[j] + 1.5 * u2[i] * v2[j];
                }
            }

            let result = lanczos_svd(&mat, k).unwrap();

            // Reconstruct: A ≈ U * diag(S) * Vt
            let mut max_err = 0.0f64;
            for i in 0..m {
                for j in 0..n {
                    let mut val = 0.0;
                    for kk in 0..k {
                        val += result.u[[i, kk]] * result.s[kk] * result.vt[[kk, j]];
                    }
                    let err = (val - mat[[i, j]]).abs();
                    max_err = max_err.max(err);
                }
            }

            assert!(
                max_err < 1e-8,
                "reconstruction max error = {max_err}, expected < 1e-8"
            );
        }
    }
}
