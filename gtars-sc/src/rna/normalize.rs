use anyhow::{Result, anyhow};

use crate::types::FeatureMatrix;

/// Log-normalize the count matrix in place.
///
/// For each cell, values become `log1p(count / total_counts * scale_factor)`.
/// Only operates on nonzero entries (preserves sparsity pattern, though values
/// will change).
pub fn log_normalize(matrix: &mut FeatureMatrix, scale_factor: f64) -> Result<()> {
    let n_cells = matrix.n_cells();

    // First pass: compute per-cell totals
    // Clone indptr to avoid borrow conflict with data_mut()
    let indptr: Vec<usize> = matrix.matrix.indptr().as_slice()
        .ok_or_else(|| anyhow!("sparse matrix has invalid indptr storage"))?
        .to_vec();

    let mut cell_totals = vec![0.0f64; n_cells];
    {
        let data = matrix.matrix.data();
        for col in 0..n_cells {
            let start = indptr[col];
            let end = indptr[col + 1];
            for &val in &data[start..end] {
                cell_totals[col] += val;
            }
        }
    }

    // Second pass: normalize in place
    let data = matrix.matrix.data_mut();
    for col in 0..n_cells {
        let total = cell_totals[col];
        if total == 0.0 {
            continue;
        }
        let start = indptr[col];
        let end = indptr[col + 1];
        for val in &mut data[start..end] {
            *val = (*val / total * scale_factor).ln_1p();
        }
    }

    Ok(())
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
    fn test_log_normalize() {
        let mut fm = read_10x(&test_data_dir()).unwrap();
        log_normalize(&mut fm, 10_000.0).unwrap();

        // Cell 0: total = 10. Gene1=5 → log1p(5/10*10000) = log1p(5000)
        let val = fm.matrix.get(0, 0).copied().unwrap_or(0.0);
        let expected = (5.0 / 10.0 * 10_000.0_f64).ln_1p();
        assert!((val - expected).abs() < 1e-10, "got {val}, expected {expected}");

        // Sparsity preserved: same nnz
        assert_eq!(fm.matrix.nnz(), 8);
    }

    #[test]
    fn test_log_normalize_preserves_zeros() {
        let mut fm = read_10x(&test_data_dir()).unwrap();
        log_normalize(&mut fm, 10_000.0).unwrap();
        // Gene2, cell 1 should still be absent (zero)
        assert_eq!(fm.matrix.get(1, 1), None);
    }
}
