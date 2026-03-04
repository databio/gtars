use anyhow::{Result, bail};
use sprs::TriMat;

use crate::types::{FeatureMatrix, RnaQcMetrics};

/// Compute RNA QC metrics for each cell.
///
/// Metrics include: number of detected features, total counts, and
/// percentage of counts from mitochondrial genes (detected by "MT-" prefix).
pub fn compute_rna_qc(matrix: &FeatureMatrix) -> Vec<RnaQcMetrics> {
    let mt_mask: Vec<bool> = matrix
        .feature_names
        .iter()
        .map(|name| name.to_uppercase().starts_with("MT-"))
        .collect();

    let n_cells = matrix.n_cells();
    let mut metrics = Vec::with_capacity(n_cells);

    for col_idx in 0..n_cells {
        let col = matrix.matrix.outer_view(col_idx);
        let (n_features, n_counts, mt_counts) = match col {
            Some(col_vec) => {
                let mut nf = 0u32;
                let mut nc = 0.0f64;
                let mut mt = 0.0f64;
                for (row_idx, &val) in col_vec.iter() {
                    if val > 0.0 {
                        nf += 1;
                        nc += val;
                        if mt_mask[row_idx] {
                            mt += val;
                        }
                    }
                }
                (nf, nc, mt)
            }
            None => (0, 0.0, 0.0),
        };

        let pct_mt = if n_counts > 0.0 {
            mt_counts / n_counts * 100.0
        } else {
            0.0
        };

        metrics.push(RnaQcMetrics {
            cell_id: matrix.cell_ids[col_idx].clone(),
            n_features,
            n_counts,
            pct_mt,
        });
    }

    metrics
}

/// Filter cells by minimum features detected and maximum mitochondrial percentage.
pub fn filter_rna_cells(
    matrix: &FeatureMatrix,
    min_features: u32,
    max_pct_mt: f64,
) -> Result<FeatureMatrix> {
    let qc = compute_rna_qc(&matrix);

    let keep: Vec<bool> = qc
        .iter()
        .map(|m| m.n_features >= min_features && m.pct_mt <= max_pct_mt)
        .collect();

    let n_keep = keep.iter().filter(|&&k| k).count();
    if n_keep == 0 {
        bail!("all cells filtered out (min_features={}, max_pct_mt={})", min_features, max_pct_mt);
    }

    subset_columns(&matrix, &keep)
}

/// Filter genes that appear in fewer than `min_cells` cells.
pub fn filter_genes(matrix: &FeatureMatrix, min_cells: u32) -> Result<FeatureMatrix> {
    let n_features = matrix.n_features();
    let mut cells_per_gene = vec![0u32; n_features];

    // Count cells with nonzero expression per gene
    for col_idx in 0..matrix.n_cells() {
        if let Some(col) = matrix.matrix.outer_view(col_idx) {
            for (row_idx, &val) in col.iter() {
                if val > 0.0 {
                    cells_per_gene[row_idx] += 1;
                }
            }
        }
    }

    let keep: Vec<bool> = cells_per_gene.iter().map(|&c| c >= min_cells).collect();
    let n_keep = keep.iter().filter(|&&k| k).count();
    if n_keep == 0 {
        bail!("all genes filtered out (min_cells={})", min_cells);
    }

    subset_rows(&matrix, &keep)
}

/// Subset columns (cells) of a FeatureMatrix. Reconstructs via TriMat.
fn subset_columns(fm: &FeatureMatrix, keep: &[bool]) -> Result<FeatureMatrix> {
    let col_map: Vec<usize> = keep
        .iter()
        .enumerate()
        .filter(|&(_, &k)| k)
        .map(|(i, _)| i)
        .collect();

    let new_n_cells = col_map.len();
    let n_features = fm.n_features();
    let mut tri = TriMat::new((n_features, new_n_cells));

    for (new_col, &old_col) in col_map.iter().enumerate() {
        if let Some(col) = fm.matrix.outer_view(old_col) {
            for (row, &val) in col.iter() {
                tri.add_triplet(row, new_col, val);
            }
        }
    }

    let new_cell_ids: Vec<String> = col_map.iter().map(|&i| fm.cell_ids[i].clone()).collect();

    Ok(FeatureMatrix {
        matrix: tri.to_csc(),
        feature_names: fm.feature_names.clone(),
        feature_ids: fm.feature_ids.clone(),
        cell_ids: new_cell_ids,
        feature_type: fm.feature_type.clone(),
    })
}

/// Subset rows (features) of a FeatureMatrix. Reconstructs via TriMat.
fn subset_rows(fm: &FeatureMatrix, keep: &[bool]) -> Result<FeatureMatrix> {
    let row_map: Vec<usize> = keep
        .iter()
        .enumerate()
        .filter(|&(_, &k)| k)
        .map(|(i, _)| i)
        .collect();

    // Build old→new row index mapping
    let mut old_to_new = vec![usize::MAX; fm.n_features()];
    for (new_row, &old_row) in row_map.iter().enumerate() {
        old_to_new[old_row] = new_row;
    }

    let new_n_features = row_map.len();
    let n_cells = fm.n_cells();
    let mut tri = TriMat::new((new_n_features, n_cells));

    for col_idx in 0..n_cells {
        if let Some(col) = fm.matrix.outer_view(col_idx) {
            for (old_row, &val) in col.iter() {
                let new_row = old_to_new[old_row];
                if new_row != usize::MAX {
                    tri.add_triplet(new_row, col_idx, val);
                }
            }
        }
    }

    let new_feature_names: Vec<String> =
        row_map.iter().map(|&i| fm.feature_names[i].clone()).collect();
    let new_feature_ids: Vec<String> =
        row_map.iter().map(|&i| fm.feature_ids[i].clone()).collect();

    Ok(FeatureMatrix {
        matrix: tri.to_csc(),
        feature_names: new_feature_names,
        feature_ids: new_feature_ids,
        cell_ids: fm.cell_ids.clone(),
        feature_type: fm.feature_type.clone(),
    })
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
    fn test_compute_qc() {
        let fm = read_10x(&test_data_dir()).unwrap();
        let qc = compute_rna_qc(&fm);

        assert_eq!(qc.len(), 3);

        // Cell 0 (AAAA-1): Gene1=5, Gene2=3, MT-CO1=2 → n_features=3, n_counts=10, pct_mt=20%
        assert_eq!(qc[0].n_features, 3);
        assert!((qc[0].n_counts - 10.0).abs() < 1e-10);
        assert!((qc[0].pct_mt - 20.0).abs() < 1e-10);

        // Cell 1 (BBBB-1): Gene4=7, Gene5=1 → n_features=2, n_counts=8, pct_mt=0%
        assert_eq!(qc[1].n_features, 2);
        assert!((qc[1].pct_mt).abs() < 1e-10);
    }

    #[test]
    fn test_filter_cells() {
        let fm = read_10x(&test_data_dir()).unwrap();
        // min_features=3 should keep AAAA-1 (3 features) and CCCC-1 (3 features), drop BBBB-1 (2 features)
        let filtered = filter_rna_cells(&fm, 3, 100.0).unwrap();
        assert_eq!(filtered.n_cells(), 2);
        assert_eq!(filtered.cell_ids, vec!["AAAA-1", "CCCC-1"]);
    }

    #[test]
    fn test_filter_genes() {
        let fm = read_10x(&test_data_dir()).unwrap();
        // min_cells=2: Gene1 (cells 0,2), MT-CO1 (cells 0,2), Gene5 (cells 1,2) → keep 3 genes
        // Gene2 (cell 0 only), Gene4 (cell 1 only) → drop
        let filtered = filter_genes(&fm, 2).unwrap();
        assert_eq!(filtered.n_features(), 3);
        assert_eq!(filtered.feature_names, vec!["Gene1", "MT-CO1", "Gene5"]);
    }
}
