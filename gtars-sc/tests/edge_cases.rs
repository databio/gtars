//! Edge-case tests for gtars-sc robustness.
//!
//! Tests degenerate inputs: single-cell matrices, constant values,
//! empty results, and boundary conditions that previously could panic
//! or produce incorrect results.

use sprs::{CsMat, TriMat};

use gtars_sc::types::{FeatureMatrix, FeatureType, MarkerConfig};
use gtars_sc::rna::hvg::find_variable_features;
use gtars_sc::rna::normalize::log_normalize;
use gtars_sc::rna::qc::{compute_rna_qc, filter_genes, filter_rna_cells};
use gtars_sc::rna::scale::scale_data;
use gtars_sc::markers::wilcoxon::find_all_markers;

/// Build a FeatureMatrix from dense data (genes × cells).
fn make_fm(data: &[Vec<f64>], names: Vec<&str>) -> FeatureMatrix {
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
        feature_names: names.iter().map(|s| s.to_string()).collect(),
        feature_ids: names.iter().map(|s| s.to_string()).collect(),
        cell_ids: (0..n_cells).map(|i| format!("cell_{i}")).collect(),
        feature_type: FeatureType::Gene,
    }
}

// --- Issue 2: scale_data with single cell ---

#[test]
fn scale_data_single_cell_errors() {
    let fm = make_fm(
        &[vec![1.0], vec![2.0], vec![3.0]],
        vec!["A", "B", "C"],
    );
    let features = vec!["A".to_string(), "B".to_string()];
    let result = scale_data(&fm, &features, None, None);
    assert!(result.is_err(), "scale_data should fail with 1 cell");
}

// --- Issue 3: find_variable_features with single cell ---

#[test]
fn hvg_single_cell_errors() {
    let fm = make_fm(
        &[vec![1.0], vec![2.0], vec![3.0]],
        vec!["A", "B", "C"],
    );
    let result = find_variable_features(&fm, 2);
    assert!(result.is_err(), "HVG should fail with 1 cell");
}

// --- Issue 5: Wilcoxon with single-element group ---

#[test]
fn wilcoxon_single_element_group() {
    // 2 genes, 4 cells, cluster 0 has 1 cell
    let fm = make_fm(
        &[
            vec![5.0, 0.1, 0.1, 0.1],
            vec![0.1, 5.0, 5.0, 5.0],
        ],
        vec!["A", "B"],
    );
    let clusters = vec![0, 1, 1, 1];
    let config = MarkerConfig {
        min_pct: 0.0,
        min_log2fc: 0.0,
        only_positive: false,
    };
    let markers = find_all_markers(&fm, &clusters, &config).unwrap();
    // Should produce valid results, no NaN or panic
    for m in &markers {
        assert!(m.pval.is_finite(), "p-value should be finite, got {}", m.pval);
        assert!(m.pval_adj.is_finite(), "adjusted p-value should be finite");
    }
}

// --- Wilcoxon with all identical values ---

#[test]
fn wilcoxon_identical_values() {
    let fm = make_fm(
        &[vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0]],
        vec!["A"],
    );
    let clusters = vec![0, 0, 0, 1, 1, 1];
    let config = MarkerConfig {
        min_pct: 0.0,
        min_log2fc: 0.0,
        only_positive: false,
    };
    let markers = find_all_markers(&fm, &clusters, &config).unwrap();
    for m in &markers {
        assert!(m.pval.is_finite(), "p-value should be finite for identical values");
    }
}

// --- Two-cell pipeline (minimum valid) ---

#[test]
fn two_cell_pipeline_no_panic() {
    let fm = make_fm(
        &[
            vec![5.0, 0.0],
            vec![0.0, 3.0],
            vec![1.0, 1.0],
        ],
        vec!["A", "B", "C"],
    );
    // HVG with 2 cells should work (n_cells=2, n_cells-1=1)
    let hvgs = find_variable_features(&fm, 2).unwrap();
    assert!(!hvgs.is_empty());

    // Scale with 2 cells should work
    let scaled = scale_data(&fm, &hvgs, None, None).unwrap();
    assert_eq!(scaled.ncols(), 2);
}

// --- filter_genes with all-zero gene ---

#[test]
fn all_zero_gene_filtered() {
    let fm = make_fm(
        &[
            vec![1.0, 2.0, 3.0],
            vec![0.0, 0.0, 0.0], // all-zero gene
            vec![4.0, 5.0, 6.0],
        ],
        vec!["A", "Zero", "C"],
    );
    let filtered = filter_genes(&fm, 1).unwrap();
    assert_eq!(filtered.n_features(), 2);
    assert!(!filtered.feature_names.contains(&"Zero".to_string()));
}

// --- all cells fail QC ---

#[test]
fn all_cells_fail_qc_errors() {
    let fm = make_fm(
        &[
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ],
        vec!["A", "B"],
    );
    // min_features=100 should filter out all cells (each has 2 features)
    let result = filter_rna_cells(&fm, 100, 100.0);
    assert!(result.is_err(), "should error when all cells filtered");
}

// --- constant gene not selected as HVG ---

#[test]
fn constant_gene_not_hvg() {
    // Gene A is constant across cells, B and C have variance
    let fm = make_fm(
        &[
            vec![5.0, 5.0, 5.0, 5.0], // constant
            vec![1.0, 10.0, 1.0, 10.0], // variable
            vec![0.0, 8.0, 0.0, 8.0],  // variable
        ],
        vec!["Constant", "Variable1", "Variable2"],
    );
    let hvgs = find_variable_features(&fm, 1).unwrap();
    assert!(!hvgs.contains(&"Constant".to_string()),
        "constant gene should not be top HVG");
}

// --- normalize with zero-count cell ---

#[test]
fn normalize_zero_count_cell() {
    let mut fm = make_fm(
        &[
            vec![5.0, 0.0, 3.0],
            vec![2.0, 0.0, 1.0],
        ],
        vec!["A", "B"],
    );
    // Cell 1 has all zeros — should not panic
    log_normalize(&mut fm, 10_000.0).unwrap();

    // Cell 1 should still have no entries
    let qc = compute_rna_qc(&fm);
    assert_eq!(qc[1].n_counts, 0.0);
}

// --- KNN with k=1 ---

#[test]
fn knn_k_equals_1() {
    use gtars_sc::reduce::neighbors::build_knn;
    use ndarray::Array2;

    let embedding = Array2::from_shape_vec(
        (4, 2),
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    ).unwrap();
    let knn = build_knn(&embedding, 1).unwrap();
    assert_eq!(knn.k, 1);
    assert_eq!(knn.indices.len(), 4);
    for neighbors in &knn.indices {
        assert_eq!(neighbors.len(), 1);
    }
}
