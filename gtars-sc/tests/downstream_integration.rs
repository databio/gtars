//! Integration tests for the downstream analysis pipeline (Phase 2).
//!
//! Uses a synthetic dataset with planted cluster structure to verify
//! KNN → SNN → Leiden → markers → silhouette works end-to-end.

use gtars_sc::cluster::leiden::leiden_clustering;
use gtars_sc::cluster::silhouette::{mean_silhouette, silhouette_scores};
use gtars_sc::markers::wilcoxon::find_all_markers;
use gtars_sc::reduce::neighbors::{build_knn, build_snn};
use gtars_sc::rna::pipeline::run_full_rna_pipeline;
use gtars_sc::types::{
    DownstreamConfig, FeatureMatrix, FeatureType, MarkerConfig, RnaPipelineConfig,
};

use ndarray::Array2;
use sprs::TriMat;

/// Build a synthetic FeatureMatrix with planted cluster structure.
///
/// Creates n_genes genes x n_cells cells with 2 clusters:
/// - First half of cells: high expression in first half of genes
/// - Second half of cells: high expression in second half of genes
fn make_planted_matrix(n_genes: usize, n_cells: usize) -> (FeatureMatrix, Vec<u32>) {
    let half_genes = n_genes / 2;
    let half_cells = n_cells / 2;

    let mut tri = TriMat::new((n_genes, n_cells));
    for gene in 0..n_genes {
        for cell in 0..n_cells {
            let is_marker_gene_for_cluster0 = gene < half_genes;
            let is_cluster0_cell = cell < half_cells;

            let val = if is_marker_gene_for_cluster0 == is_cluster0_cell {
                // High expression: marker gene matches cluster
                5.0 + (gene as f64 * 0.1) // slight variation
            } else {
                // Low background
                0.5
            };
            if val > 0.0 {
                tri.add_triplet(gene, cell, val);
            }
        }
    }

    let mat = tri.to_csc();
    let feature_names: Vec<String> = (0..n_genes).map(|i| format!("Gene{i}")).collect();
    let cell_ids: Vec<String> = (0..n_cells).map(|i| format!("Cell{i}")).collect();
    let true_clusters: Vec<u32> = (0..n_cells)
        .map(|i| if i < half_cells { 0 } else { 1 })
        .collect();

    let fm = FeatureMatrix {
        matrix: mat,
        feature_names: feature_names.clone(),
        feature_ids: feature_names,
        cell_ids,
        feature_type: FeatureType::Gene,
    };

    (fm, true_clusters)
}

#[test]
fn test_knn_snn_leiden_on_planted_data() {
    use rand::Rng;
    use rand::SeedableRng;

    // Create PCA-like embedding with planted 2-cluster structure.
    // Use random scatter around cluster centers so all cells in each
    // cluster share high KNN overlap (unlike a gradient which creates
    // sub-structure within clusters).
    let n_cells = 40;
    let n_dims = 5;
    let mut embedding = Array2::<f64>::zeros((n_cells, n_dims));
    let mut rng = rand::rngs::SmallRng::seed_from_u64(42);

    for i in 0..n_cells {
        let center = if i < 20 { 5.0 } else { -5.0 };
        for d in 0..n_dims {
            embedding[[i, d]] = center + rng.random_range(-0.3..0.3);
        }
    }

    // KNN
    let knn = build_knn(&embedding, 10).unwrap();
    assert_eq!(knn.n_cells, n_cells);
    assert_eq!(knn.k, 10);

    // All of cell 0's neighbors should be in its cluster (0-19)
    for &neighbor in &knn.indices[0] {
        assert!(neighbor < 20, "cell 0 should only have same-cluster neighbors");
    }

    // SNN
    let snn = build_snn(&knn, 1.0 / 15.0);
    let total_edges: usize = snn.edges.iter().map(|e| e.len()).sum();
    assert!(total_edges > 0, "SNN graph should have edges");

    // Leiden — should separate the two planted groups
    let result = leiden_clustering(&snn, 0.8, 10).unwrap();
    assert_eq!(result.assignments.len(), n_cells);
    assert!(result.n_clusters >= 2, "should find at least 2 clusters, got {}", result.n_clusters);

    // Verify separation: no cell from group A (0-19) shares a cluster with group B (20-39)
    let group_a_clusters: std::collections::HashSet<u32> =
        result.assignments[..20].iter().copied().collect();
    let group_b_clusters: std::collections::HashSet<u32> =
        result.assignments[20..].iter().copied().collect();
    let overlap: Vec<&u32> = group_a_clusters.intersection(&group_b_clusters).collect();
    assert!(
        overlap.is_empty(),
        "planted groups should not share clusters, but share: {:?}",
        overlap
    );
}

#[test]
fn test_markers_on_planted_data() {
    let (fm, clusters) = make_planted_matrix(20, 20);

    let config = MarkerConfig {
        min_pct: 0.0,
        min_log2fc: 0.0,
        only_positive: true,
    };

    let markers = find_all_markers(&fm, &clusters, &config).unwrap();
    assert!(!markers.is_empty(), "should find markers");

    // Check that Gene0..Gene9 are markers for cluster 0
    let c0_markers: Vec<&str> = markers
        .iter()
        .filter(|m| m.cluster == 0)
        .map(|m| m.gene.as_str())
        .collect();

    // At least some of Gene0..Gene9 should be markers for cluster 0
    let expected_c0: Vec<String> = (0..10).map(|i| format!("Gene{i}")).collect();
    let overlap = c0_markers
        .iter()
        .filter(|g| expected_c0.contains(&g.to_string()))
        .count();
    assert!(
        overlap >= 5,
        "expected at least 5 of Gene0-9 as cluster 0 markers, got {overlap}"
    );

    // All markers should have adjusted p-values
    for m in &markers {
        assert!(m.pval_adj >= 0.0 && m.pval_adj <= 1.0);
        assert!(m.pval_adj >= m.pval || (m.pval_adj - m.pval).abs() < 1e-10);
    }
}

#[test]
fn test_silhouette_on_planted_data() {
    let n_cells = 40;
    let n_dims = 5;
    let mut embedding = Array2::<f64>::zeros((n_cells, n_dims));
    let mut clusters = vec![0u32; n_cells];

    for i in 0..n_cells {
        if i < 20 {
            for d in 0..n_dims {
                embedding[[i, d]] = 5.0 + (i as f64 * 0.01);
            }
        } else {
            clusters[i] = 1;
            for d in 0..n_dims {
                embedding[[i, d]] = -5.0 + (i as f64 * 0.01);
            }
        }
    }

    let scores = silhouette_scores(&embedding, &clusters);
    let avg = mean_silhouette(&scores);

    assert_eq!(scores.len(), n_cells);
    assert!(avg > 0.8, "well-separated clusters should have high silhouette, got {avg}");
}

#[test]
fn test_full_pipeline_synthetic() {
    let (fm, _true_clusters) = make_planted_matrix(30, 20);

    let preprocess_config = RnaPipelineConfig {
        min_features: 1,
        min_cells: 1,
        max_pct_mt: 100.0,
        n_variable_features: 15,
        n_pcs: 5,
        ..Default::default()
    };

    let downstream_config = DownstreamConfig {
        k_neighbors: 5,
        resolution: 0.8,
        compute_markers: true,
        compute_silhouette: true,
        ..Default::default()
    };

    let result = run_full_rna_pipeline(fm, &preprocess_config, &downstream_config).unwrap();

    // Should produce clusters
    assert!(result.clusters.n_clusters >= 1);
    assert_eq!(
        result.clusters.assignments.len(),
        result.preprocessing.cell_metadata.len()
    );

    // Should have markers if requested
    assert!(result.markers.is_some());

    // Should have silhouette if requested
    assert!(result.silhouette_avg.is_some());

    // Cell metadata should have cluster assignments
    for meta in &result.preprocessing.cell_metadata {
        assert!(meta.cluster.is_some());
    }
}
