use std::collections::HashMap;

use anyhow::Result;
use ndarray::Array2;

use crate::reduce::lanczos::lanczos_svd;
use crate::rna::hvg::find_variable_features;
use crate::rna::normalize::log_normalize;
use crate::rna::qc::{compute_rna_qc, filter_genes, filter_rna_cells};
use crate::rna::scale::scale_data;
use crate::cluster::leiden::leiden_clustering;
use crate::cluster::silhouette::{mean_silhouette, silhouette_scores};
use crate::markers::wilcoxon::find_all_markers;
use crate::reduce::neighbors::{build_knn, build_snn};
use crate::types::{
    AnalysisResult, CellMetadata, DownstreamConfig, Embedding, EmbeddingMethod, FeatureMatrix,
    MarkerConfig, PreprocessingResult, RnaPipelineConfig,
};

/// Run the full RNA preprocessing pipeline.
///
/// Steps: filter_genes → filter_cells → QC → HVG (on counts) → log_normalize → scale → PCA
///
/// HVG selection runs on raw counts before normalization (matching Seurat's VST behavior).
pub fn run_rna_preprocessing(
    matrix: FeatureMatrix,
    config: &RnaPipelineConfig,
) -> Result<PreprocessingResult> {
    // 1. Filter genes (min_cells)
    let matrix = filter_genes(&matrix, config.min_cells)?;

    // 2. QC + filter cells
    let matrix = filter_rna_cells(&matrix, config.min_features, config.max_pct_mt)?;

    // Build cell metadata from QC (post-filter)
    let post_qc = compute_rna_qc(&matrix);
    let cell_metadata: Vec<CellMetadata> = post_qc
        .iter()
        .map(|m| CellMetadata {
            cell_id: m.cell_id.clone(),
            n_features: m.n_features,
            n_counts: m.n_counts,
            cluster: None,
            umap: None,
            custom: {
                let mut map = HashMap::new();
                map.insert("pct_mt".to_string(), m.pct_mt);
                map
            },
        })
        .collect();

    // 3. Find HVGs on raw counts (before normalization, matching Seurat's VST)
    let variable_features = find_variable_features(&matrix, config.n_variable_features)?;

    // 4. Log normalize
    let mut matrix = matrix;
    log_normalize(&mut matrix, config.scale_factor)?;

    // 5. Scale
    let scaled = scale_data(&matrix, &variable_features, config.clip_value, None)?;

    // 6. PCA via SVD
    let n_pcs = config.n_pcs.min(scaled.nrows()).min(scaled.ncols());
    let svd_result = lanczos_svd(&scaled, n_pcs)?;

    // Cell embeddings = Vt^T * diag(S) (Seurat convention: cells × components)
    let n_cells = scaled.ncols();
    let mut embeddings = Array2::<f64>::zeros((n_cells, n_pcs));
    for cell in 0..n_cells {
        for pc in 0..n_pcs {
            embeddings[[cell, pc]] = svd_result.vt[[pc, cell]] * svd_result.s[pc];
        }
    }

    let embedding = Embedding {
        values: embeddings,
        variance_explained: svd_result.variance_explained,
        method: EmbeddingMethod::Pca,
    };

    Ok(PreprocessingResult {
        matrix,
        embedding,
        variable_features,
        cell_metadata,
    })
}

/// Run downstream analysis on preprocessed data.
///
/// Steps: KNN → SNN → Leiden clustering → markers (optional) → silhouette (optional)
///
/// Updates `cell_metadata[i].cluster` with cluster assignments.
pub fn run_rna_downstream(
    result: &mut PreprocessingResult,
    config: &DownstreamConfig,
) -> Result<AnalysisResult> {
    // 1. KNN from PCA embeddings
    let knn = build_knn(&result.embedding.values, config.k_neighbors)?;

    // 2. SNN graph with Jaccard weights
    let snn = build_snn(&knn, config.prune_snn);

    // 3. Leiden clustering
    let clusters = leiden_clustering(&snn, config.resolution, config.leiden_max_iter)?;

    // Update cell metadata with cluster assignments
    for (i, meta) in result.cell_metadata.iter_mut().enumerate() {
        meta.cluster = Some(clusters.assignments[i]);
    }

    // 4. Optional: marker genes
    let markers = if config.compute_markers {
        let marker_config = MarkerConfig {
            min_pct: config.min_pct,
            min_log2fc: config.min_log2fc,
            only_positive: true,
        };
        Some(find_all_markers(
            &result.matrix,
            &clusters.assignments,
            &marker_config,
        )?)
    } else {
        None
    };

    // 5. Optional: silhouette scores
    let (sil_scores, sil_avg) = if config.compute_silhouette {
        let scores = silhouette_scores(&result.embedding.values, &clusters.assignments);
        let avg = mean_silhouette(&scores);
        (Some(scores), Some(avg))
    } else {
        (None, None)
    };

    Ok(AnalysisResult {
        preprocessing: result.clone(),
        knn,
        clusters,
        markers,
        silhouette_scores: sil_scores,
        silhouette_avg: sil_avg,
    })
}

/// Run the full RNA analysis pipeline (preprocessing + downstream).
pub fn run_full_rna_pipeline(
    matrix: FeatureMatrix,
    preprocess_config: &RnaPipelineConfig,
    downstream_config: &DownstreamConfig,
) -> Result<AnalysisResult> {
    let mut prep = run_rna_preprocessing(matrix, preprocess_config)?;
    run_rna_downstream(&mut prep, downstream_config)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::RnaPipelineConfig;
    use crate::io::read_10x;
    use std::path::PathBuf;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/small_10x")
    }

    #[test]
    fn test_pipeline_runs_on_small_data() {
        let fm = read_10x(&test_data_dir()).unwrap();
        let config = RnaPipelineConfig {
            min_features: 1, // very low for 3-cell test
            min_cells: 1,
            max_pct_mt: 100.0,
            n_variable_features: 3,
            n_pcs: 2,
            ..Default::default()
        };
        let result = run_rna_preprocessing(fm, &config).unwrap();
        assert_eq!(result.embedding.values.ncols(), 2);
        assert_eq!(result.embedding.method, EmbeddingMethod::Pca);
        assert!(!result.variable_features.is_empty());
    }
}
