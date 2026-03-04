use std::collections::HashMap;

use anyhow::Result;
use ndarray::Array2;

use crate::reduce::lanczos::lanczos_svd;
use crate::rna::hvg::find_variable_features;
use crate::rna::normalize::log_normalize;
use crate::rna::qc::{compute_rna_qc, filter_genes, filter_rna_cells};
use crate::rna::scale::scale_data;
use crate::types::{
    CellMetadata, Embedding, EmbeddingMethod, FeatureMatrix, PreprocessingResult, RnaPipelineConfig,
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
    log_normalize(&mut matrix, config.scale_factor);

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
