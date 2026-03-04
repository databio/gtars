use std::collections::HashMap;

use sprs::CsMat;

/// Sparse feature-by-cell matrix with metadata.
///
/// Stored as CSC (compressed sparse column) where each column is a cell.
/// Orientation: features (rows) × cells (columns).
#[derive(Debug, Clone)]
pub struct FeatureMatrix {
    pub matrix: CsMat<f64>,
    pub feature_names: Vec<String>,
    pub feature_ids: Vec<String>,
    pub cell_ids: Vec<String>,
    pub feature_type: FeatureType,
}

impl FeatureMatrix {
    pub fn n_features(&self) -> usize {
        self.matrix.rows()
    }

    pub fn n_cells(&self) -> usize {
        self.matrix.cols()
    }

    pub fn shape(&self) -> (usize, usize) {
        (self.n_features(), self.n_cells())
    }
}

/// Type of features stored in a FeatureMatrix.
#[derive(Debug, Clone, PartialEq)]
pub enum FeatureType {
    Gene,
    Peak,
    Custom(String),
}

/// Dense embedding matrix (cells × components) with variance explained.
#[derive(Debug, Clone)]
pub struct Embedding {
    pub values: ndarray::Array2<f64>,
    pub variance_explained: Vec<f64>,
    pub method: EmbeddingMethod,
}

/// Method used to produce an embedding.
#[derive(Debug, Clone, PartialEq)]
pub enum EmbeddingMethod {
    Pca,
    Lsi { component_1_dropped: bool },
}

/// Per-cell metadata collected during preprocessing.
#[derive(Debug, Clone)]
pub struct CellMetadata {
    pub cell_id: String,
    pub n_features: u32,
    pub n_counts: f64,
    pub cluster: Option<u32>,
    pub umap: Option<(f64, f64)>,
    pub custom: HashMap<String, f64>,
}

/// QC metrics for a single cell (RNA assay).
#[derive(Debug, Clone)]
pub struct RnaQcMetrics {
    pub cell_id: String,
    pub n_features: u32,
    pub n_counts: f64,
    pub pct_mt: f64,
}

/// K-nearest-neighbor graph (placeholder for Phase 2).
#[derive(Debug, Clone)]
pub struct KnnGraph {
    pub n_cells: usize,
    pub k: usize,
    pub indices: Vec<Vec<usize>>,
    pub distances: Vec<Vec<f64>>,
}

/// Differential expression marker result (placeholder for Phase 2).
#[derive(Debug, Clone)]
pub struct MarkerResult {
    pub gene: String,
    pub cluster: u32,
    pub avg_log2fc: f64,
    pub pval: f64,
    pub pval_adj: f64,
    pub pct_in: f64,
    pub pct_out: f64,
}

/// Result of the SVD decomposition.
#[derive(Debug, Clone)]
pub struct SvdResult {
    pub u: ndarray::Array2<f64>,
    pub s: ndarray::Array1<f64>,
    pub vt: ndarray::Array2<f64>,
    pub variance_explained: Vec<f64>,
}

/// Configuration for the RNA preprocessing pipeline.
#[derive(Debug, Clone)]
pub struct RnaPipelineConfig {
    pub min_features: u32,
    pub min_cells: u32,
    pub max_pct_mt: f64,
    pub scale_factor: f64,
    pub n_variable_features: usize,
    pub n_pcs: usize,
    pub clip_value: Option<f64>,
}

impl Default for RnaPipelineConfig {
    fn default() -> Self {
        Self {
            min_features: 200,
            min_cells: 3,
            max_pct_mt: 5.0,
            scale_factor: 10_000.0,
            n_variable_features: 2000,
            n_pcs: 50,
            clip_value: Some(10.0),
        }
    }
}

/// Full result of preprocessing pipeline.
#[derive(Debug, Clone)]
pub struct PreprocessingResult {
    pub matrix: FeatureMatrix,
    pub embedding: Embedding,
    pub variable_features: Vec<String>,
    pub cell_metadata: Vec<CellMetadata>,
}
