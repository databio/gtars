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

/// Weighted SNN (Shared Nearest Neighbor) graph for clustering.
///
/// Edge weights represent Jaccard similarity of KNN neighborhoods.
/// Edges below the prune threshold are excluded.
#[derive(Debug, Clone)]
pub struct SnnGraph {
    pub n_cells: usize,
    /// Adjacency list: for each cell, `(neighbor_index, jaccard_weight)`.
    pub edges: Vec<Vec<(usize, f64)>>,
}

/// Result of community detection (clustering).
#[derive(Debug, Clone)]
pub struct ClusterResult {
    /// Cluster assignment per cell (0-indexed).
    pub assignments: Vec<u32>,
    /// Number of distinct clusters.
    pub n_clusters: u32,
    /// Modularity / quality of the partition (if computed).
    pub quality: Option<f64>,
}

/// Configuration for the downstream analysis pipeline.
#[derive(Debug, Clone)]
pub struct DownstreamConfig {
    /// Number of nearest neighbors for KNN graph (default 20).
    pub k_neighbors: usize,
    /// SNN pruning threshold — edges with Jaccard below this are dropped (default 1/15).
    pub prune_snn: f64,
    /// Leiden resolution parameter (default 0.8).
    pub resolution: f64,
    /// Maximum Leiden iterations (default 10).
    pub leiden_max_iter: usize,
    /// Minimum fraction of cells expressing a gene to test (default 0.1).
    pub min_pct: f64,
    /// Minimum log2 fold change for marker genes (default 0.25).
    pub min_log2fc: f64,
    /// Whether to compute marker genes (default true).
    pub compute_markers: bool,
    /// Whether to compute silhouette scores (default true).
    pub compute_silhouette: bool,
}

impl Default for DownstreamConfig {
    fn default() -> Self {
        Self {
            k_neighbors: 20,
            prune_snn: 1.0 / 15.0,
            resolution: 0.8,
            leiden_max_iter: 10,
            min_pct: 0.1,
            min_log2fc: 0.25,
            compute_markers: true,
            compute_silhouette: true,
        }
    }
}

/// Configuration for marker gene detection.
#[derive(Debug, Clone)]
pub struct MarkerConfig {
    /// Minimum fraction of cells expressing gene in either group (default 0.1).
    pub min_pct: f64,
    /// Minimum absolute log2 fold change (default 0.25).
    pub min_log2fc: f64,
    /// Only return upregulated markers (default true).
    pub only_positive: bool,
}

impl Default for MarkerConfig {
    fn default() -> Self {
        Self {
            min_pct: 0.1,
            min_log2fc: 0.25,
            only_positive: true,
        }
    }
}

/// Full analysis result (preprocessing + downstream).
#[derive(Debug, Clone)]
pub struct AnalysisResult {
    pub preprocessing: PreprocessingResult,
    pub knn: KnnGraph,
    pub clusters: ClusterResult,
    pub markers: Option<Vec<MarkerResult>>,
    pub silhouette_scores: Option<Vec<f64>>,
    pub silhouette_avg: Option<f64>,
}
