use extendr_api::prelude::*;
use ndarray::Array2;

use gtars_sc::cluster::leiden::leiden_clustering;
use gtars_sc::cluster::silhouette::{mean_silhouette, silhouette_scores};
use gtars_sc::io::read_10x;
use gtars_sc::markers::wilcoxon::find_all_markers;
use gtars_sc::reduce::lanczos::lanczos_svd;
use gtars_sc::reduce::neighbors::{build_knn, build_snn};
use gtars_sc::rna::hvg::find_variable_features;
use gtars_sc::rna::normalize::log_normalize;
use gtars_sc::rna::pipeline::{run_full_rna_pipeline, run_rna_preprocessing};
use gtars_sc::rna::qc::{compute_rna_qc, filter_genes, filter_rna_cells};
use gtars_sc::rna::scale::scale_data;
use gtars_sc::types::{
    DownstreamConfig, FeatureMatrix, FeatureType, MarkerConfig, RnaPipelineConfig, SnnGraph,
};

// =========================================================================
// Helper macro
// =========================================================================

macro_rules! with_feature_matrix {
    ($ptr:expr, $fm:ident, $body:block) => {{
        let ext_ptr = <ExternalPtr<FeatureMatrix>>::try_from($ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid FeatureMatrix pointer".into()))?;
        let $fm = &*ext_ptr;
        $body
    }};
}

macro_rules! with_snn_graph {
    ($ptr:expr, $g:ident, $body:block) => {{
        let ext_ptr = <ExternalPtr<SnnGraph>>::try_from($ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid SnnGraph pointer".into()))?;
        let $g = &*ext_ptr;
        $body
    }};
}

// =========================================================================
// 1. Read 10X
// =========================================================================

/// Read a 10X Chromium directory (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
/// @export
/// @param dir Path to the 10X directory
#[extendr(r_name = "sc_read_10x")]
pub fn r_sc_read_10x(dir: &str) -> extendr_api::Result<Robj> {
    let path = std::path::PathBuf::from(dir);
    let fm = read_10x(&path)
        .map_err(|e| extendr_api::Error::Other(format!("Reading 10X: {}", e)))?;
    Ok(ExternalPtr::new(fm).into())
}

// =========================================================================
// 2. Matrix info
// =========================================================================

/// Get FeatureMatrix dimensions and metadata
/// @export
/// @param ptr External pointer to a FeatureMatrix
#[extendr(r_name = "sc_feature_matrix_info")]
pub fn r_sc_feature_matrix_info(ptr: Robj) -> extendr_api::Result<List> {
    with_feature_matrix!(ptr, fm, {
        let ft = match &fm.feature_type {
            FeatureType::Gene => "Gene".to_string(),
            FeatureType::Peak => "Peak".to_string(),
            FeatureType::Custom(s) => s.clone(),
        };
        Ok(list!(
            n_features = fm.n_features() as i32,
            n_cells = fm.n_cells() as i32,
            nnz = fm.matrix.nnz() as i32,
            feature_type = ft
        ))
    })
}

// =========================================================================
// 3. QC metrics
// =========================================================================

/// Compute per-cell RNA QC metrics (n_features, n_counts, pct_mt)
/// @export
/// @param ptr External pointer to a FeatureMatrix
#[extendr(r_name = "sc_compute_rna_qc")]
pub fn r_sc_compute_rna_qc(ptr: Robj) -> extendr_api::Result<List> {
    with_feature_matrix!(ptr, fm, {
        let qc = compute_rna_qc(fm);
        let cell_ids: Vec<String> = qc.iter().map(|m| m.cell_id.clone()).collect();
        let n_features: Vec<i32> = qc.iter().map(|m| m.n_features as i32).collect();
        let n_counts: Vec<f64> = qc.iter().map(|m| m.n_counts).collect();
        let pct_mt: Vec<f64> = qc.iter().map(|m| m.pct_mt).collect();
        Ok(list!(
            cell_id = cell_ids,
            n_features = n_features,
            n_counts = n_counts,
            pct_mt = pct_mt
        ))
    })
}

// =========================================================================
// 4. Filter genes
// =========================================================================

/// Filter genes appearing in fewer than min_cells cells
/// @export
/// @param ptr External pointer to a FeatureMatrix
/// @param min_cells Minimum number of cells a gene must appear in
#[extendr(r_name = "sc_filter_genes")]
pub fn r_sc_filter_genes(ptr: Robj, min_cells: i32) -> extendr_api::Result<Robj> {
    with_feature_matrix!(ptr, fm, {
        let filtered = filter_genes(fm, min_cells as u32)
            .map_err(|e| extendr_api::Error::Other(format!("filter_genes: {}", e)))?;
        Ok(ExternalPtr::new(filtered).into())
    })
}

// =========================================================================
// 5. Filter cells
// =========================================================================

/// Filter cells by minimum features and maximum mitochondrial percentage
/// @export
/// @param ptr External pointer to a FeatureMatrix
/// @param min_features Minimum number of features per cell
/// @param max_pct_mt Maximum mitochondrial percentage
#[extendr(r_name = "sc_filter_cells")]
pub fn r_sc_filter_cells(
    ptr: Robj,
    min_features: i32,
    max_pct_mt: f64,
) -> extendr_api::Result<Robj> {
    with_feature_matrix!(ptr, fm, {
        let filtered = filter_rna_cells(fm, min_features as u32, max_pct_mt)
            .map_err(|e| extendr_api::Error::Other(format!("filter_cells: {}", e)))?;
        Ok(ExternalPtr::new(filtered).into())
    })
}

// =========================================================================
// 6. Log normalize
// =========================================================================

/// Log-normalize count data: log1p(count / total * scale_factor)
/// @export
/// @param ptr External pointer to a FeatureMatrix
/// @param scale_factor Scale factor (default 10000)
#[extendr(r_name = "sc_log_normalize")]
pub fn r_sc_log_normalize(ptr: Robj, scale_factor: f64) -> extendr_api::Result<Robj> {
    with_feature_matrix!(ptr, fm, {
        let mut normalized = fm.clone();
        log_normalize(&mut normalized, scale_factor)
            .map_err(|e| extendr_api::Error::Other(e.to_string()))?;
        Ok(ExternalPtr::new(normalized).into())
    })
}

// =========================================================================
// 7. Find variable features (HVGs)
// =========================================================================

/// Find highly variable features using variance-stabilizing transform
/// @export
/// @param ptr External pointer to a FeatureMatrix
/// @param n_features Number of variable features to select
#[extendr(r_name = "sc_find_variable_features")]
pub fn r_sc_find_variable_features(ptr: Robj, n_features: i32) -> extendr_api::Result<Vec<String>> {
    with_feature_matrix!(ptr, fm, {
        let hvgs = find_variable_features(fm, n_features as usize)
            .map_err(|e| extendr_api::Error::Other(format!("find_variable_features: {}", e)))?;
        Ok(hvgs)
    })
}

// =========================================================================
// 8. Scale data
// =========================================================================

/// Scale data: subset to features, center, scale, optionally clip
/// @export
/// @param ptr External pointer to a FeatureMatrix
/// @param features Character vector of feature names to scale
/// @param clip_value Maximum absolute value to clip to (use NA or NULL to skip)
#[extendr(r_name = "sc_scale_data")]
pub fn r_sc_scale_data(
    ptr: Robj,
    features: Vec<String>,
    clip_value: Robj,
) -> extendr_api::Result<List> {
    with_feature_matrix!(ptr, fm, {
        let clip = if clip_value.is_na() || clip_value.is_null() {
            None
        } else {
            Some(clip_value.as_real()
                .ok_or_else(|| extendr_api::Error::Other("clip_value must be numeric or NA".into()))?)
        };

        let scaled = scale_data(fm, &features, clip, None)
            .map_err(|e| extendr_api::Error::Other(format!("scale_data: {}", e)))?;

        let (nrows, ncols) = scaled.dim();

        // Return as flat column-major vector + dimensions (wrap with matrix() in R)
        let mut col_major = Vec::with_capacity(nrows * ncols);
        for col in 0..ncols {
            for row in 0..nrows {
                col_major.push(scaled[[row, col]]);
            }
        }

        // Collect cell_ids from the FeatureMatrix
        let cell_ids: Vec<String> = fm.cell_ids.clone();

        Ok(list!(
            data = col_major,
            nrow = nrows as i32,
            ncol = ncols as i32,
            features = features,
            cells = cell_ids
        ))
    })
}

// =========================================================================
// 9. PCA
// =========================================================================

/// Run PCA via truncated SVD on a scaled matrix
/// @export
/// @param matrix_data Flat numeric vector (column-major) of scaled data
/// @param nrow Number of rows (features)
/// @param ncol Number of columns (cells)
/// @param n_pcs Number of principal components
#[extendr(r_name = "sc_run_pca")]
pub fn r_sc_run_pca(
    matrix_data: Vec<f64>,
    nrow: i32,
    ncol: i32,
    n_pcs: i32,
) -> extendr_api::Result<List> {
    let nrows = nrow as usize;
    let ncols = ncol as usize;
    let n_components = n_pcs as usize;

    // Reconstruct row-major Array2 from column-major flat vector
    let mut scaled = Array2::<f64>::zeros((nrows, ncols));
    for col in 0..ncols {
        for row in 0..nrows {
            scaled[[row, col]] = matrix_data[col * nrows + row];
        }
    }

    let n_components = n_components.min(nrows).min(ncols);
    let svd_result = lanczos_svd(&scaled, n_components)
        .map_err(|e| extendr_api::Error::Other(format!("PCA/SVD: {}", e)))?;

    // Cell embeddings = Vt^T * diag(S) (cells × components, Seurat convention)
    let mut embeddings = Vec::with_capacity(ncols * n_components);
    // Column-major: iterate PC then cell
    for pc in 0..n_components {
        for cell in 0..ncols {
            embeddings.push(svd_result.vt[[pc, cell]] * svd_result.s[pc]);
        }
    }

    let variance_explained: Vec<f64> = svd_result.variance_explained;

    Ok(list!(
        embeddings = embeddings,
        nrow = ncols as i32,
        ncol = n_components as i32,
        variance_explained = variance_explained
    ))
}

// =========================================================================
// 10. Full pipeline
// =========================================================================

/// Run the full RNA preprocessing pipeline
/// @export
/// @param ptr External pointer to a FeatureMatrix
/// @param min_features Minimum features per cell (default 200)
/// @param min_cells Minimum cells per gene (default 3)
/// @param max_pct_mt Maximum mitochondrial percentage (default 5.0)
/// @param scale_factor Normalization scale factor (default 10000)
/// @param n_variable_features Number of HVGs (default 2000)
/// @param n_pcs Number of PCs (default 50)
/// @param clip_value Clip value for scaling (default 10.0, NA to skip)
#[extendr(r_name = "sc_run_rna_pipeline")]
pub fn r_sc_run_rna_pipeline(
    ptr: Robj,
    min_features: i32,
    min_cells: i32,
    max_pct_mt: f64,
    scale_factor: f64,
    n_variable_features: i32,
    n_pcs: i32,
    clip_value: Robj,
) -> extendr_api::Result<List> {
    with_feature_matrix!(ptr, fm, {
        let clip = if clip_value.is_na() || clip_value.is_null() {
            None
        } else {
            Some(clip_value.as_real()
                .ok_or_else(|| extendr_api::Error::Other("clip_value must be numeric or NA".into()))?)
        };

        let config = RnaPipelineConfig {
            min_features: min_features as u32,
            min_cells: min_cells as u32,
            max_pct_mt,
            scale_factor,
            n_variable_features: n_variable_features as usize,
            n_pcs: n_pcs as usize,
            clip_value: clip,
        };

        let result = run_rna_preprocessing(fm.clone(), &config)
            .map_err(|e| extendr_api::Error::Other(format!("pipeline: {}", e)))?;

        // Embedding (cells × PCs) as column-major flat vector
        let (emb_nrows, emb_ncols) = result.embedding.values.dim();
        let mut emb_data = Vec::with_capacity(emb_nrows * emb_ncols);
        for col in 0..emb_ncols {
            for row in 0..emb_nrows {
                emb_data.push(result.embedding.values[[row, col]]);
            }
        }

        // Cell metadata
        let cell_ids: Vec<String> = result.cell_metadata.iter().map(|m| m.cell_id.clone()).collect();
        let n_features_vec: Vec<i32> = result.cell_metadata.iter().map(|m| m.n_features as i32).collect();
        let n_counts_vec: Vec<f64> = result.cell_metadata.iter().map(|m| m.n_counts).collect();
        let pct_mt_vec: Vec<f64> = result.cell_metadata.iter().map(|m| {
            m.custom.get("pct_mt").copied().unwrap_or(0.0)
        }).collect();

        // Matrix info
        let final_n_features = result.matrix.n_features() as i32;
        let final_n_cells = result.matrix.n_cells() as i32;

        Ok(list!(
            embeddings = emb_data,
            emb_nrow = emb_nrows as i32,
            emb_ncol = emb_ncols as i32,
            variance_explained = result.embedding.variance_explained,
            variable_features = result.variable_features,
            cell_ids = cell_ids,
            n_features = n_features_vec,
            n_counts = n_counts_vec,
            pct_mt = pct_mt_vec,
            final_n_features = final_n_features,
            final_n_cells = final_n_cells
        ))
    })
}

// =========================================================================
// 11. Find neighbors (KNN → SNN)
// =========================================================================

/// Build KNN + SNN graph from a PCA embedding matrix
/// @export
/// @param embedding_data Flat numeric vector (column-major) of cell embeddings
/// @param nrow Number of rows (cells)
/// @param ncol Number of columns (PCs)
/// @param k Number of nearest neighbors (default 20)
/// @param prune_snn Minimum Jaccard index to keep an SNN edge (default 1/15)
#[extendr(r_name = "sc_find_neighbors")]
pub fn r_sc_find_neighbors(
    embedding_data: Vec<f64>,
    nrow: i32,
    ncol: i32,
    k: i32,
    prune_snn: f64,
) -> extendr_api::Result<Robj> {
    let nrows = nrow as usize;
    let ncols = ncol as usize;

    // Reconstruct row-major Array2 from column-major flat vector
    let mut embedding = Array2::<f64>::zeros((nrows, ncols));
    for col in 0..ncols {
        for row in 0..nrows {
            embedding[[row, col]] = embedding_data[col * nrows + row];
        }
    }

    let knn = build_knn(&embedding, k as usize)
        .map_err(|e| extendr_api::Error::Other(format!("build_knn: {}", e)))?;
    let snn = build_snn(&knn, prune_snn);

    Ok(ExternalPtr::new(snn).into())
}

// =========================================================================
// 12. Find clusters (Leiden)
// =========================================================================

/// Run Leiden clustering on an SNN graph
/// @export
/// @param snn_ptr External pointer to an SnnGraph
/// @param resolution Resolution parameter (default 0.8)
/// @param max_iter Maximum iterations (default 10)
#[extendr(r_name = "sc_find_clusters")]
pub fn r_sc_find_clusters(
    snn_ptr: Robj,
    resolution: f64,
    max_iter: i32,
) -> extendr_api::Result<List> {
    with_snn_graph!(snn_ptr, snn, {
        let result = leiden_clustering(snn, resolution, max_iter as usize)
            .map_err(|e| extendr_api::Error::Other(format!("leiden: {}", e)))?;

        let assignments: Vec<i32> = result.assignments.iter().map(|&c| c as i32).collect();

        Ok(list!(
            assignments = assignments,
            n_clusters = result.n_clusters as i32,
            quality = result.quality.unwrap_or(0.0)
        ))
    })
}

// =========================================================================
// 13. Find all markers
// =========================================================================

/// Find marker genes for each cluster using Wilcoxon rank-sum test
/// @export
/// @param ptr External pointer to a FeatureMatrix (log-normalized)
/// @param clusters Integer vector of cluster assignments (0-indexed)
/// @param min_pct Minimum fraction of cells expressing the gene (default 0.1)
/// @param min_log2fc Minimum log2 fold-change (default 0.25)
#[extendr(r_name = "sc_find_all_markers")]
pub fn r_sc_find_all_markers(
    ptr: Robj,
    clusters: Vec<i32>,
    min_pct: f64,
    min_log2fc: f64,
) -> extendr_api::Result<List> {
    with_feature_matrix!(ptr, fm, {
        let clusters_u32: Vec<u32> = clusters.iter().map(|&c| c as u32).collect();
        let config = MarkerConfig {
            min_pct,
            min_log2fc,
            only_positive: true,
        };

        let markers = find_all_markers(fm, &clusters_u32, &config)
            .map_err(|e| extendr_api::Error::Other(format!("find_all_markers: {}", e)))?;

        let genes: Vec<String> = markers.iter().map(|m| m.gene.clone()).collect();
        let cluster_ids: Vec<i32> = markers.iter().map(|m| m.cluster as i32).collect();
        let avg_log2fc: Vec<f64> = markers.iter().map(|m| m.avg_log2fc).collect();
        let pval: Vec<f64> = markers.iter().map(|m| m.pval).collect();
        let pval_adj: Vec<f64> = markers.iter().map(|m| m.pval_adj).collect();
        let pct_in: Vec<f64> = markers.iter().map(|m| m.pct_in).collect();
        let pct_out: Vec<f64> = markers.iter().map(|m| m.pct_out).collect();

        Ok(list!(
            gene = genes,
            cluster = cluster_ids,
            avg_log2fc = avg_log2fc,
            pval = pval,
            pval_adj = pval_adj,
            pct_in = pct_in,
            pct_out = pct_out
        ))
    })
}

// =========================================================================
// 14. Silhouette scores
// =========================================================================

/// Compute silhouette scores for cluster quality assessment
/// @export
/// @param embedding_data Flat numeric vector (column-major) of cell embeddings
/// @param nrow Number of rows (cells)
/// @param ncol Number of columns (PCs)
/// @param clusters Integer vector of cluster assignments (0-indexed)
#[extendr(r_name = "sc_silhouette")]
pub fn r_sc_silhouette(
    embedding_data: Vec<f64>,
    nrow: i32,
    ncol: i32,
    clusters: Vec<i32>,
) -> extendr_api::Result<List> {
    let nrows = nrow as usize;
    let ncols = ncol as usize;

    let mut embedding = Array2::<f64>::zeros((nrows, ncols));
    for col in 0..ncols {
        for row in 0..nrows {
            embedding[[row, col]] = embedding_data[col * nrows + row];
        }
    }

    let clusters_u32: Vec<u32> = clusters.iter().map(|&c| c as u32).collect();
    let scores = silhouette_scores(&embedding, &clusters_u32);
    let avg = mean_silhouette(&scores);

    Ok(list!(
        scores = scores,
        avg = avg
    ))
}

// =========================================================================
// 15. Full pipeline (preprocessing + downstream)
// =========================================================================

/// Run the full RNA analysis pipeline (preprocessing + clustering + markers)
/// @export
/// @param ptr External pointer to a FeatureMatrix
/// @param min_features Minimum features per cell (default 200)
/// @param min_cells Minimum cells per gene (default 3)
/// @param max_pct_mt Maximum mitochondrial percentage (default 5.0)
/// @param scale_factor Normalization scale factor (default 10000)
/// @param n_variable_features Number of HVGs (default 2000)
/// @param n_pcs Number of PCs (default 50)
/// @param clip_value Clip value for scaling (default 10.0, NA to skip)
/// @param k_neighbors Number of nearest neighbors (default 20)
/// @param resolution Leiden resolution (default 0.8)
/// @param compute_markers Whether to compute marker genes (default TRUE)
/// @param compute_silhouette Whether to compute silhouette scores (default TRUE)
#[extendr(r_name = "sc_run_full_pipeline")]
pub fn r_sc_run_full_pipeline(
    ptr: Robj,
    min_features: i32,
    min_cells: i32,
    max_pct_mt: f64,
    scale_factor: f64,
    n_variable_features: i32,
    n_pcs: i32,
    clip_value: Robj,
    k_neighbors: i32,
    resolution: f64,
    compute_markers: bool,
    compute_silhouette: bool,
) -> extendr_api::Result<List> {
    with_feature_matrix!(ptr, fm, {
        let clip = if clip_value.is_na() || clip_value.is_null() {
            None
        } else {
            Some(clip_value.as_real()
                .ok_or_else(|| extendr_api::Error::Other("clip_value must be numeric or NA".into()))?)
        };

        let preprocess_config = RnaPipelineConfig {
            min_features: min_features as u32,
            min_cells: min_cells as u32,
            max_pct_mt,
            scale_factor,
            n_variable_features: n_variable_features as usize,
            n_pcs: n_pcs as usize,
            clip_value: clip,
        };

        let downstream_config = DownstreamConfig {
            k_neighbors: k_neighbors as usize,
            resolution,
            compute_markers,
            compute_silhouette,
            ..Default::default()
        };

        let result = run_full_rna_pipeline(fm.clone(), &preprocess_config, &downstream_config)
            .map_err(|e| extendr_api::Error::Other(format!("full_pipeline: {}", e)))?;

        // Embedding (cells × PCs) as column-major flat vector
        let (emb_nrows, emb_ncols) = result.preprocessing.embedding.values.dim();
        let mut emb_data = Vec::with_capacity(emb_nrows * emb_ncols);
        for col in 0..emb_ncols {
            for row in 0..emb_nrows {
                emb_data.push(result.preprocessing.embedding.values[[row, col]]);
            }
        }

        // Cell metadata
        let cell_ids: Vec<String> = result.preprocessing.cell_metadata.iter().map(|m| m.cell_id.clone()).collect();
        let n_features_vec: Vec<i32> = result.preprocessing.cell_metadata.iter().map(|m| m.n_features as i32).collect();
        let n_counts_vec: Vec<f64> = result.preprocessing.cell_metadata.iter().map(|m| m.n_counts).collect();
        let pct_mt_vec: Vec<f64> = result.preprocessing.cell_metadata.iter().map(|m| {
            m.custom.get("pct_mt").copied().unwrap_or(0.0)
        }).collect();
        let cluster_vec: Vec<i32> = result.preprocessing.cell_metadata.iter().map(|m| {
            m.cluster.map(|c| c as i32).unwrap_or(-1)
        }).collect();

        // Clusters
        let n_clusters = result.clusters.n_clusters as i32;
        let quality = result.clusters.quality.unwrap_or(0.0);

        // Markers (if computed)
        let (marker_genes, marker_clusters, marker_log2fc, marker_pval, marker_pval_adj, marker_pct_in, marker_pct_out) =
            if let Some(ref markers) = result.markers {
                (
                    markers.iter().map(|m| m.gene.clone()).collect::<Vec<_>>(),
                    markers.iter().map(|m| m.cluster as i32).collect::<Vec<_>>(),
                    markers.iter().map(|m| m.avg_log2fc).collect::<Vec<_>>(),
                    markers.iter().map(|m| m.pval).collect::<Vec<_>>(),
                    markers.iter().map(|m| m.pval_adj).collect::<Vec<_>>(),
                    markers.iter().map(|m| m.pct_in).collect::<Vec<_>>(),
                    markers.iter().map(|m| m.pct_out).collect::<Vec<_>>(),
                )
            } else {
                (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new())
            };

        // Silhouette
        let sil_avg = result.silhouette_avg.unwrap_or(f64::NAN);

        Ok(list!(
            embeddings = emb_data,
            emb_nrow = emb_nrows as i32,
            emb_ncol = emb_ncols as i32,
            variance_explained = result.preprocessing.embedding.variance_explained,
            variable_features = result.preprocessing.variable_features,
            cell_ids = cell_ids,
            n_features = n_features_vec,
            n_counts = n_counts_vec,
            pct_mt = pct_mt_vec,
            clusters = cluster_vec,
            n_clusters = n_clusters,
            modularity = quality,
            marker_gene = marker_genes,
            marker_cluster = marker_clusters,
            marker_avg_log2fc = marker_log2fc,
            marker_pval = marker_pval,
            marker_pval_adj = marker_pval_adj,
            marker_pct_in = marker_pct_in,
            marker_pct_out = marker_pct_out,
            silhouette_avg = sil_avg
        ))
    })
}

// =========================================================================
// Module registration
// =========================================================================

extendr_module! {
    mod sc;
    fn r_sc_read_10x;
    fn r_sc_feature_matrix_info;
    fn r_sc_compute_rna_qc;
    fn r_sc_filter_genes;
    fn r_sc_filter_cells;
    fn r_sc_log_normalize;
    fn r_sc_find_variable_features;
    fn r_sc_scale_data;
    fn r_sc_run_pca;
    fn r_sc_run_rna_pipeline;
    fn r_sc_find_neighbors;
    fn r_sc_find_clusters;
    fn r_sc_find_all_markers;
    fn r_sc_silhouette;
    fn r_sc_run_full_pipeline;
}
