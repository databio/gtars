# =========================================================================
# Single-cell RNA analysis (gtars-sc bindings)
# =========================================================================

#' Read a 10X Chromium directory
#'
#' @param dir Path to a 10X directory containing matrix.mtx.gz,
#'   features.tsv.gz, and barcodes.tsv.gz
#' @return An external pointer to a FeatureMatrix
#' @export
sc_read_10x <- function(dir) {
  stopifnot(is.character(dir), length(dir) == 1)
  .Call(wrap__r_sc_read_10x, dir)
}

#' Get FeatureMatrix info
#'
#' @param ptr External pointer to a FeatureMatrix
#' @return A named list with n_features, n_cells, nnz, feature_type
#' @export
sc_feature_matrix_info <- function(ptr) {
  .Call(wrap__r_sc_feature_matrix_info, ptr)
}

#' Compute per-cell RNA QC metrics
#'
#' @param ptr External pointer to a FeatureMatrix
#' @return A data.frame with columns: cell_id, n_features, n_counts, pct_mt
#' @export
sc_compute_rna_qc <- function(ptr) {
  res <- .Call(wrap__r_sc_compute_rna_qc, ptr)
  data.frame(
    cell_id = res$cell_id,
    n_features = res$n_features,
    n_counts = res$n_counts,
    pct_mt = res$pct_mt,
    stringsAsFactors = FALSE
  )
}

#' Filter genes by minimum cell count
#'
#' @param ptr External pointer to a FeatureMatrix
#' @param min_cells Minimum number of cells a gene must appear in (default 3)
#' @return External pointer to a new filtered FeatureMatrix
#' @export
sc_filter_genes <- function(ptr, min_cells = 3L) {
  .Call(wrap__r_sc_filter_genes, ptr, as.integer(min_cells))
}

#' Filter cells by QC criteria
#'
#' @param ptr External pointer to a FeatureMatrix
#' @param min_features Minimum features per cell (default 200)
#' @param max_pct_mt Maximum mitochondrial percentage (default 5.0)
#' @return External pointer to a new filtered FeatureMatrix
#' @export
sc_filter_cells <- function(ptr, min_features = 200L, max_pct_mt = 5.0) {
  .Call(wrap__r_sc_filter_cells, ptr, as.integer(min_features), as.double(max_pct_mt))
}

#' Log-normalize count data
#'
#' @param ptr External pointer to a FeatureMatrix
#' @param scale_factor Scale factor for normalization (default 10000)
#' @return External pointer to a new normalized FeatureMatrix
#' @export
sc_log_normalize <- function(ptr, scale_factor = 10000) {
  .Call(wrap__r_sc_log_normalize, ptr, as.double(scale_factor))
}

#' Find highly variable features
#'
#' @param ptr External pointer to a FeatureMatrix
#' @param n_features Number of variable features to select (default 2000)
#' @return Character vector of variable feature names
#' @export
sc_find_variable_features <- function(ptr, n_features = 2000L) {
  .Call(wrap__r_sc_find_variable_features, ptr, as.integer(n_features))
}

#' Scale data for selected features
#'
#' @param ptr External pointer to a FeatureMatrix
#' @param features Character vector of feature names to scale
#' @param clip_value Maximum absolute value to clip (default 10, use NA to skip)
#' @return A dense matrix (features x cells) with dimnames
#' @export
sc_scale_data <- function(ptr, features, clip_value = 10.0) {
  stopifnot(is.character(features), length(features) > 0)
  res <- .Call(wrap__r_sc_scale_data, ptr, features, clip_value)
  mat <- matrix(res$data, nrow = res$nrow, ncol = res$ncol)
  rownames(mat) <- res$features
  colnames(mat) <- res$cells
  mat
}

#' Run PCA via truncated SVD
#'
#' @param scaled_matrix A dense matrix (features x cells) from sc_scale_data
#' @param n_pcs Number of principal components (default 50)
#' @return A named list with:
#'   \item{embeddings}{Matrix (cells x PCs)}
#'   \item{variance_explained}{Numeric vector of variance explained per PC}
#' @export
sc_run_pca <- function(scaled_matrix, n_pcs = 50L) {
  stopifnot(is.matrix(scaled_matrix), is.numeric(scaled_matrix))
  nr <- nrow(scaled_matrix)
  nc <- ncol(scaled_matrix)
  # Pass as column-major flat vector (R's native storage)
  res <- .Call(wrap__r_sc_run_pca, as.double(scaled_matrix), as.integer(nr),
               as.integer(nc), as.integer(n_pcs))
  emb <- matrix(res$embeddings, nrow = res$nrow, ncol = res$ncol)
  colnames(emb) <- paste0("PC_", seq_len(res$ncol))
  if (!is.null(colnames(scaled_matrix))) {
    rownames(emb) <- colnames(scaled_matrix)
  }
  list(
    embeddings = emb,
    variance_explained = res$variance_explained
  )
}

#' Run full RNA preprocessing pipeline
#'
#' @param ptr External pointer to a FeatureMatrix
#' @param min_features Minimum features per cell (default 200)
#' @param min_cells Minimum cells per gene (default 3)
#' @param max_pct_mt Maximum mitochondrial percentage (default 5.0)
#' @param scale_factor Normalization scale factor (default 10000)
#' @param n_variable_features Number of HVGs (default 2000)
#' @param n_pcs Number of PCs (default 50)
#' @param clip_value Clip value for scaling (default 10.0, NA to skip)
#' @return A named list with embeddings, variance_explained, variable_features,
#'   and cell metadata
#' @export
sc_run_rna_pipeline <- function(ptr,
                                min_features = 200L,
                                min_cells = 3L,
                                max_pct_mt = 5.0,
                                scale_factor = 10000,
                                n_variable_features = 2000L,
                                n_pcs = 50L,
                                clip_value = 10.0) {
  res <- .Call(wrap__r_sc_run_rna_pipeline, ptr,
               as.integer(min_features), as.integer(min_cells),
               as.double(max_pct_mt), as.double(scale_factor),
               as.integer(n_variable_features), as.integer(n_pcs),
               clip_value)

  # Build embeddings matrix
  emb <- matrix(res$embeddings, nrow = res$emb_nrow, ncol = res$emb_ncol)
  colnames(emb) <- paste0("PC_", seq_len(res$emb_ncol))
  rownames(emb) <- res$cell_ids

  # Build cell metadata data.frame
  meta <- data.frame(
    cell_id = res$cell_ids,
    n_features = res$n_features,
    n_counts = res$n_counts,
    pct_mt = res$pct_mt,
    stringsAsFactors = FALSE
  )

  list(
    embeddings = emb,
    variance_explained = res$variance_explained,
    variable_features = res$variable_features,
    cell_metadata = meta,
    final_n_features = res$final_n_features,
    final_n_cells = res$final_n_cells
  )
}
