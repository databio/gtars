#!/usr/bin/env Rscript
# =========================================================================
# gtars-sc vs Seurat: step-by-step comparison on PBMC 3k
# =========================================================================
#
# Prerequisites:
#   install.packages("Seurat")
#   R CMD INSTALL path/to/gtars-r
#
# Usage:
#   Rscript tests/compare_sc_seurat.R /path/to/pbmc3k/filtered_gene_bc_matrices/hg19
#
# If no path is given, uses Seurat's built-in pbmc3k dataset (if SeuratData available)
# or downloads from 10X.
# =========================================================================

library(gtars)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

time_it <- function(expr) {
  t <- system.time(expr)
  t["elapsed"]
}

fmt_time <- function(secs) {
  if (secs < 1) {
    sprintf("%.0f ms", secs * 1000)
  } else {
    sprintf("%.2f s", secs)
  }
}

# -------------------------------------------------------------------------
# Resolve data path
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  data_dir <- args[1]
} else {
  # Try SeuratData
  if (requireNamespace("SeuratData", quietly = TRUE)) {
    SeuratData::InstallData("pbmc3k", force = FALSE)
    # SeuratData stores pre-processed; we need raw counts from the 10X files
    # Fall through to download if not already available
  }
  # Download from 10X
  data_dir <- file.path(tempdir(), "pbmc3k_10x")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
    url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
    tar_file <- file.path(tempdir(), "pbmc3k.tar.gz")
    cat("Downloading PBMC 3k dataset...\n")
    download.file(url, tar_file, quiet = TRUE)
    untar(tar_file, exdir = data_dir)
  }
  # 10X tar extracts to filtered_gene_bc_matrices/hg19/
  candidate <- file.path(data_dir, "filtered_gene_bc_matrices", "hg19")
  if (dir.exists(candidate)) {
    data_dir <- candidate
  }
}

cat("Data directory:", data_dir, "\n\n")

# =========================================================================
# Seurat pipeline
# =========================================================================

library(Seurat)

timings <- data.frame(
  step = character(),
  gtars_sec = numeric(),
  seurat_sec = numeric(),
  stringsAsFactors = FALSE
)

add_timing <- function(step, gtars_sec, seurat_sec) {
  timings[nrow(timings) + 1, ] <<- list(step, gtars_sec, seurat_sec)
}

# --- Read ---
cat("=== Step 1: Read 10X ===\n")
t_seurat <- time_it({
  seu_data <- Read10X(data.dir = data_dir)
  seu <- CreateSeuratObject(counts = seu_data, project = "pbmc3k")
})
t_gtars <- time_it({
  fm <- sc_read_10x(data_dir)
})
info <- sc_feature_matrix_info(fm)
cat(sprintf("  gtars:  %d genes x %d cells (%s nnz)\n", info$n_features, info$n_cells, format(info$nnz, big.mark = ",")))
cat(sprintf("  Seurat: %d genes x %d cells\n", nrow(seu), ncol(seu)))
add_timing("Read 10X", t_gtars, t_seurat)

# --- QC ---
cat("\n=== Step 2: QC Metrics ===\n")
t_seurat <- time_it({
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
})
t_gtars <- time_it({
  qc <- sc_compute_rna_qc(fm)
})
cat(sprintf("  gtars  median nFeatures: %.0f, median pct_mt: %.2f%%\n",
            median(qc$n_features), median(qc$pct_mt)))
cat(sprintf("  Seurat median nFeatures: %.0f, median pct_mt: %.2f%%\n",
            median(seu$nFeature_RNA), median(seu$percent.mt)))
add_timing("QC metrics", t_gtars, t_seurat)

# --- Filter genes ---
cat("\n=== Step 3: Filter Genes (min_cells=3) ===\n")
t_gtars <- time_it({
  fm_fg <- sc_filter_genes(fm, min_cells = 3L)
})
info_fg <- sc_feature_matrix_info(fm_fg)
cat(sprintf("  gtars:  %d -> %d genes\n", info$n_features, info_fg$n_features))
# Seurat doesn't separate gene filtering; it happens in CreateSeuratObject min.cells
# We re-create to match
t_seurat <- time_it({
  seu <- CreateSeuratObject(counts = seu_data, project = "pbmc3k", min.cells = 3)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
})
cat(sprintf("  Seurat: %d genes (min.cells=3)\n", nrow(seu)))
add_timing("Filter genes", t_gtars, t_seurat)

# --- Filter cells ---
cat("\n=== Step 4: Filter Cells (min_features=200, max_pct_mt=5) ===\n")
t_gtars <- time_it({
  fm_fc <- sc_filter_cells(fm_fg, min_features = 200L, max_pct_mt = 5.0)
})
info_fc <- sc_feature_matrix_info(fm_fc)
cat(sprintf("  gtars:  %d -> %d cells\n", info_fg$n_cells, info_fc$n_cells))

t_seurat <- time_it({
  seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 5)
})
cat(sprintf("  Seurat: %d cells\n", ncol(seu)))
add_timing("Filter cells", t_gtars, t_seurat)

# --- HVG (on raw counts, before normalization — matches Seurat VST) ---
cat("\n=== Step 5: Find Variable Features (n=2000, on counts) ===\n")
t_gtars <- time_it({
  hvg_gtars <- sc_find_variable_features(fm_fc, n_features = 2000L)
})
t_seurat <- time_it({
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
})
hvg_seurat <- VariableFeatures(seu)
overlap <- length(intersect(hvg_gtars, hvg_seurat))
cat(sprintf("  HVG overlap: %d / %d (%.1f%%)\n", overlap, 2000, overlap / 2000 * 100))
add_timing("Find HVGs", t_gtars, t_seurat)

# --- Normalize ---
cat("\n=== Step 6: Log Normalize ===\n")
t_gtars <- time_it({
  fm_norm <- sc_log_normalize(fm_fc, scale_factor = 10000)
})
t_seurat <- time_it({
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
})
add_timing("Normalize", t_gtars, t_seurat)

# --- Scale ---
cat("\n=== Step 7: Scale Data ===\n")
# Use gtars HVGs for gtars, Seurat HVGs for Seurat
t_gtars <- time_it({
  scaled_gtars <- sc_scale_data(fm_norm, hvg_gtars, clip_value = 10.0)
})
t_seurat <- time_it({
  seu <- ScaleData(seu, features = hvg_seurat)
})
cat(sprintf("  gtars  scaled matrix: %d x %d\n", nrow(scaled_gtars), ncol(scaled_gtars)))
cat(sprintf("  Seurat scaled matrix: %d x %d\n",
            nrow(GetAssayData(seu, layer = "scale.data")),
            ncol(GetAssayData(seu, layer = "scale.data"))))
add_timing("Scale", t_gtars, t_seurat)

# --- PCA ---
cat("\n=== Step 8: PCA (50 PCs) ===\n")
t_gtars <- time_it({
  pca_gtars <- sc_run_pca(scaled_gtars, n_pcs = 50L)
})
t_seurat <- time_it({
  seu <- RunPCA(seu, features = hvg_seurat, npcs = 50, verbose = FALSE)
})
cat(sprintf("  gtars  embeddings: %d x %d\n", nrow(pca_gtars$embeddings), ncol(pca_gtars$embeddings)))

# Compare variance explained
seu_var <- Stdev(seu, reduction = "pca")^2
seu_var_pct <- seu_var / sum(seu_var)
cat(sprintf("  Variance explained (PC1): gtars=%.4f, Seurat=%.4f\n",
            pca_gtars$variance_explained[1], seu_var_pct[1]))
add_timing("PCA", t_gtars, t_seurat)

# --- Embedding correlation ---
cat("\n=== Embedding Comparison ===\n")
# Use shared cells for comparison
shared_cells <- intersect(rownames(pca_gtars$embeddings), colnames(seu))
if (length(shared_cells) > 0) {
  gtars_emb <- pca_gtars$embeddings[shared_cells, , drop = FALSE]
  seurat_emb <- Embeddings(seu, "pca")[shared_cells, , drop = FALSE]
  n_compare <- min(ncol(gtars_emb), ncol(seurat_emb))
  cors <- sapply(seq_len(n_compare), function(i) {
    # PCA sign is arbitrary; compare absolute correlation
    abs(cor(gtars_emb[, i], seurat_emb[, i]))
  })
  cat(sprintf("  Shared cells: %d\n", length(shared_cells)))
  cat(sprintf("  Mean |cor| across %d PCs: %.4f\n", n_compare, mean(cors)))
  cat(sprintf("  PC1 |cor|: %.4f, PC2 |cor|: %.4f, PC3 |cor|: %.4f\n",
              cors[1], cors[2], cors[3]))
} else {
  cat("  No shared cells found (cell ID mismatch)\n")
}

# =========================================================================
# Timing summary
# =========================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat(sprintf("%-20s %12s %12s %10s\n", "Step", "gtars", "Seurat", "Speedup"))
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(timings))) {
  speedup <- timings$seurat_sec[i] / max(timings$gtars_sec[i], 0.001)
  cat(sprintf("%-20s %12s %12s %9.1fx\n",
              timings$step[i],
              fmt_time(timings$gtars_sec[i]),
              fmt_time(timings$seurat_sec[i]),
              speedup))
}
cat(strrep("-", 60), "\n")
total_gtars <- sum(timings$gtars_sec)
total_seurat <- sum(timings$seurat_sec)
cat(sprintf("%-20s %12s %12s %9.1fx\n",
            "TOTAL", fmt_time(total_gtars), fmt_time(total_seurat),
            total_seurat / max(total_gtars, 0.001)))
cat(strrep("=", 60), "\n")
