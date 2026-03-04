//! Smoketest against real PBMC 3k data (10X Genomics).
//!
//! Data must be downloaded first:
//!   cd tests/data/pbmc3k
//!   curl -sL https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz | tar -xz
//!
//! Run with: cargo test -p gtars-sc --test pbmc3k_smoketest

use std::path::PathBuf;

use gtars_sc::io::read_10x;
use gtars_sc::rna::hvg::find_variable_features;
use gtars_sc::rna::normalize::log_normalize;
use gtars_sc::rna::pipeline::run_rna_preprocessing;
use gtars_sc::rna::qc::{compute_rna_qc, filter_genes, filter_rna_cells};
use gtars_sc::rna::scale::scale_data;
use gtars_sc::types::{FeatureType, RnaPipelineConfig};

fn pbmc_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/pbmc3k/filtered_gene_bc_matrices/hg19")
}

fn has_pbmc_data() -> bool {
    pbmc_dir().join("matrix.mtx").exists()
}

#[test]
fn test_read_pbmc3k() {
    if !has_pbmc_data() {
        eprintln!("SKIP: PBMC 3k data not found at {}", pbmc_dir().display());
        return;
    }

    let fm = read_10x(&pbmc_dir()).unwrap();

    assert_eq!(fm.n_features(), 32738, "expected 32738 genes");
    assert_eq!(fm.n_cells(), 2700, "expected 2700 cells");
    assert_eq!(fm.feature_type, FeatureType::Gene);
    assert!(fm.matrix.nnz() > 2_000_000, "expected >2M nonzero entries");

    println!("PBMC 3k loaded: {} genes × {} cells, {} nnz",
        fm.n_features(), fm.n_cells(), fm.matrix.nnz());
}

#[test]
fn test_pbmc3k_qc() {
    if !has_pbmc_data() {
        eprintln!("SKIP: PBMC 3k data not found");
        return;
    }

    let fm = read_10x(&pbmc_dir()).unwrap();
    let qc = compute_rna_qc(&fm);

    assert_eq!(qc.len(), 2700);

    // Basic sanity: median features per cell should be ~500-2000
    let mut n_features: Vec<u32> = qc.iter().map(|m| m.n_features).collect();
    n_features.sort();
    let median_features = n_features[n_features.len() / 2];
    assert!(median_features > 200, "median features too low: {median_features}");
    assert!(median_features < 5000, "median features too high: {median_features}");

    // Some cells should have MT genes
    let mt_cells = qc.iter().filter(|m| m.pct_mt > 0.0).count();
    assert!(mt_cells > 0, "no cells with MT expression detected");

    println!("QC: median features={median_features}, cells with MT={mt_cells}/{}",
        qc.len());
}

#[test]
fn test_pbmc3k_filtering() {
    if !has_pbmc_data() {
        eprintln!("SKIP: PBMC 3k data not found");
        return;
    }

    let fm = read_10x(&pbmc_dir()).unwrap();

    // Filter genes: min 3 cells (Seurat default)
    let fm = filter_genes(&fm, 3).unwrap();
    println!("After gene filter (min_cells=3): {} genes", fm.n_features());
    assert!(fm.n_features() > 10_000, "too few genes after filter");
    assert!(fm.n_features() < 25_000, "too many genes after filter");

    // Filter cells: min 200 features, max 5% MT (Seurat defaults)
    let fm = filter_rna_cells(&fm, 200, 5.0).unwrap();
    println!("After cell filter (min_features=200, max_mt=5%): {} cells", fm.n_cells());
    assert!(fm.n_cells() > 2000, "too few cells after filter");
    assert!(fm.n_cells() < 2700, "no cells were filtered");
}

#[test]
fn test_pbmc3k_normalize_and_hvg() {
    if !has_pbmc_data() {
        eprintln!("SKIP: PBMC 3k data not found");
        return;
    }

    let fm = read_10x(&pbmc_dir()).unwrap();
    let fm = filter_genes(&fm, 3).unwrap();
    let fm = filter_rna_cells(&fm, 200, 5.0).unwrap();

    let mut fm = fm;
    log_normalize(&mut fm, 10_000.0);

    // All nonzero values should now be log-transformed (positive)
    let data = fm.matrix.data();
    for &v in data {
        assert!(v >= 0.0, "negative value after log normalize: {v}");
    }

    // Find 2000 HVGs
    let hvgs = find_variable_features(&fm, 2000).unwrap();
    assert_eq!(hvgs.len(), 2000);
    println!("Top 10 HVGs: {:?}", &hvgs[..10]);

    // HVGs should contain known highly-variable genes in PBMC data
    // (S100A8, S100A9, LYZ, CST3 are commonly in top HVGs)
    let known_hvgs = ["S100A8", "S100A9", "LYZ", "CST3", "NKG7"];
    let found: Vec<&str> = known_hvgs
        .iter()
        .filter(|&&g| hvgs.contains(&g.to_string()))
        .copied()
        .collect();
    println!("Known HVGs found in top 2000: {:?} ({}/{})", found, found.len(), known_hvgs.len());
    // We expect at least 3 of these 5 to be in the HVG list
    assert!(
        found.len() >= 3,
        "too few known HVGs found: {:?} (expected >=3 of {:?})",
        found,
        known_hvgs
    );
}

#[test]
fn test_pbmc3k_scale() {
    if !has_pbmc_data() {
        eprintln!("SKIP: PBMC 3k data not found");
        return;
    }

    let fm = read_10x(&pbmc_dir()).unwrap();
    let fm = filter_genes(&fm, 3).unwrap();
    let fm = filter_rna_cells(&fm, 200, 5.0).unwrap();

    let mut fm = fm;
    log_normalize(&mut fm, 10_000.0);

    let hvgs = find_variable_features(&fm, 2000).unwrap();
    let scaled = scale_data(&fm, &hvgs, Some(10.0), None).unwrap();

    let (n_feats, n_cells) = scaled.dim();
    assert_eq!(n_feats, 2000, "expected 2000 HVG rows");
    assert!(n_cells > 2000, "expected >2000 cells");
    println!("Scaled matrix: {} features × {} cells", n_feats, n_cells);

    // Check centering: each row mean ≈ 0
    for i in 0..n_feats.min(10) {
        let row_mean: f64 = scaled.row(i).sum() / n_cells as f64;
        assert!(row_mean.abs() < 1e-10, "row {i} mean = {row_mean}");
    }

    // Check clipping: all values in [-10, 10]
    for &v in scaled.iter() {
        assert!(v >= -10.0 - 1e-10 && v <= 10.0 + 1e-10, "value {v} outside clip range");
    }
}

#[test]
fn test_pbmc3k_full_pipeline() {
    if !has_pbmc_data() {
        eprintln!("SKIP: PBMC 3k data not found");
        return;
    }

    let fm = read_10x(&pbmc_dir()).unwrap();

    let config = RnaPipelineConfig {
        min_features: 200,
        min_cells: 3,
        max_pct_mt: 5.0,
        scale_factor: 10_000.0,
        n_variable_features: 2000,
        n_pcs: 50,
        clip_value: Some(10.0),
    };

    let result = run_rna_preprocessing(fm, &config).unwrap();

    // Cell count: should be ~2500-2700 after filtering
    let n_cells = result.embedding.values.nrows();
    println!("Pipeline result: {} cells, {} PCs", n_cells, result.embedding.values.ncols());
    assert!(n_cells > 2000, "too few cells: {n_cells}");
    assert!(n_cells <= 2700, "more cells than input: {n_cells}");

    // PCA shape
    assert_eq!(result.embedding.values.ncols(), 50, "expected 50 PCs");

    // Variance explained: first PC should capture most variance
    let ve = &result.embedding.variance_explained;
    assert_eq!(ve.len(), 50);
    assert!(ve[0] > ve[1], "PC1 should explain more variance than PC2");
    assert!(ve[0] > 0.01, "PC1 variance explained too low: {}", ve[0]);
    let total_ve: f64 = ve.iter().sum();
    println!("Variance explained: PC1={:.3}, PC2={:.3}, total={:.3}",
        ve[0], ve[1], total_ve);

    // HVGs
    assert_eq!(result.variable_features.len(), 2000);

    // Cell metadata
    assert_eq!(result.cell_metadata.len(), n_cells);
    for meta in &result.cell_metadata {
        let pct_mt = meta.custom.get("pct_mt").copied().unwrap_or(0.0);
        assert!(pct_mt <= 5.0, "cell {} has pct_mt={pct_mt}", meta.cell_id);
        assert!(meta.n_features >= 200, "cell {} has n_features={}", meta.cell_id, meta.n_features);
    }

    println!("PBMC 3k pipeline smoketest PASSED");
}
