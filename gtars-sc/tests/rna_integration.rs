use gtars_sc::io::read_10x;
use gtars_sc::rna::pipeline::run_rna_preprocessing;
use gtars_sc::types::{EmbeddingMethod, RnaPipelineConfig};

use std::path::PathBuf;

fn synthetic_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/rna_integration")
}

fn small_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/small_10x")
}

#[test]
fn test_full_pipeline_small() {
    let fm = read_10x(&small_data_dir()).unwrap();
    assert_eq!(fm.n_features(), 5);
    assert_eq!(fm.n_cells(), 3);

    let config = RnaPipelineConfig {
        min_features: 1,
        min_cells: 1,
        max_pct_mt: 100.0,
        n_variable_features: 3,
        n_pcs: 2,
        ..Default::default()
    };

    let result = run_rna_preprocessing(fm, &config).unwrap();

    // Embedding should be cells × n_pcs
    assert_eq!(result.embedding.values.nrows(), result.cell_metadata.len());
    assert_eq!(result.embedding.values.ncols(), 2);
    assert_eq!(result.embedding.method, EmbeddingMethod::Pca);

    // Variable features selected
    assert_eq!(result.variable_features.len(), 3);

    // Variance explained should be positive and sum ≤ 1
    for &ve in &result.embedding.variance_explained {
        assert!(ve >= 0.0);
    }
    let total_ve: f64 = result.embedding.variance_explained.iter().sum();
    assert!(total_ve <= 1.0 + 1e-10);

    // Cell metadata should have pct_mt
    for meta in &result.cell_metadata {
        assert!(meta.custom.contains_key("pct_mt"));
    }
}

#[test]
fn test_full_pipeline_synthetic() {
    // Generate a synthetic 50×20 matrix via MTX format
    let dir = synthetic_data_dir();
    if !dir.join("matrix.mtx").exists() {
        create_synthetic_fixture(&dir);
    }

    let fm = read_10x(&dir).unwrap();
    assert_eq!(fm.n_features(), 50);
    assert_eq!(fm.n_cells(), 20);

    let config = RnaPipelineConfig {
        min_features: 1,
        min_cells: 1,
        max_pct_mt: 100.0,
        n_variable_features: 20,
        n_pcs: 10,
        ..Default::default()
    };

    let result = run_rna_preprocessing(fm, &config).unwrap();

    assert_eq!(result.embedding.values.ncols(), 10);
    assert!(result.embedding.values.nrows() <= 20);
    assert!(result.embedding.values.nrows() > 0);
}

/// Create a synthetic 50-gene × 20-cell fixture in MTX format.
fn create_synthetic_fixture(dir: &std::path::Path) {
    use std::fs;

    fs::create_dir_all(dir).unwrap();

    let n_genes = 50;
    let n_cells = 20;

    // Features file
    let mut features = String::new();
    for i in 0..n_genes {
        let name = if i < 2 {
            format!("MT-GENE{}", i + 1)
        } else {
            format!("GENE{}", i + 1)
        };
        features.push_str(&format!("ENSG{:05}\t{}\tGene Expression\n", i + 1, name));
    }
    fs::write(dir.join("features.tsv"), &features).unwrap();

    // Barcodes file
    let mut barcodes = String::new();
    for i in 0..n_cells {
        barcodes.push_str(&format!("CELL{:04}-1\n", i + 1));
    }
    fs::write(dir.join("barcodes.tsv"), &barcodes).unwrap();

    // Matrix file - create sparse entries with some structure
    // Use a simple deterministic pattern: gene i, cell j → value if (i+j) % 3 != 0
    let mut entries = Vec::new();
    for gene in 0..n_genes {
        for cell in 0..n_cells {
            if (gene + cell) % 3 != 0 {
                let val = ((gene * 7 + cell * 13) % 50 + 1) as f64;
                entries.push((gene + 1, cell + 1, val)); // 1-indexed
            }
        }
    }

    let mut mtx = String::new();
    mtx.push_str("%%MatrixMarket matrix coordinate real general\n");
    mtx.push_str(&format!("{} {} {}\n", n_genes, n_cells, entries.len()));
    for (r, c, v) in &entries {
        mtx.push_str(&format!("{} {} {}\n", r, c, v));
    }
    fs::write(dir.join("matrix.mtx"), &mtx).unwrap();
}
