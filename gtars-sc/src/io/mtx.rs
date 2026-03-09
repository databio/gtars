use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use sprs::TriMat;

use crate::types::{FeatureMatrix, FeatureType};

/// Read a 10X Genomics-style directory into a FeatureMatrix.
///
/// Expects:
/// - `matrix.mtx` or `matrix.mtx.gz`
/// - `barcodes.tsv` or `barcodes.tsv.gz`
/// - `features.tsv` or `features.tsv.gz` (falls back to `genes.tsv[.gz]`)
pub fn read_10x(dir: &Path) -> Result<FeatureMatrix> {
    let matrix_path = find_file(dir, &["matrix.mtx.gz", "matrix.mtx"])?;
    let barcodes_path = find_file(dir, &["barcodes.tsv.gz", "barcodes.tsv"])?;
    let features_path =
        find_file(dir, &["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"])?;

    let (cell_ids, _) = read_tsv(&barcodes_path)?;
    let (feature_ids, feature_rows) = read_tsv(&features_path)?;

    // Extract feature names (col 1) and feature type (col 2 if present)
    let feature_names: Vec<String> = feature_rows
        .iter()
        .map(|cols| cols.get(1).cloned().unwrap_or_default())
        .collect();

    let feature_type = detect_feature_type(&feature_rows);

    let matrix = read_mtx(&matrix_path, feature_ids.len(), cell_ids.len())?;

    Ok(FeatureMatrix {
        matrix,
        feature_names,
        feature_ids,
        cell_ids,
        feature_type,
    })
}

/// Find the first existing file from a list of candidates in a directory.
fn find_file(dir: &Path, candidates: &[&str]) -> Result<std::path::PathBuf> {
    for name in candidates {
        let path = dir.join(name);
        if path.exists() {
            return Ok(path);
        }
    }
    bail!(
        "None of {:?} found in {}",
        candidates,
        dir.display()
    )
}

/// Open a file, decompressing with gzip if the path ends in `.gz`.
fn open_maybe_gzipped(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().is_some_and(|ext| ext == "gz") {
        Ok(Box::new(GzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// Read a TSV file, returning (first_column, all_columns_per_row).
fn read_tsv(path: &Path) -> Result<(Vec<String>, Vec<Vec<String>>)> {
    let reader = BufReader::new(open_maybe_gzipped(path)?);
    let mut first_col = Vec::new();
    let mut all_cols = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let cols: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        first_col.push(cols[0].clone());
        all_cols.push(cols);
    }

    Ok((first_col, all_cols))
}

/// Detect FeatureType from the third column of the features file.
fn detect_feature_type(rows: &[Vec<String>]) -> FeatureType {
    let first_type = rows
        .iter()
        .find_map(|cols| cols.get(2))
        .map(|s| s.as_str());

    match first_type {
        Some("Gene Expression") => FeatureType::Gene,
        Some("Peaks") => FeatureType::Peak,
        Some(other) => FeatureType::Custom(other.to_string()),
        None => FeatureType::Gene, // legacy genes.tsv has only 2 columns
    }
}

/// Read a Matrix Market file into a CSC sparse matrix.
fn read_mtx(path: &Path, n_features: usize, n_cells: usize) -> Result<sprs::CsMat<f64>> {
    let reader = BufReader::new(open_maybe_gzipped(path)?);
    let mut lines = reader.lines();

    // Skip comment lines (start with %)
    let mut header_line = String::new();
    for line in lines.by_ref() {
        let line = line?;
        if !line.starts_with('%') {
            header_line = line;
            break;
        }
    }

    // Parse header: rows cols nnz
    let header: Vec<usize> = header_line
        .split_whitespace()
        .map(|s| s.parse::<usize>())
        .collect::<std::result::Result<Vec<_>, _>>()
        .context("parsing MTX header")?;

    if header.len() != 3 {
        bail!("expected 3 values in MTX header, got {}", header.len());
    }

    let (mtx_rows, mtx_cols, mtx_nnz) = (header[0], header[1], header[2]);

    // Validate dimensions against TSV files
    if mtx_rows != n_features {
        bail!(
            "MTX row count ({}) != features count ({})",
            mtx_rows,
            n_features
        );
    }
    if mtx_cols != n_cells {
        bail!(
            "MTX column count ({}) != barcodes count ({})",
            mtx_cols,
            n_cells
        );
    }

    let mut tri = TriMat::new((n_features, n_cells));
    tri.reserve(mtx_nnz);

    for line in lines {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() != 3 {
            bail!("expected 3 values per MTX entry, got {}", parts.len());
        }
        // MTX uses 1-based indexing
        let row: usize = parts[0].parse::<usize>().context("parsing row index")? - 1;
        let col: usize = parts[1].parse::<usize>().context("parsing col index")? - 1;
        let val: f64 = parts[2].parse().context("parsing value")?;
        tri.add_triplet(row, col, val);
    }

    Ok(tri.to_csc())
}

/// Write a FeatureMatrix to a 10X-style directory (gzipped).
///
/// Creates:
/// - `matrix.mtx.gz` — Matrix Market sparse matrix
/// - `barcodes.tsv.gz` — Cell IDs
/// - `features.tsv.gz` — Feature IDs, names, types
pub fn write_10x(dir: &Path, fm: &FeatureMatrix) -> Result<()> {
    fs::create_dir_all(dir)
        .with_context(|| format!("creating output directory {}", dir.display()))?;

    write_mtx_gz(&dir.join("matrix.mtx.gz"), fm)?;
    write_barcodes_gz(&dir.join("barcodes.tsv.gz"), &fm.cell_ids)?;
    write_features_gz(
        &dir.join("features.tsv.gz"),
        &fm.feature_ids,
        &fm.feature_names,
        &fm.feature_type,
    )?;

    Ok(())
}

fn write_mtx_gz(path: &Path, fm: &FeatureMatrix) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("creating {}", path.display()))?;
    let gz = GzEncoder::new(BufWriter::new(file), Compression::default());
    let mut w = BufWriter::new(gz);

    let (n_features, n_cells) = fm.shape();
    let nnz = fm.matrix.nnz();

    writeln!(w, "%%MatrixMarket matrix coordinate real general")?;
    writeln!(w, "{} {} {}", n_features, n_cells, nnz)?;

    // Iterate CSC matrix and output in 1-based indices
    let csc = &fm.matrix;
    for col in 0..n_cells {
        let outer = csc.outer_view(col);
        if let Some(col_view) = outer {
            for (row, &val) in col_view.iter() {
                writeln!(w, "{} {} {}", row + 1, col + 1, val)?;
            }
        }
    }

    w.flush()?;
    Ok(())
}

fn write_barcodes_gz(path: &Path, cell_ids: &[String]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("creating {}", path.display()))?;
    let gz = GzEncoder::new(BufWriter::new(file), Compression::default());
    let mut w = BufWriter::new(gz);

    for id in cell_ids {
        writeln!(w, "{}", id)?;
    }
    w.flush()?;
    Ok(())
}

fn write_features_gz(
    path: &Path,
    feature_ids: &[String],
    feature_names: &[String],
    feature_type: &FeatureType,
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("creating {}", path.display()))?;
    let gz = GzEncoder::new(BufWriter::new(file), Compression::default());
    let mut w = BufWriter::new(gz);

    let type_str = match feature_type {
        FeatureType::Gene => "Gene Expression",
        FeatureType::Peak => "Peaks",
        FeatureType::Custom(s) => s.as_str(),
    };

    for (id, name) in feature_ids.iter().zip(feature_names.iter()) {
        writeln!(w, "{}\t{}\t{}", id, name, type_str)?;
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/small_10x")
    }

    fn test_data_dir_gz() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/small_10x_gz")
    }

    #[test]
    fn test_read_10x_shape() {
        let fm = read_10x(&test_data_dir()).unwrap();
        assert_eq!(fm.n_features(), 5);
        assert_eq!(fm.n_cells(), 3);
        assert_eq!(fm.shape(), (5, 3));
    }

    #[test]
    fn test_read_10x_barcodes() {
        let fm = read_10x(&test_data_dir()).unwrap();
        assert_eq!(fm.cell_ids, vec!["AAAA-1", "BBBB-1", "CCCC-1"]);
    }

    #[test]
    fn test_read_10x_features() {
        let fm = read_10x(&test_data_dir()).unwrap();
        assert_eq!(fm.feature_names, vec!["Gene1", "Gene2", "MT-CO1", "Gene4", "Gene5"]);
        assert_eq!(fm.feature_type, FeatureType::Gene);
    }

    #[test]
    fn test_read_10x_values() {
        let fm = read_10x(&test_data_dir()).unwrap();
        assert_eq!(fm.matrix.nnz(), 8);
        // Check a specific value: gene 0 (Gene1), cell 0 → 5.0
        let val = fm.matrix.get(0, 0).copied().unwrap_or(0.0);
        assert!((val - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_read_10x_gzipped() {
        let fm = read_10x(&test_data_dir_gz()).unwrap();
        assert_eq!(fm.n_features(), 5);
        assert_eq!(fm.n_cells(), 3);
        assert_eq!(fm.matrix.nnz(), 8);
    }
}
