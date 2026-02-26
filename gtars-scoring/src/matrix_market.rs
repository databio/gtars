use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::Result;
use flate2::Compression;
use flate2::write::GzEncoder;

/// Write sparse barcode counts directly to Matrix Market format
///
/// Takes the sparse HashMap from barcode_scoring_from_fragments() and writes
/// to 3 files without any intermediate dense matrix structure:
/// - {prefix}_matrix.mtx.gz: sparse triplets (row, col, value)
/// - {prefix}_barcodes.tsv.gz: cell IDs (one per line)
/// - {prefix}_features.tsv.gz: peak IDs (one per line)
///
/// IMPORTANT:
/// - Triplets are sorted by (row, col) for Matrix Market format compliance
/// - Never materializes a dense matrix - writes directly from sparse HashMap
/// - Memory efficient: only stores non-zero triplets
///
/// # Arguments
/// * `barcode_counts` - Sparse HashMap from barcode_scoring_from_fragments()
/// * `num_peaks` - Total number of peaks (from consensus.len())
/// * `output_prefix` - Prefix for output files
pub fn write_sparse_counts_to_mtx(
    barcode_counts: &HashMap<String, HashMap<usize, u32>>,
    num_peaks: usize,
    output_prefix: &str,
) -> Result<()> {
    // Sort barcodes for deterministic output
    let mut barcodes: Vec<String> = barcode_counts.keys().cloned().collect();
    barcodes.sort();

    // Collect all triplets (row, col, value)
    let mut triplets: Vec<(usize, usize, u32)> = Vec::new();
    for (row_idx, barcode) in barcodes.iter().enumerate() {
        if let Some(counts) = barcode_counts.get(barcode) {
            for (&peak_id, &count) in counts.iter() {
                triplets.push((row_idx, peak_id, count));
            }
        }
    }

    // IMPORTANT: Sort triplets by (row, col) for Matrix Market format compliance!
    // This ensures efficient loading by scipy and other sparse matrix readers.
    triplets.sort_by_key(|&(r, c, _)| (r, c));

    // 1. Write matrix.mtx.gz (sparse triplets)
    let mtx_path = format!("{}_matrix.mtx.gz", output_prefix);
    let mtx_file = File::create(&mtx_path)?;
    let mut mtx_writer = BufWriter::new(GzEncoder::new(mtx_file, Compression::default()));

    writeln!(
        mtx_writer,
        "%%MatrixMarket matrix coordinate integer general"
    )?;
    writeln!(
        mtx_writer,
        "{} {} {}",
        barcodes.len(),
        num_peaks,
        triplets.len()
    )?;

    // Write sorted triplets (1-indexed as per Matrix Market standard!)
    for (row_idx, col_idx, value) in triplets {
        writeln!(mtx_writer, "{} {} {}", row_idx + 1, col_idx + 1, value)?;
    }

    mtx_writer.flush()?;

    // 2. Write barcodes.tsv.gz (cell IDs)
    let barcodes_path = format!("{}_barcodes.tsv.gz", output_prefix);
    let barcodes_file = File::create(&barcodes_path)?;
    let mut barcodes_writer = BufWriter::new(GzEncoder::new(barcodes_file, Compression::default()));
    for barcode in &barcodes {
        writeln!(barcodes_writer, "{}", barcode)?;
    }
    barcodes_writer.flush()?;

    // 3. Write features.tsv.gz (peak IDs)
    let features_path = format!("{}_features.tsv.gz", output_prefix);
    let features_file = File::create(&features_path)?;
    let mut features_writer = BufWriter::new(GzEncoder::new(features_file, Compression::default()));
    for i in 0..num_peaks {
        writeln!(features_writer, "peak_{}", i)?;
    }
    features_writer.flush()?;

    Ok(())
}
