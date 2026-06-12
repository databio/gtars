use gtars_uniwig::bamqc::{compute_bam_qc, write_bam_qc_tsv, BamQcResult};
use std::io::Cursor;
use std::path::PathBuf;

fn test_data_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("tests/data")
}

#[test]
fn test_bamqc_chr22_small() {
    let bam_path = test_data_path().join("test_chr22_small.bam");
    if !bam_path.exists() {
        eprintln!("Skipping test - BAM file not found: {:?}", bam_path);
        return;
    }

    let result = compute_bam_qc(&bam_path).expect("Failed to compute BAM QC");

    assert!(result.total_reads > 0, "Expected reads in BAM file");
    assert!(
        result.nrf >= 0.0 && result.nrf <= 1.0,
        "NRF should be between 0 and 1, got {}",
        result.nrf
    );
    assert!(
        result.pbc1 >= 0.0 && result.pbc1 <= 1.0,
        "PBC1 should be between 0 and 1, got {}",
        result.pbc1
    );
    assert!(result.pbc2 >= 0.0, "PBC2 should be non-negative");

    assert!(
        result.distinct <= result.total_reads,
        "Distinct pairs should be <= total pairs"
    );
    assert!(
        result.m1 <= result.distinct,
        "M1 should be <= distinct pairs"
    );
}

#[test]
fn test_bamqc_tsv_output_format() {
    let result = BamQcResult {
        total_reads: 100,
        distinct: 90,
        m1: 80,
        m2: 5,
        dups: 10,
        mito_reads: 2,
        nrf: 0.8,
        pbc1: 0.888889,
        pbc2: 16.0,
    };

    let mut buffer = Cursor::new(Vec::new());
    write_bam_qc_tsv(&result, &mut buffer).expect("Failed to write TSV");

    let output = String::from_utf8(buffer.into_inner()).expect("Invalid UTF-8");
    let lines: Vec<&str> = output.lines().collect();

    assert_eq!(lines.len(), 2, "Expected header + data row");

    let headers: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(headers.len(), 10, "Expected 10 columns");
    assert_eq!(headers[0], "Total_read_pairs");
    assert_eq!(headers[7], "NRF");
    assert_eq!(headers[8], "PBC1");
    assert_eq!(headers[9], "PBC2");

    let values: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(values.len(), 10, "Expected 10 values");
    assert_eq!(values[0], "100");
    assert_eq!(values[1], "90");
    assert_eq!(values[2], "80");
    assert_eq!(values[3], "5");
}

#[test]
fn test_bamqc_metrics_calculation() {
    let result = BamQcResult {
        total_reads: 100,
        distinct: 100,
        m1: 100,
        m2: 0,
        dups: 0,
        mito_reads: 0,
        nrf: 1.0,
        pbc1: 1.0,
        pbc2: 100.0,
    };

    assert_eq!(result.dup_rate(), 0.0);
    assert_eq!(result.mito_rate(), 0.0);
}

#[test]
fn test_bamqc_with_mito() {
    let result = BamQcResult {
        total_reads: 100,
        distinct: 90,
        m1: 80,
        m2: 5,
        dups: 10,
        mito_reads: 20,
        nrf: 0.8,
        pbc1: 0.888889,
        pbc2: 16.0,
    };

    assert!((result.mito_rate() - 0.2).abs() < 0.001);
    assert!((result.dup_rate() - 0.1).abs() < 0.001);
}
