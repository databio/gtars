//! Integration tests for decode() workflow through RefgetStore
//!
//! These tests validate the end-to-end workflow of importing FASTA files,
//! storing sequences, retrieving them, and decoding them through the public API.
//!
//! Note: Unit tests for the decode() method itself (including edge cases,
//! different alphabets, encoding schemes) are in src/collection.rs

use gtars_refget::digest::{
    AlphabetType, SequenceMetadata, SequenceRecord, SequenceEncoder,
};
use gtars_refget::store::{FastaImportOptions, RefgetStore};
use std::io::Write;
use std::sync::Arc;
use tempfile::NamedTempFile;

/// Helper function to create a temporary FASTA file with test sequences
fn create_test_fasta() -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    writeln!(file, ">seq1").expect("Failed to write");
    writeln!(file, "ACGTACGTACGT").expect("Failed to write");
    writeln!(file, ">seq2").expect("Failed to write");
    writeln!(file, "TTGGCCAA").expect("Failed to write");
    writeln!(file, ">seq3").expect("Failed to write");
    writeln!(file, "NNNNAAAA").expect("Failed to write");
    file
}

#[test]
fn test_decode_workflow_encoded() {
    // Test the complete workflow: import FASTA → store (encoded) → retrieve → decode
    let fasta_file = create_test_fasta();
    let fasta_path = fasta_file.path();

    // Create an encoded (bit-packed) store (default mode)
    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(fasta_path, FastaImportOptions::new())
        .expect("Failed to import FASTA");

    // Get the collection digest
    let collections = store.list_collections(0, usize::MAX, &[]).unwrap();
    assert!(!collections.results.is_empty(), "No collections found");
    let collection_digest = collections.results[0].digest.clone();

    // Expected sequences
    let expected = vec![
        ("seq1", "ACGTACGTACGT"),
        ("seq2", "TTGGCCAA"),
        ("seq3", "NNNNAAAA"), // Tests IUPAC 'N' handling
    ];

    // Verify all sequences can be retrieved and decoded correctly
    for (name, expected_seq) in expected.iter() {
        let record = store
            .get_sequence_by_name(&collection_digest, name)
            .unwrap_or_else(|e| panic!("Failed to retrieve {}: {}", name, e));

        assert!(record.is_loaded(), "Record should have data");

        let decoded = record
            .decode()
            .unwrap_or_else(|| panic!("decode() returned None for {}", name));

        assert_eq!(
            &decoded, expected_seq,
            "Decoded sequence for {} doesn't match",
            name
        );
    }
}

#[test]
fn test_decode_workflow_raw() {
    // Test the same workflow with raw (uncompressed) storage mode
    let fasta_file = create_test_fasta();
    let fasta_path = fasta_file.path();

    // Create a raw store
    let mut store = RefgetStore::in_memory();
    store.disable_encoding(); // Switch to Raw mode
    store
        .add_sequence_collection_from_fasta(fasta_path, FastaImportOptions::new())
        .expect("Failed to import FASTA");

    // Get the collection digest
    let collections = store.list_collections(0, usize::MAX, &[]).unwrap();
    assert!(!collections.results.is_empty(), "No collections found");
    let collection_digest = collections.results[0].digest.clone();

    // Verify sequences decode correctly with raw storage
    let record = store
        .get_sequence_by_name(&collection_digest, "seq1")
        .expect("Failed to retrieve seq1");

    let decoded = record.decode().expect("decode() returned None");
    assert_eq!(
        decoded, "ACGTACGTACGT",
        "Decoded sequence from raw store doesn't match"
    );

    // Verify IUPAC sequences work with raw storage too
    let record2 = store
        .get_sequence_by_name(&collection_digest, "seq3")
        .expect("Failed to retrieve seq3");

    let decoded2 = record2.decode().expect("decode() returned None for seq3");
    assert_eq!(decoded2, "NNNNAAAA", "IUPAC sequence doesn't match");
}

// FIX 4 regression tests: raw-vs-encoded length heuristic.
// These tests directly construct Full records with bit-packed bytes to verify
// that decode() correctly identifies encoded vs raw storage for edge cases
// where old heuristic (data.len() == metadata.length) misclassified.

/// Helper: build a SequenceRecord::Full with bit-packed (encoded) bytes for the
/// given alphabet and raw sequence, using the same encoder as the import pipeline.
fn make_encoded_full_record(name: &str, raw_seq: &[u8], alphabet: AlphabetType) -> SequenceRecord {
    use gtars_refget::digest::algorithms::{md5, sha512t24u};
    let length = raw_seq.len();
    let mut encoder = SequenceEncoder::new(alphabet, length);
    encoder.update(raw_seq);
    let encoded = encoder.finalize();
    let metadata = SequenceMetadata {
        name: name.to_string(),
        description: None,
        length,
        sha512t24u: sha512t24u(raw_seq),
        md5: md5(raw_seq),
        alphabet,
        fai: None,
    };
    SequenceRecord::Full {
        metadata,
        sequence: Arc::new(encoded),
    }
}

#[test]
fn fix4_decode_dnaiupac_length1() {
    // DnaIupac uses 4 bps; for length=1 sequence "R":
    //   expected_encoded = (1 * 4).div_ceil(8) = 1
    //   data.len() = 1, metadata.length = 1
    // Old heuristic: data.len() == length → incorrectly returns raw
    // New heuristic: data.len() == length AND data.len() == expected_encoded → not raw → encoded path
    let rec = make_encoded_full_record("s1", b"R", AlphabetType::DnaIupac);
    let decoded = rec.decode().expect("decode() returned None for DnaIupac length-1");
    assert_eq!(decoded, "R", "DnaIupac length-1 'R' must decode correctly");
}

#[test]
fn fix4_decode_dna3bit_length1() {
    // Dna3bit uses 3 bps; for length=1 sequence "N":
    //   expected_encoded = (1 * 3).div_ceil(8) = 1
    //   data.len() = 1, metadata.length = 1
    // Same ambiguity — new heuristic should handle it via encoded path.
    let rec = make_encoded_full_record("s2", b"N", AlphabetType::Dna3bit);
    let decoded = rec.decode().expect("decode() returned None for Dna3bit length-1");
    assert_eq!(decoded, "N", "Dna3bit length-1 'N' must decode correctly");
}

#[test]
fn fix4_decode_dna2bit_length1() {
    // Dna2bit uses 2 bps; for length=1 sequence "A":
    //   expected_encoded = (1 * 2).div_ceil(8) = 1
    //   data.len() = 1, metadata.length = 1
    // But "A" raw is also 1 byte. The new heuristic detects this as encoded
    // (data.len() == expected_encoded → not raw → encoded path).
    let rec = make_encoded_full_record("s3", b"A", AlphabetType::Dna2bit);
    let decoded = rec.decode().expect("decode() returned None for Dna2bit length-1");
    assert_eq!(decoded, "A", "Dna2bit length-1 'A' must decode correctly");
}
