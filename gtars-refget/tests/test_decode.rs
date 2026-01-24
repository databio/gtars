//! Integration tests for decode() workflow through RefgetStore
//!
//! These tests validate the end-to-end workflow of importing FASTA files,
//! storing sequences, retrieving them, and decoding them through the public API.
//!
//! Note: Unit tests for the decode() method itself (including edge cases,
//! different alphabets, encoding schemes) are in src/collection.rs

use gtars_refget::store::RefgetStore;
use std::io::Write;
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
        .add_sequence_collection_from_fasta(fasta_path)
        .expect("Failed to import FASTA");

    // Get the collection digest
    let collections = store.list_collections();
    assert!(!collections.is_empty(), "No collections found");
    let collection_digest = collections[0].digest.clone();

    // Expected sequences
    let expected = vec![
        ("seq1", "ACGTACGTACGT"),
        ("seq2", "TTGGCCAA"),
        ("seq3", "NNNNAAAA"),  // Tests IUPAC 'N' handling
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
        .add_sequence_collection_from_fasta(fasta_path)
        .expect("Failed to import FASTA");

    // Get the collection digest
    let collections = store.list_collections();
    assert!(!collections.is_empty(), "No collections found");
    let collection_digest = collections[0].digest.clone();

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
