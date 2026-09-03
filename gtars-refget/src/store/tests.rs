//! Store integration tests -- cross-cutting concerns (import, export, persistence, lazy-loading).

use super::*;

use crate::collection::SequenceCollectionExt;
use crate::digest::types::{
    SequenceCollection, SequenceCollectionMetadata, SequenceMetadata, SequenceRecord,
};
use crate::digest::{AlphabetType, md5, sha512t24u};
use crate::hashkeyable::{DigestKey, HashKeyable};
use std::fs;
use std::path::PathBuf;
use tempfile::tempdir;

// =========================================================================
// Test helpers
// =========================================================================

/// Copy a test FASTA to a temp directory to avoid writing RGSI cache files
/// into the test data directory.
fn copy_test_fasta(temp_dir: &std::path::Path, name: &str) -> std::path::PathBuf {
    let src = format!("../tests/data/fasta/{}", name);
    let dst = temp_dir.join(name);
    std::fs::copy(&src, &dst)
        .unwrap_or_else(|e| panic!("Failed to copy {} to tempdir: {}", src, e));
    dst
}

fn calculate_test_digests(sequence: &[u8]) -> (String, String) {
    (sha512t24u(sequence), md5(sequence))
}

/// Helper: create an in-memory store with one collection from a FASTA string.
fn store_with_one_collection(fasta_content: &str) -> (RefgetStore, String) {
    let dir = tempdir().unwrap();
    let fasta = dir.path().join("test.fa");
    fs::write(&fasta, fasta_content).unwrap();

    let mut store = RefgetStore::in_memory();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    (store, meta.digest)
}

/// Creates a test store with 3 sequences for export testing.
fn setup_export_test_store(temp_path: &std::path::Path) -> (RefgetStore, DigestKey) {
    let fasta_content = ">chr1\nATGCATGCATGC\n>chr2\nGGGGAAAA\n>chr3\nTTTTCCCC\n";
    let temp_fasta_path = temp_path.join("test.fa");
    fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(&temp_fasta_path, FastaImportOptions::new())
        .unwrap();

    let collections: Vec<_> = store.collections.keys().cloned().collect();
    let collection_digest = collections[0];

    (store, collection_digest)
}

// =========================================================================
// Template and path tests (merged parametric)
// =========================================================================
// `test_expand_template`, `test_sanitize_relative_path`, and `test_mode_basics`
// are non-filesystem tests; they live in the always-compiled `nofs_tests`
// module in `store/mod.rs` so they run under `--no-default-features`.

// =========================================================================
// Mode tests
// =========================================================================

#[test]
fn test_mode_switching() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let temp_path = temp_dir.path();
    let fasta_content = ">chr1\nATGCATGCATGC\n>chr2\nGGGGAAAA\n";
    let temp_fasta_path = temp_path.join("test.fa");
    fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

    let (chr1_sha, _) = calculate_test_digests(b"ATGCATGCATGC");
    let chr1_key = chr1_sha.as_bytes().to_key();

    // Test Raw -> Encoded
    {
        let mut store = RefgetStore::in_memory();
        store.disable_encoding();
        store
            .add_sequence_collection_from_fasta(&temp_fasta_path, FastaImportOptions::new())
            .unwrap();

        if let Some(SequenceRecord::Full { sequence, .. }) = store.sequence_store.get(&chr1_key) {
            assert_eq!(&sequence[..], b"ATGCATGCATGC");
        }
        let seq_before = store.get_sequence(&chr1_sha).unwrap().decode().unwrap();

        store.set_encoding_mode(StorageMode::Encoded);

        if let Some(SequenceRecord::Full { sequence, .. }) = store.sequence_store.get(&chr1_key) {
            assert_eq!(sequence.len(), 3);
        }
        let seq_after = store.get_sequence(&chr1_sha).unwrap().decode().unwrap();
        assert_eq!(seq_before, seq_after);
    }

    // Test Encoded -> Raw
    {
        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta(&temp_fasta_path, FastaImportOptions::new())
            .unwrap();

        if let Some(SequenceRecord::Full { sequence, .. }) = store.sequence_store.get(&chr1_key) {
            assert_eq!(sequence.len(), 3);
        }
        let seq_before = store.get_sequence(&chr1_sha).unwrap().decode().unwrap();

        store.disable_encoding();

        if let Some(SequenceRecord::Full { sequence, .. }) = store.sequence_store.get(&chr1_key) {
            assert_eq!(&sequence[..], b"ATGCATGCATGC");
        }
        let seq_after = store.get_sequence(&chr1_sha).unwrap().decode().unwrap();
        assert_eq!(seq_before, seq_after);
    }
}

// =========================================================================
// Import and retrieval tests
// =========================================================================

#[test]
fn test_refget_store_retrieve_seq_and_vec() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let temp_path = temp_dir.path();

    let fasta_content = "\
>chr1
ATGCATGCATGC
>chr2
GGGGAAAA
";
    let temp_fasta_path = temp_path.join("test.fa");
    fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(&temp_fasta_path, FastaImportOptions::new())
        .unwrap();

    let collection_digest_ref: &str = "uC_UorBNf3YUu1YIDainBhI94CedlNeH";

    // Test BED-based export
    let bed_content = "\
chr1\t0\t5
chr1\t8\t12
chr2\t0\t4
";
    let temp_bed_path = temp_path.join("test.bed");
    fs::write(&temp_bed_path, bed_content).expect("Failed to write test BED file");

    let temp_output_fa_path = temp_path.join("output.fa");

    store
        .export_fasta_from_regions(
            collection_digest_ref,
            temp_bed_path.to_str().unwrap(),
            temp_output_fa_path.to_str().unwrap(),
        )
        .expect("export_fasta_from_regions failed");

    let output_fa_content =
        fs::read_to_string(&temp_output_fa_path).expect("Failed to read output FASTA file");

    let expected_fa_content = ">chr1:0-5\nATGCA\n>chr1:8-12\nATGC\n>chr2:0-4\nGGGG\n";
    assert_eq!(
        output_fa_content.trim(),
        expected_fa_content.trim(),
        "Output FASTA file content mismatch"
    );

    // Test substrings_from_regions iterator
    let vec_result: Vec<_> = store
        .substrings_from_regions(collection_digest_ref, temp_bed_path.to_str().unwrap())
        .expect("substrings_from_regions failed")
        .collect::<Result<Vec<_>, _>>()
        .expect("substrings_from_regions had errors");

    let expected_vec = vec![
        RetrievedSequence {
            sequence: "ATGCA".to_string(),
            chrom_name: "chr1".to_string(),
            start: 0,
            end: 5,
        },
        RetrievedSequence {
            sequence: "ATGC".to_string(),
            chrom_name: "chr1".to_string(),
            start: 8,
            end: 12,
        },
        RetrievedSequence {
            sequence: "GGGG".to_string(),
            chrom_name: "chr2".to_string(),
            start: 0,
            end: 4,
        },
    ];

    assert_eq!(vec_result, expected_vec, "Retrieved sequence vector mismatch");
}

#[test]
fn test_substrings_from_region_vectors() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let temp_path = temp_dir.path();

    let fasta_content = ">chr1\nATGCATGCATGC\n>chr2\nGGGGAAAA\n";
    let temp_fasta_path = temp_path.join("test.fa");
    fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(&temp_fasta_path, FastaImportOptions::new())
        .unwrap();

    let collection_digest_ref: &str = "uC_UorBNf3YUu1YIDainBhI94CedlNeH";

    // Parallel input vectors equivalent to the BED rows used in the file-based test.
    let chroms = ["chr1", "chr1", "chr2"];
    let starts = [0u32, 8, 0];
    let ends = [5u32, 12, 4];

    let vec_result: Vec<_> = store
        .substrings_from_region_vectors(collection_digest_ref, &chroms, &starts, &ends)
        .expect("substrings_from_region_vectors failed")
        .collect::<Result<Vec<_>, _>>()
        .expect("substrings_from_region_vectors had errors");

    let expected_vec = vec![
        RetrievedSequence {
            sequence: "ATGCA".to_string(),
            chrom_name: "chr1".to_string(),
            start: 0,
            end: 5,
        },
        RetrievedSequence {
            sequence: "ATGC".to_string(),
            chrom_name: "chr1".to_string(),
            start: 8,
            end: 12,
        },
        RetrievedSequence {
            sequence: "GGGG".to_string(),
            chrom_name: "chr2".to_string(),
            start: 0,
            end: 4,
        },
    ];

    assert_eq!(
        vec_result, expected_vec,
        "Vectors-based retrieval mismatch vs expected"
    );

    // Missing-sequence case: a chrom not in the collection returns an Err with the
    // "Region N" label.
    let missing_chroms = ["chrX"];
    let missing_starts = [0u32];
    let missing_ends = [4u32];
    let missing_result: Result<Vec<_>, _> = store
        .substrings_from_region_vectors(
            collection_digest_ref,
            &missing_chroms,
            &missing_starts,
            &missing_ends,
        )
        .expect("constructor should succeed for equal-length vectors")
        .collect();
    let missing_err = missing_result.unwrap_err().to_string();
    assert!(
        missing_err.contains("Region 1") && missing_err.contains("not found"),
        "unexpected missing-sequence error: {missing_err}"
    );
}

#[test]
fn test_substrings_from_region_vectors_mismatched_lengths() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let temp_path = temp_dir.path();

    let fasta_content = ">chr1\nATGCATGCATGC\n>chr2\nGGGGAAAA\n";
    let temp_fasta_path = temp_path.join("test.fa");
    fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(&temp_fasta_path, FastaImportOptions::new())
        .unwrap();

    let collection_digest_ref: &str = "uC_UorBNf3YUu1YIDainBhI94CedlNeH";

    let chroms = ["chr1", "chr2"];
    let starts = [0u32]; // intentionally shorter
    let ends = [5u32, 4];

    let result =
        store.substrings_from_region_vectors(collection_digest_ref, &chroms, &starts, &ends);
    let err_msg = match result {
        Ok(_) => panic!("mismatched lengths should error"),
        Err(e) => e.to_string(),
    };
    assert!(
        err_msg.contains("Mismatched region vector lengths"),
        "unexpected error: {err_msg}"
    );
}

#[test]
fn test_negative_bed_coordinates() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let temp_path = temp_dir.path();

    let fasta_content = ">chr1\nATGCATGCATGC\n>chr2\nGGGGAAAA\n";
    let temp_fasta_path = temp_path.join("test.fa");
    fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(&temp_fasta_path, FastaImportOptions::new())
        .unwrap();

    let collections: Vec<_> = store.collections.keys().cloned().collect();
    let collection_digest_ref: &str =
        std::str::from_utf8(&collections[0]).expect("invalid collection digest");

    // Negative start
    let bed_content = "chr1\t-5\t100\n";
    let temp_bed_path = temp_path.join("negative.bed");
    fs::write(&temp_bed_path, bed_content).expect("Failed to write test BED file");

    let result = store.substrings_from_regions(collection_digest_ref, temp_bed_path.to_str().unwrap());
    let collected: Result<Vec<_>, _> = result.unwrap().collect();
    let error_msg = collected.unwrap_err().to_string();
    assert!(error_msg.contains("invalid start or end coordinates"));

    // Negative end
    let bed_content_neg_end = "chr1\t0\t-10\n";
    let temp_bed_path_neg_end = temp_path.join("negative_end.bed");
    fs::write(&temp_bed_path_neg_end, bed_content_neg_end).expect("Failed to write test BED file");

    let result = store.substrings_from_regions(collection_digest_ref, temp_bed_path_neg_end.to_str().unwrap());
    let collected: Result<Vec<_>, _> = result.unwrap().collect();
    let error_msg = collected.unwrap_err().to_string();
    assert!(error_msg.contains("invalid start or end coordinates"));
}

#[test]
fn test_global_refget_store() {
    let sequence = b"ACGT";
    let name = "test_seq";

    let mut collection = SequenceCollection {
        metadata: SequenceCollectionMetadata {
            digest: "test_collection".to_string(),
            n_sequences: 0,
            names_digest: "test".to_string(),
            sequences_digest: "test".to_string(),
            lengths_digest: "test".to_string(),
            name_length_pairs_digest: None,
            sorted_name_length_pairs_digest: None,
            sorted_sequences_digest: None,
            file_path: None,
        },
        sequences: Vec::new(),
    };

    let seq_metadata = SequenceMetadata {
        name: name.to_string(),
        description: None,
        length: sequence.len(),
        sha512t24u: sha512t24u(sequence),
        md5: md5(sequence),
        alphabet: AlphabetType::Dna2bit,
        fai: None,
    };

    let record = SequenceRecord::Full {
        metadata: seq_metadata.clone(),
        sequence: std::sync::Arc::new(sequence.to_vec()),
    };

    collection.sequences.push(record);

    let mut store = RefgetStore::in_memory();
    store.add_sequence_collection(collection.clone()).unwrap();

    assert!(!store.sequence_store.is_empty());

    let retrieved_by_name_str = store.get_sequence_by_name(&collection.metadata.digest, name);
    assert!(retrieved_by_name_str.is_ok());
    let retrieved_record = retrieved_by_name_str.unwrap();
    assert_eq!(retrieved_record.metadata().name, name);
    assert_eq!(retrieved_record.sequence().unwrap(), sequence);

    let retrieved_by_name_key =
        store.get_sequence_by_name(collection.metadata.digest.to_key(), name);
    assert!(retrieved_by_name_key.is_ok());

    let retrieved_by_sha512_str = store.get_sequence(&seq_metadata.sha512t24u);
    assert!(retrieved_by_sha512_str.is_ok());

    let retrieved_by_sha512_key = store.get_sequence(seq_metadata.sha512t24u.to_key());
    assert!(retrieved_by_sha512_key.is_ok());
}

#[test]
fn test_get_substring_by_md5() {
    // get_substring must resolve an md5 digest via the md5_lookup map,
    // matching the behavior of get_substrings / get_sequence.
    let sequence = b"ACGTACGTAC";
    let name = "test_seq";

    let mut collection = SequenceCollection {
        metadata: SequenceCollectionMetadata {
            digest: "test_collection".to_string(),
            n_sequences: 0,
            names_digest: "test".to_string(),
            sequences_digest: "test".to_string(),
            lengths_digest: "test".to_string(),
            name_length_pairs_digest: None,
            sorted_name_length_pairs_digest: None,
            sorted_sequences_digest: None,
            file_path: None,
        },
        sequences: Vec::new(),
    };

    let seq_metadata = SequenceMetadata {
        name: name.to_string(),
        description: None,
        length: sequence.len(),
        sha512t24u: sha512t24u(sequence),
        md5: md5(sequence),
        alphabet: AlphabetType::Dna2bit,
        fai: None,
    };

    collection.sequences.push(SequenceRecord::Full {
        metadata: seq_metadata.clone(),
        sequence: std::sync::Arc::new(sequence.to_vec()),
    });

    let mut store = RefgetStore::in_memory();
    store.add_sequence_collection(collection).unwrap();

    // Baseline: resolution by sha512t24u digest.
    let by_sha = store
        .get_substring(&seq_metadata.sha512t24u, 2, 6)
        .expect("get_substring by sha512t24u should succeed");
    assert_eq!(by_sha.len(), 4);

    // The fix: resolution by md5 digest must yield identical bytes.
    let by_md5 = store
        .get_substring(&seq_metadata.md5, 2, 6)
        .expect("get_substring by md5 should succeed");
    assert_eq!(by_md5, by_sha);
}

#[test]
fn test_get_substring_zero_length_range() {
    // CONTRACT: a zero-length range (start == end) is a valid empty interval and
    // returns "" (matching wasm/JS RemoteRefgetStore.getSubstring), while an
    // inverted range (start > end) is still rejected. Verified on:
    //   1. a resident (in-memory `Full`) sequence, and
    //   2. a disk-backed `Stub` sequence (exercises get_substring_from_disk),
    // both via single-range get_substring and batch get_substrings, so all
    // retrieval paths agree on the contract.
    let temp_dir = tempdir().unwrap();
    let temp_path = temp_dir.path();
    let temp_fasta = copy_test_fasta(temp_path, "base.fa.gz");

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();

    // Pick a non-empty resident sequence.
    let (digest, length) = store
        .sequence_metadata()
        .find(|m| m.length > 1)
        .map(|m| (m.sha512t24u.clone(), m.length))
        .expect("expected at least one non-empty sequence");

    // --- Resident (Full) path ---
    // Zero-length range anywhere within bounds returns "".
    assert_eq!(
        store.get_substring(&digest, 0, 0).unwrap(),
        "",
        "resident: start == end == 0 should return empty string"
    );
    assert_eq!(
        store.get_substring(&digest, 3, 3).unwrap(),
        "",
        "resident: interior zero-length range should return empty string"
    );
    // Zero-length range at the sequence end is also valid/empty.
    assert_eq!(
        store.get_substring(&digest, length, length).unwrap(),
        "",
        "resident: zero-length range at end should return empty string"
    );
    // Inverted range (start > end) still errors.
    assert!(
        store.get_substring(&digest, 5, 2).is_err(),
        "resident: inverted range (start > end) must error"
    );
    // Batch path agrees: empty ranges -> "", interleaved with a real range.
    let batch = store
        .get_substrings(&digest, &[(0, 0), (1, 4), (3, 3)])
        .unwrap();
    assert_eq!(batch.len(), 3);
    assert_eq!(batch[0], "");
    assert_eq!(batch[1].len(), 3);
    assert_eq!(batch[2], "");
    // Batch with an inverted range still errors.
    assert!(
        store.get_substrings(&digest, &[(5, 2)]).is_err(),
        "batch: inverted range must error"
    );

    // --- Disk-backed (Stub -> get_substring_from_disk) path ---
    store
        .write_store_to_dir(temp_path, Some("sequences/%s2/%s.seq"))
        .unwrap();
    let disk_store = RefgetStore::open_local(temp_path).unwrap();
    // Sequence is a Stub on disk (not resident), so this exercises the
    // partial-read path, which must also honor the zero-length contract.
    assert!(
        !disk_store.is_sequence_loaded(&digest),
        "disk-backed sequence should be a stub, not resident"
    );
    assert_eq!(
        disk_store.get_substring(&digest, 0, 0).unwrap(),
        "",
        "disk: start == end == 0 should return empty string"
    );
    assert_eq!(
        disk_store.get_substring(&digest, 4, 4).unwrap(),
        "",
        "disk: interior zero-length range should return empty string"
    );
    assert!(
        disk_store.get_substring(&digest, 5, 2).is_err(),
        "disk: inverted range (start > end) must error"
    );
    let disk_batch = disk_store
        .get_substrings(&digest, &[(2, 2), (0, 3)])
        .unwrap();
    assert_eq!(disk_batch[0], "");
    assert_eq!(disk_batch[1].len(), 3);
}

#[test]
fn test_import_fasta() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let temp_path = temp_dir.path();

    let test_fa = "../tests/data/fasta/base.fa";
    let temp_fa = temp_path.join("base.fa");
    std::fs::copy(test_fa, &temp_fa).expect("Failed to copy test FASTA file");

    let mut store = RefgetStore::in_memory();
    store.add_sequence_collection_from_fasta(temp_fa, FastaImportOptions::new()).unwrap();

    assert!(!store.sequence_store.is_empty());

    let seq_template = "sequences/%s2/%s.seq";
    store
        .write_store_to_dir(temp_path.to_str().unwrap(), Some(seq_template))
        .unwrap();
}

// =========================================================================
// Persistence tests
// =========================================================================

#[test]
fn test_disk_persistence() {
    let temp_dir = tempdir().unwrap();
    let temp_path = temp_dir.path();
    let temp_fasta = copy_test_fasta(temp_path, "base.fa.gz");

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();

    let sequence_keys: Vec<DigestKey> = store.sequence_store.keys().cloned().collect();
    assert_eq!(sequence_keys.len(), 3);

    let sha512_key1 = sequence_keys[0];
    let sha512_key2 = sequence_keys[1];

    let original_seq1 = store.sequence_store.get(&sha512_key1).unwrap().clone();
    let original_seq2 = store.sequence_store.get(&sha512_key2).unwrap().clone();

    let seq_template = "sequences/%s2/%s.seq";
    store.write_store_to_dir(temp_path, Some(seq_template)).unwrap();

    assert!(temp_path.join("sequences").exists());
    assert!(temp_path.join("sequences").read_dir().unwrap().count() > 0);
    assert!(temp_path.join("rgstore.json").exists());
    assert!(temp_path.join("sequences.rgsi").exists());
    assert!(temp_path.join("collections.rgci").exists());
    assert!(temp_path.join("collections").exists());

    let mut loaded_store = RefgetStore::open_local(temp_path).unwrap();

    assert_eq!(loaded_store.sequence_store.len(), 3);
    assert!(loaded_store.sequence_store.contains_key(&sha512_key1));
    assert!(loaded_store.sequence_store.contains_key(&sha512_key2));

    let loaded_seq1 = loaded_store.sequence_store.get(&sha512_key1).unwrap();
    let loaded_seq2 = loaded_store.sequence_store.get(&sha512_key2).unwrap();

    assert_eq!(original_seq1.metadata().name, loaded_seq1.metadata().name);
    assert_eq!(original_seq1.metadata().length, loaded_seq1.metadata().length);
    assert_eq!(original_seq1.metadata().sha512t24u, loaded_seq1.metadata().sha512t24u);
    assert_eq!(original_seq1.metadata().md5, loaded_seq1.metadata().md5);

    assert_eq!(original_seq2.metadata().name, loaded_seq2.metadata().name);
    assert_eq!(original_seq2.metadata().length, loaded_seq2.metadata().length);
    assert_eq!(original_seq2.metadata().sha512t24u, loaded_seq2.metadata().sha512t24u);
    assert_eq!(original_seq2.metadata().md5, loaded_seq2.metadata().md5);

    assert!(!loaded_seq1.is_loaded());
    assert!(!loaded_seq2.is_loaded());

    assert_eq!(loaded_store.md5_lookup.len(), 3);
    assert_eq!(loaded_store.collections.len(), store.collections.len());

    loaded_store.load_all_collections().unwrap();
    loaded_store.load_all_sequences().unwrap();

    for (digest, original_record) in &store.sequence_store {
        let loaded_record = loaded_store.get_sequence(*digest).unwrap();
        assert_eq!(original_record.metadata().name, loaded_record.metadata().name);
        assert_eq!(original_record.metadata().length, loaded_record.metadata().length);

        if original_record.metadata().length > 0 {
            let substring_len = std::cmp::min(5, original_record.metadata().length);
            let substring = loaded_store.get_substring(digest, 0, substring_len);
            assert!(substring.is_ok());
        }
    }
}

// =========================================================================
// Export tests
// =========================================================================

#[test]
fn test_export_fasta_all_sequences() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let (store, collection_digest) = setup_export_test_store(temp_dir.path());

    let output_path = temp_dir.path().join("exported_all.fa");
    store
        .export_fasta(&collection_digest, &output_path, None, Some(80))
        .unwrap();

    let exported = fs::read_to_string(&output_path).unwrap();
    assert!(exported.contains(">chr1") && exported.contains(">chr2") && exported.contains(">chr3"));
    assert!(exported.contains("ATGCATGCATGC") && exported.contains("GGGGAAAA") && exported.contains("TTTTCCCC"));
}

#[test]
fn test_export_fasta_subset_sequences() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let (store, collection_digest) = setup_export_test_store(temp_dir.path());

    let output_path = temp_dir.path().join("exported_subset.fa");
    store
        .export_fasta(&collection_digest, &output_path, Some(vec!["chr1", "chr3"]), Some(80))
        .unwrap();

    let exported = fs::read_to_string(&output_path).unwrap();
    assert!(exported.contains(">chr1") && exported.contains(">chr3"));
    assert!(!exported.contains(">chr2") && !exported.contains("GGGGAAAA"));
}

#[test]
fn test_export_fasta_roundtrip() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let temp_path = temp_dir.path();

    let fasta_content = "\
>seq1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>seq2
GGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCC
";
    let temp_fasta_path = temp_path.join("original.fa");
    fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

    let mut store1 = RefgetStore::in_memory();
    store1
        .add_sequence_collection_from_fasta(&temp_fasta_path, FastaImportOptions::new())
        .unwrap();

    let original_digests: Vec<String> = store1
        .sequence_store
        .values()
        .map(|r| r.metadata().sha512t24u.clone())
        .collect();

    let collections: Vec<_> = store1.collections.keys().cloned().collect();
    let collection_digest = collections[0];
    let exported_path = temp_path.join("exported.fa");
    store1.export_fasta(&collection_digest, &exported_path, None, Some(60)).expect("Failed to export FASTA");

    let mut store2 = RefgetStore::in_memory();
    store2
        .add_sequence_collection_from_fasta(&exported_path, FastaImportOptions::new())
        .unwrap();

    let new_digests: Vec<String> = store2
        .sequence_store
        .values()
        .map(|r| r.metadata().sha512t24u.clone())
        .collect();

    assert_eq!(original_digests.len(), new_digests.len());
    for digest in original_digests {
        assert!(new_digests.contains(&digest), "Digest {} should be present after roundtrip", digest);
    }
}

#[test]
fn test_export_fasta_by_digests() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let (store, _) = setup_export_test_store(temp_dir.path());

    let digests: Vec<String> = store
        .sequence_store
        .values()
        .map(|r| r.metadata().sha512t24u.clone())
        .collect();
    let digest_refs: Vec<&str> = digests.iter().map(|s| s.as_str()).collect();

    let output_path = temp_dir.path().join("exported_by_digests.fa");
    store.export_fasta_by_digests(digest_refs, &output_path, Some(80)).unwrap();

    let exported = fs::read_to_string(&output_path).unwrap();
    assert!(exported.contains(">chr1") && exported.contains(">chr2") && exported.contains(">chr3"));
}

#[test]
fn test_export_fasta_error_handling() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let (store, collection_digest) = setup_export_test_store(temp_dir.path());

    let output_path = temp_dir.path().join("should_fail.fa");

    let fake_collection = b"fake_collection_digest_12345678";
    assert!(store.export_fasta(fake_collection, &output_path, None, Some(80)).is_err());

    assert!(store.export_fasta(&collection_digest, &output_path, Some(vec!["nonexistent_chr"]), Some(80)).is_err());
}

#[test]
fn test_export_fasta_after_load_local() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let temp_path = temp_dir.path();
    let store_path = temp_path.join("store");

    let fasta_content = ">chr1\nACGTACGT\n>chr2\nGGGGAAAA\n";
    let fasta_path = temp_path.join("test.fa");
    fs::write(&fasta_path, fasta_content).unwrap();

    let collection_digest: DigestKey;
    {
        let mut store = RefgetStore::on_disk(&store_path).unwrap();
        store
            .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
            .unwrap();
        let collections: Vec<_> = store.collections.keys().cloned().collect();
        assert_eq!(collections.len(), 1);
        collection_digest = collections[0];
    }

    let mut loaded_store = RefgetStore::open_local(&store_path).unwrap();

    assert!(!loaded_store.is_collection_loaded(&collection_digest));

    loaded_store.load_all_collections().unwrap();
    loaded_store.load_all_sequences().unwrap();

    let output_path = temp_path.join("exported.fa");
    loaded_store
        .export_fasta(&collection_digest, &output_path, None, Some(80))
        .expect("export_fasta should work on disk-loaded stores");

    let exported = fs::read_to_string(&output_path).unwrap();
    assert!(exported.contains(">chr1"));
    assert!(exported.contains("ACGTACGT"));
    assert!(exported.contains(">chr2"));
    assert!(exported.contains("GGGGAAAA"));
}

// =========================================================================
// FASTA header and filename tests
// =========================================================================

#[test]
fn test_sequence_names_with_spaces() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let temp_path = temp_dir.path();

    let fasta_content = "\
>JAHKSE010000016.1 unmasked:primary_assembly HG002.alt.pat.f1_v2:JAHKSE010000016.1:1:100:1
ATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGC
>JAHKSE010000012.1 unmasked:primary_assembly HG002.alt.pat.f1_v2:JAHKSE010000012.1:1:100:1
GGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGG
GGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGG
";
    let temp_fasta_path = temp_path.join("spaces_in_names.fa");
    fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(&temp_fasta_path, FastaImportOptions::new())
        .expect("Should parse FASTA headers correctly");

    assert_eq!(store.sequence_store.len(), 2);

    let name1 = "JAHKSE010000016.1";
    let name2 = "JAHKSE010000012.1";

    let collections: Vec<_> = store.collections.keys().cloned().collect();
    assert_eq!(collections.len(), 1);
    let collection_digest = collections[0];

    {
        let seq1 = store.get_sequence_by_name(&collection_digest, name1);
        assert!(seq1.is_ok());

        let seq1_meta = seq1.unwrap().metadata();
        assert_eq!(seq1_meta.name, "JAHKSE010000016.1");
        assert_eq!(
            seq1_meta.description,
            Some("unmasked:primary_assembly HG002.alt.pat.f1_v2:JAHKSE010000016.1:1:100:1".to_string())
        );
    }

    {
        let seq2 = store.get_sequence_by_name(&collection_digest, name2);
        assert!(seq2.is_ok());
    }
}

#[test]
fn test_rgsi_filename_with_dots() {
    let temp_dir = tempdir().expect("Failed to create temporary directory");
    let temp_path = temp_dir.path();

    let temp_fasta = copy_test_fasta(temp_path, "HG002.alt.pat.f1_v2.unmasked.fa");

    let _seqcol = SequenceCollection::from_path_with_cache(&temp_fasta, false, true)
        .expect("Should load FASTA");

    let correct_rgsi = temp_path.join("HG002.alt.pat.f1_v2.unmasked.rgsi");
    let wrong_rgsi = temp_path.join("HG002.rgsi");

    let files: Vec<_> = std::fs::read_dir(temp_path)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
        .collect();

    assert!(correct_rgsi.exists(), "Expected 'HG002.alt.pat.f1_v2.unmasked.rgsi' but found: {:?}", files);
    assert!(!wrong_rgsi.exists());
}

// =========================================================================
// On-disk incremental tests
// =========================================================================

#[test]
fn test_on_disk_collection_written_incrementally() {
    let temp_dir = tempdir().unwrap();
    let temp_fasta = copy_test_fasta(temp_dir.path(), "base.fa.gz");

    let cache_path = temp_dir.path().join("cache");
    let mut store = RefgetStore::on_disk(&cache_path).unwrap();

    store
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();

    let collections_dir = cache_path.join("collections");
    assert!(collections_dir.exists());

    let rgsi_files: Vec<_> = std::fs::read_dir(&collections_dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
        .collect();

    assert!(!rgsi_files.is_empty());
    assert!(rgsi_files.iter().any(|f| f.ends_with(".rgsi")));
}

#[test]
fn test_disk_size_calculation() {
    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz", FastaImportOptions::new())
        .unwrap();

    let disk_size = store.logical_sequence_bytes();
    assert!(disk_size > 0);

    let manual: usize = store
        .list_sequences()
        .iter()
        .map(|m| (m.length * m.alphabet.bits_per_symbol()).div_ceil(8))
        .sum();
    assert_eq!(disk_size, manual);
}

#[test]
fn test_incremental_index_writing() {
    let temp_dir = tempdir().unwrap();
    let temp_fasta = copy_test_fasta(temp_dir.path(), "base.fa.gz");
    let cache_path = temp_dir.path().join("store");
    let mut store = RefgetStore::on_disk(&cache_path).unwrap();

    store
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();

    assert!(cache_path.join("rgstore.json").exists());
    assert!(cache_path.join("sequences.rgsi").exists());
    assert!(cache_path.join("collections.rgci").exists());

    let _loaded = RefgetStore::on_disk(&cache_path).unwrap();
}

#[test]
fn test_write_method() {
    let temp_dir = tempdir().unwrap();
    let temp_fasta = copy_test_fasta(temp_dir.path(), "base.fa.gz");
    let cache_path = temp_dir.path().join("store");
    let mut store = RefgetStore::on_disk(&cache_path).unwrap();

    store
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();
    store.write().unwrap();

    assert!(cache_path.join("rgstore.json").exists());
}

#[test]
fn test_on_disk_smart_constructor() {
    let temp_dir = tempdir().unwrap();
    let temp_fasta = copy_test_fasta(temp_dir.path(), "base.fa.gz");
    let cache_path = temp_dir.path().join("store");

    let mut store1 = RefgetStore::on_disk(&cache_path).unwrap();
    assert_eq!(store1.mode, StorageMode::Encoded);
    store1
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();

    let store2 = RefgetStore::on_disk(&cache_path).unwrap();
    assert_eq!(store2.sequence_store.len(), store1.sequence_store.len());
    assert_eq!(store2.mode, StorageMode::Encoded);

    let cache_path_raw = temp_dir.path().join("store_raw");
    let mut store3 = RefgetStore::on_disk(&cache_path_raw).unwrap();
    store3.disable_encoding();
    assert_eq!(store3.mode, StorageMode::Raw);
    store3
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();

    let store4 = RefgetStore::on_disk(&cache_path_raw).unwrap();
    assert_eq!(store4.mode, StorageMode::Raw);

    let index_path = cache_path_raw.join("rgstore.json");
    let json = fs::read_to_string(&index_path).unwrap();
    assert!(json.contains("\"mode\":\"Raw\"") || json.contains("\"mode\": \"Raw\""));
}

// =========================================================================
// Collection metadata and loading tests
// =========================================================================

#[test]
fn test_collection_metadata_methods() {
    let temp_dir = tempdir().unwrap();
    let temp_fasta = copy_test_fasta(temp_dir.path(), "base.fa.gz");
    let cache_path = temp_dir.path().join("store");
    let mut store = RefgetStore::on_disk(&cache_path).unwrap();

    store
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();

    let collections = store.list_collections(0, usize::MAX, &[]).unwrap();
    assert_eq!(collections.results.len(), 1);
    let digest = collections.results[0].digest.clone();

    let meta = store.get_collection_metadata(&digest);
    assert!(meta.is_some());
    assert_eq!(meta.unwrap().n_sequences, 3);

    assert!(store.is_collection_loaded(&digest));

    let stats = store.stats();
    assert_eq!(stats.n_collections, 1);
    assert_eq!(stats.n_collections_in_memory, 1);
    assert_eq!(stats.n_sequences, 3);
}

#[test]
fn test_collection_explicit_loading() {
    let temp_dir = tempdir().unwrap();
    let temp_fasta = copy_test_fasta(temp_dir.path(), "base.fa.gz");
    let cache_path = temp_dir.path().join("store");

    let mut store = RefgetStore::on_disk(&cache_path).unwrap();
    store
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();
    let digest = store.list_collections(0, usize::MAX, &[]).unwrap().results[0].digest.clone();

    drop(store);
    let mut loaded_store = RefgetStore::open_local(&cache_path).unwrap();

    let meta = loaded_store.get_collection_metadata(&digest);
    assert!(meta.is_some());
    assert_eq!(meta.unwrap().n_sequences, 3);

    assert!(!loaded_store.is_collection_loaded(&digest));

    let stats_before = loaded_store.stats();
    assert_eq!(stats_before.n_collections, 1);
    assert_eq!(stats_before.n_collections_in_memory, 0);

    let seq = loaded_store.get_sequence_by_name(&digest, "chr1");
    assert!(seq.is_err());

    loaded_store.load_collection(&digest).unwrap();

    let seq = loaded_store.get_sequence_by_name(&digest, "chr1");
    assert!(seq.is_ok());
    assert_eq!(seq.unwrap().metadata().name, "chr1");

    assert!(loaded_store.is_collection_loaded(&digest));

    let stats_after = loaded_store.stats();
    assert_eq!(stats_after.n_collections_in_memory, 1);
}

#[test]
fn test_get_collection() {
    let temp_dir = tempdir().unwrap();
    let temp_fasta = copy_test_fasta(temp_dir.path(), "base.fa.gz");
    let cache_path = temp_dir.path().join("store");

    let mut store = RefgetStore::on_disk(&cache_path).unwrap();
    store
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();
    let digest = store.list_collections(0, usize::MAX, &[]).unwrap().results[0].digest.clone();
    drop(store);

    let mut loaded_store = RefgetStore::open_local(&cache_path).unwrap();

    assert!(!loaded_store.is_collection_loaded(&digest));

    let collection = loaded_store.get_collection(&digest).unwrap();
    assert!(!collection.sequences.is_empty());

    let readonly_store = RefgetStore::open_local(&cache_path).unwrap().into_readonly();
    assert!(readonly_store.get_collection(&digest).is_err());

    loaded_store.load_all_collections().unwrap();

    let collection = loaded_store.get_collection(&digest).unwrap();
    assert!(!collection.sequences.is_empty());
    assert_eq!(collection.sequences.len(), 3);

    let stats_after = loaded_store.stats();
    // RESIDENCY INVARIANT: on a disk-backed store no sequence bytes are ever
    // held in RAM by loading a collection — records are Stubs. This is 0 even
    // right after an import, which is why the field is named `_in_memory` and
    // not `_loaded`: it is not an ingest counter.
    assert_eq!(stats_after.n_sequences_in_memory, 0);
    assert_eq!(stats_after.n_collections_in_memory, 1);

    for record in loaded_store.sequence_store.values() {
        assert!(!record.is_loaded());
    }

    let seq_digest = collection.sequences[0].metadata().sha512t24u.clone();
    loaded_store.load_sequence(&seq_digest).unwrap();
    let loaded_seq = loaded_store.get_sequence(&seq_digest).unwrap();
    assert!(loaded_seq.is_loaded());
}

#[test]
fn test_get_sequence() {
    let temp_dir = tempdir().unwrap();
    let temp_fasta = copy_test_fasta(temp_dir.path(), "base.fa.gz");
    let cache_path = temp_dir.path().join("store");

    let mut store = RefgetStore::on_disk(&cache_path).unwrap();
    store
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();

    let seq_digest = store
        .sequence_store
        .values()
        .next()
        .unwrap()
        .metadata()
        .sha512t24u
        .clone();
    drop(store);

    let mut loaded_store = RefgetStore::open_local(&cache_path).unwrap();

    let seq_before = loaded_store.sequence_store.get(&seq_digest.to_key()).unwrap();
    assert!(!seq_before.is_loaded());

    let loaded_seq = loaded_store.get_sequence(&seq_digest).unwrap();
    assert!(!loaded_seq.is_loaded());

    loaded_store.load_sequence(&seq_digest).unwrap();

    let loaded_seq = loaded_store.get_sequence(&seq_digest).unwrap();
    assert!(loaded_seq.is_loaded());
    assert!(loaded_seq.sequence().is_some());
}

#[test]
fn test_get_collection_idempotent() {
    let temp_dir = tempdir().unwrap();
    let temp_fasta = copy_test_fasta(temp_dir.path(), "base.fa.gz");
    let cache_path = temp_dir.path().join("store");

    let mut store = RefgetStore::on_disk(&cache_path).unwrap();
    store
        .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
        .unwrap();
    let digest = store.list_collections(0, usize::MAX, &[]).unwrap().results[0].digest.clone();
    drop(store);

    let mut loaded_store = RefgetStore::open_local(&cache_path).unwrap();
    loaded_store.load_all_collections().unwrap();

    let result1 = loaded_store.get_collection(&digest);
    assert!(result1.is_ok());

    let result2 = loaded_store.get_collection(&digest);
    assert!(result2.is_ok());

    assert_eq!(loaded_store.stats().n_collections_in_memory, 1);
}

// =========================================================================
// Stale cache and standalone record tests
// =========================================================================

#[test]
fn test_stale_rgsi_cache_is_ignored() {
    use std::io::Write;

    let temp_dir = tempdir().unwrap();

    let fasta_path = temp_dir.path().join("test.fa");
    let mut fasta_file = fs::File::create(&fasta_path).unwrap();
    writeln!(fasta_file, ">chr1\nATGCATGC\n>chr2\nGGGGAAAA").unwrap();

    let rgsi_path = temp_dir.path().join("test.rgsi");
    let mut rgsi_file = fs::File::create(&rgsi_path).unwrap();
    writeln!(rgsi_file, "#name\tlength\talphabet\tsha512t24u\tmd5\tdescription").unwrap();

    let store_path = temp_dir.path().join("store");
    let mut store = RefgetStore::on_disk(&store_path).unwrap();

    let result = store.add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new());
    assert!(result.is_ok(), "Should handle stale cache: {:?}", result.err());

    assert_eq!(store.sequence_store.len(), 2);
}

#[test]
fn test_add_sequence_record_standalone() {
    use crate::digest::digest_sequence;

    let mut store = RefgetStore::in_memory();
    let record = digest_sequence("test", b"ACGT");
    let digest = record.metadata().sha512t24u.clone();

    store.add_sequence_record(record, false).unwrap();

    let retrieved = store.get_sequence(digest.as_bytes()).unwrap();
    assert_eq!(retrieved.metadata().length, 4);
}

#[test]
fn test_add_sequence_record_packed_bytes_in_encoded_mode() {
    // Mirrors panget's pre-packing insert pattern: digest the ASCII sequence,
    // then pack it to the alphabet's encoded byte size before inserting into
    // an Encoded-mode store. This must be accepted (not just raw ASCII).
    use crate::digest::{digest_sequence, encode_sequence, lookup_alphabet};

    let mut store = RefgetStore::in_memory();
    store.set_encoding_mode(StorageMode::Encoded);

    let record = digest_sequence("test", b"ACGTACGT");
    let digest = record.metadata().sha512t24u.clone();

    let packed_record = match record {
        SequenceRecord::Full { metadata, sequence } => {
            let alphabet = lookup_alphabet(&metadata.alphabet);
            let encoded = encode_sequence(&*sequence, alphabet);
            SequenceRecord::Full { metadata, sequence: encoded.into() }
        }
        other => other,
    };

    store.add_sequence_record(packed_record, false).unwrap();

    let substring = store.get_substring(digest.as_bytes(), 0, 8).unwrap();
    assert_eq!(substring, "ACGTACGT");
}

// =========================================================================
// Iterator error visibility tests
// =========================================================================

#[test]
fn test_iter_collections_partial_results_on_missing_rgsi() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let fasta1 = dir.path().join("a.fa");
    let fasta2 = dir.path().join("b.fa");
    fs::write(&fasta1, ">seq1\nAAAA\n").unwrap();
    fs::write(&fasta2, ">seq2\nCCCC\n").unwrap();

    let mut store = RefgetStore::in_memory();
    let (meta1, _) = store.add_sequence_collection_from_fasta(&fasta1, FastaImportOptions::new()).unwrap();
    let (meta2, _) = store.add_sequence_collection_from_fasta(&fasta2, FastaImportOptions::new()).unwrap();
    store.write_store_to_dir(&store_path, None).unwrap();

    let rgsi1 = store_path.join(format!("collections/{}.rgsi", meta1.digest));
    let rgsi2 = store_path.join(format!("collections/{}.rgsi", meta2.digest));
    assert!(rgsi1.exists());
    assert!(rgsi2.exists());

    fs::remove_file(&rgsi1).unwrap();

    let mut loaded = RefgetStore::open_local(&store_path).unwrap();
    assert_eq!(loaded.collections.len(), 2);

    let result = loaded.load_all_collections();
    assert!(result.is_err());

    loaded.load_collection(&meta2.digest).unwrap();

    let collections: Vec<_> = loaded.iter_collections().collect();
    assert_eq!(collections.len(), 1);
    assert_eq!(collections[0].metadata.digest, meta2.digest);
}

#[test]
fn test_iter_sequences_returns_stubs_for_unloaded() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let fasta = dir.path().join("test.fa");
    fs::write(&fasta, ">seq1\nATGC\n>seq2\nGGGG\n").unwrap();

    let mut store = RefgetStore::in_memory();
    store.add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new()).unwrap();
    store.write_store_to_dir(&store_path, None).unwrap();

    let seq_digests: Vec<String> = store
        .sequence_store
        .values()
        .map(|r| r.metadata().sha512t24u.clone())
        .collect();
    assert_eq!(seq_digests.len(), 2);

    let digest_to_delete = &seq_digests[0];
    let seq_file = store_path.join(format!(
        "sequences/{}/{}.seq",
        &digest_to_delete[..2],
        digest_to_delete
    ));
    assert!(seq_file.exists());
    fs::remove_file(&seq_file).unwrap();

    let mut loaded = RefgetStore::open_local(&store_path).unwrap();
    assert_eq!(loaded.sequence_store.len(), 2);

    let sequences: Vec<_> = loaded.iter_sequences().collect();
    assert_eq!(sequences.len(), 2);
    let loaded_count = sequences.iter().filter(|r| r.is_loaded()).count();
    assert_eq!(loaded_count, 0);

    let result = loaded.load_all_sequences();
    assert!(result.is_err());

    let surviving_digest = &seq_digests[1];
    loaded.load_sequence(surviving_digest).unwrap();

    let sequences: Vec<_> = loaded.iter_sequences().collect();
    let loaded_count = sequences.iter().filter(|r| r.is_loaded()).count();
    assert_eq!(loaded_count, 1);
}

// =========================================================================
// remove_collection tests
// =========================================================================

#[test]
fn test_remove_existing_collection() {
    let (mut store, digest) = store_with_one_collection(">chr1\nACGT\n>chr2\nTTTT\n");

    assert_eq!(store.list_collections(0, usize::MAX, &[]).unwrap().results.len(), 1);

    let result = store.remove_collection(&digest, false).unwrap();
    assert!(result);
    assert_eq!(store.list_collections(0, usize::MAX, &[]).unwrap().results.len(), 0);
    assert!(store.get_collection(&digest).is_err());
}

#[test]
fn test_remove_nonexistent_collection() {
    let (mut store, _digest) = store_with_one_collection(">chr1\nACGT\n");

    let result = store.remove_collection("nonexistent_digest_value", false).unwrap();
    assert!(!result);
    assert_eq!(store.list_collections(0, usize::MAX, &[]).unwrap().results.len(), 1);
}

#[test]
fn test_remove_without_orphan_cleanup_keeps_sequences() {
    let (mut store, digest) = store_with_one_collection(">chr1\nACGT\n>chr2\nTTTT\n");

    assert_eq!(store.list_sequences().len(), 2);
    store.remove_collection(&digest, false).unwrap();
    assert_eq!(store.list_sequences().len(), 2);
}

#[test]
fn test_remove_with_orphan_cleanup_removes_sequences() {
    let (mut store, digest) = store_with_one_collection(">chr1\nACGT\n>chr2\nTTTT\n");

    assert_eq!(store.list_sequences().len(), 2);
    store.remove_collection(&digest, true).unwrap();
    assert_eq!(store.list_sequences().len(), 0);
}

#[test]
fn test_remove_with_orphan_cleanup_retains_shared_sequences() {
    let dir = tempdir().unwrap();
    let fasta1 = dir.path().join("a.fa");
    let fasta2 = dir.path().join("b.fa");
    fs::write(&fasta1, ">chr1\nACGT\n>chr2\nTTTT\n").unwrap();
    fs::write(&fasta2, ">chr1\nACGT\n>chr3\nGGGG\n").unwrap();

    let mut store = RefgetStore::in_memory();
    let (meta1, _) = store
        .add_sequence_collection_from_fasta(&fasta1, FastaImportOptions::new())
        .unwrap();
    let (meta2, _) = store
        .add_sequence_collection_from_fasta(&fasta2, FastaImportOptions::new())
        .unwrap();

    assert_eq!(store.list_collections(0, usize::MAX, &[]).unwrap().results.len(), 2);
    assert_eq!(store.list_sequences().len(), 3);

    store.remove_collection(&meta1.digest, true).unwrap();

    assert_eq!(store.list_collections(0, usize::MAX, &[]).unwrap().results.len(), 1);
    assert_eq!(store.list_sequences().len(), 2);

    let coll = store.get_collection(&meta2.digest).unwrap();
    assert_eq!(coll.sequences.len(), 2);
}

#[test]
fn test_remove_with_orphan_cleanup_clears_md5_lookup() {
    let (mut store, digest) = store_with_one_collection(">chr1\nACGT\n>chr2\nTTTT\n");

    assert_eq!(store.list_sequences().len(), 2);
    assert_eq!(store.md5_lookup.len(), 2);

    store.remove_collection(&digest, true).unwrap();

    assert_eq!(store.list_sequences().len(), 0);
    // Regression: md5_lookup must be emptied along with the sequences.
    assert_eq!(store.md5_lookup.len(), 0);
}

#[test]
fn test_remove_with_orphan_cleanup_md5_lookup_retains_shared_sequences() {
    let dir = tempdir().unwrap();
    let fasta1 = dir.path().join("a.fa");
    let fasta2 = dir.path().join("b.fa");
    // chr1/ACGT is shared; TTTT is unique to coll1; GGGG is unique to coll2.
    fs::write(&fasta1, ">chr1\nACGT\n>chr2\nTTTT\n").unwrap();
    fs::write(&fasta2, ">chr1\nACGT\n>chr3\nGGGG\n").unwrap();

    let mut store = RefgetStore::in_memory();
    let (meta1, _) = store
        .add_sequence_collection_from_fasta(&fasta1, FastaImportOptions::new())
        .unwrap();
    let (meta2, _) = store
        .add_sequence_collection_from_fasta(&fasta2, FastaImportOptions::new())
        .unwrap();

    let md5_shared = md5(b"ACGT").to_key(); // retained (shared)
    let md5_orphan = md5(b"TTTT").to_key(); // reclaimed (unique to coll1)
    let md5_other = md5(b"GGGG").to_key(); // retained (unique to coll2)

    assert_eq!(store.md5_lookup.len(), 3);

    store.remove_collection(&meta1.digest, true).unwrap();

    // Reclaimed sequence's md5 entry is gone...
    assert!(!store.md5_lookup.contains_key(&md5_orphan));
    // ...but retained sequences' md5 entries survive.
    assert!(store.md5_lookup.contains_key(&md5_shared));
    assert!(store.md5_lookup.contains_key(&md5_other));
    assert_eq!(store.md5_lookup.len(), 2);

    // And the surviving collection is still fully readable.
    let coll = store.get_collection(&meta2.digest).unwrap();
    assert_eq!(coll.sequences.len(), 2);
    assert_eq!(store.list_sequences().len(), 2);
}

#[test]
fn test_remove_collection_on_disk() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let fasta = dir.path().join("test.fa");
    fs::write(&fasta, ">chr1\nACGT\n>chr2\nTTTT\n").unwrap();

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();

    let digest = meta.digest.clone();

    let rgsi_path = store_path.join(format!("collections/{}.rgsi", digest));
    assert!(rgsi_path.exists());
    let rgci_path = store_path.join("collections.rgci");
    assert!(rgci_path.exists());

    let result = store.remove_collection(&digest, false).unwrap();
    assert!(result);

    assert!(!rgsi_path.exists());

    let rgci_content = fs::read_to_string(&rgci_path).unwrap();
    assert!(!rgci_content.contains(&digest));

    assert!(store_path.join("sequences.rgsi").exists());
    assert!(store_path.join("rgstore.json").exists());
}

#[test]
fn test_remove_collection_on_disk_with_orphan_sequences() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let fasta = dir.path().join("test.fa");
    fs::write(&fasta, ">chr1\nACGT\n").unwrap();

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();

    let digest = meta.digest.clone();

    let seq_digest = sha512t24u(b"ACGT");
    let seq_file = store_path.join(format!("sequences/{}/{}.seq", &seq_digest[..2], seq_digest));
    assert!(seq_file.exists());

    store.remove_collection(&digest, true).unwrap();

    assert!(!seq_file.exists());

    let rgsi_content = fs::read_to_string(store_path.join("sequences.rgsi")).unwrap();
    let non_comment_lines: Vec<_> = rgsi_content.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(non_comment_lines.is_empty());
}

#[test]
fn test_remove_collection_on_disk_cleans_up_empty_shard_dirs() {
    // The rmdir of the two-char shard directories is hoisted OUT of the
    // per-file unlink loop (see `remove_collection` in readonly.rs). This test
    // pins that the hoisted pass actually runs: every orphan `.seq` file is
    // unlinked AND every shard dir it emptied is removed.
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let fasta = dir.path().join("test.fa");
    fs::write(&fasta, ">chr1\nACGT\n>chr2\nTTTT\n>chr3\nGGGG\n").unwrap();

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    let digest = meta.digest.clone();

    let seq_files: Vec<PathBuf> = [b"ACGT".as_slice(), b"TTTT".as_slice(), b"GGGG".as_slice()]
        .iter()
        .map(|seq| {
            let d = sha512t24u(seq);
            store_path.join(format!("sequences/{}/{}.seq", &d[..2], d))
        })
        .collect();

    for f in &seq_files {
        assert!(f.exists(), "expected seq file to exist before removal: {:?}", f);
    }
    // Distinct shard dirs (dedup via the set) that should be emptied by removal.
    let shard_dirs: std::collections::HashSet<PathBuf> = seq_files
        .iter()
        .map(|f| f.parent().unwrap().to_path_buf())
        .collect();
    for d in &shard_dirs {
        assert!(d.is_dir(), "expected shard dir to exist before removal: {:?}", d);
    }

    store.remove_collection(&digest, true).unwrap();

    for f in &seq_files {
        assert!(!f.exists(), "orphan seq file was not removed: {:?}", f);
    }
    for d in &shard_dirs {
        assert!(
            !d.exists(),
            "emptied shard dir was not cleaned up (hoisted rmdir pass did not run): {:?}",
            d
        );
    }
    // The `sequences/` root itself is part of the store layout and stays put.
    assert!(store_path.join("sequences").is_dir());
}

// =========================================================================
// Misc tests
// =========================================================================

#[test]
fn test_store_exists() {
    let dir = tempdir().unwrap();
    let path = dir.path();

    assert!(!RefgetStore::store_exists(path));

    let mut store = RefgetStore::on_disk(path).unwrap();
    let fasta_path = dir.path().join("test.fa");
    fs::write(&fasta_path, ">seq1\nACGT\n").unwrap();
    store
        .add_sequence_collection_from_fasta(fasta_path.to_str().unwrap(), FastaImportOptions::new())
        .unwrap();
    store.write().unwrap();

    assert!(RefgetStore::store_exists(path));
    assert!(!RefgetStore::store_exists("/nonexistent/path/to/store"));
}

#[test]
fn test_clear() {
    let dir = tempdir().unwrap();
    let fasta_path = dir.path().join("test.fa");
    fs::write(&fasta_path, ">seq1\nACGT\n>seq2\nTGCA\n").unwrap();

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(fasta_path.to_str().unwrap(), FastaImportOptions::new())
        .unwrap();

    assert_eq!(store.sequence_digests().count(), 2);
    assert_eq!(store.collections.len(), 1);
    let name_count = store.name_lookup.len();

    store.clear();

    assert_eq!(store.sequence_digests().count(), 0);
    assert_eq!(store.collections.len(), 1);
    assert_eq!(store.name_lookup.len(), name_count);

    let fasta2 = dir.path().join("test2.fa");
    fs::write(&fasta2, ">seq3\nGGGG\n").unwrap();
    store
        .add_sequence_collection_from_fasta(fasta2.to_str().unwrap(), FastaImportOptions::new())
        .unwrap();
    assert_eq!(store.sequence_digests().count(), 1);
    assert_eq!(store.collections.len(), 2);
}

// =========================================================================
// State digest tests
// =========================================================================

#[test]
fn test_rgstore_json_contains_state_digests() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    store
        .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
        .unwrap();

    let json = fs::read_to_string(store_path.join("rgstore.json")).unwrap();
    let metadata: serde_json::Value = serde_json::from_str(&json).unwrap();

    assert!(metadata.get("modified").is_some(), "modified field should be present");
    assert!(metadata.get("collections_digest").is_some(), "collections_digest should be present");
    assert!(metadata.get("sequences_digest").is_some(), "sequences_digest should be present");
}

#[test]
fn test_state_digests_change_on_add() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    store
        .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
        .unwrap();

    let meta1 = store.store_metadata().unwrap();
    let coll_digest1 = meta1.get("collections_digest").cloned().unwrap();

    // Add another collection
    let fasta2 = dir.path().join("test2.fa");
    fs::write(&fasta2, ">seq_new\nTTTTAAAA\n").unwrap();
    store
        .add_sequence_collection_from_fasta(fasta2.to_str().unwrap(), FastaImportOptions::new())
        .unwrap();

    let meta2 = store.store_metadata().unwrap();
    let coll_digest2 = meta2.get("collections_digest").cloned().unwrap();

    assert_ne!(coll_digest1, coll_digest2, "collections_digest should change when a collection is added");
}

#[test]
fn test_aliases_digest_changes_on_alias_add() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    let (meta, _) = store
        .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
        .unwrap();

    let meta1 = store.store_metadata().unwrap();
    assert!(meta1.get("aliases_digest").is_none(), "no aliases yet");

    // Add an alias, then trigger index re-write by adding another collection
    store.add_sequence_alias("test_ns", "my_alias", &meta.sequences_digest).unwrap();

    let fasta2 = dir.path().join("test2.fa");
    fs::write(&fasta2, ">alias_test_seq\nAAAACCCC\n").unwrap();
    store
        .add_sequence_collection_from_fasta(fasta2.to_str().unwrap(), FastaImportOptions::new())
        .unwrap();

    let meta2 = store.store_metadata().unwrap();
    assert!(meta2.get("aliases_digest").is_some(), "aliases_digest should appear after alias + index rewrite");
}

#[test]
fn test_logical_sequence_bytes_persisted_in_manifest() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    store
        .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
        .unwrap();

    let expected = store.logical_sequence_bytes() as u64;
    assert!(expected > 0);

    // Read the raw rgstore.json and confirm the cached value matches.
    let json = fs::read_to_string(store_path.join("rgstore.json")).unwrap();
    let manifest: serde_json::Value = serde_json::from_str(&json).unwrap();
    assert_eq!(
        manifest.get("logical_sequence_bytes").and_then(|v| v.as_u64()),
        Some(expected),
        "manifest should cache logical_sequence_bytes matching the store"
    );

    // store_metadata() should surface it too.
    let meta = store.store_metadata().unwrap();
    assert_eq!(
        meta.get("logical_sequence_bytes").map(|s| s.as_str()),
        Some(expected.to_string().as_str())
    );
}

#[test]
fn test_logical_sequence_bytes_survives_alias_refresh() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    let (meta, _) = store
        .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
        .unwrap();

    let expected = store.logical_sequence_bytes() as u64;
    assert!(expected > 0);

    // Alias-only mutation goes through the cheap manifest-refresh path, which must
    // preserve the cached logical_sequence_bytes untouched.
    store.add_sequence_alias("test_ns", "my_alias", &meta.sequences_digest).unwrap();

    let json = fs::read_to_string(store_path.join("rgstore.json")).unwrap();
    let manifest: serde_json::Value = serde_json::from_str(&json).unwrap();
    assert_eq!(
        manifest.get("logical_sequence_bytes").and_then(|v| v.as_u64()),
        Some(expected),
        "logical_sequence_bytes should survive an alias-only refresh"
    );
}

#[test]
fn test_old_rgstore_json_without_state_digests_loads() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");
    fs::create_dir_all(&store_path).unwrap();

    // Write a minimal rgstore.json without the new fields
    let old_metadata = serde_json::json!({
        "version": 1,
        "seqdata_path_template": "sequences/%s2/%s.seq",
        "collections_path_template": "collections/%s.rgsi",
        "sequence_index": "sequences.rgsi",
        "collection_index": "collections.rgci",
        "mode": "Raw",
        "created_at": "2026-01-01T00:00:00Z",
        "ancillary_digests": true,
        "attribute_index": false
    });
    fs::write(store_path.join("rgstore.json"), serde_json::to_string_pretty(&old_metadata).unwrap()).unwrap();

    // Write empty index files
    fs::write(store_path.join("sequences.rgsi"), "#name\tlength\talphabet\tsha512t24u\tmd5\tdescription\n").unwrap();
    fs::write(store_path.join("collections.rgci"), "#digest\tn_sequences\tnames_digest\tsequences_digest\tlengths_digest\tname_length_pairs_digest\tsorted_name_length_pairs_digest\tsorted_sequences_digest\n").unwrap();

    // Should load without error
    let store = RefgetStore::open_local(&store_path).unwrap();
    assert_eq!(store.stats().n_sequences, 0);
}

// =========================================================================
// Order preservation tests
// =========================================================================

#[test]
fn test_collection_order_preserved_after_roundtrip() {
    // Build a FASTA with many sequences in a specific order
    let fasta_content = "\
>chr1\nACGTACGT\n\
>chr2\nGGGGAAAA\n\
>chr3\nTTTTCCCC\n\
>chr10\nAAAAAAAA\n\
>chr11\nCCCCCCCC\n\
>chr22\nGGGGGGGG\n\
>chrX\nTTTTTTTT\n\
>chrY\nACACACCA\n\
";
    let expected_names: Vec<&str> = vec!["chr1","chr2","chr3","chr10","chr11","chr22","chrX","chrY"];

    // Write FASTA to a temp dir and add to a disk-backed store
    let dir = tempdir().unwrap();
    let fasta_path = dir.path().join("order_test.fa");
    fs::write(&fasta_path, fasta_content).unwrap();

    let store_path = dir.path().join("store");
    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
        .unwrap();
    let digest = meta.digest.clone();
    store.write().unwrap();

    // Load the collection before reopening and record the names order
    store.load_all_collections().unwrap();
    let original_collection = store.get_collection(&digest).unwrap();
    let original_names: Vec<String> = original_collection
        .sequences
        .iter()
        .map(|s| s.metadata().name.clone())
        .collect();
    assert_eq!(
        original_names,
        expected_names,
        "Names should match FASTA order before roundtrip"
    );

    // Drop and reopen from disk
    drop(store);
    let mut reloaded = RefgetStore::open_local(&store_path).unwrap();
    reloaded.load_all_collections().unwrap();
    let reloaded_collection = reloaded.get_collection(&digest).unwrap();
    let reloaded_names: Vec<String> = reloaded_collection
        .sequences
        .iter()
        .map(|s| s.metadata().name.clone())
        .collect();

    assert_eq!(
        reloaded_names, original_names,
        "Sequence order must be identical after save/reload roundtrip"
    );

    // Also verify lengths and digests are in corresponding order
    let original_lengths: Vec<usize> = original_collection
        .sequences
        .iter()
        .map(|s| s.metadata().length)
        .collect();
    let reloaded_lengths: Vec<usize> = reloaded_collection
        .sequences
        .iter()
        .map(|s| s.metadata().length)
        .collect();
    assert_eq!(
        reloaded_lengths, original_lengths,
        "Lengths must be in same order after roundtrip"
    );

    let original_digests: Vec<String> = original_collection
        .sequences
        .iter()
        .map(|s| s.metadata().sha512t24u.clone())
        .collect();
    let reloaded_digests: Vec<String> = reloaded_collection
        .sequences
        .iter()
        .map(|s| s.metadata().sha512t24u.clone())
        .collect();
    assert_eq!(
        reloaded_digests, original_digests,
        "SHA512t24u digests must be in same order after roundtrip"
    );
}

/// Test that multiple collections sharing sequences under different names and different orderings
/// all preserve their correct per-collection names and element orderings across a disk roundtrip.
///
/// This covers the intersection of two previously-fixed bugs:
/// 1. HashMap ordering (fixed: inner map now IndexMap)
/// 2. Global name leakage (fixed: get_collection() overrides meta.name from name_lookup)
#[test]
fn test_shared_sequences_order_preserved_after_disk_roundtrip() {
    // FASTA A: base ordering — chrX first, then chr1, then chr2
    let fasta_a = ">chrX\nTTGGGGAA\n>chr1\nGGAA\n>chr2\nGCGC\n";
    // FASTA B: different order — chr1 first, same sequences as A
    let fasta_b = ">chr1\nGGAA\n>chr2\nGCGC\n>chrX\nTTGGGGAA\n";
    // FASTA C: name swap — chr2 has GGAA, chr1 has GCGC (opposite of A/B)
    let fasta_c = ">chrX\nTTGGGGAA\n>chr2\nGGAA\n>chr1\nGCGC\n";

    let dir = tempdir().unwrap();
    let fasta_a_path = dir.path().join("a.fa");
    let fasta_b_path = dir.path().join("b.fa");
    let fasta_c_path = dir.path().join("c.fa");
    fs::write(&fasta_a_path, fasta_a).unwrap();
    fs::write(&fasta_b_path, fasta_b).unwrap();
    fs::write(&fasta_c_path, fasta_c).unwrap();

    let store_path = dir.path().join("store");
    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    store.set_quiet(true);

    let (meta_a, _) = store.add_sequence_collection_from_fasta(&fasta_a_path, FastaImportOptions::new()).unwrap();
    let (meta_b, _) = store.add_sequence_collection_from_fasta(&fasta_b_path, FastaImportOptions::new()).unwrap();
    let (meta_c, _) = store.add_sequence_collection_from_fasta(&fasta_c_path, FastaImportOptions::new()).unwrap();

    let digest_a = meta_a.digest.clone();
    let digest_b = meta_b.digest.clone();
    let digest_c = meta_c.digest.clone();

    // Load collections before write and record level2 output
    store.load_all_collections().unwrap();
    let pre_a = store.get_collection_level2(&digest_a).unwrap();
    let pre_b = store.get_collection_level2(&digest_b).unwrap();
    let pre_c = store.get_collection_level2(&digest_c).unwrap();

    // Verify pre-write ordering for FASTA A: chrX, chr1, chr2
    assert_eq!(pre_a.names, vec!["chrX", "chr1", "chr2"], "A: names before roundtrip");
    // Verify pre-write ordering for FASTA B: chr1, chr2, chrX
    assert_eq!(pre_b.names, vec!["chr1", "chr2", "chrX"], "B: names before roundtrip");
    // Verify pre-write ordering for FASTA C: chrX, chr2, chr1 (name swap)
    assert_eq!(pre_c.names, vec!["chrX", "chr2", "chr1"], "C: names before roundtrip");

    store.write().unwrap();

    // Drop and reopen from disk
    drop(store);
    let mut reloaded = RefgetStore::open_local(&store_path).unwrap();
    reloaded.load_all_collections().unwrap();

    let post_a = reloaded.get_collection_level2(&digest_a).unwrap();
    let post_b = reloaded.get_collection_level2(&digest_b).unwrap();
    let post_c = reloaded.get_collection_level2(&digest_c).unwrap();

    // Names must match exactly (order-sensitive) after roundtrip
    assert_eq!(post_a.names, pre_a.names, "A: names after roundtrip");
    assert_eq!(post_b.names, pre_b.names, "B: names after roundtrip");
    assert_eq!(post_c.names, pre_c.names, "C: names after roundtrip");

    // Lengths must match exactly after roundtrip
    assert_eq!(post_a.lengths, pre_a.lengths, "A: lengths after roundtrip");
    assert_eq!(post_b.lengths, pre_b.lengths, "B: lengths after roundtrip");
    assert_eq!(post_c.lengths, pre_c.lengths, "C: lengths after roundtrip");

    // Sequence digests must match exactly after roundtrip
    assert_eq!(post_a.sequences, pre_a.sequences, "A: sequences after roundtrip");
    assert_eq!(post_b.sequences, pre_b.sequences, "B: sequences after roundtrip");
    assert_eq!(post_c.sequences, pre_c.sequences, "C: sequences after roundtrip");

    // Cross-check: FASTA C has chr2=GGAA and chr1=GCGC (opposite of A's chr1=GGAA, chr2=GCGC)
    // The sequence digest for chr2 in C should equal chr1 in A
    assert_eq!(
        post_c.sequences[1], post_a.sequences[1],
        "C.chr2 and A.chr1 share GGAA bytes, should have same sequence digest"
    );
    assert_eq!(
        post_c.sequences[2], post_a.sequences[2],
        "C.chr1 and A.chr2 share GCGC bytes, should have same sequence digest"
    );
}

// =========================================================================
// Name source tests
// =========================================================================

/// Test that get_collection() returns per-collection names, not global last-written names.
///
/// When the same sequence bytes appear in two collections under different names
/// (e.g., "chr1" in one, "chr2" in another), get_collection() must return the
/// correct name for each collection — not the name from whichever was loaded last.
#[test]
fn test_shared_sequence_different_names() {
    let dir = tempdir().unwrap();

    // subset.fa: contains chr1=GGAA (among others)
    let fasta1_content = ">chrX\nTTGGGGAA\n>chr1\nGGAA\n";
    let fasta1_path = dir.path().join("subset.fa");
    fs::write(&fasta1_path, fasta1_content).unwrap();

    // swap_wo_coords.fa: contains chr2=GGAA (same bytes, different name)
    let fasta2_content = ">chrX\nTTGGGGAA\n>chr2\nGGAA\n>chr1\nGCGC\n";
    let fasta2_path = dir.path().join("swap_wo_coords.fa");
    fs::write(&fasta2_path, fasta2_content).unwrap();

    let mut store = RefgetStore::in_memory();

    let (meta1, _) = store
        .add_sequence_collection_from_fasta(&fasta1_path, FastaImportOptions::new())
        .unwrap();

    let (meta2, _) = store
        .add_sequence_collection_from_fasta(&fasta2_path, FastaImportOptions::new())
        .unwrap();

    // Retrieve the collection from the first FASTA (subset.fa)
    let coll1 = store.get_collection(&meta1.digest).unwrap();
    let names1: Vec<&str> = coll1.sequences.iter().map(|s| s.metadata().name.as_str()).collect();

    // The first FASTA has chrX and chr1 — "chr1" must appear, not "chr2"
    assert!(
        names1.contains(&"chr1"),
        "Collection 1 (subset.fa) should contain 'chr1', got: {:?}",
        names1
    );
    assert!(
        !names1.contains(&"chr2"),
        "Collection 1 (subset.fa) must NOT contain 'chr2', got: {:?}",
        names1
    );

    // Retrieve the collection from the second FASTA (swap_wo_coords.fa)
    let coll2 = store.get_collection(&meta2.digest).unwrap();
    let names2: Vec<&str> = coll2.sequences.iter().map(|s| s.metadata().name.as_str()).collect();

    // The second FASTA has chrX, chr2, and chr1 — the sequence GGAA is "chr2" here
    assert!(
        names2.contains(&"chr2"),
        "Collection 2 (swap_wo_coords.fa) should contain 'chr2', got: {:?}",
        names2
    );
}

// =========================================================================
// Import collection tests
// =========================================================================

/// Helper: create a disk-backed store with one collection from a FASTA string.
fn disk_store_with_one_collection(fasta_content: &str) -> (RefgetStore, String, tempfile::TempDir, tempfile::TempDir) {
    let store_dir = tempdir().unwrap();
    let fasta_dir = tempdir().unwrap();
    let fasta = fasta_dir.path().join("test.fa");
    fs::write(&fasta, fasta_content).unwrap();

    let mut store = RefgetStore::on_disk(store_dir.path()).unwrap();
    store.set_quiet(true);
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    let digest = meta.digest.clone();
    // Load collection so name_lookup is populated
    store.load_all_collections().unwrap();
    (store, digest, store_dir, fasta_dir)
}

#[test]
fn test_import_collection_basic() {
    let (mut source, digest, _src_dir, _fasta_dir) =
        disk_store_with_one_collection(">chr1\nATGC\n>chr2\nGGGG\n");
    let target_dir = tempdir().unwrap();
    let mut target = RefgetStore::on_disk(target_dir.path()).unwrap();

    target.import_collection(&mut source, &digest).unwrap();

    // Target should have the collection
    let coll = target.get_collection(&digest).unwrap();
    assert_eq!(coll.sequences.len(), 2);
}

#[test]
fn test_import_collection_copies_sequence_aliases() {
    let (mut source, digest, _src_dir, _fasta_dir) =
        disk_store_with_one_collection(">chr1\nATGC\n>chr2\nGGGG\n");

    // Add sequence aliases in source
    let coll = source.get_collection(&digest).unwrap();
    let seq0_digest = coll.sequences[0].metadata().sha512t24u.clone();
    let seq1_digest = coll.sequences[1].metadata().sha512t24u.clone();
    source.add_sequence_alias("ncbi", "NC_000001.1", &seq0_digest).unwrap();
    source.add_sequence_alias("ucsc", "chr1", &seq0_digest).unwrap();
    source.add_sequence_alias("ncbi", "NC_000002.1", &seq1_digest).unwrap();

    let target_dir = tempdir().unwrap();
    let mut target = RefgetStore::on_disk(target_dir.path()).unwrap();
    target.import_collection(&mut source, &digest).unwrap();

    // Target should have the sequence aliases
    let ns = target.list_sequence_alias_namespaces();
    assert!(ns.contains(&"ncbi".to_string()), "Missing ncbi namespace: {:?}", ns);
    assert!(ns.contains(&"ucsc".to_string()), "Missing ucsc namespace: {:?}", ns);

    // Verify forward lookup
    let resolved = target.get_sequence_metadata_by_alias("ncbi", "NC_000001.1");
    assert!(resolved.is_some(), "ncbi alias NC_000001.1 not found in target");
    assert_eq!(resolved.unwrap().sha512t24u, seq0_digest);
}

#[test]
fn test_import_collection_copies_collection_aliases() {
    let (mut source, digest, _src_dir, _fasta_dir) =
        disk_store_with_one_collection(">chr1\nATGC\n>chr2\nGGGG\n");

    source.add_collection_alias("insdc", "GCA_000001.1", &digest).unwrap();
    source.add_collection_alias("refseq", "GCF_000001.1", &digest).unwrap();

    let target_dir = tempdir().unwrap();
    let mut target = RefgetStore::on_disk(target_dir.path()).unwrap();
    target.import_collection(&mut source, &digest).unwrap();

    let ns = target.list_collection_alias_namespaces();
    assert!(ns.contains(&"insdc".to_string()), "Missing insdc namespace: {:?}", ns);
    assert!(ns.contains(&"refseq".to_string()), "Missing refseq namespace: {:?}", ns);

    let aliases = target.get_aliases_for_collection(&digest);
    assert_eq!(aliases.len(), 2, "Expected 2 collection aliases, got {:?}", aliases);
}

#[test]
fn test_import_collection_copies_fhr_metadata() {
    use super::fhr_metadata::FhrMetadata;

    let (mut source, digest, _src_dir, _fasta_dir) =
        disk_store_with_one_collection(">chr1\nATGC\n");

    let fhr = FhrMetadata {
        genome: Some("Homo sapiens".to_string()),
        ..Default::default()
    };
    source.set_fhr_metadata(&digest, fhr).unwrap();

    let target_dir = tempdir().unwrap();
    let mut target = RefgetStore::on_disk(target_dir.path()).unwrap();
    target.import_collection(&mut source, &digest).unwrap();

    let fhr = target.get_fhr_metadata(&digest);
    assert!(fhr.is_some(), "FHR metadata not copied");
    assert_eq!(fhr.unwrap().genome.as_deref(), Some("Homo sapiens"));
}

#[test]
fn test_import_collection_disk_roundtrip_aliases() {
    // This test catches the bug where aliases loaded from disk aren't available
    // for reverse lookup during import.
    let source_dir = tempdir().unwrap();
    let target_dir = tempdir().unwrap();
    let fasta_dir = tempdir().unwrap();
    let fasta = fasta_dir.path().join("test.fa");
    fs::write(&fasta, ">chr1\nATGC\n>chr2\nGGGG\n").unwrap();

    // Create source store on disk, add aliases, then drop it
    let digest;
    let seq0_digest;
    {
        let mut source = RefgetStore::on_disk(source_dir.path()).unwrap();
        source.set_quiet(true);
        let (meta, _) = source
            .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
            .unwrap();
        digest = meta.digest.clone();

        let coll = source.get_collection(&digest).unwrap();
        seq0_digest = coll.sequences[0].metadata().sha512t24u.clone();

        source.add_sequence_alias("ncbi", "NC_000001.1", &seq0_digest).unwrap();
        source.add_collection_alias("insdc", "GCA_000001.1", &digest).unwrap();
    }
    // Source is dropped here

    // Reopen from disk (this tests that aliases are loaded from disk)
    let mut source = RefgetStore::on_disk(source_dir.path()).unwrap();
    source.load_all_collections().unwrap();

    // Verify aliases were loaded from disk
    let seq_ns = source.list_sequence_alias_namespaces();
    assert!(seq_ns.contains(&"ncbi".to_string()), "Source lost seq aliases after reopen: {:?}", seq_ns);
    let coll_ns = source.list_collection_alias_namespaces();
    assert!(coll_ns.contains(&"insdc".to_string()), "Source lost coll aliases after reopen: {:?}", coll_ns);

    // Import into a new disk-backed target
    let mut target = RefgetStore::on_disk(target_dir.path()).unwrap();
    target.import_collection(&mut source, &digest).unwrap();

    // Verify target has the aliases
    let ns = target.list_sequence_alias_namespaces();
    assert!(ns.contains(&"ncbi".to_string()), "Target missing ncbi seq alias: {:?}", ns);

    let resolved = target.get_sequence_metadata_by_alias("ncbi", "NC_000001.1");
    assert!(resolved.is_some(), "Alias NC_000001.1 not found in target after disk roundtrip");
    assert_eq!(resolved.unwrap().sha512t24u, seq0_digest);

    let coll_aliases = target.get_aliases_for_collection(&digest);
    assert_eq!(coll_aliases.len(), 1, "Expected 1 collection alias in target: {:?}", coll_aliases);
}

#[test]
fn test_import_collection_file_copy_roundtrip() {
    // Verify RGSI and .seq files are byte-for-byte identical after import
    // when ancillary digests match between source and dest.
    let (mut source, digest, src_dir, _fasta_dir) =
        disk_store_with_one_collection(">chr1\nATGC\n>chr2\nGGGG\n");
    let target_dir = tempdir().unwrap();
    let mut target = RefgetStore::on_disk(target_dir.path()).unwrap();

    target.import_collection(&mut source, &digest).unwrap();

    // Verify RGSI file is byte-for-byte identical
    let src_rgsi = fs::read(
        src_dir.path().join(format!("collections/{}.rgsi", digest)),
    ).unwrap();
    let dst_rgsi = fs::read(
        target_dir.path().join(format!("collections/{}.rgsi", digest)),
    ).unwrap();
    assert_eq!(src_rgsi, dst_rgsi, "RGSI files should be byte-identical");

    // Verify .seq files are byte-for-byte identical
    let coll = target.get_collection(&digest).unwrap();
    for seq in &coll.sequences {
        let seq_digest = &seq.metadata().sha512t24u;
        let src_seq_path = source.sequence_file_path(seq_digest).unwrap();
        let dst_seq_path = target.sequence_file_path(seq_digest).unwrap();
        let src_data = fs::read(&src_seq_path).unwrap();
        let dst_data = fs::read(&dst_seq_path).unwrap();
        assert_eq!(src_data, dst_data, "Sequence file for {} should be byte-identical", seq_digest);
    }

    // Verify in-memory metadata matches
    let src_coll = source.get_collection(&digest).unwrap();
    assert_eq!(coll.metadata.digest, src_coll.metadata.digest);
    assert_eq!(coll.sequences.len(), src_coll.sequences.len());
    for (src_seq, dst_seq) in src_coll.sequences.iter().zip(coll.sequences.iter()) {
        assert_eq!(src_seq.metadata().sha512t24u, dst_seq.metadata().sha512t24u);
        assert_eq!(src_seq.metadata().name, dst_seq.metadata().name);
        assert_eq!(src_seq.metadata().length, dst_seq.metadata().length);
    }
}

#[test]
fn test_import_collection_ancillary_digest_enrichment() {
    // Source store has ancillary_digests: false, destination has ancillary_digests: true.
    // The destination RGSI should contain ancillary digest headers that the source lacks.
    let source_dir = tempdir().unwrap();
    let fasta_dir = tempdir().unwrap();
    let fasta = fasta_dir.path().join("test.fa");
    fs::write(&fasta, ">chr1\nATGC\n>chr2\nGGGG\n").unwrap();

    // Create source store with ancillary_digests: false
    let mut source = RefgetStore::on_disk(source_dir.path()).unwrap();
    source.set_quiet(true);
    source.disable_ancillary_digests();
    let (meta, _) = source
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    let digest = meta.digest.clone();
    source.load_all_collections().unwrap();

    // Verify source RGSI lacks ancillary digests
    let src_rgsi_content = fs::read_to_string(
        source_dir.path().join(format!("collections/{}.rgsi", digest)),
    ).unwrap();
    assert!(
        !src_rgsi_content.contains("name_length_pairs_digest"),
        "Source should NOT have ancillary digests",
    );

    // Create destination store with ancillary_digests: true (default)
    let target_dir = tempdir().unwrap();
    let mut target = RefgetStore::on_disk(target_dir.path()).unwrap();
    target.enable_ancillary_digests();
    target.import_collection(&mut source, &digest).unwrap();

    // Verify destination RGSI has ancillary digest headers
    let dst_rgsi_content = fs::read_to_string(
        target_dir.path().join(format!("collections/{}.rgsi", digest)),
    ).unwrap();
    assert!(
        dst_rgsi_content.contains("name_length_pairs_digest"),
        "Destination should have ancillary digests. RGSI:\n{}",
        dst_rgsi_content,
    );
    assert!(
        dst_rgsi_content.contains("sorted_name_length_pairs_digest"),
        "Destination should have sorted_name_length_pairs_digest",
    );
    assert!(
        dst_rgsi_content.contains("sorted_sequences_digest"),
        "Destination should have sorted_sequences_digest",
    );

    // Verify in-memory metadata has non-None ancillary fields
    let coll_meta = target.get_collection_metadata(&digest).unwrap();
    assert!(coll_meta.name_length_pairs_digest.is_some(), "name_length_pairs_digest should be Some");
    assert!(coll_meta.sorted_name_length_pairs_digest.is_some(), "sorted_name_length_pairs_digest should be Some");
    assert!(coll_meta.sorted_sequences_digest.is_some(), "sorted_sequences_digest should be Some");
}

#[test]
fn test_import_collection_mode_mismatch_error() {
    // Source with Raw mode, destination with Encoded mode should fail.
    let source_dir = tempdir().unwrap();
    let fasta_dir = tempdir().unwrap();
    let fasta = fasta_dir.path().join("test.fa");
    fs::write(&fasta, ">chr1\nATGC\n").unwrap();

    // Create source store in Raw mode
    let mut source = RefgetStore::on_disk(source_dir.path()).unwrap();
    source.set_quiet(true);
    source.disable_encoding();
    let (meta, _) = source
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    let digest = meta.digest.clone();
    source.load_all_collections().unwrap();

    // Create destination store in Encoded mode (default)
    let target_dir = tempdir().unwrap();
    let mut target = RefgetStore::on_disk(target_dir.path()).unwrap();
    // target uses Encoded mode by default

    let result = target.import_collection(&mut source, &digest);
    assert!(result.is_err(), "Should fail with mode mismatch");
    let err_msg = result.unwrap_err().to_string();
    assert!(
        err_msg.contains("matching storage modes"),
        "Error should mention storage modes: {}",
        err_msg,
    );
}

// =========================================================================
// Parallel encode parity / determinism tests
// =========================================================================

/// Recursively collect every `*.seq` file under `dir`, keyed by path relative
/// to `dir`, with its bytes. Used to compare two on-disk stores.
fn collect_seq_files(dir: &std::path::Path) -> std::collections::BTreeMap<String, Vec<u8>> {
    fn walk(
        base: &std::path::Path,
        cur: &std::path::Path,
        out: &mut std::collections::BTreeMap<String, Vec<u8>>,
    ) {
        for entry in std::fs::read_dir(cur).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.is_dir() {
                walk(base, &path, out);
            } else if path.extension().and_then(|e| e.to_str()) == Some("seq") {
                let rel = path.strip_prefix(base).unwrap().to_string_lossy().to_string();
                out.insert(rel, std::fs::read(&path).unwrap());
            }
        }
    }
    let mut out = std::collections::BTreeMap::new();
    let seqdir = dir.join("sequences");
    if seqdir.exists() {
        walk(dir, &seqdir, &mut out);
    }
    out
}

/// Multi-record FASTA with names deliberately NOT in sorted order so that a
/// digest-sorted index would differ from FASTA (name_lookup) order, exercising
/// the order-preservation logic.
const PARITY_FASTA: &str =
    ">chrX\nTTGGGGAACCCCTTTT\n>chr1\nGGAATTCCGGAATTCC\n>chr2\nACGTACGTACGTACGT\n\
     >chrM\nGGGGCCCCAAAATTTT\n>chr10\nTACGTACGTACGTACG\n";

fn build_on_disk_store(dir: &std::path::Path, fasta: &std::path::Path, jobs: usize) -> String {
    let mut store = RefgetStore::on_disk(dir).unwrap();
    store.set_quiet(true);
    let (meta, _) = store
        .add_sequence_collection_from_fasta(fasta, FastaImportOptions::new().jobs(jobs))
        .unwrap();
    meta.digest
}

#[test]
fn test_parallel_equals_serial_on_disk() {
    let work = tempdir().unwrap();
    let fasta = work.path().join("parity.fa");
    fs::write(&fasta, PARITY_FASTA).unwrap();

    let dir_serial = tempdir().unwrap();
    let dir_parallel = tempdir().unwrap();

    let digest_serial = build_on_disk_store(dir_serial.path(), &fasta, 1);
    let digest_parallel = build_on_disk_store(dir_parallel.path(), &fasta, 8);

    // Collection digest identical.
    assert_eq!(digest_serial, digest_parallel, "collection digest must match");

    // sequences.rgsi byte-identical (digest-sorted, so order-independent).
    let rgsi_serial = fs::read(dir_serial.path().join("sequences.rgsi")).unwrap();
    let rgsi_parallel = fs::read(dir_parallel.path().join("sequences.rgsi")).unwrap();
    assert_eq!(rgsi_serial, rgsi_parallel, "sequences.rgsi must be byte-identical");

    // collections.rgci byte-identical.
    let rgci_serial = fs::read(dir_serial.path().join("collections.rgci")).unwrap();
    let rgci_parallel = fs::read(dir_parallel.path().join("collections.rgci")).unwrap();
    assert_eq!(rgci_serial, rgci_parallel, "collections.rgci must be byte-identical");

    // Every .seq file present with identical bytes (both directions).
    let seqs_serial = collect_seq_files(dir_serial.path());
    let seqs_parallel = collect_seq_files(dir_parallel.path());
    assert_eq!(
        seqs_serial, seqs_parallel,
        ".seq files must match between serial and parallel builds"
    );
    assert!(!seqs_serial.is_empty(), "expected some .seq files");
}

#[test]
fn test_parallel_name_lookup_matches_fasta_order() {
    let work = tempdir().unwrap();
    let fasta = work.path().join("parity.fa");
    fs::write(&fasta, PARITY_FASTA).unwrap();

    // Raw FASTA header order.
    let fasta_order = vec!["chrX", "chr1", "chr2", "chrM", "chr10"];

    let dir_serial = tempdir().unwrap();
    let dir_parallel = tempdir().unwrap();
    let digest = build_on_disk_store(dir_serial.path(), &fasta, 1);
    let _ = build_on_disk_store(dir_parallel.path(), &fasta, 8);

    let mut serial = RefgetStore::open_local(dir_serial.path()).unwrap();
    let mut parallel = RefgetStore::open_local(dir_parallel.path()).unwrap();
    serial.load_all_collections().unwrap();
    parallel.load_all_collections().unwrap();

    let key = digest.to_key();
    let serial_names: Vec<String> = serial
        .name_lookup
        .get(&key)
        .unwrap()
        .keys()
        .cloned()
        .collect();
    let parallel_names: Vec<String> = parallel
        .name_lookup
        .get(&key)
        .unwrap()
        .keys()
        .cloned()
        .collect();

    assert_eq!(serial_names, fasta_order, "serial name_lookup must follow FASTA order");
    assert_eq!(
        parallel_names, fasta_order,
        "parallel name_lookup must follow FASTA order"
    );
}

#[test]
fn test_parallel_equals_serial_in_memory() {
    let work = tempdir().unwrap();
    let fasta = work.path().join("parity.fa");
    fs::write(&fasta, PARITY_FASTA).unwrap();

    let mut serial = RefgetStore::in_memory();
    serial.set_quiet(true);
    let (meta_s, _) = serial
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new().jobs(1))
        .unwrap();

    let mut parallel = RefgetStore::in_memory();
    parallel.set_quiet(true);
    let (meta_p, _) = parallel
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new().jobs(8))
        .unwrap();

    assert_eq!(meta_s.digest, meta_p.digest);

    let key = meta_s.digest.to_key();

    // Identical name_lookup order and contents.
    let names_s: Vec<(String, DigestKey)> = serial
        .name_lookup
        .get(&key)
        .unwrap()
        .iter()
        .map(|(n, d)| (n.clone(), *d))
        .collect();
    let names_p: Vec<(String, DigestKey)> = parallel
        .name_lookup
        .get(&key)
        .unwrap()
        .iter()
        .map(|(n, d)| (n.clone(), *d))
        .collect();
    assert_eq!(names_s, names_p, "name_lookup order+contents must match");

    // Identical set of digests.
    let digs_s: std::collections::BTreeSet<DigestKey> = serial.sequence_digests().collect();
    let digs_p: std::collections::BTreeSet<DigestKey> = parallel.sequence_digests().collect();
    assert_eq!(digs_s, digs_p, "digest sets must match");

    // Per-digest metadata + decoded bytes parity.
    for name in ["chrX", "chr1", "chr2", "chrM", "chr10"] {
        let rec_s = serial.get_sequence_by_name(&meta_s.digest, name).unwrap();
        let rec_p = parallel.get_sequence_by_name(&meta_p.digest, name).unwrap();
        let m_s = rec_s.metadata();
        let m_p = rec_p.metadata();
        assert_eq!(m_s.sha512t24u, m_p.sha512t24u);
        assert_eq!(m_s.md5, m_p.md5);
        assert_eq!(m_s.alphabet, m_p.alphabet);
        assert_eq!(m_s.length, m_p.length);

        let sub_s = serial.get_substring(&m_s.sha512t24u, 0, m_s.length).unwrap();
        let sub_p = parallel.get_substring(&m_p.sha512t24u, 0, m_p.length).unwrap();
        assert_eq!(sub_s, sub_p, "decoded bytes must match for {}", name);
    }
}

#[test]
fn test_parallel_cached_metadata_parity() {
    // First import (no cache) writes the .rgsi cache; second import hits the
    // cached-metadata fast path. Build serial+parallel on-disk stores via the
    // cached path and assert full equality.

    // Build the rgsi cache by importing once into a throwaway store. The cache
    // is written next to the FASTA, so use a per-arm FASTA copy.
    let build_cached = |jobs: usize| -> (tempfile::TempDir, String) {
        let src_dir = tempdir().unwrap();
        let fasta = src_dir.path().join("cached.fa");
        fs::write(&fasta, PARITY_FASTA).unwrap();

        // Prime the cache (writes cached.rgsi next to the FASTA).
        let prime_dir = tempdir().unwrap();
        let mut prime = RefgetStore::on_disk(prime_dir.path()).unwrap();
        prime.set_quiet(true);
        prime
            .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new().jobs(1))
            .unwrap();
        assert!(
            src_dir.path().join("cached.rgsi").exists(),
            "rgsi cache should have been written"
        );

        // Now build a fresh store via the cached fast path.
        let dir = tempdir().unwrap();
        let mut store = RefgetStore::on_disk(dir.path()).unwrap();
        store.set_quiet(true);
        let (meta, _) = store
            .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new().jobs(jobs))
            .unwrap();
        (dir, meta.digest)
    };

    let (dir_serial, digest_serial) = build_cached(1);
    let (dir_parallel, digest_parallel) = build_cached(8);

    assert_eq!(digest_serial, digest_parallel);

    let rgsi_s = fs::read(dir_serial.path().join("sequences.rgsi")).unwrap();
    let rgsi_p = fs::read(dir_parallel.path().join("sequences.rgsi")).unwrap();
    assert_eq!(rgsi_s, rgsi_p, "cached-path sequences.rgsi must match");

    let seqs_s = collect_seq_files(dir_serial.path());
    let seqs_p = collect_seq_files(dir_parallel.path());
    assert_eq!(seqs_s, seqs_p, "cached-path .seq files must match");

    // name_lookup order should follow FASTA order in both.
    let mut serial = RefgetStore::open_local(dir_serial.path()).unwrap();
    let mut parallel = RefgetStore::open_local(dir_parallel.path()).unwrap();
    serial.load_all_collections().unwrap();
    parallel.load_all_collections().unwrap();
    let key = digest_serial.to_key();
    let names_s: Vec<String> = serial.name_lookup.get(&key).unwrap().keys().cloned().collect();
    let names_p: Vec<String> = parallel.name_lookup.get(&key).unwrap().keys().cloned().collect();
    assert_eq!(names_s, vec!["chrX", "chr1", "chr2", "chrM", "chr10"]);
    assert_eq!(names_s, names_p);
}

#[test]
fn test_parallel_duplicate_contig_dedup() {
    // Two records with identical sequence but different names. Must dedup to a
    // single .seq file / single sequence_store entry; name_lookup keeps both
    // names in FASTA order mapping to the same digest.
    let fasta_content = ">dupA\nACGTACGTACGT\n>uniq\nGGGGCCCCAAAA\n>dupB\nACGTACGTACGT\n";
    let work = tempdir().unwrap();
    let fasta = work.path().join("dup.fa");
    fs::write(&fasta, fasta_content).unwrap();

    let dir = tempdir().unwrap();
    let mut store = RefgetStore::on_disk(dir.path()).unwrap();
    store.set_quiet(true);
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new().jobs(8))
        .unwrap();

    // Two distinct digests (the duplicate collapses).
    let digs: std::collections::BTreeSet<DigestKey> = store.sequence_digests().collect();
    assert_eq!(digs.len(), 2, "duplicate sequence must dedup to 2 entries");

    let seqs = collect_seq_files(dir.path());
    assert_eq!(seqs.len(), 2, "duplicate sequence must dedup to 2 .seq files");

    let key = meta.digest.to_key();
    let nl = store.name_lookup.get(&key).unwrap();
    let names: Vec<String> = nl.keys().cloned().collect();
    assert_eq!(names, vec!["dupA", "uniq", "dupB"], "name_lookup keeps FASTA order");
    // dupA and dupB map to same digest.
    assert_eq!(nl.get("dupA"), nl.get("dupB"));
    assert_ne!(nl.get("dupA"), nl.get("uniq"));
}

#[test]
fn test_parallel_determinism_across_runs() {
    let work = tempdir().unwrap();
    let fasta = work.path().join("parity.fa");
    fs::write(&fasta, PARITY_FASTA).unwrap();

    let mut prev: Option<Vec<u8>> = None;
    for _ in 0..3 {
        let dir = tempdir().unwrap();
        build_on_disk_store(dir.path(), &fasta, 8);
        let rgsi = fs::read(dir.path().join("sequences.rgsi")).unwrap();
        if let Some(p) = &prev {
            assert_eq!(p, &rgsi, "sequences.rgsi must be identical across parallel runs");
        }
        prev = Some(rgsi);
    }
}

// =========================================================================
// Multi-file (file-level) parallelism determinism tests
// =========================================================================

/// Collect every data/index file under a store dir EXCEPT the `rgstore.json`
/// manifest (which embeds wall-clock timestamps). Keyed by path relative to
/// `dir`, with bytes. Used to assert byte-identical stores.
fn collect_store_data_files(
    dir: &std::path::Path,
) -> std::collections::BTreeMap<String, Vec<u8>> {
    fn walk(
        base: &std::path::Path,
        cur: &std::path::Path,
        out: &mut std::collections::BTreeMap<String, Vec<u8>>,
    ) {
        for entry in std::fs::read_dir(cur).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.is_dir() {
                walk(base, &path, out);
            } else {
                let rel = path.strip_prefix(base).unwrap().to_string_lossy().to_string();
                // rgstore.json embeds created_at/modified timestamps and the
                // hashes of those files; skip it for byte-identity comparison.
                if rel == "rgstore.json" {
                    continue;
                }
                out.insert(rel, std::fs::read(&path).unwrap());
            }
        }
    }
    let mut out = std::collections::BTreeMap::new();
    walk(dir, dir, &mut out);
    out
}

/// Write `content` gzip-compressed to `path`.
fn write_gzipped(path: &std::path::Path, content: &str) {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;
    let f = fs::File::create(path).unwrap();
    let mut enc = GzEncoder::new(f, Compression::default());
    enc.write_all(content.as_bytes()).unwrap();
    enc.finish().unwrap();
}

/// Build an on-disk store from multiple FASTA files with the given `jobs`
/// (files-in-flight) setting. Returns the per-file collection digests in input
/// order.
fn build_multi(
    dir: &std::path::Path,
    files: &[std::path::PathBuf],
    jobs: usize,
) -> Vec<String> {
    let mut store = RefgetStore::on_disk(dir).unwrap();
    store.set_quiet(true);
    let opts = FastaImportOptions::new().jobs(jobs);
    store
        .add_sequence_collections_from_fastas(files, opts)
        .unwrap()
        .collections
        .into_iter()
        .map(|(m, _)| m.digest)
        .collect()
}

/// The per-run `ImportReport` counters must reflect what an import ACTUALLY
/// did, unlike `stats()`, which is a RAM-residency snapshot.
#[test]
fn test_import_report_per_run_counters() {
    let work = tempdir().unwrap();

    // Same fixture shape as the parallel-equals-serial test: file "c" shares
    // the chr1 sequence content with file "a", so 6 records are seen but only
    // 5 distinct digests exist.
    let fa_a = work.path().join("a.fa");
    let fa_b = work.path().join("b.fa");
    let fa_c = work.path().join("c.fa");
    fs::write(&fa_a, ">chr1\nGGAATTCCGGAATTCC\n>chr2\nACGTACGTACGTACGT\n").unwrap();
    fs::write(&fa_b, ">chrX\nTTGGGGAACCCCTTTT\n>chrM\nGGGGCCCCAAAATTTT\n").unwrap();
    fs::write(&fa_c, ">altchr1\nGGAATTCCGGAATTCC\n>chr9\nTACGTACGTACGTACG\n").unwrap();
    let files = vec![fa_a, fa_b, fa_c];

    let store_dir = tempdir().unwrap();
    let mut store = RefgetStore::on_disk(store_dir.path()).unwrap();
    store.set_quiet(true);

    let report = store
        .add_sequence_collections_from_fastas(&files, FastaImportOptions::new().jobs(4))
        .unwrap();

    assert_eq!(report.collections.len(), 3);
    assert_eq!(report.n_collections_new, 3);
    // Each distinct digest is written exactly once; the shared chr1 content is
    // seen twice, so exactly one record is deduped.
    assert_eq!(report.n_sequences_written, 5);
    assert_eq!(report.n_sequences_deduped, 1);
    // Every record seen across all processed files is accounted for.
    assert_eq!(
        report.n_sequences_written + report.n_sequences_deduped,
        6,
        "written + deduped must equal the records seen"
    );
    assert_eq!(collect_seq_files(store_dir.path()).len(), 5);

    // Re-importing the same files into the same store adds nothing. Whether the
    // records are counted as deduped or not seen at all depends on the `.rgsi`
    // sidecar short-circuit (which skips a present collection before its FASTA
    // is ever opened), so only the "nothing new" half is pinned here.
    let report2 = store
        .add_sequence_collections_from_fastas(&files, FastaImportOptions::new().jobs(4))
        .unwrap();
    assert_eq!(report2.collections.len(), 3);
    assert_eq!(report2.n_collections_new, 0);
    assert_eq!(
        report2.n_sequences_written, 0,
        "a re-import must write no sequence bytes"
    );
    assert!(report2.collections.iter().all(|(_, was_new)| !*was_new));

    // Contrast: the residency gauge says nothing about either run.
    assert_eq!(store.stats().n_sequences_in_memory, 0);

    // An in-memory store dispatches no disk writes at all, so it is classified
    // by a separate branch; it must produce the same counts.
    let mut mem_store = RefgetStore::in_memory();
    mem_store.set_quiet(true);
    let mem_report = mem_store
        .add_sequence_collections_from_fastas(&files, FastaImportOptions::new().jobs(4))
        .unwrap();
    assert_eq!(mem_report.n_collections_new, 3);
    assert_eq!(mem_report.n_sequences_written, 5);
    assert_eq!(mem_report.n_sequences_deduped, 1);
}

#[test]
fn test_multifile_parallel_equals_serial() {
    let work = tempdir().unwrap();

    // Three distinct collections. File "c" deliberately SHARES a sequence
    // (the chr1 record) with file "a" to exercise cross-collection dedup.
    let fa_a = work.path().join("a.fa");
    let fa_b = work.path().join("b.fa");
    let fa_c = work.path().join("c.fa");
    fs::write(&fa_a, ">chr1\nGGAATTCCGGAATTCC\n>chr2\nACGTACGTACGTACGT\n").unwrap();
    fs::write(&fa_b, ">chrX\nTTGGGGAACCCCTTTT\n>chrM\nGGGGCCCCAAAATTTT\n").unwrap();
    // Shares the chr1 sequence content with file a (different collection).
    fs::write(&fa_c, ">altchr1\nGGAATTCCGGAATTCC\n>chr9\nTACGTACGTACGTACG\n").unwrap();
    let files = vec![fa_a.clone(), fa_b.clone(), fa_c.clone()];

    let dir_serial = tempdir().unwrap();
    let dir_parallel = tempdir().unwrap();

    let digests_serial = build_multi(dir_serial.path(), &files, 1);
    let digests_parallel = build_multi(dir_parallel.path(), &files, 8);

    // Per-file collection digests identical and in the same (input) order.
    assert_eq!(
        digests_serial, digests_parallel,
        "per-file collection digests must match in input order"
    );

    // The whole store (sans rgstore.json timestamp manifest) is byte-identical.
    let serial = collect_store_data_files(dir_serial.path());
    let parallel = collect_store_data_files(dir_parallel.path());
    assert_eq!(
        serial, parallel,
        "multi-file parallel store must be byte-identical to serial store"
    );

    // Cross-collection dedup: the chr1 sequence is shared by file a and file c,
    // so the total number of distinct .seq files must be 5 (chr1, chr2, chrX,
    // chrM, chr9) -- NOT 6.
    let seqs = collect_seq_files(dir_serial.path());
    assert_eq!(seqs.len(), 5, "shared chr1 sequence must be stored once");
}

#[test]
fn test_multifile_parallel_gzipped() {
    let work = tempdir().unwrap();

    // Same content as the plain-text case, but gzipped to exercise K concurrent
    // MultiGzDecoders.
    let fa_a = work.path().join("a.fa.gz");
    let fa_b = work.path().join("b.fa.gz");
    let fa_c = work.path().join("c.fa.gz");
    write_gzipped(&fa_a, ">chr1\nGGAATTCCGGAATTCC\n>chr2\nACGTACGTACGTACGT\n");
    write_gzipped(&fa_b, ">chrX\nTTGGGGAACCCCTTTT\n>chrM\nGGGGCCCCAAAATTTT\n");
    write_gzipped(&fa_c, ">altchr1\nGGAATTCCGGAATTCC\n>chr9\nTACGTACGTACGTACG\n");
    let files = vec![fa_a, fa_b, fa_c];

    let dir_serial = tempdir().unwrap();
    let dir_parallel = tempdir().unwrap();

    let digests_serial = build_multi(dir_serial.path(), &files, 1);
    let digests_parallel = build_multi(dir_parallel.path(), &files, 8);

    assert_eq!(digests_serial, digests_parallel, "gzipped digests must match");

    let serial = collect_store_data_files(dir_serial.path());
    let parallel = collect_store_data_files(dir_parallel.path());
    assert_eq!(
        serial, parallel,
        "gzipped multi-file parallel store must be byte-identical to serial"
    );
}

/// Regression test for the parallel-import "resume re-decode" bug.
///
/// On RESUME (re-importing a FASTA whose collection is already in the store and
/// whose sibling `.rgsi` cache exists), the import MUST short-circuit BEFORE
/// opening/decoding the FASTA. The old serial path did this; a refactor moved
/// the "already exists" check to `End`, AFTER the whole file had been
/// re-decoded + re-encoded -- wasting enormous CPU on resume.
///
/// This test makes the regression deterministic WITHOUT any production hooks:
/// after the first import we CORRUPT the FASTA bytes while keeping the `.rgsi`
/// sidecar intact. A correct (short-circuiting) implementation never opens the
/// corrupt FASTA, so the resume succeeds, reports `was_new = false`, and leaves
/// the store unchanged. A regressed implementation would try to decode the
/// garbage FASTA and either error or mis-build -- failing the test.
#[test]
fn test_resume_skips_present_collection_before_decoding_fasta() {
    let temp_dir = tempdir().unwrap();
    let temp_path = temp_dir.path();

    // FASTA whose collection we will import, then "resume".
    let fasta_path = temp_path.join("genome.fa");
    fs::write(
        &fasta_path,
        ">chr1\nACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAA\n",
    )
    .unwrap();

    // Disk-backed store so RGSI caching is enabled (use_cache = local_path set).
    let cache_path = temp_path.join("store");
    let mut store = RefgetStore::on_disk(&cache_path).unwrap();

    // First import: builds the collection and writes the sibling genome.rgsi.
    let (meta_first, was_new_first) = store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new().jobs(4))
        .unwrap();
    assert!(was_new_first, "first import must add a new collection");

    let rgsi_sidecar = temp_path.join("genome.rgsi");
    assert!(
        rgsi_sidecar.exists(),
        "first import must write the sibling .rgsi cache"
    );

    // --force must NOT skip: with the FASTA still intact and the collection
    // already present, forcing must RE-PROCESS it (was_new = true), proving the
    // short-circuit respects --force. (Done before corrupting so the forced
    // decode has a valid FASTA to read.)
    let (_meta_forced, was_new_forced) = store
        .add_sequence_collection_from_fasta(
            &fasta_path,
            FastaImportOptions::new().jobs(4).force(true),
        )
        .expect("forced re-import of an intact FASTA must succeed");
    assert!(
        was_new_forced,
        "--force must bypass the skip and re-process the present collection (was_new = true)"
    );

    // Snapshot store contents to prove the resume changes nothing.
    let collections_before: Vec<DigestKey> = {
        let mut v: Vec<DigestKey> = store.collections.keys().cloned().collect();
        v.sort();
        v
    };
    let sequences_before: Vec<DigestKey> = {
        let mut v: Vec<DigestKey> = store.sequence_store.keys().cloned().collect();
        v.sort();
        v
    };

    // CORRUPT the FASTA: replace it with a record whose sequence NAME is absent
    // from the cached metadata. Any code path that actually opens + decodes this
    // FASTA fails ("Sequence '...' not found in cached metadata"). The .rgsi
    // sidecar stays intact, so the collection digest is still discoverable for
    // free -- a correct short-circuit never reaches the decode.
    fs::write(&fasta_path, b">totally_bogus_name_not_in_cache\nNNNNNNNN\n").unwrap();

    // Resume: re-import the SAME (now-corrupt) FASTA in parallel (jobs > 1).
    let (meta_resume, was_new_resume) = store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new().jobs(4))
        .expect(
            "resume must succeed: a present collection with an intact .rgsi must be skipped \
             BEFORE the (now-corrupt) FASTA is opened/decoded",
        );

    assert!(
        !was_new_resume,
        "resume of an already-present collection must report was_new = false"
    );
    assert_eq!(
        meta_resume.digest, meta_first.digest,
        "resume must report the same collection digest from the cached .rgsi"
    );

    // Store must be byte-for-byte unchanged by the resume.
    let collections_after: Vec<DigestKey> = {
        let mut v: Vec<DigestKey> = store.collections.keys().cloned().collect();
        v.sort();
        v
    };
    let sequences_after: Vec<DigestKey> = {
        let mut v: Vec<DigestKey> = store.sequence_store.keys().cloned().collect();
        v.sort();
        v
    };
    assert_eq!(
        collections_before, collections_after,
        "resume must not change the set of collections"
    );
    assert_eq!(
        sequences_before, sequences_after,
        "resume must not change the set of sequences"
    );
}

// =========================================================================
// Many-tiny-sequence throughput (writer-pool) tests
// =========================================================================

/// Generate a FASTA with `n` short, DISTINCT sequences. Each record encodes its
/// index in base-4 as DNA so every sequence has a unique content digest (no
/// dedup), mimicking a transcriptome with hundreds of thousands of tiny records
/// -- the workload the writer pool is meant to accelerate (the old single
/// inserter did a File::create + create_dir_all per record).
fn write_many_tiny_fasta(path: &std::path::Path, n: usize) {
    use std::io::Write;
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let f = fs::File::create(path).unwrap();
    let mut w = std::io::BufWriter::new(f);
    for i in 0..n {
        // 16-base sequence uniquely encoding `i` (4^16 ≈ 4.3e9 >> n), padded so
        // every record is distinct and a constant small size.
        let mut seq = [b'A'; 16];
        let mut v = i;
        let mut p = 15usize;
        while v > 0 {
            seq[p] = BASES[v & 0b11];
            v >>= 2;
            p = p.wrapping_sub(1);
        }
        writeln!(w, ">seq{}", i).unwrap();
        w.write_all(&seq).unwrap();
        w.write_all(b"\n").unwrap();
    }
    w.flush().unwrap();
}

/// Throughput regression test for the parallel `.seq` writer pool.
///
/// Builds an on-disk store from a synthetic FASTA with MANY tiny sequences. The
/// pre-writer-pool code funneled every per-sequence `File::create` +
/// `create_dir_all` through the single inserter thread; this test asserts the
/// build finishes well under a generous bound (so a regression that re-serializes
/// the writes -- or reintroduces a per-sequence `create_dir_all` -- shows up),
/// while staying fast and bounded (no benchmark loop, no giant dataset).
#[test]
fn test_many_tiny_sequences_throughput() {
    let work = tempdir().unwrap();
    let fasta = work.path().join("tiny.fa");
    // ~100k distinct tiny records: enough to expose a single-thread file-create
    // bottleneck, small enough to build in a few seconds.
    let n = 100_000usize;
    write_many_tiny_fasta(&fasta, n);

    let dir = tempdir().unwrap();
    let start = std::time::Instant::now();
    let mut store = RefgetStore::on_disk(dir.path()).unwrap();
    store.set_quiet(true);
    let (meta, was_new) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new().jobs(4))
        .unwrap();
    let elapsed = start.elapsed();

    assert!(was_new, "collection must be new");
    assert_eq!(
        meta.n_sequences, n,
        "all {} distinct sequences must be registered",
        n
    );

    // Every distinct sequence must have produced a .seq file on disk (the writer
    // pool drained at the end-of-run barrier before indexes were written).
    let seqs = collect_seq_files(dir.path());
    assert_eq!(
        seqs.len(),
        n,
        "every distinct tiny sequence must have a .seq file on disk after the barrier"
    );

    // Generous wall-clock bound: this is NOT a benchmark, just a guard against a
    // re-serialized (single-thread file-create) regression. With the writer pool
    // this completes in a couple seconds on CI; 60s leaves ample slack.
    assert!(
        elapsed.as_secs() < 60,
        "building {} tiny sequences took {:?}; expected well under the single-inserter bound",
        n,
        elapsed
    );

    eprintln!(
        "throughput: built {} tiny seqs in {:?} ({:.0} seqs/sec)",
        n,
        elapsed,
        n as f64 / elapsed.as_secs_f64()
    );
}

// =========================================================================
// stream_sequence tests
// =========================================================================

fn first_seq_digest(store: &RefgetStore) -> String {
    store
        .sequence_store
        .values()
        .next()
        .unwrap()
        .metadata()
        .sha512t24u
        .clone()
}

fn first_seq_length(store: &RefgetStore) -> usize {
    store
        .sequence_store
        .values()
        .next()
        .unwrap()
        .metadata()
        .length
}

fn build_on_disk_store_streaming(mode: StorageMode) -> (tempfile::TempDir, RefgetStore, String, usize) {
    let dir = tempdir().unwrap();
    let fasta = dir.path().join("test.fa");
    fs::write(&fasta, ">chr1\nACGTACGTACGTACGTACGT\n").unwrap();
    let store_path = dir.path().join("store");
    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    store.set_encoding_mode(mode);
    store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    let digest = first_seq_digest(&store);
    let length = first_seq_length(&store);
    (dir, store, digest, length)
}

#[test]
fn test_stream_sequence_local_full() {
    use std::io::Read;
    let (_dir, store, digest, length) = build_on_disk_store_streaming(StorageMode::Encoded);
    let mut reader = store.stream_sequence(&digest, None, None).unwrap();
    let mut out = Vec::new();
    reader.read_to_end(&mut out).unwrap();
    assert_eq!(out.len(), length);
    assert_eq!(&out[..], b"ACGTACGTACGTACGTACGT");
}

#[test]
fn test_stream_sequence_local_substring() {
    use std::io::Read;
    let (_dir, mut store, digest, length) = build_on_disk_store_streaming(StorageMode::Encoded);
    store.load_sequence(&digest).unwrap();
    let ranges: Vec<(u64, u64)> = vec![
        (0, 4),
        (1, 5),
        (2, 10),
        (5, 6),
        (0, length as u64),
        (length as u64 - 1, length as u64),
    ];
    for (s, e) in ranges {
        let mut reader = store.stream_sequence(&digest, Some(s), Some(e)).unwrap();
        let mut streamed = Vec::new();
        reader.read_to_end(&mut streamed).unwrap();

        let expected = store.get_substring(&digest, s as usize, e as usize).unwrap();
        assert_eq!(
            String::from_utf8(streamed).unwrap(),
            expected,
            "mismatch for range {}..{}",
            s,
            e
        );
    }
}

// =========================================================================
// Partial-read fd-cache eviction test
// =========================================================================

/// Tiny deterministic xorshift RNG so the test is reproducible without adding
/// a `rand` dependency.
struct XorShift(u64);
impl XorShift {
    fn new(seed: u64) -> Self {
        XorShift(seed | 1)
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn range(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

/// Build a FASTA with several distinct, multi-base sequences so byte windows
/// are non-trivial and span multiple packed bytes in Encoded mode.
fn multi_seq_fasta() -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut out = String::new();
    for (i, len) in [400usize, 333, 512, 271, 600].iter().enumerate() {
        out.push_str(&format!(">seq{}\n", i));
        let mut s = String::with_capacity(*len);
        // deterministic but distinct per-sequence content
        for j in 0..*len {
            s.push(bases[(i * 7 + j * 3 + (j / 5)) % 4] as char);
        }
        out.push_str(&s);
        out.push('\n');
    }
    out
}

/// With a TINY fd-cache cap (2) and reads interleaved across >2 sequences,
/// every partial-read result must be byte-identical to the resident decode
/// path. Interleaving across 5 sequences with cap 2 forces constant
/// eviction+reopen, exercising the LRU close/reopen logic. Run for both
/// Encoded and Raw storage modes.
fn run_fd_cache_eviction_for_mode(raw: bool) {
    let temp_dir = tempdir().unwrap();
    let fasta_content = multi_seq_fasta();
    let fasta = temp_dir.path().join("multi.fa");
    fs::write(&fasta, &fasta_content).unwrap();
    let store_path = temp_dir.path().join("store");

    // Build the store on disk.
    let mut builder = RefgetStore::on_disk(&store_path).unwrap();
    if raw {
        builder.disable_encoding();
    }
    let (coll_meta, _) = builder
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    let coll_digest = coll_meta.digest.clone();
    drop(builder);

    // Reference: a resident store with everything loaded (whole-sequence decode).
    let mut resident = RefgetStore::open_local(&store_path).unwrap();
    resident.set_quiet(true);
    resident.load_all_collections().unwrap();
    resident.load_all_sequences().unwrap();

    // Subject: a stub-only store driving the partial-read path with cap=2.
    let mut partial = RefgetStore::open_local(&store_path).unwrap();
    partial.set_quiet(true);
    partial.load_all_collections().unwrap();
    partial.set_seq_fd_cache_cap_for_test(2);

    // Collect (digest, length) for each of the >2 sequences in the collection.
    let names: Vec<String> = (0..5).map(|i| format!("seq{}", i)).collect();
    let mut seqs: Vec<(String, usize)> = Vec::new();
    for name in &names {
        let rec = resident.get_sequence_by_name(&coll_digest, name).unwrap();
        seqs.push((rec.metadata().sha512t24u.clone(), rec.metadata().length));
    }
    assert!(seqs.len() > 2, "need >2 sequences to force eviction with cap 2");

    let mut rng = XorShift::new(0xC0FFEE_u64 ^ (raw as u64));
    for _ in 0..2000 {
        let (digest, len) = &seqs[rng.range(seqs.len())];
        let len = *len;
        let a = rng.range(len);
        let b = rng.range(len);
        let (start, end) = if a == b {
            (a.min(len - 1), a.min(len - 1) + 1)
        } else {
            (a.min(b), a.max(b))
        };
        let expected = resident.get_substring(digest, start, end).unwrap();
        let got = partial.get_substring(digest, start, end).unwrap();
        assert_eq!(
            got, expected,
            "partial-read mismatch (raw={}) digest={} [{}, {})",
            raw, digest, start, end
        );
    }
}

#[test]
fn test_stream_sequence_raw_mode() {
    use std::io::Read;
    let (_dir, store, digest, length) = build_on_disk_store_streaming(StorageMode::Raw);
    let mut reader = store.stream_sequence(&digest, Some(2), Some(8)).unwrap();
    let mut out = Vec::new();
    reader.read_to_end(&mut out).unwrap();
    assert_eq!(&out[..], b"GTACGT");
    assert!(length >= 8);
}

#[test]
fn test_stream_sequence_u64_bounds_compile_and_match() {
    // Compile-level guarantee: stream_sequence accepts u64 bounds.
    // Sequences > u32::MAX (>4 Gb) cannot be fixtured in unit tests, but
    // the arithmetic path (start_bit/end_bit/byte_start/byte_end/bases_to_emit)
    // is u64 throughout, so this test exercises the widened API surface.
    use std::io::Read;
    let (_dir, store, digest, _length) = build_on_disk_store_streaming(StorageMode::Encoded);
    let s: u64 = 0;
    let e: u64 = 8;
    let mut reader = store.stream_sequence(&digest, Some(s), Some(e)).unwrap();
    let mut out = Vec::new();
    reader.read_to_end(&mut out).unwrap();
    assert_eq!(out.len(), 8);
}

#[test]
fn test_stream_sequence_zero_length() {
    use std::io::Read;
    let (_dir, store, digest, _) = build_on_disk_store_streaming(StorageMode::Encoded);
    let mut reader = store.stream_sequence(&digest, Some(5), Some(5)).unwrap();
    let mut out = Vec::new();
    reader.read_to_end(&mut out).unwrap();
    assert!(out.is_empty());
}

#[test]
fn test_stream_sequence_invalid_range() {
    let (_dir, store, digest, length) = build_on_disk_store_streaming(StorageMode::Encoded);
    assert!(store.stream_sequence(&digest, Some(5), Some(3)).is_err());
    assert!(
        store
            .stream_sequence(&digest, Some(0), Some(length as u64 + 1))
            .is_err()
    );
}

#[test]
fn test_stream_sequence_bounded_memory() {
    use std::io::Read;
    let dir = tempdir().unwrap();
    let fasta = dir.path().join("big.fa");
    let mut content = String::from(">chr1\n");
    let bases = [b'A', b'C', b'G', b'T'];
    let mut seq = Vec::with_capacity(1_000_000);
    for i in 0..1_000_000 {
        seq.push(bases[i % 4]);
    }
    content.push_str(std::str::from_utf8(&seq).unwrap());
    content.push('\n');
    fs::write(&fasta, &content).unwrap();

    let store_path = dir.path().join("store");
    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    let digest = first_seq_digest(&store);

    let mut reader = store.stream_sequence(&digest, None, None).unwrap();
    let mut buf = [0u8; 4096];
    let mut collected = Vec::with_capacity(1_000_000);
    loop {
        let n = reader.read(&mut buf).unwrap();
        if n == 0 {
            break;
        }
        collected.extend_from_slice(&buf[..n]);
    }
    assert_eq!(collected, seq);
}

// =========================================================================
// Bounded-memory tests for stream_sequence
// =========================================================================
//
// These tests verify the O(1) memory claim of `stream_sequence`: regardless
// of sequence size, neither the Rust decoder nor the napi binding should
// materialize the full (sub)sequence in memory.

/// Install peak_alloc as the global allocator (tests-only) so we can measure
/// allocation deltas.
#[global_allocator]
static PEAK_ALLOC: peak_alloc::PeakAlloc = peak_alloc::PeakAlloc;

/// Serialize allocator-sensitive tests so concurrent tests don't pollute the
/// process-wide peak counter.
static ALLOC_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

// NOTE: `#[ignore]` because this measures the *process-global* `PEAK_ALLOC`
// high-water mark, which is polluted by any other test allocating concurrently.
// `ALLOC_TEST_LOCK` only serializes the two allocator tests against each other,
// not against the rest of the suite, so under default parallel `cargo test` the
// peak is unreliable. Run the allocation-profiling guards explicitly and serially:
//   cargo test -p gtars-refget -- --ignored --test-threads=1
#[test]
#[ignore = "allocation profiling: process-global PeakAlloc; run with --ignored --test-threads=1"]
fn test_stream_sequence_bounded_memory_full_record() {
    use std::io::Read;

    let _guard = ALLOC_TEST_LOCK.lock().unwrap_or_else(|p| p.into_inner());

    // Build a large FASTA (~1M bases) and import into an in-memory store in
    // Encoded mode. The in-memory store holds the encoded sequence as a
    // `SequenceRecord::Full { sequence: Vec<u8>, .. }` — for a 2-bit alphabet
    // that's ~250 KB encoded, decoding to ~1 MB raw bases.
    const SEQ_LEN: usize = 1_000_000;
    let dir = tempdir().unwrap();
    let fasta = dir.path().join("big.fa");
    let bases = [b'A', b'C', b'G', b'T'];
    let mut seq = Vec::with_capacity(SEQ_LEN);
    for i in 0..SEQ_LEN {
        seq.push(bases[i % 4]);
    }
    let mut content = String::from(">chr1\n");
    content.push_str(std::str::from_utf8(&seq).unwrap());
    content.push('\n');
    fs::write(&fasta, &content).unwrap();

    let mut store = RefgetStore::in_memory();
    // Use Raw mode so the Full record holds SEQ_LEN bytes verbatim. The old
    // buggy code path cloned that entire buffer during streaming; with the
    // Arc-backed reader it must not.
    store.set_encoding_mode(StorageMode::Raw);
    store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    let digest = first_seq_digest(&store);

    // Confirm the record is Full (in-memory resident).
    {
        use crate::hashkeyable::HashKeyable;
        let key = digest.to_key();
        let rec = store.sequence_store.get(&key).unwrap();
        assert!(rec.is_loaded(), "expected Full record for in-memory store");
    }

    // Measure peak allocation across stream_sequence + drain. The whole
    // buffer is already resident; streaming should not allocate another
    // copy of it.
    PEAK_ALLOC.reset_peak_usage();
    let baseline = PEAK_ALLOC.current_usage();

    {
        let mut reader = store.stream_sequence(&digest, None, None).unwrap();
        let mut chunk = [0u8; 4096];
        let mut total = 0usize;
        loop {
            let n = reader.read(&mut chunk).unwrap();
            if n == 0 {
                break;
            }
            total += n;
        }
        assert_eq!(total, SEQ_LEN);
    }

    let peak = PEAK_ALLOC.peak_usage();
    let peak_delta = peak.saturating_sub(baseline);

    // Streaming must not allocate anywhere near the sequence size. The
    // previous bug cloned the entire byte range (1 MB for this test); with
    // the fix, allocations are bounded by a few reader frames. The bound
    // below is loose enough to tolerate concurrent test allocator noise
    // but tight enough to catch a full-buffer clone (which would be
    // >= SEQ_LEN bytes).
    let bound: usize = SEQ_LEN / 4; // 250 KiB; a clone would cost >= 1 MiB.
    assert!(
        peak_delta < bound,
        "stream_sequence peak allocation delta {} bytes exceeds bound {} bytes \
         (SEQ_LEN={} bytes raw). \
         Streaming is cloning the buffer instead of sharing ownership.",
        peak_delta,
        bound,
        SEQ_LEN,
    );
}

// See `test_stream_sequence_bounded_memory_full_record` for why this is
// `#[ignore]`d (process-global PeakAlloc). Run via --ignored --test-threads=1.
#[test]
#[ignore = "allocation profiling: process-global PeakAlloc; run with --ignored --test-threads=1"]
fn test_stream_sequence_bounded_memory_stub_record() {
    // Regression guard for the production path: the Node proxy opens a store
    // via `open_local` / `open_remote`, so records are `Stub` (metadata only,
    // never loaded into memory) and `stream_sequence` must fall back to the
    // on-disk `.seq` file (or HTTP Range GET) without buffering the whole
    // file. This test directly measures peak allocation on that Stub-backed
    // path — analogous to `test_stream_sequence_bounded_memory_full_record`
    // but for the disk-streaming branch.
    use std::io::Read;
    use crate::hashkeyable::HashKeyable;

    let _guard = ALLOC_TEST_LOCK.lock().unwrap_or_else(|p| p.into_inner());

    // Measure peak at two very different sizes; the bound must not depend
    // on SEQ_LEN. (A full-file buffer would scale linearly; a fixed-size
    // BufReader would stay constant.) We drain into a fixed-size scratch
    // buffer — not `read_to_end` — because `read_to_end`'s geometric Vec
    // growth introduces transient realloc peaks of ~0.5 * SEQ_LEN that
    // swamp the actual streaming memory and falsely suggest scaling.
    for &SEQ_LEN in &[1_000_000usize, 16_000_000] {
    let dir = tempdir().unwrap();
    let fasta = dir.path().join("big.fa");
    let bases = [b'A', b'C', b'G', b'T'];
    let mut seq = Vec::with_capacity(SEQ_LEN);
    for i in 0..SEQ_LEN {
        seq.push(bases[i % 4]);
    }
    let mut content = String::from(">chr1\n");
    content.push_str(std::str::from_utf8(&seq).unwrap());
    content.push('\n');
    fs::write(&fasta, &content).unwrap();

    let store_path = dir.path().join("store");
    let digest;
    {
        let mut builder = RefgetStore::on_disk(&store_path).unwrap();
        builder.set_encoding_mode(StorageMode::Raw);
        builder
            .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
            .unwrap();
        digest = first_seq_digest(&builder);
        // Drop `builder` here — `add_sequence_collection_from_fasta` populates
        // an in-memory Full record; we want a fresh open so records are Stub.
    }

    // Reopen the store fresh — records should be Stub (metadata only).
    let store = RefgetStore::open_local(&store_path).unwrap();

    // Confirm the record is a Stub going in.
    {
        let key = digest.to_key();
        let rec = store.inner.sequence_store.get(&key).unwrap();
        assert!(
            !rec.is_loaded(),
            "expected Stub record after open_local, got Full"
        );
    }

    // Reset peak AFTER opening the store so we don't count store-open
    // allocations.
    PEAK_ALLOC.reset_peak_usage();
    let baseline = PEAK_ALLOC.current_usage();

    let mut reader = store.stream_sequence(&digest, None, None).unwrap();
    // Drain into a fixed-size buffer, checksumming as we go.
    // Avoids `read_to_end`'s geometric Vec growth, whose realloc
    // doubles transiently dominate peak (giving a false ~0.5N signal).
    let mut scratch = [0u8; 8192];
    let mut total = 0usize;
    let mut checksum: u64 = 0;
    loop {
        let n = reader.read(&mut scratch).unwrap();
        if n == 0 { break; }
        for &b in &scratch[..n] { checksum = checksum.wrapping_add(b as u64); }
        total += n;
    }

    let peak = PEAK_ALLOC.peak_usage();
    let peak_delta = peak.saturating_sub(baseline);

    assert_eq!(total, SEQ_LEN);
    let expected_checksum: u64 = seq.iter().map(|&b| b as u64).sum();
    assert_eq!(checksum, expected_checksum);

    // A full-file buffer in the Stub fallback would give peak_delta
    // proportional to SEQ_LEN. A BufReader-based stream is constant
    // (~8 KiB default). 64 KiB leaves headroom for allocator noise.
    assert!(
        peak_delta < 64 * 1024,
        "stream_sequence (Stub path) peak allocation {} bytes exceeds \
         64 KiB at SEQ_LEN={} — streaming is not O(1) in sequence size",
        peak_delta, SEQ_LEN,
    );
    }
}

#[test]
fn test_stream_sequence_no_preload_for_stub_record() {
    // Structural assertion for fix B: streaming must work when the record is
    // a Stub (not in memory), pulling bytes directly from the backing store,
    // without populating the in-memory cache. This lets the napi layer drop
    // its pre-`load_sequence` call.
    use std::io::Read;
    use crate::hashkeyable::HashKeyable;

    let (_dir, store, digest, length) = build_on_disk_store_streaming(StorageMode::Encoded);

    // Drop the in-memory sequence buffer back to Stub, preserving metadata.
    let mut store = store;
    {
        let key = digest.to_key();
        let rec = store.inner.sequence_store.get(&key).unwrap().clone();
        let meta = rec.metadata().clone();
        store
            .inner
            .sequence_store
            .insert(key, SequenceRecord::Stub(meta));
    }

    // Confirm the record is a Stub going in.
    {
        let key = digest.to_key();
        let rec = store.inner.sequence_store.get(&key).unwrap();
        assert!(!rec.is_loaded(), "expected Stub before stream_sequence");
    }

    // Stream without calling load_sequence() first.
    let mut reader = store.stream_sequence(&digest, None, None).unwrap();
    let mut out = Vec::new();
    reader.read_to_end(&mut out).unwrap();
    assert_eq!(out.len(), length);
    assert_eq!(&out[..], b"ACGTACGTACGTACGTACGT");

    // And the store should still be a Stub (no pre-loading side effect).
    let key = digest.to_key();
    let rec = store.inner.sequence_store.get(&key).unwrap();
    assert!(
        !rec.is_loaded(),
        "stream_sequence must not populate the in-memory Full record"
    );
}

// =========================================================================
// Lazy sequence index loading tests
// =========================================================================

#[test]
fn test_list_collections_without_sequence_index() {
    // Simulate a remote store where the sequence index is not yet loaded.
    // list_collections should work without triggering the sequence index download.
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    // Create a store with a collection
    let fasta_content = ">chr1\nATGCATGC\n>chr2\nGGGGAAAA\n";
    let fasta_path = dir.path().join("test.fa");
    fs::write(&fasta_path, fasta_content).unwrap();

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    store.add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new()).unwrap();

    // Re-open and clear the sequence store + mark index as not loaded
    // (simulating what open_remote now does)
    let mut store = RefgetStore::open_local(&store_path).unwrap();
    store.inner.sequence_store.clear();
    store.inner.md5_lookup.clear();
    store.inner.sequence_index_loaded = false;
    store.inner.sequence_index_path = Some("sequences.rgsi".to_string());

    // list_collections should work — it only touches collection metadata
    let result = store.list_collections(0, 100, &[]).unwrap();
    assert_eq!(result.results.len(), 1, "Should find 1 collection");
    assert!(!store.inner.sequence_index_loaded, "Sequence index should NOT be loaded after list_collections");
}

#[test]
fn test_load_all_sequences_triggers_index_load() {
    // When calling load_all_sequences on a store with deferred index,
    // it should lazily load the sequence index first.
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");

    let fasta_content = ">chr1\nATGCATGC\n";
    let fasta_path = dir.path().join("test.fa");
    fs::write(&fasta_path, fasta_content).unwrap();

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    store.add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new()).unwrap();

    // Simulate deferred state
    let mut store = RefgetStore::open_local(&store_path).unwrap();
    store.inner.sequence_store.clear();
    store.inner.md5_lookup.clear();
    store.inner.sequence_index_loaded = false;
    store.inner.sequence_index_path = Some("sequences.rgsi".to_string());

    assert!(!store.inner.sequence_index_loaded);

    // load_all_sequences should trigger the index load
    store.load_all_sequences().unwrap();
    assert!(store.inner.sequence_index_loaded, "Sequence index should be loaded after load_all_sequences");
    assert!(!store.inner.sequence_store.is_empty(), "Sequences should be populated");
}

#[test]
fn test_fd_cache_eviction_matches_resident_encoded() {
    run_fd_cache_eviction_for_mode(false);
}

#[test]
fn test_fd_cache_eviction_matches_resident_raw() {
    run_fd_cache_eviction_for_mode(true);
}

/// `get_substrings` must produce output byte-identical to (a) a loop of
/// `get_substring` on the same disk-backed store and (b) the resident decode,
/// for BOTH branches of the adaptive heuristic:
///   - sparse small ranges (total coverage < WHOLE_DECODE_FRACTION) -> partial branch
///   - wide/overlapping ranges (total coverage >= WHOLE_DECODE_FRACTION) -> whole-decode branch
/// Run for both Encoded and Raw storage modes.
fn run_get_substrings_matches_for_mode(raw: bool) {
    let temp_dir = tempdir().unwrap();
    let fasta_content = multi_seq_fasta();
    let fasta = temp_dir.path().join("multi.fa");
    fs::write(&fasta, &fasta_content).unwrap();
    let store_path = temp_dir.path().join("store");

    let mut builder = RefgetStore::on_disk(&store_path).unwrap();
    if raw {
        builder.disable_encoding();
    }
    let (coll_meta, _) = builder
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    let coll_digest = coll_meta.digest.clone();
    drop(builder);

    // Reference: resident store (whole-sequence decode held in memory).
    let mut resident = RefgetStore::open_local(&store_path).unwrap();
    resident.set_quiet(true);
    resident.load_all_collections().unwrap();
    resident.load_all_sequences().unwrap();

    // Subject: stub-only disk-backed store.
    let mut disk = RefgetStore::open_local(&store_path).unwrap();
    disk.set_quiet(true);
    disk.load_all_collections().unwrap();

    // Use seq0 (length 400) for the per-sequence range sets.
    let rec = resident.get_sequence_by_name(&coll_digest, "seq0").unwrap();
    let digest = rec.metadata().sha512t24u.clone();
    let len = rec.metadata().length;
    assert_eq!(len, 400);

    // BRANCH A (partial): many sparse small ranges; total coverage tiny relative
    // to len (well below the 0.125 threshold).
    let sparse: Vec<(usize, usize)> = vec![
        (0, 3),
        (10, 14),
        (50, 51),
        (100, 105),
        (200, 202),
        (399, 400),
    ];
    let sparse_total: usize = sparse.iter().map(|&(s, e)| e - s).sum();
    assert!((sparse_total as f64) < 0.125 * (len as f64), "branch A must be partial");

    // BRANCH B (whole-decode): wide and overlapping ranges whose total coverage
    // far exceeds the threshold.
    let wide: Vec<(usize, usize)> = vec![
        (0, 200),
        (100, 350),
        (50, 400),
        (10, 390),
    ];
    let wide_total: usize = wide.iter().map(|&(s, e)| e - s).sum();
    assert!((wide_total as f64) >= 0.125 * (len as f64), "branch B must be whole-decode");

    for ranges in [&sparse, &wide] {
        // Expected from resident decode.
        let expected_resident: Vec<String> = ranges
            .iter()
            .map(|&(s, e)| resident.get_substring(&digest, s, e).unwrap())
            .collect();
        // Loop of disk get_substring.
        let expected_disk_loop: Vec<String> = ranges
            .iter()
            .map(|&(s, e)| disk.get_substring(&digest, s, e).unwrap())
            .collect();
        // Batch.
        let batch = disk.get_substrings(&digest, ranges).unwrap();

        assert_eq!(
            batch, expected_resident,
            "batch != resident (raw={}) ranges={:?}",
            raw, ranges
        );
        assert_eq!(
            batch, expected_disk_loop,
            "batch != disk loop (raw={}) ranges={:?}",
            raw, ranges
        );
    }

    // Resident store batch path (Full arm) must also match.
    let batch_resident = resident.get_substrings(&digest, &wide).unwrap();
    let expected: Vec<String> = wide
        .iter()
        .map(|&(s, e)| resident.get_substring(&digest, s, e).unwrap())
        .collect();
    assert_eq!(batch_resident, expected, "resident Full-arm batch mismatch (raw={})", raw);

    // Invalid range must error.
    assert!(disk.get_substrings(&digest, &[(0, len + 1)]).is_err());
}

#[test]
fn test_get_substrings_matches_encoded() {
    run_get_substrings_matches_for_mode(false);
}

#[test]
fn test_get_substrings_matches_raw() {
    run_get_substrings_matches_for_mode(true);
}

// =========================================================================
// Import-time collection alias tests
// =========================================================================

/// REQUIRED REGRESSION GUARD: omitting `.collection_alias(..)` must change
/// nothing. An import with plain options registers no collection alias at all.
#[test]
fn test_import_without_collection_alias_registers_none() {
    let dir = tempdir().unwrap();
    let fasta = copy_test_fasta(dir.path(), "base.fa");

    let mut store = RefgetStore::in_memory();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();

    assert!(
        store.list_collection_alias_namespaces().is_empty(),
        "plain import must not create any collection alias namespace"
    );
    assert!(
        store.get_aliases_for_collection(&meta.digest).is_empty(),
        "plain import must not name the collection"
    );
}

/// `namespaces` feeds the SEQUENCE index only; it must never name the
/// collection. Proves the two alias indexes stay decoupled.
#[test]
fn test_import_without_collection_alias_but_with_namespaces() {
    let dir = tempdir().unwrap();
    let fasta = dir.path().join("ns.fa");
    fs::write(&fasta, ">chr1 ucsc:chr1\nAAAACCCC\n>chr2 ucsc:chr2\nGGGGTTTT\n").unwrap();

    let mut store = RefgetStore::in_memory();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new().namespaces(&["ucsc"]))
        .unwrap();

    let seq_aliases = store
        .list_sequence_aliases("ucsc")
        .expect("ucsc sequence alias namespace should exist");
    assert!(!seq_aliases.is_empty(), "header aliases should be registered");

    assert!(
        store.list_collection_alias_namespaces().is_empty(),
        "namespaces must not leak into the collection alias index"
    );
    assert!(store.get_aliases_for_collection(&meta.digest).is_empty());
}

#[test]
fn test_import_with_collection_alias_registers_collection_alias() {
    let dir = tempdir().unwrap();
    let fasta = copy_test_fasta(dir.path(), "base.fa");

    let mut store = RefgetStore::in_memory();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(
            &fasta,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .unwrap();

    let resolved = store
        .get_collection_metadata_by_alias("ucsc", "hg38")
        .expect("ucsc:hg38 should resolve");
    assert_eq!(resolved.digest, meta.digest);

    assert!(
        store.list_sequence_alias_namespaces().is_empty(),
        "collection alias must not leak into the sequence alias index"
    );
}

/// The collection alias works with `namespaces` unset (decision 1: the option is
/// explicit, not derived from the namespaces list).
#[test]
fn test_collection_alias_without_namespaces() {
    let dir = tempdir().unwrap();
    let fasta = dir.path().join("bare.fa");
    // Bare headers: no `ns:value` tokens anywhere in the file.
    fs::write(&fasta, ">chr1\nAAAACCCC\n>chr2\nGGGGTTTT\n").unwrap();

    let mut store = RefgetStore::in_memory();
    let opts = FastaImportOptions::new().collection_alias("ucsc", "hg38");
    assert!(opts.namespaces.is_empty(), "namespaces deliberately unset");
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, opts)
        .unwrap();

    let resolved = store
        .get_collection_metadata_by_alias("ucsc", "hg38")
        .expect("ucsc:hg38 should resolve without namespaces set");
    assert_eq!(resolved.digest, meta.digest);
}

#[test]
fn test_collection_alias_persists_to_disk() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");
    let fasta = copy_test_fasta(dir.path(), "base.fa");

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(
            &fasta,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .unwrap();
    store.write().unwrap();

    let tsv = store_path.join("aliases").join("collections").join("ucsc.tsv");
    assert!(tsv.exists(), "collection alias TSV should be written");
    let contents = fs::read_to_string(&tsv).unwrap();
    assert!(
        contents.contains(&format!("hg38\t{}", meta.digest)),
        "TSV should map hg38 to the collection digest, got: {}",
        contents
    );

    // `store_metadata()` only surfaces the state digests, so read rgstore.json
    // directly for the advertised namespace list.
    let manifest_json: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(store_path.join("rgstore.json")).unwrap()).unwrap();
    let coll_ns = manifest_json
        .get("collection_alias_namespaces")
        .and_then(|v| v.as_array())
        .expect("rgstore.json should carry collection_alias_namespaces");
    assert!(
        coll_ns.iter().any(|v| v.as_str() == Some("ucsc")),
        "rgstore.json should advertise the ucsc collection alias namespace, got {:?}",
        coll_ns
    );
    assert!(
        manifest_json
            .get("sequence_alias_namespaces")
            .and_then(|v| v.as_array())
            .map(|a| a.is_empty())
            .unwrap_or(true),
        "no sequence alias namespaces should be advertised"
    );

    let manifest = store.store_metadata().unwrap();
    assert!(
        manifest.get("aliases_digest").is_some(),
        "rgstore.json should carry an aliases_digest"
    );

    // Round-trip through a fresh open.
    let reopened = RefgetStore::open_local(&store_path).unwrap();
    let resolved = reopened
        .get_collection_metadata_by_alias("ucsc", "hg38")
        .expect("alias should survive a reopen");
    assert_eq!(resolved.digest, meta.digest);
}

/// Re-importing the same FASTA with the same alias is a no-op, on both the
/// in-run duplicate path and the build-side `Skip` path (fresh store from disk).
#[test]
fn test_collection_alias_idempotent_on_reimport() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");
    let fasta = copy_test_fasta(dir.path(), "base.fa");

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    let (meta1, was_new1) = store
        .add_sequence_collection_from_fasta(
            &fasta,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .unwrap();
    assert!(was_new1);

    // In-run duplicate path: same store, same file again.
    let (meta2, was_new2) = store
        .add_sequence_collection_from_fasta(
            &fasta,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .unwrap();
    assert!(!was_new2, "second import of the same FASTA is not new");
    assert_eq!(meta2.digest, meta1.digest);
    assert_eq!(
        store
            .get_collection_metadata_by_alias("ucsc", "hg38")
            .unwrap()
            .digest,
        meta1.digest
    );
    store.write().unwrap();

    // Build-side Skip path: reopen from disk so the collection is in
    // `present_collections` before any decode happens.
    let mut reopened = RefgetStore::open_local(&store_path).unwrap();
    let (meta3, was_new3) = reopened
        .add_sequence_collection_from_fasta(
            &fasta,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .unwrap();
    assert!(!was_new3, "already-present collection should be skipped");
    assert_eq!(meta3.digest, meta1.digest);
    assert_eq!(
        reopened
            .get_collection_metadata_by_alias("ucsc", "hg38")
            .unwrap()
            .digest,
        meta1.digest
    );
}

/// A previously-imported collection that was NOT named still gets its name on a
/// later re-import via the build-side Skip path. This is the silent-drop bug.
#[test]
fn test_collection_alias_registered_on_build_side_skip() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");
    let fasta = copy_test_fasta(dir.path(), "base.fa");

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
        .unwrap();
    store.write().unwrap();
    assert!(store.list_collection_alias_namespaces().is_empty());

    // Reopen: the collection is already present, so the builder emits Skip and
    // finalize_collection is never reached. The alias must still land.
    let mut reopened = RefgetStore::open_local(&store_path).unwrap();
    let (_, was_new) = reopened
        .add_sequence_collection_from_fasta(
            &fasta,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .unwrap();
    assert!(!was_new, "collection was already present");
    assert_eq!(
        reopened
            .get_collection_metadata_by_alias("ucsc", "hg38")
            .expect("alias must be registered even on the skip path")
            .digest,
        meta.digest
    );
}

#[test]
fn test_collection_alias_conflict_errors() {
    let dir = tempdir().unwrap();
    let fasta_a = dir.path().join("a.fa");
    fs::write(&fasta_a, ">chr1\nAAAACCCC\n").unwrap();
    let fasta_b = dir.path().join("b.fa");
    fs::write(&fasta_b, ">chr1\nGGGGTTTT\n").unwrap();

    let mut store = RefgetStore::in_memory();
    let (meta_a, _) = store
        .add_sequence_collection_from_fasta(
            &fasta_a,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .unwrap();

    let err = store
        .add_sequence_collection_from_fasta(
            &fasta_b,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .expect_err("conflicting alias must error");
    let msg = format!("{}", err);
    assert!(msg.contains("ucsc"), "error should name the namespace: {}", msg);
    assert!(msg.contains("hg38"), "error should name the alias: {}", msg);

    // The failed import must not have remapped the alias.
    assert_eq!(
        store
            .get_collection_metadata_by_alias("ucsc", "hg38")
            .unwrap()
            .digest,
        meta_a.digest
    );
}

/// A failed conflicting-alias import must not leave the collection behind.
///
/// The alias conflict is validated up front, before any part of the import is
/// committed, so a rejected import is a true no-op: the store must hold exactly
/// the one collection that succeeded. Previously the collection was inserted
/// first and the alias checked afterwards, so a "failed" import still added a
/// second collection to `store.collections`.
#[test]
fn test_collection_alias_conflict_does_not_leak_collection() {
    let dir = tempdir().unwrap();
    let fasta_a = dir.path().join("a.fa");
    fs::write(&fasta_a, ">chr1\nAAAACCCC\n").unwrap();
    let fasta_b = dir.path().join("b.fa");
    fs::write(&fasta_b, ">chr1\nGGGGTTTT\n").unwrap();

    let mut store = RefgetStore::in_memory();
    let (meta_a, _) = store
        .add_sequence_collection_from_fasta(
            &fasta_a,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .unwrap();
    assert_eq!(store.collections.len(), 1);

    store
        .add_sequence_collection_from_fasta(
            &fasta_b,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .expect_err("conflicting alias must error");

    // The rejected import must not have been committed.
    assert_eq!(
        store.collections.len(),
        1,
        "a failed import must not leave its collection in the store"
    );
    assert!(
        store.collections.contains_key(&meta_a.digest.to_key()),
        "the surviving collection should be the one that imported successfully"
    );
    assert!(
        store.name_lookup.len() <= 1,
        "no name_lookup entry should be installed for the rejected collection"
    );
}

/// Disk-backed counterpart: a rejected conflicting-alias import must not leave
/// an orphaned `collections/<digest>.rgsi` that the top-level index never
/// references. The per-collection `.rgsi` is written immediately, but
/// `collections.rgci` is only written at the end of a successful import, so
/// erroring after the collection write used to strand a file on disk.
#[test]
fn test_collection_alias_conflict_leaves_no_orphan_on_disk() {
    let dir = tempdir().unwrap();
    let store_path = dir.path().join("store");
    let fasta_a = dir.path().join("a.fa");
    fs::write(&fasta_a, ">chr1\nAAAACCCC\n").unwrap();
    let fasta_b = dir.path().join("b.fa");
    fs::write(&fasta_b, ">chr1\nGGGGTTTT\n").unwrap();

    let mut store = RefgetStore::on_disk(&store_path).unwrap();
    let (meta_a, _) = store
        .add_sequence_collection_from_fasta(
            &fasta_a,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .unwrap();
    store.write().unwrap();

    store
        .add_sequence_collection_from_fasta(
            &fasta_b,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .expect_err("conflicting alias must error");

    // Exactly one per-collection index file, and it belongs to the collection
    // that actually imported.
    let coll_dir = store_path.join("collections");
    let mut rgsi: Vec<String> = fs::read_dir(&coll_dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".rgsi"))
        .collect();
    rgsi.sort();
    assert_eq!(
        rgsi,
        vec![format!("{}.rgsi", meta_a.digest)],
        "failed import must not leave an orphaned collection .rgsi on disk"
    );

    assert_eq!(store.collections.len(), 1);

    // The store on disk is still coherent: reopening sees exactly the one
    // collection, still reachable through the alias.
    let mut reopened = RefgetStore::open_local(&store_path).unwrap();
    assert_eq!(reopened.collections.len(), 1);
    assert_eq!(
        reopened
            .get_collection_metadata_by_alias("ucsc", "hg38")
            .expect("alias should still resolve after the failed import")
            .digest,
        meta_a.digest
    );
    // And the surviving collection is fully readable, not a dangling index row.
    reopened
        .get_collection(&meta_a.digest)
        .expect("the committed collection must still load");
}

#[test]
fn test_collection_alias_conflict_force_overwrites() {
    let dir = tempdir().unwrap();
    let fasta_a = dir.path().join("a.fa");
    fs::write(&fasta_a, ">chr1\nAAAACCCC\n").unwrap();
    let fasta_b = dir.path().join("b.fa");
    fs::write(&fasta_b, ">chr1\nGGGGTTTT\n").unwrap();

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(
            &fasta_a,
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .unwrap();

    let (meta_b, _) = store
        .add_sequence_collection_from_fasta(
            &fasta_b,
            FastaImportOptions::new()
                .collection_alias("ucsc", "hg38")
                .force(true),
        )
        .unwrap();

    assert_eq!(
        store
            .get_collection_metadata_by_alias("ucsc", "hg38")
            .unwrap()
            .digest,
        meta_b.digest,
        "force should remap the alias to the new collection"
    );
}

#[test]
fn test_collection_alias_rejected_for_multiple_files() {
    let dir = tempdir().unwrap();
    let fasta_a = dir.path().join("a.fa");
    fs::write(&fasta_a, ">chr1\nAAAACCCC\n").unwrap();
    let fasta_b = dir.path().join("b.fa");
    fs::write(&fasta_b, ">chr1\nGGGGTTTT\n").unwrap();

    let mut store = RefgetStore::in_memory();
    let err = store
        .add_sequence_collections_from_fastas(
            &[fasta_a, fasta_b],
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .expect_err("collection_alias with multiple files must error");
    let msg = format!("{}", err);
    assert!(
        msg.contains("add_collection_alias"),
        "error should teach the per-file workaround: {}",
        msg
    );

    // The guard runs before any thread spawns, so the store must be untouched.
    assert!(store.sequence_store.is_empty(), "no sequences imported");
    assert!(store.collections.is_empty(), "no collections imported");
    assert!(store.list_collection_alias_namespaces().is_empty());
    assert!(store.list_sequence_alias_namespaces().is_empty());
}

/// The multi-file guard must be a `> 1` check, not a `!= 1` check -- the
/// single-file wrapper goes through this same entry point with a 1-element slice.
#[test]
fn test_collection_alias_single_element_slice_ok() {
    let dir = tempdir().unwrap();
    let fasta = dir.path().join("one.fa");
    fs::write(&fasta, ">chr1\nAAAACCCC\n").unwrap();

    let mut store = RefgetStore::in_memory();
    let report = store
        .add_sequence_collections_from_fastas(
            &[fasta],
            FastaImportOptions::new().collection_alias("ucsc", "hg38"),
        )
        .expect("one-element slice with an alias must succeed");
    assert_eq!(report.collections.len(), 1);
    assert_eq!(
        store
            .get_collection_metadata_by_alias("ucsc", "hg38")
            .unwrap()
            .digest,
        report.collections[0].0.digest
    );
}

// =========================================================================
// Concurrent-writer safety
//
// These are the regression tests for the 2026-07-23 incident, where four
// concurrent `genome_init` jobs wrote one store directory and one genome's
// collection was written to disk but dropped from both indexes -- its alias
// stopped resolving and the nightly build failed. Neither writer errored.
// =========================================================================

/// Build a one-collection store at `path` from an inline FASTA.
fn add_collection_to_store(store: &mut RefgetStore, dir: &std::path::Path, name: &str, fasta: &str) -> String {
    let fasta_path = dir.join(format!("{}.fa", name));
    fs::write(&fasta_path, fasta).unwrap();
    let (meta, _) = store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
        .unwrap();
    meta.digest
}

/// THE incident, reproduced: two handles open the same store, each adds a
/// different collection, and they commit in an interleaved order. Before the
/// merge, whoever committed last rewrote the index from its own open-time
/// snapshot and the other's collection vanished.
#[test]
fn test_interleaved_writers_both_survive() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    // Seed the store so both writers have something to load at open time.
    let mut seed = RefgetStore::on_disk(&store_dir).unwrap();
    let seed_digest = add_collection_to_store(&mut seed, work.path(), "seed", ">chrS\nAAAACCCC\n");
    seed.write().unwrap();
    drop(seed);

    // Both writers snapshot the store at the same moment.
    let mut writer_a = RefgetStore::open_local(&store_dir).unwrap();
    let mut writer_b = RefgetStore::open_local(&store_dir).unwrap();

    // A stages its collection first...
    let digest_a = add_collection_to_store(&mut writer_a, work.path(), "a", ">chrA\nGGGGTTTT\n");
    // ...B stages and commits while A is still holding a stale snapshot...
    let digest_b = add_collection_to_store(&mut writer_b, work.path(), "b", ">chrB\nTTTTGGGG\n");
    writer_b.write().unwrap();
    // ...and only then does A commit, from the snapshot that predates B.
    writer_a.write().unwrap();
    drop(writer_a);
    drop(writer_b);

    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    for (label, digest) in [("seed", &seed_digest), ("A", &digest_a), ("B", &digest_b)] {
        assert!(
            reopened.get_collection_metadata(digest).is_some(),
            "collection {} ({}) was dropped from the index by the other writer",
            label,
            digest
        );
    }
}

/// The same shape at the sequence level: rows another writer published must not
/// be dropped by a commit from a stale snapshot.
#[test]
fn test_interleaved_writers_preserve_sequence_rows() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut seed = RefgetStore::on_disk(&store_dir).unwrap();
    add_collection_to_store(&mut seed, work.path(), "seed", ">chrS\nAAAACCCC\n");
    seed.write().unwrap();
    drop(seed);

    let mut writer_a = RefgetStore::open_local(&store_dir).unwrap();
    let mut writer_b = RefgetStore::open_local(&store_dir).unwrap();

    add_collection_to_store(&mut writer_a, work.path(), "a", ">chrA\nGGGGTTTT\n");
    add_collection_to_store(&mut writer_b, work.path(), "b", ">chrB\nTTTTGGGG\n");
    writer_b.write().unwrap();
    writer_a.write().unwrap();
    drop(writer_a);
    drop(writer_b);

    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    let names: std::collections::HashSet<String> = reopened
        .list_sequences()
        .iter()
        .map(|m| m.name.clone())
        .collect();
    for expected in ["chrS", "chrA", "chrB"] {
        assert!(names.contains(expected), "sequence {} was dropped; have {:?}", expected, names);
    }
}

/// An alias namespace one writer creates must not be dropped from the manifest
/// (or deleted from disk) by another writer that never loaded it. The manifest
/// is the only discovery mechanism over HTTP, so an unadvertised TSV is
/// unreachable.
#[test]
fn test_alias_namespace_from_other_writer_survives_commit() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut seed = RefgetStore::on_disk(&store_dir).unwrap();
    let seed_digest = add_collection_to_store(&mut seed, work.path(), "seed", ">chrS\nAAAACCCC\n");
    seed.write().unwrap();
    drop(seed);

    // A opens the store BEFORE the namespace exists, so it will never load it.
    let mut writer_a = RefgetStore::open_local(&store_dir).unwrap();

    let mut writer_b = RefgetStore::open_local(&store_dir).unwrap();
    writer_b
        .add_collection_alias("refgenie", "athaliana", &seed_digest)
        .unwrap();
    drop(writer_b);

    // A commits from its stale snapshot.
    add_collection_to_store(&mut writer_a, work.path(), "a", ">chrA\nGGGGTTTT\n");
    writer_a.write().unwrap();
    drop(writer_a);

    assert!(
        store_dir.join("aliases/collections/refgenie.tsv").exists(),
        "the other writer's alias TSV was deleted"
    );
    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    assert!(
        reopened
            .get_collection_metadata_by_alias("refgenie", "athaliana")
            .is_some(),
        "alias stopped resolving after another writer committed"
    );
}

/// Removal must still work: a blind union would resurrect whatever
/// `remove_collection` deleted.
#[test]
fn test_removal_is_not_resurrected_by_the_merge() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    let keep = add_collection_to_store(&mut store, work.path(), "keep", ">chrK\nAAAACCCC\n");
    let drop_me = add_collection_to_store(&mut store, work.path(), "drop", ">chrD\nGGGGTTTT\n");
    store.write().unwrap();

    assert!(store.remove_collection(&drop_me, true).unwrap());
    drop(store);

    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    assert!(reopened.get_collection_metadata(&keep).is_some());
    assert!(
        reopened.get_collection_metadata(&drop_me).is_none(),
        "removed collection came back through the merge"
    );
}

/// Orphan GC on a freshly-opened store. `open_local` loads collection STUBS and
/// never populates `name_lookup`, so the old `name_lookup`-derived live set was
/// empty here and unlinked sequences the surviving collection still needs.
/// This test fails against the pre-merge implementation.
#[test]
fn test_orphan_gc_keeps_sequences_shared_with_another_collection() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    // Two collections that SHARE chrShared, plus one sequence unique to each.
    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    let coll_a = add_collection_to_store(
        &mut store,
        work.path(),
        "a",
        ">chrShared\nAAAACCCCGGGGTTTT\n>chrOnlyA\nACACACAC\n",
    );
    let coll_b = add_collection_to_store(
        &mut store,
        work.path(),
        "b",
        ">chrShared\nAAAACCCCGGGGTTTT\n>chrOnlyB\nGTGTGTGT\n",
    );
    store.write().unwrap();
    drop(store);

    // Reopen fresh: stub load, so name_lookup is empty.
    let mut store = RefgetStore::open_local(&store_dir).unwrap();
    assert!(
        store.name_lookup.is_empty(),
        "precondition: a fresh open must not populate name_lookup"
    );

    // The dry-run must agree with what actually happens.
    let planned = store.plan_orphan_removal(&coll_a).unwrap();
    assert_eq!(
        planned.len(),
        1,
        "only chrOnlyA is an orphan; planned {:?}",
        planned
    );

    assert!(store.remove_collection(&coll_a, true).unwrap());
    drop(store);

    let mut reopened = RefgetStore::open_local(&store_dir).unwrap();
    reopened.load_collection(&coll_b).unwrap();
    let collection = reopened.get_collection(&coll_b).unwrap();
    assert_eq!(collection.sequences.len(), 2);
    for seq in &collection.sequences {
        let digest = &seq.metadata().sha512t24u;
        assert!(
            reopened.get_sequence_metadata(digest).is_some(),
            "sequence {} ({}) survived in the index but not in the store",
            seq.metadata().name,
            digest
        );
        assert!(
            reopened.get_sequence(digest).is_ok(),
            "sequence {} ({}) was unlinked from disk by the orphan GC",
            seq.metadata().name,
            digest
        );
    }
}

/// Orphan GC must fail closed. If a collection listed in the index has no
/// readable `.rgsi`, the live set is unknowable and nothing may be deleted --
/// treating a missing input as "references nothing" is how live data gets
/// unlinked.
#[test]
fn test_orphan_gc_refuses_when_a_collection_rgsi_is_missing() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    let coll_a = add_collection_to_store(&mut store, work.path(), "a", ">chrA\nAAAACCCC\n");
    let coll_b = add_collection_to_store(&mut store, work.path(), "b", ">chrB\nGGGGTTTT\n");
    store.write().unwrap();
    drop(store);

    fs::remove_file(store_dir.join(format!("collections/{}.rgsi", coll_b))).unwrap();

    let mut store = RefgetStore::open_local(&store_dir).unwrap();
    let err = store.remove_collection(&coll_a, true).unwrap_err();
    assert!(
        err.to_string().contains("refusing to remove orphan sequences"),
        "expected a fail-closed error, got: {}",
        err
    );
}

/// Conflicting alias bindings are a semantic conflict, not a merge detail.
/// Error by default; `--force-alias` takes the in-memory value.
#[test]
fn test_conflicting_alias_errors_unless_forced() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    let coll_a = add_collection_to_store(&mut store, work.path(), "a", ">chrA\nAAAACCCC\n");
    let coll_b = add_collection_to_store(&mut store, work.path(), "b", ">chrB\nGGGGTTTT\n");
    store.add_collection_alias("ucsc", "hg38", &coll_a).unwrap();
    store.write().unwrap();
    drop(store);

    // Another writer rebinds the same alias to a different collection.
    let mut other = RefgetStore::open_local(&store_dir).unwrap();
    let err = other
        .add_collection_alias("ucsc", "hg38", &coll_b)
        .expect_err("a conflicting alias binding must not be silently resolved");
    assert!(err.to_string().contains("alias conflict"), "got: {}", err);
    // The published binding is untouched by the failed commit.
    assert_eq!(
        fs::read_to_string(store_dir.join("aliases/collections/ucsc.tsv")).unwrap(),
        format!("hg38\t{}\n", coll_a)
    );

    other.set_force_alias(true);
    other
        .add_collection_alias("ucsc", "hg38", &coll_b)
        .expect("force-alias must overwrite the published binding");
    drop(other);

    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    assert_eq!(
        reopened
            .get_collection_metadata_by_alias("ucsc", "hg38")
            .unwrap()
            .digest,
        coll_b
    );
}

/// `rgstore.json` is published LAST, so the digests it advertises always
/// describe files that are already on disk.
#[test]
fn test_manifest_digests_match_published_indexes() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    add_collection_to_store(&mut store, work.path(), "a", ">chrA\nAAAACCCC\n");
    store.write().unwrap();

    let metadata = store.store_metadata().unwrap();
    let sha = |name: &str| {
        use sha2::{Digest, Sha256};
        format!("{:x}", Sha256::digest(fs::read(store_dir.join(name)).unwrap()))
    };
    assert_eq!(metadata.get("sequences_digest").unwrap(), &sha("sequences.rgsi"));
    assert_eq!(metadata.get("collections_digest").unwrap(), &sha("collections.rgci"));
}

/// `created_at` now means what it says. It used to be stamped with `now()` on
/// every commit, making it a duplicate of `modified`.
#[test]
fn test_created_at_survives_subsequent_commits() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    add_collection_to_store(&mut store, work.path(), "a", ">chrA\nAAAACCCC\n");
    store.write().unwrap();

    let read_created_at = || -> String {
        let json = fs::read_to_string(store_dir.join("rgstore.json")).unwrap();
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        v["created_at"].as_str().unwrap().to_string()
    };
    let first = read_created_at();

    add_collection_to_store(&mut store, work.path(), "b", ">chrB\nGGGGTTTT\n");
    store.write().unwrap();

    assert_eq!(read_created_at(), first, "created_at was overwritten by a later commit");
}

/// The lock serializes writers, and a store that already holds a batch lock
/// must not deadlock against the commit each mutation triggers.
#[test]
fn test_batch_lock_is_reentrant_and_exclusive() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    add_collection_to_store(&mut store, work.path(), "a", ">chrA\nAAAACCCC\n");
    store.write().unwrap();

    store.lock_for_batch("test-batch").unwrap();
    assert!(store.holds_batch_lock());

    // Nested commits must not block on the lock this store already holds.
    add_collection_to_store(&mut store, work.path(), "b", ">chrB\nGGGGTTTT\n");
    store.write().unwrap();

    // Meanwhile another process would be locked out.
    assert!(
        super::lock_status(&store_dir).unwrap().is_some(),
        "the batch lock should be visible on disk"
    );

    store.release_batch_lock();
    assert!(super::lock_status(&store_dir).unwrap().is_none());
}

/// Transient lock/temp files must be recognizable so anything mirroring a store
/// directory (`aws s3 sync`, integrity checkers) can skip them.
#[test]
fn test_commit_leaves_no_transient_files_behind() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    add_collection_to_store(&mut store, work.path(), "a", ">chrA\nAAAACCCC\n");
    store.write().unwrap();
    drop(store);

    let leftovers: Vec<String> = fs::read_dir(&store_dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter_map(|e| e.file_name().to_str().map(String::from))
        .filter(|n| super::is_transient_store_file(n))
        .collect();
    assert!(leftovers.is_empty(), "transient files left behind: {:?}", leftovers);
}

// =========================================================================
// Delta-at-commit: a stale handle must not resurrect another writer's removal
//
// These are the regression tests for the 2026-07-28 finding. Merge-at-commit
// wrote back EVERYTHING a handle had in memory, and `open_local` loads a stub
// for every collection in the index purely so it can be read. A handle that
// opened before a removal therefore re-added the removed rows on its next
// commit -- after the `.seq` files were already unlinked. Not a lock race: the
// removal completes entirely before the offending commit begins.
// =========================================================================

/// Sequence digests of one collection, read from its on-disk `.rgsi`.
fn sequence_digests_of(store: &mut RefgetStore, digest: &str) -> Vec<String> {
    store.load_collection(digest).unwrap();
    store
        .get_collection(digest)
        .unwrap()
        .sequences
        .iter()
        .map(|s| s.metadata().sha512t24u.clone())
        .collect()
}

fn index_contains(store_dir: &std::path::Path, file: &str, needle: &str) -> bool {
    fs::read_to_string(store_dir.join(file))
        .unwrap()
        .lines()
        .any(|l| l.contains(needle))
}

/// THE resurrection bug. W opens the store, R removes a collection and its
/// orphan sequences and commits, then W commits an unrelated addition. W's
/// commit must not bring the removed collection, its sequences, or its alias
/// back.
///
/// Fails against merge-at-commit with dangling index rows: the `.seq` files are
/// gone but the rows point at them.
#[test]
fn test_stale_handle_does_not_resurrect_a_removed_collection() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut seed = RefgetStore::on_disk(&store_dir).unwrap();
    let keep = add_collection_to_store(&mut seed, work.path(), "keep", ">chrK\nAAAACCCC\n");
    let victim = add_collection_to_store(
        &mut seed,
        work.path(),
        "victim",
        ">chrV1\nGGGGTTTT\n>chrV2\nTTTTGGGG\n",
    );
    seed.add_collection_alias("refgenie", "oluc", &victim).unwrap();
    seed.write().unwrap();
    let victim_sequences = sequence_digests_of(&mut seed, &victim);
    drop(seed);

    // W opens BEFORE the removal. Its snapshot contains the victim.
    let mut writer = RefgetStore::open_local(&store_dir).unwrap();
    assert!(writer.get_collection_metadata(&victim).is_some());

    // R removes the victim and commits. Disk is correct at this point.
    let mut remover = RefgetStore::open_local(&store_dir).unwrap();
    assert!(remover.remove_collection(&victim, true).unwrap());
    drop(remover);
    assert!(!index_contains(&store_dir, "collections.rgci", &victim));

    // W now commits an unrelated addition from its stale snapshot.
    let added = add_collection_to_store(&mut writer, work.path(), "added", ">chrN\nACACACAC\n");
    writer.write().unwrap();
    drop(writer);

    assert!(
        !index_contains(&store_dir, "collections.rgci", &victim),
        "the removed collection came back into collections.rgci"
    );
    for seq in &victim_sequences {
        assert!(
            !index_contains(&store_dir, "sequences.rgsi", seq),
            "orphan sequence {} came back into sequences.rgsi -- its .seq file is gone",
            seq
        );
    }
    assert!(
        !store_dir.join("aliases/collections/refgenie.tsv").exists(),
        "the removed collection's alias namespace came back"
    );

    // ...and the concurrent addition is not collateral damage.
    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    assert!(reopened.get_collection_metadata(&keep).is_some());
    assert!(reopened.get_collection_metadata(&added).is_some());
    assert!(reopened.get_collection_metadata(&victim).is_none());
    assert!(
        reopened
            .get_collection_metadata_by_alias("refgenie", "oluc")
            .is_none(),
        "the alias still resolves to a collection that no longer exists"
    );
}

/// Every index row a delta commit publishes must have its `.seq` file on disk.
/// The observable symptom of the resurrection bug was not a missing row but a
/// present one that could not be read.
#[test]
fn test_no_dangling_sequence_rows_after_a_stale_commit() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut seed = RefgetStore::on_disk(&store_dir).unwrap();
    add_collection_to_store(&mut seed, work.path(), "keep", ">chrK\nAAAACCCC\n");
    let victim = add_collection_to_store(&mut seed, work.path(), "victim", ">chrV\nGGGGTTTT\n");
    seed.write().unwrap();
    drop(seed);

    let mut writer = RefgetStore::open_local(&store_dir).unwrap();

    let mut remover = RefgetStore::open_local(&store_dir).unwrap();
    remover.remove_collection(&victim, true).unwrap();
    drop(remover);

    add_collection_to_store(&mut writer, work.path(), "added", ">chrN\nACACACAC\n");
    writer.write().unwrap();
    drop(writer);

    let mut reopened = RefgetStore::open_local(&store_dir).unwrap();
    let digests: Vec<String> = reopened
        .list_sequences()
        .iter()
        .map(|m| m.sha512t24u.clone())
        .collect();
    for digest in &digests {
        // `load_sequence` reads the `.seq` file, so this fails on a row whose
        // bytes were unlinked. `get_sequence` would NOT: it happily returns the
        // Stub built from the index row, which is why the corruption was
        // invisible until someone actually asked for bases.
        assert!(
            reopened.load_sequence(digest).is_ok(),
            "sequences.rgsi lists {} but its .seq file is gone",
            digest
        );
    }
}

/// The lost-update case merge-at-commit existed to prevent must not regress:
/// two handles open together, one commits, the other commits after, and both
/// additions survive.
#[test]
fn test_delta_commit_still_preserves_a_concurrent_addition() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut seed = RefgetStore::on_disk(&store_dir).unwrap();
    let seed_digest = add_collection_to_store(&mut seed, work.path(), "seed", ">chrS\nAAAACCCC\n");
    seed.write().unwrap();
    drop(seed);

    let mut writer_a = RefgetStore::open_local(&store_dir).unwrap();
    let mut writer_b = RefgetStore::open_local(&store_dir).unwrap();

    let digest_a = add_collection_to_store(&mut writer_a, work.path(), "a", ">chrA\nGGGGTTTT\n");
    let digest_b = add_collection_to_store(&mut writer_b, work.path(), "b", ">chrB\nTTTTGGGG\n");
    writer_b.write().unwrap();
    writer_a.write().unwrap();
    drop(writer_a);
    drop(writer_b);

    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    for (label, digest) in [("seed", &seed_digest), ("A", &digest_a), ("B", &digest_b)] {
        assert!(
            reopened.get_collection_metadata(digest).is_some(),
            "collection {} ({}) was dropped by the other writer's commit",
            label,
            digest
        );
    }
}

/// Removing one of two collections that share sequences reclaims exactly the
/// unshared digests, and the shared ones stay readable byte-for-byte.
#[test]
fn test_removal_reclaims_only_unshared_sequences() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    const SHARED: &str = "AAAACCCCGGGGTTTT";
    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    let coll_a = add_collection_to_store(
        &mut store,
        work.path(),
        "a",
        &format!(">chrShared\n{}\n>chrOnlyA\nACACACAC\n", SHARED),
    );
    let coll_b = add_collection_to_store(
        &mut store,
        work.path(),
        "b",
        &format!(">chrShared\n{}\n>chrOnlyB\nGTGTGTGT\n", SHARED),
    );
    store.write().unwrap();
    let a_sequences = sequence_digests_of(&mut store, &coll_a);
    let b_sequences = sequence_digests_of(&mut store, &coll_b);
    drop(store);

    let shared: Vec<&String> = a_sequences.iter().filter(|d| b_sequences.contains(d)).collect();
    assert_eq!(shared.len(), 1, "precondition: exactly one shared sequence");
    let only_a: Vec<&String> = a_sequences.iter().filter(|d| !b_sequences.contains(d)).collect();
    assert_eq!(only_a.len(), 1);

    let mut store = RefgetStore::open_local(&store_dir).unwrap();
    let planned = store.plan_orphan_removal(&coll_a).unwrap();
    assert_eq!(planned, vec![only_a[0].clone()], "the dry-run must name exactly chrOnlyA");
    assert!(store.remove_collection(&coll_a, true).unwrap());
    drop(store);

    let mut reopened = RefgetStore::open_local(&store_dir).unwrap();
    assert!(
        reopened.get_sequence_metadata(only_a[0]).is_none(),
        "the unshared sequence was not reclaimed"
    );
    reopened.load_sequence(shared[0]).unwrap();
    assert_eq!(
        reopened.get_substring(shared[0], 0, SHARED.len()).unwrap(),
        SHARED,
        "the shared sequence did not survive intact"
    );
}

/// A namespace this handle empties must still not be deleted when another writer
/// has added an alias to it. Intent is not the deciding factor; the merged
/// result is.
#[test]
fn test_emptying_a_namespace_spares_another_writers_alias() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut seed = RefgetStore::on_disk(&store_dir).unwrap();
    let coll_a = add_collection_to_store(&mut seed, work.path(), "a", ">chrA\nAAAACCCC\n");
    let coll_b = add_collection_to_store(&mut seed, work.path(), "b", ">chrB\nGGGGTTTT\n");
    seed.add_collection_alias("refgenie", "mine", &coll_a).unwrap();
    seed.write().unwrap();
    drop(seed);

    // A opens with only `mine` in the namespace.
    let mut writer_a = RefgetStore::open_local(&store_dir).unwrap();

    // B adds a second alias to the same namespace and commits.
    let mut writer_b = RefgetStore::open_local(&store_dir).unwrap();
    writer_b.add_collection_alias("refgenie", "theirs", &coll_b).unwrap();
    drop(writer_b);

    // A removes its only alias -- emptying the namespace as far as A can see.
    assert!(writer_a.remove_collection_alias("refgenie", "mine").unwrap());
    drop(writer_a);

    let reopened = RefgetStore::open_local(&store_dir).unwrap();
    assert!(
        reopened
            .get_collection_metadata_by_alias("refgenie", "theirs")
            .is_some(),
        "the other writer's alias was deleted with the namespace"
    );
    assert!(
        reopened
            .get_collection_metadata_by_alias("refgenie", "mine")
            .is_none(),
        "the removed alias came back"
    );
}

/// Uncommitted work is visible to callers, so a caller that must not lose a
/// removal can check before dropping the handle.
#[test]
fn test_has_uncommitted_changes_tracks_the_commit() {
    let work = tempdir().unwrap();
    let store_dir = work.path().join("store");

    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    // on_disk commits after each add, so the store starts clean.
    add_collection_to_store(&mut store, work.path(), "a", ">chrA\nAAAACCCC\n");
    store.write().unwrap();
    assert!(!store.has_uncommitted_changes());

    store.lock_for_batch("test-batch").unwrap();
    let digest = add_collection_to_store(&mut store, work.path(), "b", ">chrB\nGGGGTTTT\n");
    assert!(store.get_collection_metadata(&digest).is_some());
    store.write().unwrap();
    assert!(
        !store.has_uncommitted_changes(),
        "a successful commit must clear the pending set"
    );
    store.release_batch_lock();
}
