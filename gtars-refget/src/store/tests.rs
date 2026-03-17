//! Store integration tests -- cross-cutting concerns (import, export, persistence, lazy-loading).

use super::*;

use crate::collection::{
    SequenceCollection, SequenceCollectionExt, SequenceCollectionMetadata, SequenceMetadata, SequenceRecord,
};
use crate::digest::{AlphabetType, md5, sha512t24u};
use crate::hashkeyable::HashKeyable;
use std::fs;
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
fn setup_export_test_store(temp_path: &std::path::Path) -> (RefgetStore, [u8; 32]) {
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

#[test]
fn test_expand_template() {
    let digest = "ABCDEFghijklmnop";

    let result = ReadonlyRefgetStore::expand_template(digest, "sequences/%s2/%s.seq");
    assert_eq!(result, std::path::PathBuf::from("sequences/AB/ABCDEFghijklmnop.seq"));

    let result = ReadonlyRefgetStore::expand_template(digest, "sequences/%s2/%s4/%s.seq");
    assert_eq!(result, std::path::PathBuf::from("sequences/AB/ABCD/ABCDEFghijklmnop.seq"));

    let result = ReadonlyRefgetStore::expand_template(digest, "sequences/%s.seq");
    assert_eq!(result, std::path::PathBuf::from("sequences/ABCDEFghijklmnop.seq"));
}

#[test]
fn test_sanitize_relative_path() {
    // Rejects traversal
    assert!(ReadonlyRefgetStore::sanitize_relative_path("../etc/passwd").is_err());
    assert!(ReadonlyRefgetStore::sanitize_relative_path("foo/../bar").is_err());
    assert!(ReadonlyRefgetStore::sanitize_relative_path("foo/../../bar").is_err());
    assert!(ReadonlyRefgetStore::sanitize_relative_path("..").is_err());

    // Rejects absolute
    assert!(ReadonlyRefgetStore::sanitize_relative_path("/etc/passwd").is_err());
    assert!(ReadonlyRefgetStore::sanitize_relative_path("\\windows\\system32").is_err());

    // Accepts valid
    assert!(ReadonlyRefgetStore::sanitize_relative_path("sequences/ab/abc123.seq").is_ok());
    assert!(ReadonlyRefgetStore::sanitize_relative_path("collections/xyz.rgsi").is_ok());
    assert!(ReadonlyRefgetStore::sanitize_relative_path("rgstore.json").is_ok());
    assert!(ReadonlyRefgetStore::sanitize_relative_path("sequences/%s2/%s.seq").is_ok());
}

// =========================================================================
// Mode tests
// =========================================================================

#[test]
fn test_mode_basics() {
    let mut store = RefgetStore::in_memory();

    assert_eq!(store.mode, StorageMode::Encoded);

    store.disable_encoding();
    assert_eq!(store.mode, StorageMode::Raw);
    store.enable_encoding();
    assert_eq!(store.mode, StorageMode::Encoded);

    store.set_encoding_mode(StorageMode::Raw);
    assert_eq!(store.mode, StorageMode::Raw);
    store.set_encoding_mode(StorageMode::Encoded);
    assert_eq!(store.mode, StorageMode::Encoded);
}

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
            assert_eq!(sequence, b"ATGCATGCATGC");
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
            assert_eq!(sequence, b"ATGCATGCATGC");
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

    let (chr1_sha, chr1_md5) = calculate_test_digests(b"ATGCATGCATGC");
    let (chr2_sha, chr2_md5) = calculate_test_digests(b"GGGGAAAA");

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

    let expected_fa_content = format!(
        ">chr1 12 dna2bit {} {}\nATGCAATGC\n>chr2 8 dna2bit {} {}\nGGGG\n",
        chr1_sha, chr1_md5, chr2_sha, chr2_md5
    );
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
        sequence: sequence.to_vec(),
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

    let sequence_keys: Vec<[u8; 32]> = store.sequence_store.keys().cloned().collect();
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

    let collection_digest: [u8; 32];
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

    let disk_size = store.total_disk_size();
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
    assert_eq!(stats.n_collections_loaded, 1);
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
    assert_eq!(stats_before.n_collections_loaded, 0);

    let seq = loaded_store.get_sequence_by_name(&digest, "chr1");
    assert!(seq.is_err());

    loaded_store.load_collection(&digest).unwrap();

    let seq = loaded_store.get_sequence_by_name(&digest, "chr1");
    assert!(seq.is_ok());
    assert_eq!(seq.unwrap().metadata().name, "chr1");

    assert!(loaded_store.is_collection_loaded(&digest));

    let stats_after = loaded_store.stats();
    assert_eq!(stats_after.n_collections_loaded, 1);
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
    assert_eq!(stats_after.n_sequences_loaded, 0);
    assert_eq!(stats_after.n_collections_loaded, 1);

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

    assert_eq!(loaded_store.stats().n_collections_loaded, 1);
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

// =========================================================================
// Decoded cache tests
// =========================================================================

#[test]
fn test_ensure_decoded_then_sequence_bytes() {
    use crate::digest::digest_sequence;

    let mut store = RefgetStore::in_memory();
    let record = digest_sequence("chr1", b"ACGTACGT");
    let digest = record.metadata().sha512t24u.clone();

    store.add_sequence_record(record, false).unwrap();
    store.ensure_decoded(digest.as_bytes()).unwrap();

    let bytes = store.sequence_bytes(digest.as_bytes()).unwrap();
    assert_eq!(bytes, b"ACGTACGT");
}

#[test]
fn test_sequence_bytes_returns_none_before_ensure_decoded() {
    use crate::digest::digest_sequence;

    let mut store = RefgetStore::in_memory();
    let record = digest_sequence("chr1", b"ACGT");
    let digest = record.metadata().sha512t24u.clone();

    store.add_sequence_record(record, false).unwrap();
    assert!(store.sequence_bytes(digest.as_bytes()).is_none());
}

#[test]
fn test_ensure_decoded_is_idempotent() {
    use crate::digest::digest_sequence;

    let mut store = RefgetStore::in_memory();
    let record = digest_sequence("chr1", b"GATTACA");
    let digest = record.metadata().sha512t24u.clone();

    store.add_sequence_record(record, false).unwrap();

    store.ensure_decoded(digest.as_bytes()).unwrap();
    store.ensure_decoded(digest.as_bytes()).unwrap();
    store.ensure_decoded(digest.as_bytes()).unwrap();

    let bytes = store.sequence_bytes(digest.as_bytes()).unwrap();
    assert_eq!(bytes, b"GATTACA");
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
