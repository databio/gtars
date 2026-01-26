//! # Rust implementation of GA4GH Refget sequence collection functions

//! This module provides functions managing and retrieving sequences from a sequence collection.
//!
//! # Module Structure
//!
//! The library is organized into two main parts:
//!
//! ## Core (WASM-compatible)
//!
//! The `digest` module contains all WASM-compatible code:
//! - `digest::algorithms` - Hash functions (sha512t24u, md5, canonicalize_json)
//! - `digest::alphabet` - Sequence alphabets and encoding tables
//! - `digest::encoder` - Sequence bit-packing
//! - `digest::types` - Core data structures (SequenceRecord, SequenceCollection)
//! - `digest::fasta` - Bytes-based FASTA parsing
//! - `digest::stream` - Streaming FASTA hasher for chunk-by-chunk processing
//!
//! ## Filesystem (requires `filesystem` feature)
//!
//! - `fasta` - File-based FASTA parsing (wraps digest::fasta with file I/O)
//! - `collection` - Extended SequenceCollection with filesystem operations
//! - `store` - RefgetStore for persistent sequence storage
//!
//! # Feature Flags
//!
//! - `filesystem` (default): Enables file-based operations
//! - Without `filesystem`: Only WASM-compatible code in `digest` module

// ============================================================================
// Core WASM-compatible module
// ============================================================================

/// Core digest and encoding functionality - WASM-safe.
/// All code in this module works without filesystem access.
pub mod digest;

// Re-export commonly used items from digest at crate root for convenience
pub use digest::{
    // Algorithms
    sha512t24u, md5, canonicalize_json,
    // Alphabet
    Alphabet, AlphabetType, AlphabetGuesser,
    lookup_alphabet, guess_alphabet,
    ASCII_ALPHABET, DNA_2BIT_ALPHABET, DNA_3BIT_ALPHABET, DNA_IUPAC_ALPHABET, PROTEIN_ALPHABET,
    // Encoder
    encode_sequence, decode_string_from_bytes, decode_substring_from_bytes, SequenceEncoder,
    // Types
    SequenceRecord, SequenceMetadata, SequenceCollection, SequenceCollectionMetadata,
    SequenceCollectionRecord, SeqColDigestLvl1, FaiMetadata,
    digest_sequence, digest_sequence_with_description, parse_rgsi_line,
    // Fasta (bytes-based, WASM-compatible)
    digest_fasta_bytes, load_fasta_bytes, parse_fasta_header, ParseOptions,
    // Streaming
    FastaStreamHasher,
};

// ============================================================================
// Filesystem-dependent modules (require `filesystem` feature)
// ============================================================================

/// File-based FASTA operations.
/// Wraps the WASM-compatible digest::fasta with filesystem I/O.
#[cfg(feature = "filesystem")]
pub mod fasta;

/// Extended SequenceCollection with filesystem operations.
/// Adds methods for RGSI file I/O, caching, and file-based construction.
#[cfg(feature = "filesystem")]
pub mod collection;

/// Persistent sequence storage (RefgetStore).
#[cfg(feature = "filesystem")]
pub mod store;

// Internal modules for filesystem operations
#[cfg(feature = "filesystem")]
mod hashkeyable;
#[cfg(feature = "filesystem")]
mod utils;

// Re-export filesystem functions at crate root for backward compatibility
#[cfg(feature = "filesystem")]
pub use fasta::{compute_fai, digest_fasta, load_fasta, FaiRecord};
#[cfg(feature = "filesystem")]
pub use collection::{
    read_rgsi_file,
    SequenceCollectionExt, SequenceCollectionRecordExt, SequenceMetadataExt, SequenceRecordExt,
};

// ============================================================================
// Tests
// ============================================================================

#[cfg(all(test, feature = "filesystem"))]
mod tests {
    use super::*;

    use std::time::Instant;
    use store::RefgetStore;
    use tempfile::tempdir;

    #[test]
    #[ignore]
    fn test_loading_large_fasta_file() {
        // Path to a large FASTA file
        let fasta_path =
            std::env::var("FASTA_PATH").expect("FASTA_PATH environment variable not set");
        println!("Loading large FASTA file: {}", &fasta_path);

        // Create a new sequence store, and dd sequences to the store
        println!("Adding sequences from FASTA file...");
        let start = Instant::now();
        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta(&fasta_path)
            .unwrap();
        let duration = start.elapsed();
        println!("Time taken to load: {:.2?}", duration);

        let mut store2 = RefgetStore::in_memory();
        store2.disable_encoding(); // Switch to Raw mode
        store2
            .add_sequence_collection_from_fasta(&fasta_path)
            .unwrap();

        // Get list of sequences
        let sequences: Vec<_> = store.sequence_digests().collect();
        assert!(!sequences.is_empty(), "No sequences found in the store");

        // Look up the first sequence by digest
        println!("Look up a sequence by digest...");
        let digest = &sequences[0];
        let digest_str = String::from_utf8(digest.to_vec()).expect("Invalid ASCII data");

        // Test retrieval of a substring
        println!("Retrieving a substring of sequence named: {:?}", digest_str);
        let start_basic = 0;
        let end_basic = 3;
        let substring = store.get_substring(digest, start_basic, end_basic);
        assert!(
            substring.is_ok(),
            "Failed to retrieve substring with name: {:?}",
            digest_str
        );
        println!("Retrieved substring: {:?}", substring.unwrap());

        // Retrieve substring via digest
        let start = 148 * 70;
        let end = 148 * 70 + 70;
        let substring2 = store.get_substring(digest, start, end);
        assert!(
            substring2.is_ok(),
            "Failed to retrieve substring with name: {:?}",
            digest_str
        );

        let substring3 = store2.get_substring(digest, start, end);
        assert_eq!(substring2.as_ref().unwrap(), substring3.as_ref().unwrap());
        println!("Retrieved substring: {:?}", substring2.unwrap());
        println!("Retrieved substring: {:?}", substring3.unwrap());
    }

    #[test]
    fn test_get_sequence_encoded() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();
        // Create a new sequence store
        let mut store = RefgetStore::in_memory();
        let fasta_path = "../tests/data/fasta/base.fa.gz";
        let temp_fasta = temp_path.join("base.fa.gz");
        std::fs::copy(fasta_path, &temp_fasta).expect("Failed to copy base.fa.gz to tempdir");

        // Add sequences to the store
        store
            .add_sequence_collection_from_fasta(temp_fasta)
            .unwrap();
        println!("Listing sequences in the store...");
        let digest = "iYtREV555dUFKg2_agSJW6suquUyPpMw"; // from base.fa.gz
        let digest_str =
            String::from_utf8(digest.as_bytes().to_vec()).expect("Invalid ASCII data");

        // Test retrieval of a substring
        println!(
            "Retrieving a substring of sequence named: {:?}",
            digest_str
        );
        let start = 2;
        let end = start + 5;
        let substring = store.get_substring(digest, start, end);
        assert!(
            substring.is_ok(),
            "Failed to retrieve substring with name: {:?}",
            digest_str
        );
        println!(
            "Retrieved substring: {:?}",
            substring.as_ref().unwrap()
        );
        assert_eq!(substring.unwrap(), "GGGGA");

        println!(
            "Retrieving a substring of sequence named: {:?}",
            digest_str
        );
        let start = 3;
        let end = start + 2;
        let substring = store.get_substring(digest, start, end);
        assert!(
            substring.is_ok(),
            "Failed to retrieve substring with name: {:?}",
            digest_str
        );
        println!(
            "Retrieved substring: {:?}",
            substring.as_ref().unwrap()
        );
        assert_eq!(substring.unwrap(), "GG");
    }
}
