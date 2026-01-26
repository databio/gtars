//! # Rust implementation of GA4GH Refget sequence collection functions

//! This module provides functions managing and retrieving sequences from a sequence collection.
//!
//! # Functions
//!
//! The module includes the following main components:
//!
//! * `alphabet.rs` - Defines various sequence alphabets (e.g., DNA, protein, ASCII).
//! * `collection.rs` - Contains the `SequenceCollection` struct and methods for managing sequence collections.
//! * `digest.rs` - Implements functions for calculating and verifying sha512t24u and other digests.
//! * `encoder.rs` - Contains functions for encoding sequences into compact representations.
//! * `fasta.rs` - Provides functions for reading and writing FASTA files.
//! * `store.rs` - Implements a sequence store that allows for efficient storage and retrieval of sequences indexed by sha512t24u digest.
pub mod alphabet;
pub mod collection;
pub mod digest;
pub mod encoder;
pub mod fasta;
pub mod store;

// Used internally to make it easy to convert types to a 32-byte key for hash tables
mod hashkeyable;
mod utils;

#[cfg(test)]
mod tests {
    use super::*;

    use std::time::Instant;
    use store::RefgetStore;
    use tempfile::tempdir;
    #[test]
    #[ignore]
    fn test_loading_large_fasta_file() {
        // Path to a large FASTA file
        // let fasta_path = "GRCh38_full_analysis_set_plus_decoy_hla.fa";
        let fasta_path =
            std::env::var("FASTA_PATH").expect("FASTA_PATH environment variable not set");
        // let fasta_path = "../tests/data/subset.fa.gz";
        // let fasta_path = "../tests/data/fasta/base.fa.gz";
        println!("Loading large FASTA file: {}", &fasta_path);

        // Create a new sequence store, and dd sequences to the store
        println!("Adding sequences from FASTA file...");
        let start = Instant::now();
        let mut store = RefgetStore::in_memory();
        store.add_sequence_collection_from_fasta(&fasta_path).unwrap();
        let duration = start.elapsed();
        println!("⏱️  Time taken to load: {:.2?}", duration);

        let mut store2 = RefgetStore::in_memory();
        store2.disable_encoding();  // Switch to Raw mode
        store2.add_sequence_collection_from_fasta(&fasta_path).unwrap();

        // Get list of sequences
        let sequences: Vec<_> = store.sequence_digests().collect();
        assert!(!sequences.is_empty(), "No sequences found in the store");

        // Look up the first sequence by digest
        println!("Look up a sequence by digest...");
        let digest = &sequences[0];
        let digest_str = String::from_utf8(digest.to_vec()).expect("Invalid ASCII data");
        // let seq = store.get_sequence(name);
        // assert!(seq.is_some(), "Failed to retrieve sequence with name: {}", name);
        // println!("Retrieved sequence: {:?}", seq.unwrap());

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
        // let fasta_path = "../tests/data/subset.fa.gz";
        let fasta_path = "../tests/data/fasta/base.fa.gz";
        let temp_fasta = temp_path.join("base.fa.gz");
        std::fs::copy(fasta_path, &temp_fasta).expect("Failed to copy base.fa.gz to tempdir");

        // Add sequences to the store
        store.add_sequence_collection_from_fasta(temp_fasta).unwrap();
        println!("Listing sequences in the store...");
        // let sequences = store.sequence_digests();
        // let digest = &sequences[0];
        // let digest_str = String::from_utf8(digest.to_vec()).expect("Invalid ASCII data");
        // let digest = "Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO";  // from subset.fa.gz
        let digest = "iYtREV555dUFKg2_agSJW6suquUyPpMw"; // from base.fa.gz
        let digest_str = String::from_utf8(digest.as_bytes().to_vec()).expect("Invalid ASCII data");

        // Test retrieval of a substring
        println!("Retrieving a substring of sequence named: {:?}", digest_str);
        let start = 2;
        let end = start + 5;
        let substring = store.get_substring(digest, start, end);
        assert!(
            substring.is_ok(),
            "Failed to retrieve substring with name: {:?}",
            digest_str
        );
        println!("Retrieved substring: {:?}", substring.as_ref().unwrap());
        assert_eq!(substring.unwrap(), "GGGGA");
        // assert!(substring.unwrap() == "CCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCC");

        println!("Retrieving a substring of sequence named: {:?}", digest_str);
        let start = 3;
        let end = start + 2;
        let substring = store.get_substring(digest, start, end);
        assert!(
            substring.is_ok(),
            "Failed to retrieve substring with name: {:?}",
            digest_str
        );
        println!("Retrieved substring: {:?}", substring.as_ref().unwrap());
        assert_eq!(substring.unwrap(), "GG");
        // assert!(substring.unwrap() == "TCTGACCTGAGGAGAACTGTGCTCCGCCTTCAGAGTACCACCGAAATCTGTGCAGAGGACAACGCAGCTC");
    }
}
