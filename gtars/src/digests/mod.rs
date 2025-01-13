//! # Fast digest computations for genomic sequences
//!
//! This module provides functions for computing digests of strings.
//!
//! # Functions
//!
//! The following functions are available:
//!
//! * `sha512t24u` - Processes a given string to compute its GA4GH sha512t24 checksum.
//! * `md5` - Processes a given string to compute its MD5 checksum.
//! * `digest_fasta` - Processes a FASTA file to compute the digests of each sequence in the file.
//!
//! # Usage
//!
//! The `sha512t24u` function can be used to compute the GA4GH sha512t24 checksum of a string.
//!
//! ```rust
//! use gtars::digests::sha512t24u;
//!
//! let digest = sha512t24u("hello world");
//! ```
use std::fs::File;
use std::io;
use std::io::prelude::{Read, Write};
use std::path::Path;

use anyhow::Result;
use md5::Md5;
use seq_io::fasta::{Reader, Record, RefRecord};
use sha2::{Digest, Sha512};

use crate::common::utils::get_dynamic_reader;

/// A struct representing the digest of a given string.
#[derive(Debug)]
pub struct DigestResult {
    pub id: String,
    pub length: usize,
    pub sha512t24u: String,
    pub md5: String,
}

/// Processes a given string to compute its GA4GH sha512t24u digest.
///
/// # Arguments
///
/// * `string` - The input string to be processed.
///
/// # Returns
///
/// A string SHA-512 digest of the input string.
pub fn sha512t24u(string: &str) -> String {
    let mut sha512_hasher_box = Box::new(Sha512::new());
    for s in string.as_bytes().chunks(800) {
        sha512_hasher_box.as_mut().update(s);
    }
    base64_url::encode(&sha512_hasher_box.as_mut().finalize_reset()[0..24])
}

/// Process a string to compute its md5 digest
///
/// # Arguments
///
/// * `string` - The input string to be processed.
///
/// # Returns
///
/// A string MD5 digest of the input string.
pub fn md5(string: &str) -> String {
    let mut hasher = Md5::new();
    for s in string.as_bytes().chunks(800) {
        hasher.update(s);
    }
    let result = hasher.finalize();
    format!("{:x}", result)
}

/// Processes a FASTA file to compute the digests of each sequence in the file.
///
/// This function reads a FASTA file, computes the SHA-512 and MD5 digests for each sequence,
/// and returns a vector of `DigestResult` structs containing the results. It can also handle
// gzipped FASTA files (ending in `.gz`).
///
/// # Arguments
///
/// * `file_path` - A string slice that holds the path to the FASTA file to be processed.
///
/// # Returns
///
/// A vector of `DigestResult` structs, each containing the length, SHA-512 digest, and MD5 digest
/// of a sequence in the FASTA file.
///
/// # Panics
///
/// This function will panic if the file cannot be opened or if there is an error reading the file.
///
/// # Examples
///
///
pub fn digest_fasta(file_path: &str) -> Result<Vec<DigestResult>> {
    let path = Path::new(&file_path);
    let file_reader = get_dynamic_reader(&path)?;
    let mut fasta_reader = Reader::new(file_reader);
    let mut results = Vec::new();
    while let Some(record) = fasta_reader.next() {
        // returns a RefRecord object
        let record = record.expect("Error found when retrieving next record.");
        let id = record.id().expect("No ID found for the FASTA record");
        let mut sha512_hasher = Sha512::new();
        let mut md5_hasher = Md5::new();
        let mut length = 0;
        // let result = process_sequence(record, verbose);
        for seq_line in record.seq_lines() {
            // let seq_line = seq_line.expect("Error found when retrieving next sequence line.");
            sha512_hasher.update(seq_line.to_ascii_uppercase());
            md5_hasher.update(seq_line.to_ascii_uppercase());
            length += seq_line.len();
        }
        // let result = sha512_hasher.finalize();
        let sha512 = base64_url::encode(&sha512_hasher.finalize_reset()[0..24]);
        let md5 = format!("{:x}", md5_hasher.finalize_reset());
        results.push(DigestResult {
            id: id.to_string(),
            length: length,
            sha512t24u: sha512,
            md5: md5,
        });
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sha512t24u() {
        let digest = sha512t24u("hello world");
        assert_eq!(digest, "MJ7MSJwS1utMxA9QyQLytNDtd-5RGnx6");
    }

    #[test]
    fn test_md5() {
        let digest = md5("hello world");
        assert_eq!(digest, "5eb63bbbe01eeed093cb22bb8f5acdc3");
    }

    #[test]
    fn test_digest_fasta() {
        let results = digest_fasta("tests/data/base.fa").expect("Can't open test fasta file");
        println!("{:?}", results);
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].length, 8);
        assert_eq!(results[0].sha512t24u, "iYtREV555dUFKg2_agSJW6suquUyPpMw");
        assert_eq!(results[0].md5, "5f63cfaa3ef61f88c9635fb9d18ec945");
        assert_eq!(results[1].length, 4);
        assert_eq!(results[1].sha512t24u, "YBbVX0dLKG1ieEDCiMmkrTZFt_Z5Vdaj");
        assert_eq!(results[1].md5, "31fc6ca291a32fb9df82b85e5f077e31");
        assert_eq!(results[2].length, 4);
        assert_eq!(results[2].sha512t24u, "AcLxtBuKEPk_7PGE_H4dGElwZHCujwH6");
        assert_eq!(results[2].md5, "92c6a56c9e9459d8a42b96f7884710bc");
    }

    #[test]
    fn test_digest_gzipped_fasta() {
        let results = digest_fasta("tests/data/base.fa.gz").expect("Can't open test fasta file");
        println!("{:?}", results);
        assert_eq!(results[0].length, 8);
        assert_eq!(results[0].sha512t24u, "iYtREV555dUFKg2_agSJW6suquUyPpMw");
        assert_eq!(results[0].md5, "5f63cfaa3ef61f88c9635fb9d18ec945");
    }

    #[test]
    fn bogus_fasta_file() {
        let result = digest_fasta("tests/data/bogus.fa");
        assert!(result.is_err(), "Expected an error for a bogus fasta file");
    }
}
