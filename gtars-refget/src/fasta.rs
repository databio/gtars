use anyhow::{Context, Result};
use std::io::BufRead;
use std::path::Path;

use super::alphabet::{AlphabetGuesser, AlphabetType};
use super::collection::{SeqColDigestLvl1, SequenceCollection, SequenceMetadata, SequenceRecord};

use gtars_core::utils::get_dynamic_reader;

use md5::Md5;
use seq_io::fasta::{Reader, Record};
use sha2::{Digest, Sha512};

/// Processes a FASTA file to compute the digests of each sequence in the file.
///
/// This function reads a FASTA file, computes the SHA-512 and MD5 digests for each sequence,
/// and returns a vector of `SequenceDigest` structs containing the results. It can also handle
// gzipped FASTA files (ending in `.gz`).
///
/// # Arguments
///
/// * `file_path` - A string slice that holds the path to the FASTA file to be processed.
///
/// # Returns
///
/// A vector of `SequenceDigest` structs, each containing the length, SHA-512 digest, and MD5 digest
/// of a sequence in the FASTA file.
///
/// # Errors
///
/// This function will return an error if:
/// - The file cannot be opened or read
/// - The FASTA format is invalid or corrupted
/// - A sequence record is missing an ID
///
/// # Examples
///
///
pub fn digest_fasta<T: AsRef<Path>>(file_path: T) -> Result<SequenceCollection> {
    let file_reader = get_dynamic_reader(file_path.as_ref())?;
    let mut fasta_reader = Reader::new(file_reader);
    let mut results = Vec::new();
    println!("Processing FASTA file: {}", file_path.as_ref().display());
    while let Some(record) = fasta_reader.next() {
        // returns a RefRecord object
        let record = record.with_context(|| {
            format!(
                "Failed to read FASTA record from file: {}",
                file_path.as_ref().display()
            )
        })?;
        let id = record.id().with_context(|| {
            format!(
                "FASTA record #{} is missing a sequence ID",
                results.len() + 1
            )
        })?;
        let mut sha512_hasher = Sha512::new();
        let mut md5_hasher = Md5::new();
        let mut length = 0;
        let mut sequence = Vec::new();
        let mut alphabet_guesser = AlphabetGuesser::new();
        for seq_line in record.seq_lines() {
            let seq_line = seq_line.to_ascii_uppercase();
            sha512_hasher.update(&seq_line);
            md5_hasher.update(&seq_line);
            length += seq_line.len();
            sequence.extend_from_slice(&seq_line);
            alphabet_guesser.update(&seq_line);
        }
        let sha512 = base64_url::encode(&sha512_hasher.finalize_reset()[0..24]);
        let md5 = format!("{:x}", md5_hasher.finalize_reset());
        let alphabet = alphabet_guesser.guess();

        // Create SequenceRecord
        let metadata = SequenceMetadata {
            name: id.to_string(),
            length,
            sha512t24u: sha512,
            md5,
            alphabet,
        };

        results.push(SequenceRecord {
            metadata,
            data: None,
        });
    }

    // Compute lvl1 digests from the sequence records
    let metadata_refs: Vec<&SequenceMetadata> = results.iter().map(|r| &r.metadata).collect();
    let lvl1 = SeqColDigestLvl1::from_metadata(&metadata_refs);

    // Compute collection digest from lvl1 digests
    let collection_digest = lvl1.to_digest();

    Ok(SequenceCollection {
        sequences: results,
        digest: collection_digest,
        lvl1,
        file_path: Some(file_path.as_ref().to_path_buf()),
        has_data: true,
    })
}

/// Read a FARG file and return a SequenceCollection struct with all metadata.
///
/// This will not read the actual sequence data, only the metadata.
/// # Arguments
/// * `file_path` - The path to the FARG file to be read.
pub fn read_fasta_refget_file<T: AsRef<Path>>(file_path: T) -> Result<SequenceCollection> {
    let file = std::fs::File::open(&file_path)?;
    let reader = std::io::BufReader::new(file);
    let mut results = Vec::new();

    // Variables to store header metadata
    let mut seqcol_digest = String::new();
    let mut names_digest = String::new();
    let mut sequences_digest = String::new();
    let mut lengths_digest = String::new();

    for line in reader.lines() {
        let line = line?;

        // Parse header metadata lines
        if line.starts_with("##") {
            if let Some(stripped) = line.strip_prefix("##")
                && let Some((key, value)) = stripped.split_once('=')
            {
                match key {
                    "seqcol_digest" => seqcol_digest = value.to_string(),
                    "names_digest" => names_digest = value.to_string(),
                    "sequences_digest" => sequences_digest = value.to_string(),
                    "lengths_digest" => lengths_digest = value.to_string(),
                    _ => {} // Ignore unknown metadata keys
                }
            }
            continue;
        }

        // Skip other comment lines
        if line.starts_with('#') {
            continue;
        }

        // Parse sequence data lines
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 5 {
            continue; // Skip lines that don't have exactly 5 columns
        }

        let result = SequenceMetadata {
            name: parts[0].to_string(),
            length: parts[1].parse().unwrap_or(0),
            alphabet: parts[2].parse().unwrap_or(AlphabetType::Unknown),
            sha512t24u: parts[3].to_string(),
            md5: parts[4].to_string(),
        };

        let record = SequenceRecord {
            metadata: result,
            data: None,
        };
        results.push(record);
    }
    // If the digests were not found in the file, compute them
    #[allow(clippy::needless_late_init)]
    let lvl1: SeqColDigestLvl1;
    if (sequences_digest.is_empty() || names_digest.is_empty() || lengths_digest.is_empty())
        && !results.is_empty()
    {
        let metadata_vec: Vec<&SequenceMetadata> = results.iter().map(|r| &r.metadata).collect();
        lvl1 = SeqColDigestLvl1::from_metadata(&metadata_vec);
    } else {
        lvl1 = SeqColDigestLvl1 {
            sequences_digest,
            names_digest,
            lengths_digest,
        };
    }

    if seqcol_digest.is_empty() {
        // If seqcol_digest is not provided, compute it from lvl1
        seqcol_digest = lvl1.to_digest();
    }

    // Return the complete SequenceCollection
    Ok(SequenceCollection {
        sequences: results,
        digest: seqcol_digest,
        lvl1,
        file_path: Some(file_path.as_ref().to_path_buf()),
        has_data: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn digests_digest_fasta() {
        let seqcol =
            digest_fasta("../tests/data/fasta/base.fa").expect("Can't open test fasta file");
        let results = seqcol.sequences;
        println!("{:?}", results);
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].metadata.length, 8);
        assert_eq!(
            results[0].metadata.sha512t24u,
            "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        );
        assert_eq!(results[0].metadata.md5, "5f63cfaa3ef61f88c9635fb9d18ec945");
        assert_eq!(results[0].metadata.alphabet, AlphabetType::Dna2bit);
        assert_eq!(results[1].metadata.length, 4);
        assert_eq!(
            results[1].metadata.sha512t24u,
            "YBbVX0dLKG1ieEDCiMmkrTZFt_Z5Vdaj"
        );
        assert_eq!(results[1].metadata.md5, "31fc6ca291a32fb9df82b85e5f077e31");
        assert_eq!(results[1].metadata.alphabet, AlphabetType::Dna2bit);
        assert_eq!(results[2].metadata.length, 4);
        assert_eq!(
            results[2].metadata.sha512t24u,
            "AcLxtBuKEPk_7PGE_H4dGElwZHCujwH6"
        );
        assert_eq!(results[2].metadata.md5, "92c6a56c9e9459d8a42b96f7884710bc");
        assert_eq!(results[2].metadata.alphabet, AlphabetType::Dna2bit);
    }

    #[test]
    fn digests_digest_gzipped_fasta() {
        let seqcol =
            digest_fasta("../tests/data/fasta/base.fa.gz").expect("Can't open test fasta file");
        let results = seqcol.sequences;
        println!("{:?}", results);
        assert_eq!(results[0].metadata.length, 8);
        assert_eq!(
            results[0].metadata.sha512t24u,
            "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        );
        assert_eq!(results[0].metadata.md5, "5f63cfaa3ef61f88c9635fb9d18ec945");
        assert_eq!(results[0].metadata.alphabet, AlphabetType::Dna2bit);
    }

    #[test]
    fn digests_bogus_fasta_file() {
        let result = digest_fasta("../tests/data/bogus.fa");
        assert!(result.is_err(), "Expected an error for a bogus fasta file");
    }

    #[test]
    fn digests_fa_to_farg() {
        let seqcol = SequenceCollection::from_path_no_cache("../tests/data/fasta/base.fa")
            .expect("Failed to create SequenceCollection from FASTA file");
        seqcol.to_farg().expect("Failed to write farg file");

        let loaded_seqcol = read_fasta_refget_file("../tests/data/fasta/base.farg")
            .expect("Failed to read refget file");
        println!("Original SequenceCollection: {}", seqcol);
        println!("Loaded SequenceCollection: {}", loaded_seqcol);
        // Test round-trip integrity
        for (original, read) in seqcol.sequences.iter().zip(loaded_seqcol.sequences.iter()) {
            assert_eq!(original.metadata.name, read.metadata.name);
            assert_eq!(original.metadata.length, read.metadata.length);
            assert_eq!(original.metadata.sha512t24u, read.metadata.sha512t24u);
            assert_eq!(original.metadata.md5, read.metadata.md5);
            assert_eq!(original.metadata.alphabet, read.metadata.alphabet);
        }
    }

    #[test]
    fn digests_seqcol_from_fasta() {
        let seqcol = SequenceCollection::from_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to create SequenceCollection from FASTA file");
        println!("SeqCol digest: {:?}", seqcol.digest);
        println!(
            "SeqCol sequences_digest: {:?}",
            seqcol.lvl1.sequences_digest
        );
        println!("SequenceCollection: {:?}", seqcol);
        assert_eq!(
            seqcol.lvl1.sequences_digest,
            "0uDQVLuHaOZi1u76LjV__yrVUIz9Bwhr"
        );
        assert_eq!(seqcol.lvl1.names_digest, "Fw1r9eRxfOZD98KKrhlYQNEdSRHoVxAG");
        assert_eq!(
            seqcol.lvl1.lengths_digest,
            "cGRMZIb3AVgkcAfNv39RN7hnT5Chk7RX"
        );
    }
}
