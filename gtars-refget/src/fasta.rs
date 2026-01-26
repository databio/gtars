//! File-based FASTA operations.
//!
//! This module provides FASTA parsing functions that work with file paths.
//! It wraps the WASM-compatible digest::fasta module with filesystem I/O.
//!
//! For in-memory/WASM-compatible parsing, use `crate::digest::fasta` directly.

use anyhow::Result;
use std::path::Path;

use crate::digest;
use crate::digest::types::{
    FaiMetadata, SequenceCollection, SequenceCollectionMetadata, SequenceMetadata, SequenceRecord,
};

// Re-export types from digest::fasta
pub use digest::fasta::{parse_fasta_header, ParseOptions};

/// A lightweight record containing only FAI (FASTA index) metadata for a sequence.
/// Returned by `compute_fai()` for fast FAI-only computation without digest overhead.
#[derive(Clone, Debug)]
pub struct FaiRecord {
    pub name: String,
    pub length: usize,
    pub fai: Option<FaiMetadata>,
}

/// Configuration for file-based FASTA parsing
struct FileParseOptions {
    /// Compute SHA512, MD5 digests and alphabet
    compute_digests: bool,
    /// Store sequence data in memory
    store_sequence: bool,
    /// Track FAI metadata (disabled for gzipped files)
    fai_enabled: bool,
}

impl FileParseOptions {
    const FAI_ONLY: Self = Self {
        compute_digests: false,
        store_sequence: false,
        fai_enabled: true,
    };
    const DIGEST_ONLY: Self = Self {
        compute_digests: true,
        store_sequence: false,
        fai_enabled: true, // Will be disabled for gzipped files
    };
    const FULL: Self = Self {
        compute_digests: true,
        store_sequence: true,
        fai_enabled: true, // Will be disabled for gzipped files
    };
}

/// Result of internal file-based FASTA parsing
enum FileParseResult {
    Fai(Vec<FaiRecord>),
    Sequences(Vec<SequenceRecord>),
}

/// Internal file-based FASTA parser.
fn parse_fasta_file<P: AsRef<Path>>(file_path: P, mut opts: FileParseOptions) -> Result<FileParseResult> {
    use gtars_core::utils::get_dynamic_reader;
    use md5::Md5;
    use sha2::{Digest, Sha512};
    use std::io::BufRead;

    // Detect if file is gzipped - skip FAI data for gzipped files
    let is_gzipped = file_path
        .as_ref()
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s == "gz")
        .unwrap_or(false);

    // Disable FAI for gzipped files (can't seek into compressed files)
    if is_gzipped {
        opts.fai_enabled = false;
    }

    let file_reader = get_dynamic_reader(file_path.as_ref())?;
    let mut reader = std::io::BufReader::new(file_reader);

    // State
    let mut byte_position: u64 = 0;
    let mut line = String::new();

    let mut current_id: Option<String> = None;
    let mut current_description: Option<String> = None;
    let mut current_offset: u64 = 0;
    let mut current_line_bases: Option<u32> = None;
    let mut current_line_bytes: Option<u32> = None;
    let mut length = 0;

    // Digest state
    let mut sha512_hasher = if opts.compute_digests {
        Some(Sha512::new())
    } else {
        None
    };
    let mut md5_hasher = if opts.compute_digests {
        Some(Md5::new())
    } else {
        None
    };
    let mut alphabet_guesser = if opts.compute_digests {
        Some(digest::AlphabetGuesser::new())
    } else {
        None
    };

    // Sequence data
    let mut sequence_data: Vec<u8> = Vec::new();

    // Results
    let mut fai_results: Vec<FaiRecord> = Vec::new();
    let mut seq_results: Vec<SequenceRecord> = Vec::new();

    loop {
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            // EOF - finalize the last sequence
            if let Some(id) = current_id.take() {
                let fai = make_fai(
                    opts.fai_enabled,
                    current_offset,
                    current_line_bases,
                    current_line_bytes,
                );

                if opts.compute_digests {
                    let sha512 =
                        base64_url::encode(&sha512_hasher.as_mut().unwrap().finalize_reset()[0..24]);
                    let md5 = format!("{:x}", md5_hasher.as_mut().unwrap().finalize_reset());
                    let alphabet = alphabet_guesser.as_mut().unwrap().guess();

                    let metadata = SequenceMetadata {
                        name: id,
                        description: current_description.take(),
                        length,
                        sha512t24u: sha512,
                        md5,
                        alphabet,
                        fai,
                    };

                    if opts.store_sequence {
                        seq_results.push(SequenceRecord::Full {
                            metadata,
                            sequence: std::mem::take(&mut sequence_data),
                        });
                    } else {
                        seq_results.push(SequenceRecord::Stub(metadata));
                    }
                } else {
                    // FAI-only mode
                    fai_results.push(FaiRecord {
                        name: id,
                        length,
                        fai,
                    });
                }
            }
            break;
        }

        if line.starts_with('>') {
            // Save previous sequence
            if let Some(id) = current_id.take() {
                let fai = make_fai(
                    opts.fai_enabled,
                    current_offset,
                    current_line_bases,
                    current_line_bytes,
                );

                if opts.compute_digests {
                    let sha512 =
                        base64_url::encode(&sha512_hasher.as_mut().unwrap().finalize_reset()[0..24]);
                    let md5 = format!("{:x}", md5_hasher.as_mut().unwrap().finalize_reset());
                    let alphabet = alphabet_guesser.as_mut().unwrap().guess();

                    let metadata = SequenceMetadata {
                        name: id,
                        description: current_description.take(),
                        length,
                        sha512t24u: sha512,
                        md5,
                        alphabet,
                        fai,
                    };

                    if opts.store_sequence {
                        seq_results.push(SequenceRecord::Full {
                            metadata,
                            sequence: std::mem::take(&mut sequence_data),
                        });
                    } else {
                        seq_results.push(SequenceRecord::Stub(metadata));
                    }
                } else {
                    fai_results.push(FaiRecord {
                        name: id,
                        length,
                        fai,
                    });
                }
            }

            // Parse header
            let (name, description) = parse_fasta_header(&line[1..]);
            current_id = Some(name);
            current_description = description;

            // Track position for FAI
            if opts.fai_enabled {
                byte_position += bytes_read as u64;
                current_offset = byte_position;
            }
            current_line_bases = None;
            current_line_bytes = None;
            length = 0;

            // Reset digest state
            if opts.compute_digests {
                *sha512_hasher.as_mut().unwrap() = Sha512::new();
                *md5_hasher.as_mut().unwrap() = Md5::new();
                *alphabet_guesser.as_mut().unwrap() = digest::AlphabetGuesser::new();
            }
        } else if current_id.is_some() && !line.trim().is_empty() {
            // Sequence line
            let trimmed = line.trim_end();
            let line_len_bytes = bytes_read as u32;
            let line_len_bases = trimmed.len() as u32;

            if current_line_bases.is_none() {
                current_line_bases = Some(line_len_bases);
                current_line_bytes = Some(line_len_bytes);
            }

            length += trimmed.len();

            if opts.compute_digests || opts.store_sequence {
                let seq_upper = trimmed.to_ascii_uppercase();

                if opts.compute_digests {
                    sha512_hasher.as_mut().unwrap().update(seq_upper.as_bytes());
                    md5_hasher.as_mut().unwrap().update(seq_upper.as_bytes());
                    alphabet_guesser
                        .as_mut()
                        .unwrap()
                        .update(seq_upper.as_bytes());
                }

                if opts.store_sequence {
                    sequence_data.extend_from_slice(seq_upper.as_bytes());
                }
            }

            if opts.fai_enabled {
                byte_position += bytes_read as u64;
            }
        } else if opts.fai_enabled {
            byte_position += bytes_read as u64;
        }

        line.clear();
    }

    if !opts.compute_digests {
        Ok(FileParseResult::Fai(fai_results))
    } else {
        Ok(FileParseResult::Sequences(seq_results))
    }
}

/// Helper to construct FAI metadata
fn make_fai(
    fai_enabled: bool,
    offset: u64,
    line_bases: Option<u32>,
    line_bytes: Option<u32>,
) -> Option<FaiMetadata> {
    if !fai_enabled {
        return None;
    }
    match (line_bases, line_bytes) {
        (Some(lb), Some(lby)) => Some(FaiMetadata {
            offset,
            line_bases: lb,
            line_bytes: lby,
        }),
        _ => None,
    }
}

/// Processes a FASTA file to compute the digests of each sequence.
///
/// This function reads a FASTA file, computes the SHA-512 and MD5 digests for each sequence,
/// and returns a SequenceCollection. It handles both plain and gzipped FASTA files.
///
/// # Arguments
///
/// * `file_path` - Path to the FASTA file to be processed
///
/// # Returns
///
/// A SequenceCollection containing sequence metadata (Stub records, no sequence data).
///
/// # Examples
///
/// ```
/// use gtars_refget::fasta::digest_fasta;
///
/// let seqcol = digest_fasta("../tests/data/fasta/base.fa")
///     .expect("Failed to digest FASTA file");
/// println!("Digest: {}", seqcol.metadata.digest);
/// ```
pub fn digest_fasta<T: AsRef<Path>>(file_path: T) -> Result<SequenceCollection> {
    let results = match parse_fasta_file(file_path.as_ref(), FileParseOptions::DIGEST_ONLY)? {
        FileParseResult::Sequences(seqs) => seqs,
        _ => unreachable!("DIGEST_ONLY mode returns Sequences"),
    };

    let collection_metadata =
        SequenceCollectionMetadata::from_sequences(&results, Some(file_path.as_ref().to_path_buf()));

    Ok(SequenceCollection {
        metadata: collection_metadata,
        sequences: results,
    })
}

/// Computes FAI (FASTA index) metadata for sequences in a FASTA file without computing digests.
///
/// This is a lightweight alternative to `digest_fasta()` for cases where only FAI metadata
/// is needed. It skips the expensive digest computation.
///
/// # Arguments
///
/// * `file_path` - Path to the FASTA file to index
///
/// # Returns
///
/// A vector of `FaiRecord` structs containing sequence name, length, and FAI metadata.
/// For gzipped files, the `fai` field will be `None`.
///
/// # Examples
///
/// ```
/// use gtars_refget::fasta::compute_fai;
///
/// let fai_records = compute_fai("../tests/data/fasta/base.fa")
///     .expect("Failed to compute FAI");
/// for record in fai_records {
///     println!("{}: {} bp", record.name, record.length);
/// }
/// ```
pub fn compute_fai<T: AsRef<Path>>(file_path: T) -> Result<Vec<FaiRecord>> {
    match parse_fasta_file(file_path, FileParseOptions::FAI_ONLY)? {
        FileParseResult::Fai(records) => Ok(records),
        _ => unreachable!("FAI_ONLY mode returns Fai"),
    }
}

/// Loads a FASTA file with sequence data in memory.
///
/// This function reads a FASTA file and returns a SequenceCollection where each
/// SequenceRecord contains both metadata AND the actual sequence data.
///
/// # Arguments
///
/// * `file_path` - Path to the FASTA file to be loaded
///
/// # Returns
///
/// A SequenceCollection with Full records (metadata + sequence data).
///
/// # Examples
///
/// ```
/// use gtars_refget::fasta::load_fasta;
///
/// let seqcol = load_fasta("../tests/data/fasta/base.fa")
///     .expect("Failed to load FASTA file");
/// for seq in &seqcol.sequences {
///     println!("{}: {} bp", seq.metadata().name, seq.metadata().length);
/// }
/// ```
pub fn load_fasta<P: AsRef<Path>>(file_path: P) -> Result<SequenceCollection> {
    let results = match parse_fasta_file(file_path.as_ref(), FileParseOptions::FULL)? {
        FileParseResult::Sequences(seqs) => seqs,
        _ => unreachable!("FULL mode returns Sequences"),
    };

    let collection_metadata =
        SequenceCollectionMetadata::from_sequences(&results, Some(file_path.as_ref().to_path_buf()));

    Ok(SequenceCollection {
        metadata: collection_metadata,
        sequences: results,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::AlphabetType;

    #[test]
    fn test_digest_fasta() {
        let seqcol = digest_fasta("../tests/data/fasta/base.fa").expect("Can't open test fasta file");
        let results = &seqcol.sequences;
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].metadata().length, 8);
        assert_eq!(
            results[0].metadata().sha512t24u,
            "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        );
        assert_eq!(results[0].metadata().md5, "5f63cfaa3ef61f88c9635fb9d18ec945");
        assert_eq!(results[0].metadata().alphabet, AlphabetType::Dna2bit);
    }

    #[test]
    fn test_digest_fasta_gzipped() {
        let seqcol =
            digest_fasta("../tests/data/fasta/base.fa.gz").expect("Can't open test fasta file");
        let results = &seqcol.sequences;
        assert_eq!(results[0].metadata().length, 8);
        assert_eq!(
            results[0].metadata().sha512t24u,
            "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        );
        // Gzipped files should have no FAI data
        assert!(results[0].metadata().fai.is_none());
    }

    #[test]
    fn test_load_fasta() {
        let seqcol = load_fasta("../tests/data/fasta/base.fa").expect("Failed to load FASTA");
        assert_eq!(seqcol.sequences.len(), 3);
        assert!(seqcol.sequences[0].is_loaded());
        assert_eq!(
            seqcol.sequences[0].sequence().unwrap(),
            b"TTGGGGAA"
        );
    }

    #[test]
    fn test_compute_fai() {
        let fai_records = compute_fai("../tests/data/fasta/base.fa").expect("Failed to compute FAI");
        assert_eq!(fai_records.len(), 3);
        assert_eq!(fai_records[0].name, "chrX");
        assert_eq!(fai_records[0].length, 8);
        assert!(fai_records[0].fai.is_some());
    }

    #[test]
    fn test_fai_matches_samtools() {
        // Parse the samtools-generated FAI file
        fn parse_fai_file(fai_path: &str) -> Vec<(String, u64, u64, u32, u32)> {
            use std::fs::File;
            use std::io::{BufRead, BufReader};

            let file = File::open(fai_path).expect("Failed to open FAI file");
            let reader = BufReader::new(file);
            let mut results = Vec::new();

            for line in reader.lines() {
                let line = line.expect("Failed to read line");
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 5 {
                    let name = parts[0].to_string();
                    let length = parts[1].parse::<u64>().expect("Failed to parse length");
                    let offset = parts[2].parse::<u64>().expect("Failed to parse offset");
                    let line_bases = parts[3].parse::<u32>().expect("Failed to parse line_bases");
                    let line_bytes = parts[4].parse::<u32>().expect("Failed to parse line_bytes");
                    results.push((name, length, offset, line_bases, line_bytes));
                }
            }
            results
        }

        let seqcol = digest_fasta("../tests/data/fasta/base.fa").expect("Failed to digest FASTA file");
        let fai_data = parse_fai_file("../tests/data/fasta/base.fa.fai");

        assert_eq!(seqcol.sequences.len(), fai_data.len());

        for (seq, (name, length, offset, line_bases, line_bytes)) in
            seqcol.sequences.iter().zip(fai_data.iter())
        {
            assert_eq!(seq.metadata().name, *name);
            assert_eq!(seq.metadata().length, *length as usize);
            assert!(seq.metadata().fai.is_some(), "FAI data should be present");
            let fai = seq.metadata().fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(
                fai.line_bases, *line_bases,
                "Line bases mismatch for {}",
                name
            );
            assert_eq!(
                fai.line_bytes, *line_bytes,
                "Line bytes mismatch for {}",
                name
            );
        }
    }
}
