use anyhow::Result;
use std::io::BufRead;
use std::path::Path;

use super::alphabet::AlphabetGuesser;
use super::collection::{FaiMetadata, SequenceCollection, SequenceCollectionMetadata, SequenceMetadata, SequenceRecord};

use md5::Md5;
use sha2::{Digest, Sha512};

/// Helper to construct FAI metadata when enabled and line info is available.
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

/// Configuration for FASTA parsing to avoid unnecessary work.
#[derive(Clone, Copy)]
struct ParseOptions {
    /// Compute SHA512, MD5 digests and alphabet (requires uppercase conversion)
    compute_digests: bool,
    /// Store sequence data in memory
    store_sequence: bool,
}

impl ParseOptions {
    const FAI_ONLY: Self = Self { compute_digests: false, store_sequence: false };
    const DIGEST_ONLY: Self = Self { compute_digests: true, store_sequence: false };
    const FULL: Self = Self { compute_digests: true, store_sequence: true };
}

/// Result of internal FASTA parsing - either FaiRecords or SequenceRecords
enum ParseResult {
    Fai(Vec<FaiRecord>),
    Sequences(Vec<SequenceRecord>),
}

/// Unified internal FASTA parser. Conditionally performs work based on options.
fn parse_fasta_internal<P: AsRef<Path>>(file_path: P, opts: ParseOptions) -> Result<ParseResult> {
    use gtars_core::utils::get_dynamic_reader;

    // Detect if file is gzipped
    let is_gzipped = file_path.as_ref().extension()
        .and_then(|s| s.to_str())
        .map(|s| s == "gz")
        .unwrap_or(false);

    // Skip FAI data for gzipped files (can't seek into compressed files)
    let fai_enabled = !is_gzipped;

    let mut byte_position: u64 = 0;

    let file_reader = get_dynamic_reader(file_path.as_ref())?;
    let mut reader = std::io::BufReader::new(file_reader);
    let mut line = String::new();

    // Common state
    let mut current_id: Option<String> = None;
    let mut current_description: Option<String> = None;
    let mut current_offset: u64 = 0;
    let mut current_line_bases: Option<u32> = None;
    let mut current_line_bytes: Option<u32> = None;
    let mut length = 0;

    // Digest state (only initialized if needed)
    let mut sha512_hasher = if opts.compute_digests { Some(Sha512::new()) } else { None };
    let mut md5_hasher = if opts.compute_digests { Some(Md5::new()) } else { None };
    let mut alphabet_guesser = if opts.compute_digests { Some(AlphabetGuesser::new()) } else { None };

    // Sequence data (only if storing)
    let mut sequence_data: Vec<u8> = Vec::new();

    // Results
    let mut fai_results: Vec<FaiRecord> = Vec::new();
    let mut seq_results: Vec<SequenceRecord> = Vec::new();

    loop {
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            // EOF - finalize the last sequence if any
            if let Some(id) = current_id.take() {
                let fai = make_fai(fai_enabled, current_offset, current_line_bases, current_line_bytes);

                if opts.compute_digests {
                    let sha512 = base64_url::encode(&sha512_hasher.as_mut().unwrap().finalize_reset()[0..24]);
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
                    fai_results.push(FaiRecord { name: id, length, fai });
                }
            }
            break;
        }

        if line.starts_with('>') {
            // Save previous sequence if any
            if let Some(id) = current_id.take() {
                let fai = make_fai(fai_enabled, current_offset, current_line_bases, current_line_bytes);

                if opts.compute_digests {
                    let sha512 = base64_url::encode(&sha512_hasher.as_mut().unwrap().finalize_reset()[0..24]);
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
                    fai_results.push(FaiRecord { name: id, length, fai });
                }
            }

            // Start new sequence - parse header
            let (name, description) = parse_fasta_header(&line[1..]);
            current_id = Some(name);
            current_description = description;

            // Track position for FAI
            if fai_enabled {
                byte_position += bytes_read as u64;
                current_offset = byte_position;
            }
            current_line_bases = None;
            current_line_bytes = None;
            length = 0;

            // Reset digest state if needed
            if opts.compute_digests {
                *sha512_hasher.as_mut().unwrap() = Sha512::new();
                *md5_hasher.as_mut().unwrap() = Md5::new();
                *alphabet_guesser.as_mut().unwrap() = AlphabetGuesser::new();
            }
        } else if current_id.is_some() && !line.trim().is_empty() {
            // Sequence line
            let trimmed = line.trim_end();
            let line_len_bytes = bytes_read as u32;
            let line_len_bases = trimmed.len() as u32;

            // Record line dimensions from first sequence line
            if current_line_bases.is_none() {
                current_line_bases = Some(line_len_bases);
                current_line_bytes = Some(line_len_bytes);
            }

            length += trimmed.len();

            // Only do expensive work if needed
            if opts.compute_digests || opts.store_sequence {
                let seq_upper = trimmed.to_ascii_uppercase();

                if opts.compute_digests {
                    sha512_hasher.as_mut().unwrap().update(seq_upper.as_bytes());
                    md5_hasher.as_mut().unwrap().update(seq_upper.as_bytes());
                    alphabet_guesser.as_mut().unwrap().update(seq_upper.as_bytes());
                }

                if opts.store_sequence {
                    sequence_data.extend_from_slice(seq_upper.as_bytes());
                }
            }

            // Track position for FAI
            if fai_enabled {
                byte_position += bytes_read as u64;
            }
        } else {
            // Track position for empty lines or other content
            if fai_enabled {
                byte_position += bytes_read as u64;
            }
        }

        line.clear();
    }

    if opts.compute_digests {
        Ok(ParseResult::Sequences(seq_results))
    } else {
        Ok(ParseResult::Fai(fai_results))
    }
}

/// Parse a FASTA header line (without the leading '>') into name and description.
///
/// Following FASTA standard: the sequence ID is the first word (up to first whitespace),
/// and everything after is the description.
///
/// # Arguments
/// * `header` - The header text (everything after '>')
///
/// # Returns
/// Tuple of (name, description) where description is None if no whitespace in header.
///
/// # Examples
/// ```ignore
/// let (name, desc) = parse_fasta_header("chr1 some description here");
/// assert_eq!(name, "chr1");
/// assert_eq!(desc, Some("some description here".to_string()));
///
/// let (name, desc) = parse_fasta_header("chr1");
/// assert_eq!(name, "chr1");
/// assert_eq!(desc, None);
/// ```
pub fn parse_fasta_header(header: &str) -> (String, Option<String>) {
    let header = header.trim();
    match header.split_once(char::is_whitespace) {
        Some((id, desc)) => (id.to_string(), Some(desc.trim().to_string())),
        None => (header.to_string(), None),
    }
}

/// A lightweight record containing only FAI (FASTA index) metadata for a sequence.
/// Returned by `compute_fai()` for fast FAI-only computation without digest overhead.
#[derive(Clone, Debug)]
pub struct FaiRecord {
    pub name: String,
    pub length: usize,
    pub fai: Option<FaiMetadata>,
}

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
    let results = match parse_fasta_internal(file_path.as_ref(), ParseOptions::DIGEST_ONLY)? {
        ParseResult::Sequences(seqs) => seqs,
        ParseResult::Fai(_) => unreachable!("DIGEST_ONLY mode returns Sequences"),
    };

    let collection_metadata = SequenceCollectionMetadata::from_sequences(
        &results,
        Some(file_path.as_ref().to_path_buf())
    );

    Ok(SequenceCollection {
        metadata: collection_metadata,
        sequences: results,
    })
}

/// Computes FAI (FASTA index) metadata for sequences in a FASTA file without computing digests.
///
/// This is a lightweight alternative to `digest_fasta()` for cases where only FAI metadata
/// is needed. It skips the expensive digest computation (SHA512, MD5) and alphabet detection,
/// making it much faster for FAI-only indexing (similar to `samtools faidx`).
///
/// # Arguments
///
/// * `file_path` - Path to the FASTA file to index
///
/// # Returns
///
/// A vector of `FaiRecord` structs containing sequence name, length, and FAI metadata.
/// For gzipped files, the `fai` field will be `None` since byte offsets are not meaningful
/// for compressed data.
///
/// # Examples
///
/// ```
/// use gtars_refget::fasta::compute_fai;
///
/// // Using a test file from the repository
/// let fai_records = compute_fai("../tests/data/fasta/base.fa").expect("Failed to compute FAI");
/// assert_eq!(fai_records.len(), 3);
///
/// for record in fai_records {
///     println!("{}: {} bp", record.name, record.length);
///     if let Some(fai) = record.fai {
///         println!("  offset: {}, line_bases: {}, line_bytes: {}",
///                  fai.offset, fai.line_bases, fai.line_bytes);
///     }
/// }
/// ```
pub fn compute_fai<T: AsRef<Path>>(file_path: T) -> Result<Vec<FaiRecord>> {
    match parse_fasta_internal(file_path, ParseOptions::FAI_ONLY)? {
        ParseResult::Fai(records) => Ok(records),
        ParseResult::Sequences(_) => unreachable!("FAI_ONLY mode returns Fai"),
    }
}

/// Loads a FASTA file with sequence data in memory (not just metadata).
///
/// This function reads a FASTA file and returns a SequenceCollection where each
/// SequenceRecord contains both metadata AND the actual sequence data. This is useful
/// when you need direct access to sequences without going through a RefgetStore.
///
/// # Arguments
///
/// * `file_path` - A path to the FASTA file to be loaded.
///
/// # Returns
///
/// A SequenceCollection with:
/// - `metadata`: Complete SequenceMetadata (name, length, digests, alphabet)
/// - `fai`: FAI metadata if not gzipped
/// - `data`: Some(Vec<u8>) containing the raw sequence bytes (uppercased)
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
/// ```
/// use gtars_refget::fasta::load_fasta;
///
/// let seqcol = load_fasta("../tests/data/fasta/base.fa").expect("Failed to load FASTA");
/// for seq_record in &seqcol.sequences {
///     println!("Sequence: {}", seq_record.metadata().name);
///     if let Some(data) = seq_record.sequence() {
///         println!("  Data: {} bytes", data.len());
///     }
/// }
/// ```
pub fn load_fasta<P: AsRef<Path>>(file_path: P) -> Result<SequenceCollection> {
    let results = match parse_fasta_internal(file_path.as_ref(), ParseOptions::FULL)? {
        ParseResult::Sequences(seqs) => seqs,
        ParseResult::Fai(_) => unreachable!("FULL mode returns Sequences"),
    };

    let collection_metadata = SequenceCollectionMetadata::from_sequences(
        &results,
        Some(file_path.as_ref().to_path_buf())
    );

    Ok(SequenceCollection {
        metadata: collection_metadata,
        sequences: results,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::AlphabetType;
    use crate::collection::read_rgsi_file;

    #[test]
    fn digests_digest_fasta() {
        let seqcol =
            digest_fasta("../tests/data/fasta/base.fa").expect("Can't open test fasta file");
        let results = seqcol.sequences;
        println!("{:?}", results);
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].metadata().length, 8);
        assert_eq!(
            results[0].metadata().sha512t24u,
            "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        );
        assert_eq!(results[0].metadata().md5, "5f63cfaa3ef61f88c9635fb9d18ec945");
        assert_eq!(results[0].metadata().alphabet, AlphabetType::Dna2bit);
        assert_eq!(results[1].metadata().length, 4);
        assert_eq!(
            results[1].metadata().sha512t24u,
            "YBbVX0dLKG1ieEDCiMmkrTZFt_Z5Vdaj"
        );
        assert_eq!(results[1].metadata().md5, "31fc6ca291a32fb9df82b85e5f077e31");
        assert_eq!(results[1].metadata().alphabet, AlphabetType::Dna2bit);
        assert_eq!(results[2].metadata().length, 4);
        assert_eq!(
            results[2].metadata().sha512t24u,
            "AcLxtBuKEPk_7PGE_H4dGElwZHCujwH6"
        );
        assert_eq!(results[2].metadata().md5, "92c6a56c9e9459d8a42b96f7884710bc");
        assert_eq!(results[2].metadata().alphabet, AlphabetType::Dna2bit);

        // Test FAI metadata
        assert!(results[0].metadata().fai.is_some());
        let fai0 = results[0].metadata().fai.as_ref().unwrap();
        assert_eq!(fai0.offset, 6);  // After ">chrX\n"
        assert_eq!(fai0.line_bases, 8);
        assert_eq!(fai0.line_bytes, 9);  // Includes \n

        assert!(results[1].metadata().fai.is_some());
        let fai1 = results[1].metadata().fai.as_ref().unwrap();
        assert_eq!(fai1.offset, 21);  // After ">chr1\n"
        assert_eq!(fai1.line_bases, 4);
        assert_eq!(fai1.line_bytes, 5);

        assert!(results[2].metadata().fai.is_some());
        let fai2 = results[2].metadata().fai.as_ref().unwrap();
        assert_eq!(fai2.offset, 32);  // After ">chr2\n"
        assert_eq!(fai2.line_bases, 4);
        assert_eq!(fai2.line_bytes, 5);
    }

    #[test]
    fn digests_digest_gzipped_fasta() {
        let seqcol =
            digest_fasta("../tests/data/fasta/base.fa.gz").expect("Can't open test fasta file");
        let results = seqcol.sequences;
        println!("{:?}", results);
        assert_eq!(results[0].metadata().length, 8);
        assert_eq!(
            results[0].metadata().sha512t24u,
            "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        );
        assert_eq!(results[0].metadata().md5, "5f63cfaa3ef61f88c9635fb9d18ec945");
        assert_eq!(results[0].metadata().alphabet, AlphabetType::Dna2bit);
    }

    #[test]
    fn digests_bogus_fasta_file() {
        let result = digest_fasta("../tests/data/bogus.fa");
        assert!(result.is_err(), "Expected an error for a bogus fasta file");
    }

    #[test]
    fn digests_fa_to_rgsi() {
        let seqcol = SequenceCollection::from_path_no_cache("../tests/data/fasta/base.fa")
            .expect("Failed to create SequenceCollection from FASTA file");
        seqcol.write_rgsi().expect("Failed to write rgsi file");

        let loaded_seqcol = read_rgsi_file("../tests/data/fasta/base.rgsi")
            .expect("Failed to read RGSI file");
        println!("Original SequenceCollection: {}", seqcol);
        println!("Loaded SequenceCollection: {}", loaded_seqcol);
        // Test round-trip integrity
        for (original, read) in seqcol.sequences.iter().zip(loaded_seqcol.sequences.iter()) {
            assert_eq!(original.metadata().name, read.metadata().name);
            assert_eq!(original.metadata().length, read.metadata().length);
            assert_eq!(original.metadata().sha512t24u, read.metadata().sha512t24u);
            assert_eq!(original.metadata().md5, read.metadata().md5);
            assert_eq!(original.metadata().alphabet, read.metadata().alphabet);
        }
    }

    #[test]
    fn digests_seqcol_from_fasta() {
        let seqcol = SequenceCollection::from_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to create SequenceCollection from FASTA file");
        println!("SeqCol digest: {:?}", seqcol.metadata.digest);
        println!(
            "SeqCol sequences_digest: {:?}",
            seqcol.metadata.sequences_digest
        );
        println!("SequenceCollection: {:?}", seqcol);
        assert_eq!(
            seqcol.metadata.sequences_digest,
            "0uDQVLuHaOZi1u76LjV__yrVUIz9Bwhr"
        );
        assert_eq!(seqcol.metadata.names_digest, "Fw1r9eRxfOZD98KKrhlYQNEdSRHoVxAG");
        assert_eq!(
            seqcol.metadata.lengths_digest,
            "cGRMZIb3AVgkcAfNv39RN7hnT5Chk7RX"
        );
    }

    // Helper function to parse samtools FAI file
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

    #[test]
    fn fai_wrapped_40() {
        let seqcol = digest_fasta("../tests/data/fasta/wrapped_40.fa")
            .expect("Failed to digest FASTA file");
        let fai_data = parse_fai_file("../tests/data/fasta/wrapped_40.fa.fai");

        assert_eq!(seqcol.sequences.len(), fai_data.len());

        for (seq, (name, length, offset, line_bases, line_bytes)) in
            seqcol.sequences.iter().zip(fai_data.iter())
        {
            assert_eq!(seq.metadata().name, *name);
            assert_eq!(seq.metadata().length, *length as usize);
            assert!(seq.metadata().fai.is_some(), "FAI data should be present");
            let fai = seq.metadata().fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(fai.line_bases, *line_bases, "Line bases mismatch for {}", name);
            assert_eq!(fai.line_bytes, *line_bytes, "Line bytes mismatch for {}", name);
        }
    }

    #[test]
    fn fai_wrapped_20() {
        let seqcol = digest_fasta("../tests/data/fasta/wrapped_20.fa")
            .expect("Failed to digest FASTA file");
        let fai_data = parse_fai_file("../tests/data/fasta/wrapped_20.fa.fai");

        assert_eq!(seqcol.sequences.len(), fai_data.len());

        for (seq, (name, length, offset, line_bases, line_bytes)) in
            seqcol.sequences.iter().zip(fai_data.iter())
        {
            assert_eq!(seq.metadata().name, *name);
            assert_eq!(seq.metadata().length, *length as usize);
            assert!(seq.metadata().fai.is_some(), "FAI data should be present");
            let fai = seq.metadata().fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(fai.line_bases, *line_bases, "Line bases mismatch for {}", name);
            assert_eq!(fai.line_bytes, *line_bytes, "Line bytes mismatch for {}", name);
        }
    }

    #[test]
    fn fai_unwrapped() {
        let seqcol = digest_fasta("../tests/data/fasta/unwrapped.fa")
            .expect("Failed to digest FASTA file");
        let fai_data = parse_fai_file("../tests/data/fasta/unwrapped.fa.fai");

        assert_eq!(seqcol.sequences.len(), fai_data.len());

        for (seq, (name, length, offset, line_bases, line_bytes)) in
            seqcol.sequences.iter().zip(fai_data.iter())
        {
            assert_eq!(seq.metadata().name, *name);
            assert_eq!(seq.metadata().length, *length as usize);
            assert!(seq.metadata().fai.is_some(), "FAI data should be present");
            let fai = seq.metadata().fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(fai.line_bases, *line_bases, "Line bases mismatch for {}", name);
            assert_eq!(fai.line_bytes, *line_bytes, "Line bytes mismatch for {}", name);
        }
    }

    #[test]
    fn fai_base_mixed_lengths() {
        // base.fa has mixed sequence lengths (8bp, 4bp, 4bp) on single lines
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to digest FASTA file");
        let fai_data = parse_fai_file("../tests/data/fasta/base.fa.fai");

        assert_eq!(seqcol.sequences.len(), fai_data.len());
        assert_eq!(seqcol.sequences.len(), 3, "base.fa should have 3 sequences");

        for (seq, (name, length, offset, line_bases, line_bytes)) in
            seqcol.sequences.iter().zip(fai_data.iter())
        {
            assert_eq!(seq.metadata().name, *name);
            assert_eq!(seq.metadata().length, *length as usize);
            assert!(seq.metadata().fai.is_some(), "FAI data should be present");
            let fai = seq.metadata().fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(fai.line_bases, *line_bases, "Line bases mismatch for {}", name);
            assert_eq!(fai.line_bytes, *line_bytes, "Line bytes mismatch for {}", name);
        }
    }

    #[test]
    fn fai_gzipped_returns_none() {
        let seqcol = digest_fasta("../tests/data/fasta/base.fa.gz")
            .expect("Failed to digest gzipped FASTA file");

        // All sequences should have None for FAI data
        for seq in seqcol.sequences.iter() {
            assert!(seq.metadata().fai.is_none(),
                "Gzipped files should not have FAI data for sequence {}",
                seq.metadata().name);
        }
    }

    #[test]
    fn fai_crlf_endings() {
        let seqcol = digest_fasta("../tests/data/fasta/crlf_endings.fa")
            .expect("Failed to digest FASTA file with CRLF endings");
        let fai_data = parse_fai_file("../tests/data/fasta/crlf_endings.fa.fai");

        assert_eq!(seqcol.sequences.len(), fai_data.len());

        for (seq, (name, length, offset, line_bases, line_bytes)) in
            seqcol.sequences.iter().zip(fai_data.iter())
        {
            assert_eq!(seq.metadata().name, *name);
            assert_eq!(seq.metadata().length, *length as usize);
            assert!(seq.metadata().fai.is_some(), "FAI data should be present");
            let fai = seq.metadata().fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(fai.line_bases, *line_bases, "Line bases mismatch for {}", name);
            assert_eq!(fai.line_bytes, *line_bytes, "Line bytes mismatch for {} (should be +2 for CRLF)", name);
            // Verify CRLF: line_bytes should be line_bases + 2
            assert_eq!(fai.line_bytes, fai.line_bases + 2,
                "CRLF files should have line_bytes = line_bases + 2 for {}",
                name);
        }
    }

    #[test]
    fn fai_same_sequence_different_wrapping() {
        // Test that wrapped_40, wrapped_20, and unwrapped all contain identical sequences
        // despite different line wrapping styles
        let seqcol_40 = digest_fasta("../tests/data/fasta/wrapped_40.fa")
            .expect("Failed to digest wrapped_40.fa");
        let seqcol_20 = digest_fasta("../tests/data/fasta/wrapped_20.fa")
            .expect("Failed to digest wrapped_20.fa");
        let seqcol_unwrapped = digest_fasta("../tests/data/fasta/unwrapped.fa")
            .expect("Failed to digest unwrapped.fa");

        // All should have same number of sequences
        assert_eq!(seqcol_40.sequences.len(), seqcol_20.sequences.len());
        assert_eq!(seqcol_40.sequences.len(), seqcol_unwrapped.sequences.len());

        // Compare each sequence across all three files
        for ((seq40, seq20), seq_unwrap) in seqcol_40.sequences.iter()
            .zip(seqcol_20.sequences.iter())
            .zip(seqcol_unwrapped.sequences.iter())
        {
            // Same name
            assert_eq!(seq40.metadata().name, seq20.metadata().name);
            assert_eq!(seq40.metadata().name, seq_unwrap.metadata().name);

            // Same length
            assert_eq!(seq40.metadata().length, seq20.metadata().length,
                "Length mismatch for {}", seq40.metadata().name);
            assert_eq!(seq40.metadata().length, seq_unwrap.metadata().length,
                "Length mismatch for {}", seq40.metadata().name);

            // Same digests (proves sequences are identical)
            assert_eq!(seq40.metadata().sha512t24u, seq20.metadata().sha512t24u,
                "SHA512 digest mismatch for {} - sequences not identical",
                seq40.metadata().name);
            assert_eq!(seq40.metadata().sha512t24u, seq_unwrap.metadata().sha512t24u,
                "SHA512 digest mismatch for {} - sequences not identical",
                seq40.metadata().name);
            assert_eq!(seq40.metadata().md5, seq20.metadata().md5,
                "MD5 digest mismatch for {} - sequences not identical",
                seq40.metadata().name);
            assert_eq!(seq40.metadata().md5, seq_unwrap.metadata().md5,
                "MD5 digest mismatch for {} - sequences not identical",
                seq40.metadata().name);

            // But different FAI data (different wrapping)
            assert!(seq40.metadata().fai.is_some() && seq20.metadata().fai.is_some() && seq_unwrap.metadata().fai.is_some());
            let fai40 = seq40.metadata().fai.as_ref().unwrap();
            let fai20 = seq20.metadata().fai.as_ref().unwrap();
            let fai_unwrap = seq_unwrap.metadata().fai.as_ref().unwrap();

            // Different line_bases due to different wrapping (for multi-line sequences)
            // Note: short sequences like chr3 that fit on one line may have same line_bases
            // across different wrapping styles, which is correct behavior
            if seq40.metadata().length > 40 {
                // Long enough to require wrapping in wrapped_40
                assert_ne!(fai40.line_bases, fai20.line_bases,
                    "wrapped_40 and wrapped_20 should have different line_bases for {}",
                    seq40.metadata().name);
            }
            if seq40.metadata().length > 20 {
                // Long enough to require wrapping in wrapped_20
                assert_ne!(fai20.line_bases, fai_unwrap.line_bases,
                    "wrapped_20 and unwrapped should have different line_bases for {}",
                    seq40.metadata().name);
            }
        }
    }

    // Tests for compute_fai() function

    #[test]
    fn test_compute_fai_wrapped_40() {
        let fai_records = compute_fai("../tests/data/fasta/wrapped_40.fa")
            .expect("Failed to compute FAI");
        let fai_data = parse_fai_file("../tests/data/fasta/wrapped_40.fa.fai");

        assert_eq!(fai_records.len(), fai_data.len());

        for (record, (name, length, offset, line_bases, line_bytes)) in
            fai_records.iter().zip(fai_data.iter())
        {
            assert_eq!(record.name, *name);
            assert_eq!(record.length, *length as usize);
            assert!(record.fai.is_some(), "FAI data should be present");
            let fai = record.fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(fai.line_bases, *line_bases, "Line bases mismatch for {}", name);
            assert_eq!(fai.line_bytes, *line_bytes, "Line bytes mismatch for {}", name);
        }
    }

    #[test]
    fn test_compute_fai_wrapped_20() {
        let fai_records = compute_fai("../tests/data/fasta/wrapped_20.fa")
            .expect("Failed to compute FAI");
        let fai_data = parse_fai_file("../tests/data/fasta/wrapped_20.fa.fai");

        assert_eq!(fai_records.len(), fai_data.len());

        for (record, (name, length, offset, line_bases, line_bytes)) in
            fai_records.iter().zip(fai_data.iter())
        {
            assert_eq!(record.name, *name);
            assert_eq!(record.length, *length as usize);
            assert!(record.fai.is_some(), "FAI data should be present");
            let fai = record.fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(fai.line_bases, *line_bases, "Line bases mismatch for {}", name);
            assert_eq!(fai.line_bytes, *line_bytes, "Line bytes mismatch for {}", name);
        }
    }

    #[test]
    fn test_compute_fai_unwrapped() {
        let fai_records = compute_fai("../tests/data/fasta/unwrapped.fa")
            .expect("Failed to compute FAI");
        let fai_data = parse_fai_file("../tests/data/fasta/unwrapped.fa.fai");

        assert_eq!(fai_records.len(), fai_data.len());

        for (record, (name, length, offset, line_bases, line_bytes)) in
            fai_records.iter().zip(fai_data.iter())
        {
            assert_eq!(record.name, *name);
            assert_eq!(record.length, *length as usize);
            assert!(record.fai.is_some(), "FAI data should be present");
            let fai = record.fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(fai.line_bases, *line_bases, "Line bases mismatch for {}", name);
            assert_eq!(fai.line_bytes, *line_bytes, "Line bytes mismatch for {}", name);
        }
    }

    #[test]
    fn test_compute_fai_base() {
        let fai_records = compute_fai("../tests/data/fasta/base.fa")
            .expect("Failed to compute FAI");
        let fai_data = parse_fai_file("../tests/data/fasta/base.fa.fai");

        assert_eq!(fai_records.len(), fai_data.len());
        assert_eq!(fai_records.len(), 3, "base.fa should have 3 sequences");

        for (record, (name, length, offset, line_bases, line_bytes)) in
            fai_records.iter().zip(fai_data.iter())
        {
            assert_eq!(record.name, *name);
            assert_eq!(record.length, *length as usize);
            assert!(record.fai.is_some(), "FAI data should be present");
            let fai = record.fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(fai.line_bases, *line_bases, "Line bases mismatch for {}", name);
            assert_eq!(fai.line_bytes, *line_bytes, "Line bytes mismatch for {}", name);
        }
    }

    #[test]
    fn test_compute_fai_crlf() {
        let fai_records = compute_fai("../tests/data/fasta/crlf_endings.fa")
            .expect("Failed to compute FAI");
        let fai_data = parse_fai_file("../tests/data/fasta/crlf_endings.fa.fai");

        assert_eq!(fai_records.len(), fai_data.len());

        for (record, (name, length, offset, line_bases, line_bytes)) in
            fai_records.iter().zip(fai_data.iter())
        {
            assert_eq!(record.name, *name);
            assert_eq!(record.length, *length as usize);
            assert!(record.fai.is_some(), "FAI data should be present");
            let fai = record.fai.as_ref().unwrap();
            assert_eq!(fai.offset, *offset, "Offset mismatch for {}", name);
            assert_eq!(fai.line_bases, *line_bases, "Line bases mismatch for {}", name);
            assert_eq!(fai.line_bytes, *line_bytes, "Line bytes mismatch for {} (should be +2 for CRLF)", name);
            // Verify CRLF: line_bytes should be line_bases + 2
            assert_eq!(fai.line_bytes, fai.line_bases + 2,
                "CRLF files should have line_bytes = line_bases + 2 for {}",
                name);
        }
    }

    #[test]
    fn test_compute_fai_gzipped() {
        let fai_records = compute_fai("../tests/data/fasta/base.fa.gz")
            .expect("Failed to compute FAI for gzipped file");

        // All sequences should have None for FAI data
        for record in fai_records.iter() {
            assert!(record.fai.is_none(),
                "Gzipped files should not have FAI data for sequence {}",
                record.name);
        }
    }

    #[test]
    fn test_load_fasta_with_data() {
        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load FASTA with data");

        assert_eq!(seqcol.sequences.len(), 3);

        // Check first sequence (chrX)
        let seq0 = &seqcol.sequences[0];
        assert_eq!(seq0.metadata().name, "chrX");
        assert_eq!(seq0.metadata().length, 8);
        assert!(seq0.is_loaded(), "Sequence should have data");

        if let Some(data) = seq0.sequence() {
            assert_eq!(data.len(), 8);
            assert_eq!(std::str::from_utf8(data).unwrap(), "TTGGGGAA");
        }

        // Check second sequence (chr1)
        let seq1 = &seqcol.sequences[1];
        assert_eq!(seq1.metadata().name, "chr1");
        assert_eq!(seq1.metadata().length, 4);
        assert!(seq1.is_loaded(), "Sequence should have data");

        if let Some(data) = seq1.sequence() {
            assert_eq!(data.len(), 4);
            assert_eq!(std::str::from_utf8(data).unwrap(), "GGAA");
        }

        // Check third sequence (chr2)
        let seq2 = &seqcol.sequences[2];
        assert_eq!(seq2.metadata().name, "chr2");
        assert_eq!(seq2.metadata().length, 4);
        assert!(seq2.is_loaded(), "Sequence should have data");

        if let Some(data) = seq2.sequence() {
            assert_eq!(data.len(), 4);
            assert_eq!(std::str::from_utf8(data).unwrap(), "GCGC");
        }

        // Verify digests match digest_fasta() output
        let digest_seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to digest FASTA");

        for (loaded, digested) in seqcol.sequences.iter().zip(digest_seqcol.sequences.iter()) {
            assert_eq!(loaded.metadata().name, digested.metadata().name);
            assert_eq!(loaded.metadata().length, digested.metadata().length);
            assert_eq!(loaded.metadata().sha512t24u, digested.metadata().sha512t24u);
            assert_eq!(loaded.metadata().md5, digested.metadata().md5);
            assert_eq!(loaded.metadata().alphabet, digested.metadata().alphabet);
        }
    }

    #[test]
    fn test_load_fasta_wrapped() {
        let seqcol = load_fasta("../tests/data/fasta/wrapped_40.fa")
            .expect("Failed to load wrapped FASTA");

        // Verify all sequences have data
        for seq in &seqcol.sequences {
            assert!(seq.is_loaded(), "Sequence {} should have data", seq.metadata().name);

            if let Some(data) = seq.sequence() {
                assert_eq!(data.len(), seq.metadata().length,
                    "Data length mismatch for {}", seq.metadata().name);

                // Verify data is uppercase
                let data_str = std::str::from_utf8(data).unwrap();
                assert_eq!(data_str, data_str.to_uppercase(),
                    "Sequence data should be uppercased");
            }
        }
    }

    #[test]
    fn test_load_fasta_gzipped() {
        let seqcol = load_fasta("../tests/data/fasta/base.fa.gz")
            .expect("Failed to load gzipped FASTA");

        // Should still load data from gzipped files
        assert_eq!(seqcol.sequences.len(), 3);
        assert!(seqcol.sequences[0].is_loaded(), "Should have sequence data");

        for seq in &seqcol.sequences {
            assert!(seq.is_loaded(), "Gzipped file should still have sequence data");
            assert!(seq.metadata().fai.is_none(), "Gzipped file should not have FAI metadata");
        }
    }
}
