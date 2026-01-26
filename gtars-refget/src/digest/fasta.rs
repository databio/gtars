//! Bytes-based FASTA parsing - WASM-safe.
//!
//! This module provides FASTA parsing functions that work with in-memory data
//! (byte slices) rather than file paths. This makes them WASM-compatible.

use anyhow::Result;
use md5::Md5;
use sha2::{Digest, Sha512};
use std::io::BufRead;

use super::alphabet::AlphabetGuesser;
use super::types::{
    FaiMetadata, SequenceCollection, SequenceCollectionMetadata, SequenceMetadata, SequenceRecord,
};

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
/// ```
/// use gtars_refget::digest::fasta::parse_fasta_header;
///
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
pub struct ParseOptions {
    /// Compute SHA512, MD5 digests and alphabet (requires uppercase conversion)
    pub compute_digests: bool,
    /// Store sequence data in memory
    pub store_sequence: bool,
    /// Track FAI metadata (only valid for uncompressed file reading)
    pub fai_enabled: bool,
}

impl ParseOptions {
    pub const DIGEST_ONLY: Self = Self {
        compute_digests: true,
        store_sequence: false,
        fai_enabled: false,
    };
    pub const FULL: Self = Self {
        compute_digests: true,
        store_sequence: true,
        fai_enabled: false,
    };
}

/// Core FASTA parser that works with any `BufRead` implementation.
/// This is WASM-compatible as it doesn't require filesystem access.
pub(crate) fn parse_fasta_reader<R: BufRead>(mut reader: R, opts: ParseOptions) -> Result<Vec<SequenceRecord>> {
    let fai_enabled = opts.fai_enabled;
    let mut byte_position: u64 = 0;
    let mut line = String::new();

    // Common state
    let mut current_id: Option<String> = None;
    let mut current_description: Option<String> = None;
    let mut current_offset: u64 = 0;
    let mut current_line_bases: Option<u32> = None;
    let mut current_line_bytes: Option<u32> = None;
    let mut length = 0;

    // Digest state (only initialized if needed)
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
        Some(AlphabetGuesser::new())
    } else {
        None
    };

    // Sequence data (only if storing)
    let mut sequence_data: Vec<u8> = Vec::new();

    // Results
    let mut seq_results: Vec<SequenceRecord> = Vec::new();

    loop {
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            // EOF - finalize the last sequence if any
            if let Some(id) = current_id.take() {
                let fai = make_fai(fai_enabled, current_offset, current_line_bases, current_line_bytes);

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
                }
            }
            break;
        }

        if line.starts_with('>') {
            // Save previous sequence if any
            if let Some(id) = current_id.take() {
                let fai = make_fai(fai_enabled, current_offset, current_line_bases, current_line_bytes);

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

    Ok(seq_results)
}

/// Compute seqcol digest from FASTA content bytes.
///
/// This is the WASM-friendly version of `digest_fasta()` that works
/// with in-memory data instead of file paths. It automatically detects
/// and decompresses gzip-compressed content.
///
/// # Arguments
///
/// * `content` - The FASTA file content as bytes
///
/// # Returns
///
/// A `SequenceCollection` containing all sequences with their computed digests.
///
/// # Examples
///
/// ```
/// use gtars_refget::digest::fasta::digest_fasta_bytes;
///
/// let fasta_content = b">chr1\nACGT\n>chr2\nTGCA\n";
/// let collection = digest_fasta_bytes(fasta_content).expect("Failed to digest");
/// assert_eq!(collection.sequences.len(), 2);
/// ```
pub fn digest_fasta_bytes(content: &[u8]) -> Result<SequenceCollection> {
    // Detect gzip by magic bytes [0x1f, 0x8b]
    let is_gzipped = content.len() >= 2 && content[0] == 0x1f && content[1] == 0x8b;

    let results = if is_gzipped {
        // Decompress gzip content
        use flate2::read::GzDecoder;
        let decoder = GzDecoder::new(content);
        let reader = std::io::BufReader::new(decoder);
        parse_fasta_reader(reader, ParseOptions::DIGEST_ONLY)?
    } else {
        // Parse directly from bytes
        let reader = std::io::BufReader::new(content);
        parse_fasta_reader(reader, ParseOptions::DIGEST_ONLY)?
    };

    let collection_metadata = SequenceCollectionMetadata::from_sequences(&results, None);

    Ok(SequenceCollection {
        metadata: collection_metadata,
        sequences: results,
    })
}

/// Load FASTA content with sequence data in memory.
///
/// This is the WASM-friendly version of `load_fasta()` that works
/// with in-memory data instead of file paths. It automatically detects
/// and decompresses gzip-compressed content.
///
/// # Arguments
///
/// * `content` - The FASTA file content as bytes
///
/// # Returns
///
/// A `SequenceCollection` containing all sequences with both metadata
/// and actual sequence data loaded.
///
/// # Examples
///
/// ```
/// use gtars_refget::digest::fasta::load_fasta_bytes;
///
/// let fasta_content = b">chr1\nACGT\n>chr2\nTGCA\n";
/// let collection = load_fasta_bytes(fasta_content).expect("Failed to load");
/// for seq in &collection.sequences {
///     assert!(seq.is_loaded());
/// }
/// ```
pub fn load_fasta_bytes(content: &[u8]) -> Result<SequenceCollection> {
    // Detect gzip by magic bytes [0x1f, 0x8b]
    let is_gzipped = content.len() >= 2 && content[0] == 0x1f && content[1] == 0x8b;

    let results = if is_gzipped {
        // Decompress gzip content
        use flate2::read::GzDecoder;
        let decoder = GzDecoder::new(content);
        let reader = std::io::BufReader::new(decoder);
        parse_fasta_reader(reader, ParseOptions::FULL)?
    } else {
        // Parse directly from bytes
        let reader = std::io::BufReader::new(content);
        parse_fasta_reader(reader, ParseOptions::FULL)?
    };

    let collection_metadata = SequenceCollectionMetadata::from_sequences(&results, None);

    Ok(SequenceCollection {
        metadata: collection_metadata,
        sequences: results,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::alphabet::AlphabetType;

    #[test]
    fn test_digest_fasta_bytes_basic() {
        let fasta = b">chr1\nACGT\n>chr2\nTGCA\n";
        let collection = digest_fasta_bytes(fasta).expect("Failed to digest");
        assert_eq!(collection.sequences.len(), 2);
        assert_eq!(collection.sequences[0].metadata().name, "chr1");
        assert_eq!(collection.sequences[0].metadata().length, 4);
        assert_eq!(collection.sequences[1].metadata().name, "chr2");
        assert_eq!(collection.sequences[1].metadata().length, 4);
    }

    #[test]
    fn test_digest_fasta_bytes_with_description() {
        let fasta = b">chr1 some description here\nACGT\n";
        let collection = digest_fasta_bytes(fasta).expect("Failed to digest");
        assert_eq!(collection.sequences[0].metadata().name, "chr1");
        assert_eq!(
            collection.sequences[0].metadata().description,
            Some("some description here".to_string())
        );
    }

    #[test]
    fn test_load_fasta_bytes_basic() {
        let fasta = b">chr1\nACGT\n>chr2\nTGCA\n";
        let collection = load_fasta_bytes(fasta).expect("Failed to load");
        assert_eq!(collection.sequences.len(), 2);

        // Verify data is loaded
        assert!(collection.sequences[0].is_loaded());
        assert!(collection.sequences[1].is_loaded());

        // Verify sequence data
        assert_eq!(collection.sequences[0].sequence().unwrap(), b"ACGT");
        assert_eq!(collection.sequences[1].sequence().unwrap(), b"TGCA");
    }

    #[test]
    fn test_digest_fasta_bytes_empty() {
        let fasta = b"";
        let collection = digest_fasta_bytes(fasta).expect("Failed to digest empty");
        assert_eq!(collection.sequences.len(), 0);
    }

    #[test]
    fn test_digest_fasta_bytes_multiline() {
        // Test with wrapped sequence lines
        let fasta = b">chr1\nACGT\nTGCA\nAAAA\n";
        let collection = digest_fasta_bytes(fasta).expect("Failed to digest");
        assert_eq!(collection.sequences.len(), 1);
        assert_eq!(collection.sequences[0].metadata().length, 12);
    }

    #[test]
    fn test_digest_fasta_bytes_uppercase() {
        // Verify lowercase is uppercased for digest computation
        let fasta_lower = b">chr1\nacgt\n";
        let fasta_upper = b">chr1\nACGT\n";
        let lower_collection = digest_fasta_bytes(fasta_lower).expect("Failed to digest lower");
        let upper_collection = digest_fasta_bytes(fasta_upper).expect("Failed to digest upper");

        // Digests should be identical
        assert_eq!(
            lower_collection.sequences[0].metadata().sha512t24u,
            upper_collection.sequences[0].metadata().sha512t24u
        );
        assert_eq!(
            lower_collection.sequences[0].metadata().md5,
            upper_collection.sequences[0].metadata().md5
        );
    }

    #[test]
    fn test_digest_fasta_bytes_known_value() {
        // Test against known digest value from base.fa (chrX: TTGGGGAA)
        let fasta = b">chrX\nTTGGGGAA\n";
        let collection = digest_fasta_bytes(fasta).expect("Failed to digest");
        assert_eq!(
            collection.sequences[0].metadata().sha512t24u,
            "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        );
        assert_eq!(
            collection.sequences[0].metadata().md5,
            "5f63cfaa3ef61f88c9635fb9d18ec945"
        );
        assert_eq!(
            collection.sequences[0].metadata().alphabet,
            AlphabetType::Dna2bit
        );
    }

    #[test]
    fn test_load_fasta_bytes_uppercase() {
        // Verify sequence data is uppercased
        let fasta = b">chr1\nacgt\n";
        let collection = load_fasta_bytes(fasta).expect("Failed to load");
        assert_eq!(collection.sequences[0].sequence().unwrap(), b"ACGT");
    }

    #[test]
    fn test_digest_fasta_bytes_gzipped() {
        // Test with gzip content (pre-compressed base.fa content)
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let fasta = b">chr1\nACGT\n";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(fasta).expect("Failed to compress");
        let compressed = encoder.finish().expect("Failed to finish compression");

        let collection = digest_fasta_bytes(&compressed).expect("Failed to digest gzipped");
        assert_eq!(collection.sequences.len(), 1);
        assert_eq!(collection.sequences[0].metadata().name, "chr1");
        assert_eq!(collection.sequences[0].metadata().length, 4);
    }

    #[test]
    fn test_parse_fasta_header() {
        let (name, desc) = parse_fasta_header("chr1 some description here");
        assert_eq!(name, "chr1");
        assert_eq!(desc, Some("some description here".to_string()));

        let (name, desc) = parse_fasta_header("chr1");
        assert_eq!(name, "chr1");
        assert_eq!(desc, None);

        let (name, desc) = parse_fasta_header("chr1\tdescription with tab");
        assert_eq!(name, "chr1");
        assert_eq!(desc, Some("description with tab".to_string()));
    }
}
