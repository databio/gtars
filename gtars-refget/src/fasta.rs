use anyhow::Result;
use std::io::BufRead;
use std::path::Path;

use super::alphabet::{AlphabetGuesser, AlphabetType};
use super::collection::{FaiMetadata, SeqColDigestLvl1, SequenceCollection, SequenceMetadata, SequenceRecord};

use md5::Md5;
use sha2::{Digest, Sha512};

/// A lightweight record containing only FAI (FASTA index) metadata for a sequence.
/// Used by `compute_fai()` for fast FAI-only computation without digest overhead.
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
    use gtars_core::utils::get_dynamic_reader;

    println!("Processing FASTA file: {}", file_path.as_ref().display());

    // Detect if file is gzipped
    let is_gzipped = file_path.as_ref().extension()
        .and_then(|s| s.to_str())
        .map(|s| s == "gz")
        .unwrap_or(false);

    // For non-gzipped files, compute FAI data
    // For gzipped files, FAI data is not meaningful (offsets are for compressed data)
    let fai_enabled = !is_gzipped;

    let mut byte_position: u64 = 0;

    let file_reader = get_dynamic_reader(file_path.as_ref())?;
    let mut reader = std::io::BufReader::new(file_reader);
    let mut results = Vec::new();
    let mut line = String::new();

    let mut current_id: Option<String> = None;
    let mut current_offset: u64 = 0;
    let mut current_line_bases: Option<u32> = None;
    let mut current_line_bytes: Option<u32> = None;
    let mut sha512_hasher = Sha512::new();
    let mut md5_hasher = Md5::new();
    let mut length = 0;
    let mut alphabet_guesser = AlphabetGuesser::new();

    loop {
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            // EOF - finalize the last sequence if any
            if let Some(id) = current_id.take() {
                let sha512 = base64_url::encode(&sha512_hasher.finalize_reset()[0..24]);
                let md5 = format!("{:x}", md5_hasher.finalize_reset());
                let alphabet = alphabet_guesser.guess();

                let metadata = SequenceMetadata {
                    name: id.to_string(),
                    length,
                    sha512t24u: sha512,
                    md5,
                    alphabet,
                };

                let fai = if fai_enabled {
                    if let (Some(lb), Some(lby)) = (current_line_bases, current_line_bytes) {
                        Some(super::collection::FaiMetadata {
                            offset: current_offset,
                            line_bases: lb,
                            line_bytes: lby,
                        })
                    } else {
                        None
                    }
                } else {
                    None
                };

                results.push(SequenceRecord {
                    metadata,
                    fai,
                    data: None,
                });
            }
            break;
        }

        if line.starts_with('>') {
            // Save previous sequence if any
            if let Some(id) = current_id.take() {
                let sha512 = base64_url::encode(&sha512_hasher.finalize_reset()[0..24]);
                let md5 = format!("{:x}", md5_hasher.finalize_reset());
                let alphabet = alphabet_guesser.guess();

                let metadata = SequenceMetadata {
                    name: id.to_string(),
                    length,
                    sha512t24u: sha512,
                    md5,
                    alphabet,
                };

                let fai = if fai_enabled {
                    if let (Some(lb), Some(lby)) = (current_line_bases, current_line_bytes) {
                        Some(super::collection::FaiMetadata {
                            offset: current_offset,
                            line_bases: lb,
                            line_bytes: lby,
                        })
                    } else {
                        None
                    }
                } else {
                    None
                };

                results.push(SequenceRecord {
                    metadata,
                    fai,
                    data: None,
                });
            }

            // Start new sequence
            let id = line[1..].trim().to_string();
            current_id = Some(id);

            // Track position for FAI
            if fai_enabled {
                byte_position += bytes_read as u64;
                current_offset = byte_position;
            }
            current_line_bases = None;
            current_line_bytes = None;
            sha512_hasher = Sha512::new();
            md5_hasher = Md5::new();
            length = 0;
            alphabet_guesser = AlphabetGuesser::new();
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

            // Update digests and length
            let seq_upper = trimmed.to_ascii_uppercase();
            sha512_hasher.update(seq_upper.as_bytes());
            md5_hasher.update(seq_upper.as_bytes());
            length += trimmed.len();
            alphabet_guesser.update(seq_upper.as_bytes());

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
    use gtars_core::utils::get_dynamic_reader;

    // Detect if file is gzipped
    let is_gzipped = file_path.as_ref().extension()
        .and_then(|s| s.to_str())
        .map(|s| s == "gz")
        .unwrap_or(false);

    // For gzipped files, FAI data is not meaningful (offsets are for compressed data)
    let fai_enabled = !is_gzipped;

    let mut byte_position: u64 = 0;

    let file_reader = get_dynamic_reader(file_path.as_ref())?;
    let mut reader = std::io::BufReader::new(file_reader);
    let mut results = Vec::new();
    let mut line = String::new();

    let mut current_id: Option<String> = None;
    let mut current_offset: u64 = 0;
    let mut current_line_bases: Option<u32> = None;
    let mut current_line_bytes: Option<u32> = None;
    let mut length = 0;

    loop {
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            // EOF - finalize the last sequence if any
            if let Some(id) = current_id.take() {
                let fai = if fai_enabled {
                    if let (Some(lb), Some(lby)) = (current_line_bases, current_line_bytes) {
                        Some(FaiMetadata {
                            offset: current_offset,
                            line_bases: lb,
                            line_bytes: lby,
                        })
                    } else {
                        None
                    }
                } else {
                    None
                };

                results.push(FaiRecord {
                    name: id,
                    length,
                    fai,
                });
            }
            break;
        }

        if line.starts_with('>') {
            // Save previous sequence if any
            if let Some(id) = current_id.take() {
                let fai = if fai_enabled {
                    if let (Some(lb), Some(lby)) = (current_line_bases, current_line_bytes) {
                        Some(FaiMetadata {
                            offset: current_offset,
                            line_bases: lb,
                            line_bytes: lby,
                        })
                    } else {
                        None
                    }
                } else {
                    None
                };

                results.push(FaiRecord {
                    name: id,
                    length,
                    fai,
                });
            }

            // Start new sequence
            let id = line[1..].trim().to_string();
            current_id = Some(id);

            // Track position for FAI
            if fai_enabled {
                byte_position += bytes_read as u64;
                current_offset = byte_position;
            }
            current_line_bases = None;
            current_line_bytes = None;
            length = 0;
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

            // Update length (no need for uppercase or hashing)
            length += trimmed.len();

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

    Ok(results)
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
///     println!("Sequence: {}", seq_record.metadata.name);
///     if let Some(data) = &seq_record.data {
///         println!("  Data: {} bytes", data.len());
///     }
/// }
/// ```
pub fn load_fasta<P: AsRef<Path>>(file_path: P) -> Result<SequenceCollection> {
    use gtars_core::utils::get_dynamic_reader;

    println!("Loading FASTA file with data: {}", file_path.as_ref().display());

    // Detect if file is gzipped
    let is_gzipped = file_path.as_ref().extension()
        .and_then(|s| s.to_str())
        .map(|s| s == "gz")
        .unwrap_or(false);

    // For non-gzipped files, compute FAI data
    // For gzipped files, FAI data is not meaningful (offsets are for compressed data)
    let fai_enabled = !is_gzipped;

    let mut byte_position: u64 = 0;

    let file_reader = get_dynamic_reader(file_path.as_ref())?;
    let mut reader = std::io::BufReader::new(file_reader);
    let mut results = Vec::new();
    let mut line = String::new();

    let mut current_id: Option<String> = None;
    let mut current_offset: u64 = 0;
    let mut current_line_bases: Option<u32> = None;
    let mut current_line_bytes: Option<u32> = None;
    let mut sha512_hasher = Sha512::new();
    let mut md5_hasher = Md5::new();
    let mut length = 0;
    let mut alphabet_guesser = AlphabetGuesser::new();
    let mut sequence_data = Vec::new(); // Accumulate sequence data

    loop {
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            // EOF - finalize the last sequence if any
            if let Some(id) = current_id.take() {
                let sha512 = base64_url::encode(&sha512_hasher.finalize_reset()[0..24]);
                let md5 = format!("{:x}", md5_hasher.finalize_reset());
                let alphabet = alphabet_guesser.guess();

                let metadata = SequenceMetadata {
                    name: id.to_string(),
                    length,
                    sha512t24u: sha512,
                    md5,
                    alphabet,
                };

                let fai = if fai_enabled {
                    if let (Some(lb), Some(lby)) = (current_line_bases, current_line_bytes) {
                        Some(super::collection::FaiMetadata {
                            offset: current_offset,
                            line_bases: lb,
                            line_bytes: lby,
                        })
                    } else {
                        None
                    }
                } else {
                    None
                };

                results.push(SequenceRecord {
                    metadata,
                    fai,
                    data: Some(sequence_data.clone()),
                });
            }
            break;
        }

        if line.starts_with('>') {
            // Save previous sequence if any
            if let Some(id) = current_id.take() {
                let sha512 = base64_url::encode(&sha512_hasher.finalize_reset()[0..24]);
                let md5 = format!("{:x}", md5_hasher.finalize_reset());
                let alphabet = alphabet_guesser.guess();

                let metadata = SequenceMetadata {
                    name: id.to_string(),
                    length,
                    sha512t24u: sha512,
                    md5,
                    alphabet,
                };

                let fai = if fai_enabled {
                    if let (Some(lb), Some(lby)) = (current_line_bases, current_line_bytes) {
                        Some(super::collection::FaiMetadata {
                            offset: current_offset,
                            line_bases: lb,
                            line_bytes: lby,
                        })
                    } else {
                        None
                    }
                } else {
                    None
                };

                results.push(SequenceRecord {
                    metadata,
                    fai,
                    data: Some(sequence_data.clone()),
                });
            }

            // Start new sequence
            let id = line[1..].trim().to_string();
            current_id = Some(id);

            // Track position for FAI
            if fai_enabled {
                byte_position += bytes_read as u64;
                current_offset = byte_position;
            }
            current_line_bases = None;
            current_line_bytes = None;
            sha512_hasher = Sha512::new();
            md5_hasher = Md5::new();
            length = 0;
            alphabet_guesser = AlphabetGuesser::new();
            sequence_data = Vec::new(); // Reset for new sequence
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

            // Update digests and length
            let seq_upper = trimmed.to_ascii_uppercase();
            sha512_hasher.update(seq_upper.as_bytes());
            md5_hasher.update(seq_upper.as_bytes());
            length += trimmed.len();
            alphabet_guesser.update(seq_upper.as_bytes());

            // Accumulate sequence data
            sequence_data.extend_from_slice(seq_upper.as_bytes());

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
            if let Some(stripped) = line.strip_prefix("##") {
                if let Some((key, value)) = stripped.split_once('=') {
                    match key {
                        "seqcol_digest" => seqcol_digest = value.to_string(),
                        "names_digest" => names_digest = value.to_string(),
                        "sequences_digest" => sequences_digest = value.to_string(),
                        "lengths_digest" => lengths_digest = value.to_string(),
                        _ => {} // Ignore unknown metadata keys
                    }
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
            fai: None,  // FARG files don't store FAI data
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

        // Test FAI metadata
        assert!(results[0].fai.is_some());
        let fai0 = results[0].fai.as_ref().unwrap();
        assert_eq!(fai0.offset, 6);  // After ">chrX\n"
        assert_eq!(fai0.line_bases, 8);
        assert_eq!(fai0.line_bytes, 9);  // Includes \n

        assert!(results[1].fai.is_some());
        let fai1 = results[1].fai.as_ref().unwrap();
        assert_eq!(fai1.offset, 21);  // After ">chr1\n"
        assert_eq!(fai1.line_bases, 4);
        assert_eq!(fai1.line_bytes, 5);

        assert!(results[2].fai.is_some());
        let fai2 = results[2].fai.as_ref().unwrap();
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
        seqcol.write_farg().expect("Failed to write farg file");

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
            assert_eq!(seq.metadata.name, *name);
            assert_eq!(seq.metadata.length, *length as usize);
            assert!(seq.fai.is_some(), "FAI data should be present");
            let fai = seq.fai.as_ref().unwrap();
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
            assert_eq!(seq.metadata.name, *name);
            assert_eq!(seq.metadata.length, *length as usize);
            assert!(seq.fai.is_some(), "FAI data should be present");
            let fai = seq.fai.as_ref().unwrap();
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
            assert_eq!(seq.metadata.name, *name);
            assert_eq!(seq.metadata.length, *length as usize);
            assert!(seq.fai.is_some(), "FAI data should be present");
            let fai = seq.fai.as_ref().unwrap();
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
            assert_eq!(seq.metadata.name, *name);
            assert_eq!(seq.metadata.length, *length as usize);
            assert!(seq.fai.is_some(), "FAI data should be present");
            let fai = seq.fai.as_ref().unwrap();
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
            assert!(seq.fai.is_none(),
                "Gzipped files should not have FAI data for sequence {}",
                seq.metadata.name);
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
            assert_eq!(seq.metadata.name, *name);
            assert_eq!(seq.metadata.length, *length as usize);
            assert!(seq.fai.is_some(), "FAI data should be present");
            let fai = seq.fai.as_ref().unwrap();
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
            assert_eq!(seq40.metadata.name, seq20.metadata.name);
            assert_eq!(seq40.metadata.name, seq_unwrap.metadata.name);

            // Same length
            assert_eq!(seq40.metadata.length, seq20.metadata.length,
                "Length mismatch for {}", seq40.metadata.name);
            assert_eq!(seq40.metadata.length, seq_unwrap.metadata.length,
                "Length mismatch for {}", seq40.metadata.name);

            // Same digests (proves sequences are identical)
            assert_eq!(seq40.metadata.sha512t24u, seq20.metadata.sha512t24u,
                "SHA512 digest mismatch for {} - sequences not identical",
                seq40.metadata.name);
            assert_eq!(seq40.metadata.sha512t24u, seq_unwrap.metadata.sha512t24u,
                "SHA512 digest mismatch for {} - sequences not identical",
                seq40.metadata.name);
            assert_eq!(seq40.metadata.md5, seq20.metadata.md5,
                "MD5 digest mismatch for {} - sequences not identical",
                seq40.metadata.name);
            assert_eq!(seq40.metadata.md5, seq_unwrap.metadata.md5,
                "MD5 digest mismatch for {} - sequences not identical",
                seq40.metadata.name);

            // But different FAI data (different wrapping)
            assert!(seq40.fai.is_some() && seq20.fai.is_some() && seq_unwrap.fai.is_some());
            let fai40 = seq40.fai.as_ref().unwrap();
            let fai20 = seq20.fai.as_ref().unwrap();
            let fai_unwrap = seq_unwrap.fai.as_ref().unwrap();

            // Different line_bases due to different wrapping (for multi-line sequences)
            // Note: short sequences like chr3 that fit on one line may have same line_bases
            // across different wrapping styles, which is correct behavior
            if seq40.metadata.length > 40 {
                // Long enough to require wrapping in wrapped_40
                assert_ne!(fai40.line_bases, fai20.line_bases,
                    "wrapped_40 and wrapped_20 should have different line_bases for {}",
                    seq40.metadata.name);
            }
            if seq40.metadata.length > 20 {
                // Long enough to require wrapping in wrapped_20
                assert_ne!(fai20.line_bases, fai_unwrap.line_bases,
                    "wrapped_20 and unwrapped should have different line_bases for {}",
                    seq40.metadata.name);
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
        assert_eq!(seq0.metadata.name, "chrX");
        assert_eq!(seq0.metadata.length, 8);
        assert!(seq0.data.is_some(), "Sequence should have data");

        if let Some(ref data) = seq0.data {
            assert_eq!(data.len(), 8);
            assert_eq!(std::str::from_utf8(data).unwrap(), "TTGGGGAA");
        }

        // Check second sequence (chr1)
        let seq1 = &seqcol.sequences[1];
        assert_eq!(seq1.metadata.name, "chr1");
        assert_eq!(seq1.metadata.length, 4);
        assert!(seq1.data.is_some(), "Sequence should have data");

        if let Some(ref data) = seq1.data {
            assert_eq!(data.len(), 4);
            assert_eq!(std::str::from_utf8(data).unwrap(), "GGAA");
        }

        // Check third sequence (chr2)
        let seq2 = &seqcol.sequences[2];
        assert_eq!(seq2.metadata.name, "chr2");
        assert_eq!(seq2.metadata.length, 4);
        assert!(seq2.data.is_some(), "Sequence should have data");

        if let Some(ref data) = seq2.data {
            assert_eq!(data.len(), 4);
            assert_eq!(std::str::from_utf8(data).unwrap(), "GCGC");
        }

        // Verify digests match digest_fasta() output
        let digest_seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to digest FASTA");

        for (loaded, digested) in seqcol.sequences.iter().zip(digest_seqcol.sequences.iter()) {
            assert_eq!(loaded.metadata.name, digested.metadata.name);
            assert_eq!(loaded.metadata.length, digested.metadata.length);
            assert_eq!(loaded.metadata.sha512t24u, digested.metadata.sha512t24u);
            assert_eq!(loaded.metadata.md5, digested.metadata.md5);
            assert_eq!(loaded.metadata.alphabet, digested.metadata.alphabet);
        }
    }

    #[test]
    fn test_load_fasta_wrapped() {
        let seqcol = load_fasta("../tests/data/fasta/wrapped_40.fa")
            .expect("Failed to load wrapped FASTA");

        // Verify all sequences have data
        for seq in &seqcol.sequences {
            assert!(seq.data.is_some(), "Sequence {} should have data", seq.metadata.name);

            if let Some(ref data) = seq.data {
                assert_eq!(data.len(), seq.metadata.length,
                    "Data length mismatch for {}", seq.metadata.name);

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
        assert!(seqcol.sequences[0].data.is_some(), "Should have sequence data");

        for seq in &seqcol.sequences {
            assert!(seq.data.is_some(), "Gzipped file should still have sequence data");
            assert!(seq.fai.is_none(), "Gzipped file should not have FAI metadata");
        }
    }
}
