//! Core types for sequence collections - WASM-safe.
//!
//! This module contains the fundamental data structures for representing sequences
//! and sequence collections. All types here are WASM-compatible and don't require
//! filesystem access.

use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::path::PathBuf;

use super::algorithms::{canonicalize_json, md5, sha512t24u};
use super::alphabet::{AlphabetType, guess_alphabet};

/// Metadata for a single sequence, including its name, length, digests, and alphabet type.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SequenceMetadata {
    pub name: String,
    /// Description from FASTA header (text after first whitespace).
    #[serde(default)]
    pub description: Option<String>,
    pub length: usize,
    pub sha512t24u: String,
    pub md5: String,
    pub alphabet: AlphabetType,
    pub fai: Option<FaiMetadata>,
}

impl Default for SequenceMetadata {
    fn default() -> Self {
        Self {
            name: String::new(),
            description: None,
            length: 0,
            sha512t24u: String::new(),
            md5: String::new(),
            alphabet: AlphabetType::Ascii,
            fai: None,
        }
    }
}

/// FASTA index (FAI) metadata for a sequence.
/// This data is only present when a sequence was loaded from a FASTA file.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FaiMetadata {
    pub offset: u64,     // byte offset to first base of sequence data
    pub line_bases: u32, // number of bases per line
    pub line_bytes: u32, // number of bytes per line (including newline chars)
}

/// A representation of a single sequence that includes metadata and optionally data.
/// Combines sequence metadata with optional raw/encoded data.
///
/// This enum has two variants:
/// - `Stub`: Contains only metadata, no sequence data loaded
/// - `Full`: Contains both metadata and the actual sequence data
#[derive(Clone, Debug)]
pub enum SequenceRecord {
    /// A sequence record with only metadata, no sequence data
    Stub(SequenceMetadata),
    /// A sequence record with both metadata and sequence data
    Full {
        metadata: SequenceMetadata,
        sequence: Vec<u8>,
    },
}

impl SequenceRecord {
    /// Get metadata regardless of variant
    pub fn metadata(&self) -> &SequenceMetadata {
        match self {
            SequenceRecord::Stub(meta) => meta,
            SequenceRecord::Full { metadata, .. } => metadata,
        }
    }

    /// Get sequence data if present
    pub fn sequence(&self) -> Option<&[u8]> {
        match self {
            SequenceRecord::Stub(_) => None,
            SequenceRecord::Full { sequence, .. } => Some(sequence),
        }
    }

    /// Check if sequence data is loaded (Full) or just metadata (Stub).
    pub fn is_loaded(&self) -> bool {
        matches!(self, SequenceRecord::Full { .. })
    }

    /// Load data into a Stub record, or replace data in a Full record (takes ownership)
    pub fn with_data(self, sequence: Vec<u8>) -> Self {
        let metadata = match self {
            SequenceRecord::Stub(m) => m,
            SequenceRecord::Full { metadata, .. } => metadata,
        };
        SequenceRecord::Full { metadata, sequence }
    }

    /// Load data into a Stub record in-place, converting it to Full.
    /// If already Full, replaces the existing sequence data.
    ///
    /// This is more efficient than `with_data()` when you have a mutable reference,
    /// as it avoids cloning the metadata.
    pub fn load_data(&mut self, sequence: Vec<u8>) {
        match self {
            SequenceRecord::Stub(metadata) => {
                // Take ownership of metadata without cloning
                let metadata = std::mem::take(metadata);
                *self = SequenceRecord::Full { metadata, sequence };
            }
            SequenceRecord::Full {
                sequence: existing, ..
            } => {
                // Just replace the sequence data
                *existing = sequence;
            }
        }
    }

    /// Decodes the sequence data to a string.
    ///
    /// This method attempts to decode the sequence data stored in this record.
    /// It handles both raw (uncompressed UTF-8) and encoded (bit-packed) data.
    /// The decoding strategy depends on the alphabet type:
    /// - For ASCII alphabet: data is already in raw form, just convert to string
    /// - For other alphabets: attempt encoded decoding first, fall back to raw
    ///
    /// # Returns
    ///
    /// * `Some(String)` - The decoded sequence if data is loaded
    /// * `None` - If no data is loaded in this record
    pub fn decode(&self) -> Option<String> {
        use super::alphabet::lookup_alphabet;
        use super::encoder::decode_substring_from_bytes;

        let (metadata, data) = match self {
            SequenceRecord::Stub(_) => return None,
            SequenceRecord::Full { metadata, sequence } => (metadata, sequence),
        };

        // For ASCII alphabet (8 bits per symbol), the data is always stored raw
        if metadata.alphabet == AlphabetType::Ascii {
            return String::from_utf8(data.clone()).ok();
        }

        // Try to detect if data is raw or encoded
        // Heuristic: for encoded data, the size should be approximately length * bits_per_symbol / 8
        // For raw data, the size should be approximately equal to length
        let alphabet = lookup_alphabet(&metadata.alphabet);

        // If data size matches the expected length (not the encoded size), it's probably raw
        if data.len() == metadata.length {
            // Try to decode as UTF-8
            if let Ok(raw_string) = String::from_utf8(data.clone()) {
                // Data appears to be raw UTF-8
                return Some(raw_string);
            }
        }

        // Data is probably encoded (size matches expected encoded size), try to decode it
        let decoded_bytes = decode_substring_from_bytes(data, 0, metadata.length, alphabet);

        // Convert to string
        String::from_utf8(decoded_bytes).ok()
    }
}

impl Display for SequenceRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SequenceRecord: {} (length: {}, alphabet: {}, ga4gh: {:02x?}, md5: {:02x?})",
            &self.metadata().name,
            &self.metadata().length,
            &self.metadata().alphabet,
            &self.metadata().sha512t24u,
            &self.metadata().md5
        )?;
        Ok(())
    }
}

/// A struct representing the first level of digests for a refget sequence collection.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SeqColDigestLvl1 {
    pub sequences_digest: String,
    pub names_digest: String,
    pub lengths_digest: String,
}

impl SeqColDigestLvl1 {
    /// Compute collection digest from lvl1 digests
    pub fn to_digest(&self) -> String {
        // Create JSON object with the lvl1 digest strings
        let mut lvl1_object = serde_json::Map::new();
        lvl1_object.insert(
            "names".to_string(),
            serde_json::Value::String(self.names_digest.clone()),
        );
        lvl1_object.insert(
            "sequences".to_string(),
            serde_json::Value::String(self.sequences_digest.clone()),
        );

        let lvl1_json = serde_json::Value::Object(lvl1_object);

        // Canonicalize the JSON object and compute collection digest
        let lvl1_canonical = canonicalize_json(&lvl1_json);
        sha512t24u(lvl1_canonical.as_bytes())
    }

    /// Compute lvl1 digests from a collection of SequenceMetadata
    pub fn from_metadata(metadata_vec: &[&SequenceMetadata]) -> Self {
        use serde_json::Value;

        // Extract arrays for each field
        let sequences: Vec<String> = metadata_vec
            .iter()
            .map(|md| format!("SQ.{}", md.sha512t24u))
            .collect();
        let names: Vec<&str> = metadata_vec.iter().map(|md| md.name.as_str()).collect();
        let lengths: Vec<usize> = metadata_vec.iter().map(|md| md.length).collect();

        // Convert to JSON Values and canonicalize
        let sequences_json = Value::Array(
            sequences
                .iter()
                .map(|s| Value::String(s.to_string()))
                .collect(),
        );
        let names_json = Value::Array(names.iter().map(|s| Value::String(s.to_string())).collect());
        let lengths_json = Value::Array(
            lengths
                .iter()
                .map(|l| Value::Number(serde_json::Number::from(*l)))
                .collect(),
        );

        // Canonicalize to JCS format
        let sequences_canonical = canonicalize_json(&sequences_json);
        let names_canonical = canonicalize_json(&names_json);
        let lengths_canonical = canonicalize_json(&lengths_json);

        // Hash the canonicalized arrays
        SeqColDigestLvl1 {
            sequences_digest: sha512t24u(sequences_canonical.as_bytes()),
            names_digest: sha512t24u(names_canonical.as_bytes()),
            lengths_digest: sha512t24u(lengths_canonical.as_bytes()),
        }
    }
}

/// Metadata for a sequence collection (parallel to SequenceMetadata).
/// Contains the collection digest and level 1 digests for names, sequences, and lengths.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SequenceCollectionMetadata {
    /// Top-level seqcol digest
    pub digest: String,
    /// Number of sequences in the collection
    pub n_sequences: usize,
    /// Level 1 digest of names array
    pub names_digest: String,
    /// Level 1 digest of sequences array
    pub sequences_digest: String,
    /// Level 1 digest of lengths array
    pub lengths_digest: String,
    /// Optional path to the source file
    pub file_path: Option<PathBuf>,
}

impl SequenceCollectionMetadata {
    /// Compute metadata from sequence records
    pub fn from_sequences(sequences: &[SequenceRecord], file_path: Option<PathBuf>) -> Self {
        // Extract metadata refs
        let metadata_refs: Vec<&SequenceMetadata> =
            sequences.iter().map(|r| r.metadata()).collect();

        // Compute level 1 digests
        let lvl1 = SeqColDigestLvl1::from_metadata(&metadata_refs);

        // Compute top-level digest from level 1 digests
        let digest = lvl1.to_digest();

        Self {
            digest,
            n_sequences: sequences.len(),
            names_digest: lvl1.names_digest,
            sequences_digest: lvl1.sequences_digest,
            lengths_digest: lvl1.lengths_digest,
            file_path,
        }
    }

    /// Create from an existing SequenceCollection
    pub fn from_collection(collection: &SequenceCollection) -> Self {
        collection.metadata.clone()
    }

    /// Convert to SeqColDigestLvl1 for compatibility
    pub fn to_lvl1(&self) -> SeqColDigestLvl1 {
        SeqColDigestLvl1 {
            sequences_digest: self.sequences_digest.clone(),
            names_digest: self.names_digest.clone(),
            lengths_digest: self.lengths_digest.clone(),
        }
    }
}

/// A single Sequence Collection, which may or may not hold data.
#[derive(Clone, Debug)]
pub struct SequenceCollection {
    /// Collection metadata (digest, level 1 digests, n_sequences, file_path)
    pub metadata: SequenceCollectionMetadata,

    /// Vector of SequenceRecords, which contain metadata (name, length, digests, alphabet)
    /// and optionally the actual sequence data.
    pub sequences: Vec<SequenceRecord>,
}

impl SequenceCollection {
    /// Create a SequenceCollection from a vector of SequenceRecords.
    pub fn from_records(records: Vec<SequenceRecord>) -> Self {
        // Compute metadata from the sequence records
        let metadata = SequenceCollectionMetadata::from_sequences(&records, None);

        SequenceCollection {
            metadata,
            sequences: records,
        }
    }
}

impl Display for SequenceCollection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SequenceCollection with {} sequences, digest: {}",
            self.sequences.len(),
            self.metadata.digest
        )?;
        write!(f, "\nFirst 3 sequences:")?;
        for seqrec in self.sequences.iter().take(3) {
            write!(f, "\n- {}", seqrec)?;
        }
        Ok(())
    }
}

// Iterator implementations for SequenceCollection
// Allows: for seq in &collection { ... }
impl<'a> IntoIterator for &'a SequenceCollection {
    type Item = &'a SequenceRecord;
    type IntoIter = std::slice::Iter<'a, SequenceRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.sequences.iter()
    }
}

// Consuming iterator
// Allows: for seq in collection { ... } (consumes the collection)
impl IntoIterator for SequenceCollection {
    type Item = SequenceRecord;
    type IntoIter = std::vec::IntoIter<SequenceRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.sequences.into_iter()
    }
}

/// A collection record that may or may not have its sequence list loaded.
/// Parallel to SequenceRecord.
#[derive(Clone, Debug)]
pub enum SequenceCollectionRecord {
    /// Collection with only metadata, sequence list not loaded
    Stub(SequenceCollectionMetadata),
    /// Collection with metadata and the actual sequence list
    Full {
        metadata: SequenceCollectionMetadata,
        sequences: Vec<SequenceRecord>,
    },
}

impl SequenceCollectionRecord {
    /// Get metadata regardless of variant
    pub fn metadata(&self) -> &SequenceCollectionMetadata {
        match self {
            SequenceCollectionRecord::Stub(meta) => meta,
            SequenceCollectionRecord::Full { metadata, .. } => metadata,
        }
    }

    /// Get sequences if loaded
    pub fn sequences(&self) -> Option<&[SequenceRecord]> {
        match self {
            SequenceCollectionRecord::Stub(_) => None,
            SequenceCollectionRecord::Full { sequences, .. } => Some(sequences),
        }
    }

    /// Check if sequences are loaded
    pub fn has_sequences(&self) -> bool {
        matches!(self, SequenceCollectionRecord::Full { .. })
    }

    /// Load sequences into a Stub record, converting to Full
    pub fn with_sequences(self, sequences: Vec<SequenceRecord>) -> Self {
        let metadata = match self {
            SequenceCollectionRecord::Stub(m) => m,
            SequenceCollectionRecord::Full { metadata, .. } => metadata,
        };
        SequenceCollectionRecord::Full {
            metadata,
            sequences,
        }
    }

    /// Convert to a SequenceCollection (requires Full variant or empty collection for Stub)
    pub fn to_collection(&self) -> SequenceCollection {
        match self {
            SequenceCollectionRecord::Stub(meta) => {
                // Create empty collection with metadata
                SequenceCollection {
                    metadata: meta.clone(),
                    sequences: Vec::new(),
                }
            }
            SequenceCollectionRecord::Full {
                metadata,
                sequences,
            } => SequenceCollection {
                metadata: metadata.clone(),
                sequences: sequences.clone(),
            },
        }
    }
}

impl From<SequenceCollection> for SequenceCollectionRecord {
    fn from(collection: SequenceCollection) -> Self {
        SequenceCollectionRecord::Full {
            metadata: collection.metadata,
            sequences: collection.sequences,
        }
    }
}

// ============================================================================
// Pure computation functions
// ============================================================================

/// Create a SequenceRecord from raw data, computing all metadata.
///
/// This is the sequence-level parallel to `digest_fasta()` for collections.
/// It computes the GA4GH sha512t24u digest, MD5 digest, detects the alphabet,
/// and packages everything into a SequenceRecord with Full variant.
///
/// # Arguments
/// * `name` - The sequence name (e.g., "chr1")
/// * `data` - The raw sequence bytes (e.g., b"ACGTACGT")
///
/// # Returns
/// A SequenceRecord::Full with computed metadata and the original data
///
/// # Example
/// ```
/// use gtars_refget::digest::types::digest_sequence;
///
/// let seq = digest_sequence("chr1", b"ACGTACGT");
/// assert_eq!(seq.metadata().name, "chr1");
/// assert_eq!(seq.metadata().length, 8);
/// assert!(!seq.metadata().sha512t24u.is_empty());
/// ```
pub fn digest_sequence(name: &str, data: &[u8]) -> SequenceRecord {
    // Uppercase the data for consistent digest computation (matches FASTA processing)
    let uppercased: Vec<u8> = data.iter().map(|b| b.to_ascii_uppercase()).collect();

    let metadata = SequenceMetadata {
        name: name.to_string(),
        description: None,
        length: data.len(),
        sha512t24u: sha512t24u(&uppercased),
        md5: md5(&uppercased),
        alphabet: guess_alphabet(&uppercased),
        fai: None, // No FAI data for programmatically created sequences
    };
    SequenceRecord::Full {
        metadata,
        sequence: uppercased,
    }
}

/// Create a SequenceRecord with a description field.
///
/// Same as `digest_sequence()` but allows specifying an optional description.
///
/// # Arguments
/// * `name` - The sequence name (e.g., "chr1")
/// * `description` - Optional description text
/// * `data` - The raw sequence bytes (e.g., b"ACGTACGT")
///
/// # Returns
/// A SequenceRecord::Full with computed metadata and the original data
pub fn digest_sequence_with_description(
    name: &str,
    description: Option<&str>,
    data: &[u8],
) -> SequenceRecord {
    let mut seq = digest_sequence(name, data);
    if let SequenceRecord::Full {
        ref mut metadata, ..
    } = seq
    {
        metadata.description = description.map(String::from);
    }
    seq
}

/// Parse a single RGSI line into SequenceMetadata.
///
/// Supports two formats:
/// - 5-column (no description): `name\tlength\talphabet\tsha512t24u\tmd5`
/// - 6-column (with description): `name\tlength\talphabet\tsha512t24u\tmd5\tdescription`
///
/// Returns None if the line is a comment, empty, or has wrong column count.
pub fn parse_rgsi_line(line: &str) -> Option<SequenceMetadata> {
    // Skip empty lines
    if line.trim().is_empty() {
        return None;
    }

    let parts: Vec<&str> = line.split('\t').collect();

    match parts.len() {
        // 5-column format: no description
        5 => Some(SequenceMetadata {
            name: parts[0].to_string(),
            description: None,
            length: parts[1].parse().ok()?,
            alphabet: parts[2].parse().unwrap_or(AlphabetType::Unknown),
            sha512t24u: parts[3].to_string(),
            md5: parts[4].to_string(),
            fai: None,
        }),
        // 6-column format: description at end
        6 => Some(SequenceMetadata {
            name: parts[0].to_string(),
            description: if parts[5].is_empty() {
                None
            } else {
                Some(parts[5].to_string())
            },
            length: parts[1].parse().ok()?,
            alphabet: parts[2].parse().unwrap_or(AlphabetType::Unknown),
            sha512t24u: parts[3].to_string(),
            md5: parts[4].to_string(),
            fai: None,
        }),
        _ => None,
    }
}
