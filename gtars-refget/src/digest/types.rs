//! Core types for sequence collections - WASM-safe.
//!
//! This module contains the fundamental data structures for representing sequences
//! and sequence collections. All types here are WASM-compatible and don't require
//! filesystem access.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
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

    /// Compute name_length_pairs digest.
    ///
    /// Algorithm: for each sequence, create {"length": L, "name": "N"},
    /// canonicalize each to JSON, digest each, collect into array,
    /// canonicalize array, digest array.
    pub fn compute_name_length_pairs_digest(metadata: &[&SequenceMetadata]) -> String {
        use serde_json::Value;

        // Build array of {"length": N, "name": "X"} pair objects
        let pairs: Vec<Value> = metadata
            .iter()
            .map(|md| {
                let mut obj = serde_json::Map::new();
                obj.insert(
                    "length".to_string(),
                    Value::Number(serde_json::Number::from(md.length)),
                );
                obj.insert("name".to_string(), Value::String(md.name.clone()));
                Value::Object(obj)
            })
            .collect();

        // Canonicalize the entire array of pair objects, then digest
        let canonical = canonicalize_json(&Value::Array(pairs));
        sha512t24u(canonical.as_bytes())
    }

    /// Compute sorted_name_length_pairs digest (order-invariant coordinate system).
    ///
    /// Algorithm: same as name_length_pairs but sort the individual pair digests
    /// lexicographically before digesting the array.
    pub fn compute_sorted_name_length_pairs_digest(metadata: &[&SequenceMetadata]) -> String {
        use serde_json::Value;

        let mut pair_digests: Vec<String> = metadata
            .iter()
            .map(|md| {
                let mut obj = serde_json::Map::new();
                obj.insert(
                    "length".to_string(),
                    Value::Number(serde_json::Number::from(md.length)),
                );
                obj.insert("name".to_string(), Value::String(md.name.clone()));
                let canonical = canonicalize_json(&Value::Object(obj));
                sha512t24u(canonical.as_bytes())
            })
            .collect();

        pair_digests.sort();

        let array_json = Value::Array(
            pair_digests
                .iter()
                .map(|d| Value::String(d.clone()))
                .collect(),
        );
        let canonical = canonicalize_json(&array_json);
        sha512t24u(canonical.as_bytes())
    }

    /// Compute sorted_sequences digest.
    ///
    /// Algorithm: take sequences array (with SQ. prefix), sort lexicographically,
    /// canonicalize, digest.
    pub fn compute_sorted_sequences_digest(metadata: &[&SequenceMetadata]) -> String {
        use serde_json::Value;

        let mut sequences: Vec<String> = metadata
            .iter()
            .map(|md| format!("SQ.{}", md.sha512t24u))
            .collect();

        sequences.sort();

        let array_json = Value::Array(
            sequences
                .iter()
                .map(|s| Value::String(s.clone()))
                .collect(),
        );
        let canonical = canonicalize_json(&array_json);
        sha512t24u(canonical.as_bytes())
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
    /// Ancillary: digest of name_length_pairs array
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name_length_pairs_digest: Option<String>,
    /// Ancillary: digest of sorted_name_length_pairs array (order-invariant coordinate system)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sorted_name_length_pairs_digest: Option<String>,
    /// Ancillary: digest of sorted sequences array
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sorted_sequences_digest: Option<String>,
    /// Optional path to the source file
    pub file_path: Option<PathBuf>,
}

impl SequenceCollectionMetadata {
    /// Compute metadata from sequence records (core digests only).
    pub fn from_sequences(
        sequences: &[SequenceRecord],
        file_path: Option<PathBuf>,
    ) -> Self {
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
            name_length_pairs_digest: None,
            sorted_name_length_pairs_digest: None,
            sorted_sequences_digest: None,
            file_path,
        }
    }

    /// Compute ancillary digests (name_length_pairs, sorted_name_length_pairs,
    /// sorted_sequences) from sequence records. No-op if already computed.
    pub fn compute_ancillary_digests(&mut self, sequences: &[SequenceRecord]) {
        if self.name_length_pairs_digest.is_some() {
            return;
        }
        let metadata_refs: Vec<&SequenceMetadata> =
            sequences.iter().map(|r| r.metadata()).collect();
        self.name_length_pairs_digest =
            Some(SeqColDigestLvl1::compute_name_length_pairs_digest(&metadata_refs));
        self.sorted_name_length_pairs_digest =
            Some(SeqColDigestLvl1::compute_sorted_name_length_pairs_digest(&metadata_refs));
        self.sorted_sequences_digest =
            Some(SeqColDigestLvl1::compute_sorted_sequences_digest(&metadata_refs));
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

    /// Return level 1 representation (attribute digests with spec-compliant field names).
    pub fn to_level1(&self) -> CollectionLevel1 {
        CollectionLevel1 {
            names: self.names_digest.clone(),
            lengths: self.lengths_digest.clone(),
            sequences: self.sequences_digest.clone(),
            name_length_pairs: self.name_length_pairs_digest.clone(),
            sorted_name_length_pairs: self.sorted_name_length_pairs_digest.clone(),
            sorted_sequences: self.sorted_sequences_digest.clone(),
        }
    }
}

/// Level 1 representation: attribute digests with spec-compliant JSON field names.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CollectionLevel1 {
    pub names: String,
    pub lengths: String,
    pub sequences: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub name_length_pairs: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sorted_name_length_pairs: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sorted_sequences: Option<String>,
}

/// Level 2 representation: full arrays with spec-compliant JSON field names.
/// Sequences include SQ. prefix per spec.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CollectionLevel2 {
    pub names: Vec<String>,
    pub lengths: Vec<usize>,
    /// Sequence digests with SQ. prefix per spec
    pub sequences: Vec<String>,
}

/// Result of comparing two sequence collections.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeqColComparison {
    pub digests: ComparisonDigests,
    pub attributes: AttributeComparison,
    pub array_elements: ArrayElementComparison,
}

/// The digests of the two compared collections.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComparisonDigests {
    pub a: String,
    pub b: String,
}

/// Which attributes (array names) are in A only, B only, or both.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AttributeComparison {
    pub a_only: Vec<String>,
    pub b_only: Vec<String>,
    pub a_and_b: Vec<String>,
}

/// Element-level comparison for each shared attribute.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArrayElementComparison {
    pub a_count: HashMap<String, usize>,
    pub b_count: HashMap<String, usize>,
    pub a_and_b_count: HashMap<String, usize>,
    pub a_and_b_same_order: HashMap<String, Option<bool>>,
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
        // Compute metadata from the sequence records (with ancillary digests)
        let metadata = SequenceCollectionMetadata::from_sequences(&records, None);

        SequenceCollection {
            metadata,
            sequences: records,
        }
    }

    /// Return level 2 representation (full arrays, spec format).
    /// Transposes Vec<SequenceRecord> into parallel arrays.
    pub fn to_level2(&self) -> CollectionLevel2 {
        let names: Vec<String> = self
            .sequences
            .iter()
            .map(|r| r.metadata().name.clone())
            .collect();
        let lengths: Vec<usize> = self.sequences.iter().map(|r| r.metadata().length).collect();
        let sequences: Vec<String> = self
            .sequences
            .iter()
            .map(|r| format!("SQ.{}", r.metadata().sha512t24u))
            .collect();

        CollectionLevel2 {
            names,
            lengths,
            sequences,
        }
    }

    /// Compare this collection with another, following the seqcol spec comparison algorithm.
    ///
    /// Dynamically includes ancillary attributes (name_length_pairs, sorted_name_length_pairs,
    /// sorted_sequences) when present in each collection's metadata.
    pub fn compare(&self, other: &SequenceCollection) -> SeqColComparison {
        // Build string-array maps for each collection, including ancillary when present
        let arrays_a = self.to_comparison_arrays();
        let arrays_b = other.to_comparison_arrays();

        let a_keys: std::collections::BTreeSet<&str> = arrays_a.keys().map(|s| s.as_str()).collect();
        let b_keys: std::collections::BTreeSet<&str> = arrays_b.keys().map(|s| s.as_str()).collect();

        let mut a_only = Vec::new();
        let mut b_only = Vec::new();
        let mut a_and_b = Vec::new();

        let mut all_keys: Vec<&str> = a_keys.union(&b_keys).copied().collect();
        all_keys.sort();

        for key in all_keys {
            let in_a = a_keys.contains(key);
            let in_b = b_keys.contains(key);
            match (in_a, in_b) {
                (true, true) => a_and_b.push(key.to_string()),
                (true, false) => a_only.push(key.to_string()),
                (false, true) => b_only.push(key.to_string()),
                (false, false) => unreachable!(),
            }
        }

        let mut a_count = HashMap::new();
        let mut b_count = HashMap::new();
        let mut a_and_b_count = HashMap::new();
        let mut a_and_b_same_order = HashMap::new();

        // Counts for all attributes present in each collection
        for (k, v) in &arrays_a {
            a_count.insert(k.clone(), v.len());
        }
        for (k, v) in &arrays_b {
            b_count.insert(k.clone(), v.len());
        }

        // Element comparison only for shared attributes
        for attr in &a_and_b {
            let arr_a = arrays_a.get(attr).unwrap();
            let arr_b = arrays_b.get(attr).unwrap();
            let (overlap, same_order) = compare_elements(arr_a, arr_b);
            a_and_b_count.insert(attr.clone(), overlap);
            a_and_b_same_order.insert(attr.clone(), same_order);
        }

        SeqColComparison {
            digests: ComparisonDigests {
                a: self.metadata.digest.clone(),
                b: other.metadata.digest.clone(),
            },
            attributes: AttributeComparison {
                a_only,
                b_only,
                a_and_b,
            },
            array_elements: ArrayElementComparison {
                a_count,
                b_count,
                a_and_b_count,
                a_and_b_same_order,
            },
        }
    }

    /// Build string arrays for comparison, including ancillary attributes when present.
    fn to_comparison_arrays(&self) -> HashMap<String, Vec<String>> {
        let mut map = HashMap::new();

        // Core 3 are always present
        map.insert(
            "names".to_string(),
            self.sequences.iter().map(|r| r.metadata().name.clone()).collect(),
        );
        map.insert(
            "lengths".to_string(),
            self.sequences.iter().map(|r| r.metadata().length.to_string()).collect(),
        );
        map.insert(
            "sequences".to_string(),
            self.sequences.iter().map(|r| format!("SQ.{}", r.metadata().sha512t24u)).collect(),
        );

        // Ancillary: only include if this collection has them computed
        if self.metadata.sorted_sequences_digest.is_some() {
            // sorted_sequences: sort the SQ.-prefixed digests
            let mut sorted_seqs: Vec<String> = self.sequences
                .iter()
                .map(|r| format!("SQ.{}", r.metadata().sha512t24u))
                .collect();
            sorted_seqs.sort();
            map.insert("sorted_sequences".to_string(), sorted_seqs);
        }

        if self.metadata.name_length_pairs_digest.is_some() {
            // name_length_pairs: canonical JSON of each {length, name} pair
            let nlp: Vec<String> = self.sequences
                .iter()
                .map(|r| {
                    let md = r.metadata();
                    let mut obj = serde_json::Map::new();
                    obj.insert("length".to_string(), serde_json::Value::Number(serde_json::Number::from(md.length)));
                    obj.insert("name".to_string(), serde_json::Value::String(md.name.clone()));
                    canonicalize_json(&serde_json::Value::Object(obj))
                })
                .collect();
            map.insert("name_length_pairs".to_string(), nlp);
        }

        if self.metadata.sorted_name_length_pairs_digest.is_some() {
            // sorted_name_length_pairs: digest each pair, sort the digests
            let mut snlp: Vec<String> = self.sequences
                .iter()
                .map(|r| {
                    let md = r.metadata();
                    let mut obj = serde_json::Map::new();
                    obj.insert("length".to_string(), serde_json::Value::Number(serde_json::Number::from(md.length)));
                    obj.insert("name".to_string(), serde_json::Value::String(md.name.clone()));
                    sha512t24u(canonicalize_json(&serde_json::Value::Object(obj)).as_bytes())
                })
                .collect();
            snlp.sort();
            map.insert("sorted_name_length_pairs".to_string(), snlp);
        }

        map
    }
}

/// Compare two arrays element-wise following the seqcol spec.
/// Returns (overlap_count, same_order).
///
/// Port of Python `_compare_elements()`:
/// - Filter A to elements present in B, filter B to elements present in A
/// - overlap = min of filtered lengths
/// - same_order: None if fewer than 2 overlapping elements or unbalanced duplicates,
///   otherwise filtered_a == filtered_b
fn compare_elements(a: &[String], b: &[String]) -> (usize, Option<bool>) {
    use std::collections::HashSet;

    let set_a: HashSet<&str> = a.iter().map(|s| s.as_str()).collect();
    let set_b: HashSet<&str> = b.iter().map(|s| s.as_str()).collect();

    // Filter each to elements present in the other
    let filtered_a: Vec<&str> = a.iter().filter(|x| set_b.contains(x.as_str())).map(|s| s.as_str()).collect();
    let filtered_b: Vec<&str> = b.iter().filter(|x| set_a.contains(x.as_str())).map(|s| s.as_str()).collect();

    let overlap = filtered_a.len().min(filtered_b.len());

    let same_order = if overlap < 2 {
        None
    } else if filtered_a.len() != filtered_b.len() || filtered_a.len() != overlap {
        // Unbalanced duplicates
        None
    } else {
        Some(filtered_a == filtered_b)
    };

    (overlap, same_order)
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

/// Parse a single line from an RGCI (collection index) file.
///
/// RGCI format is tab-separated with 5+ columns:
/// digest, n_sequences, names_digest, sequences_digest, lengths_digest,
/// [name_length_pairs_digest, sorted_name_length_pairs_digest, sorted_sequences_digest]
///
/// Lines starting with '#' are treated as comments and return None.
/// Lines with fewer than 5 columns return None.
/// Columns 5-7 are optional ancillary digests (empty string = None).
pub fn parse_rgci_line(line: &str) -> Option<SequenceCollectionMetadata> {
    if line.starts_with('#') {
        return None;
    }
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 5 {
        return None;
    }
    // Parse optional ancillary digest columns (empty string -> None)
    let opt_col = |i: usize| -> Option<String> {
        parts.get(i).and_then(|s| {
            if s.is_empty() { None } else { Some(s.to_string()) }
        })
    };
    Some(SequenceCollectionMetadata {
        digest: parts[0].to_string(),
        n_sequences: parts[1].parse().ok()?,
        names_digest: parts[2].to_string(),
        sequences_digest: parts[3].to_string(),
        lengths_digest: parts[4].to_string(),
        name_length_pairs_digest: opt_col(5),
        sorted_name_length_pairs_digest: opt_col(6),
        sorted_sequences_digest: opt_col(7),
        file_path: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_metadata() -> Vec<SequenceMetadata> {
        vec![
            SequenceMetadata {
                name: "chrX".to_string(),
                description: None,
                length: 8,
                sha512t24u: "abc123".to_string(),
                md5: "md5abc".to_string(),
                alphabet: AlphabetType::Dna2bit,
                fai: None,
            },
            SequenceMetadata {
                name: "chr1".to_string(),
                description: None,
                length: 4,
                sha512t24u: "def456".to_string(),
                md5: "md5def".to_string(),
                alphabet: AlphabetType::Dna2bit,
                fai: None,
            },
        ]
    }

    #[test]
    fn test_ancillary_digest_nlp() {
        let metadata = test_metadata();
        let refs: Vec<&SequenceMetadata> = metadata.iter().collect();
        let nlp = SeqColDigestLvl1::compute_name_length_pairs_digest(&refs);
        assert!(!nlp.is_empty());
        assert_eq!(nlp.len(), 32); // SHA-512/24u is 32 chars base64url
    }

    #[test]
    fn test_ancillary_digest_snlp() {
        let metadata = test_metadata();
        let refs: Vec<&SequenceMetadata> = metadata.iter().collect();
        let snlp = SeqColDigestLvl1::compute_sorted_name_length_pairs_digest(&refs);
        assert!(!snlp.is_empty());
        assert_eq!(snlp.len(), 32);
    }

    #[test]
    fn test_ancillary_digest_sorted_sequences() {
        let metadata = test_metadata();
        let refs: Vec<&SequenceMetadata> = metadata.iter().collect();
        let ss = SeqColDigestLvl1::compute_sorted_sequences_digest(&refs);
        assert!(!ss.is_empty());
        assert_eq!(ss.len(), 32);
    }

    #[test]
    fn test_nlp_and_snlp_both_valid() {
        let metadata = test_metadata();
        let refs: Vec<&SequenceMetadata> = metadata.iter().collect();
        let nlp = SeqColDigestLvl1::compute_name_length_pairs_digest(&refs);
        let snlp = SeqColDigestLvl1::compute_sorted_name_length_pairs_digest(&refs);
        // Both should be valid 32-char digests
        assert_eq!(nlp.len(), 32);
        assert_eq!(snlp.len(), 32);
        // NLP and SNLP may or may not be equal depending on whether pair digests
        // are already sorted. The key property is that SNLP is order-invariant
        // (tested in test_snlp_order_invariant).
    }

    #[test]
    fn test_snlp_order_invariant() {
        let metadata = test_metadata();
        // Reverse the order
        let reversed: Vec<SequenceMetadata> = metadata.iter().rev().cloned().collect();

        let refs1: Vec<&SequenceMetadata> = metadata.iter().collect();
        let refs2: Vec<&SequenceMetadata> = reversed.iter().collect();

        let snlp1 = SeqColDigestLvl1::compute_sorted_name_length_pairs_digest(&refs1);
        let snlp2 = SeqColDigestLvl1::compute_sorted_name_length_pairs_digest(&refs2);

        // Sorted NLP should be the same regardless of input order
        assert_eq!(snlp1, snlp2);
    }

    #[test]
    fn test_sorted_sequences_order_invariant() {
        let metadata = test_metadata();
        let reversed: Vec<SequenceMetadata> = metadata.iter().rev().cloned().collect();

        let refs1: Vec<&SequenceMetadata> = metadata.iter().collect();
        let refs2: Vec<&SequenceMetadata> = reversed.iter().collect();

        let ss1 = SeqColDigestLvl1::compute_sorted_sequences_digest(&refs1);
        let ss2 = SeqColDigestLvl1::compute_sorted_sequences_digest(&refs2);

        // Sorted sequences should be the same regardless of input order
        assert_eq!(ss1, ss2);
    }

    #[test]
    fn test_nlp_order_sensitive() {
        let metadata = test_metadata();
        let reversed: Vec<SequenceMetadata> = metadata.iter().rev().cloned().collect();

        let refs1: Vec<&SequenceMetadata> = metadata.iter().collect();
        let refs2: Vec<&SequenceMetadata> = reversed.iter().collect();

        let nlp1 = SeqColDigestLvl1::compute_name_length_pairs_digest(&refs1);
        let nlp2 = SeqColDigestLvl1::compute_name_length_pairs_digest(&refs2);

        // Unsorted NLP should differ when order changes
        assert_ne!(nlp1, nlp2);
    }

    #[test]
    fn test_to_level1() {
        let metadata = test_metadata();
        let records: Vec<_> = metadata
            .iter()
            .map(|m| SequenceRecord::Stub(m.clone()))
            .collect();
        let mut coll_meta = SequenceCollectionMetadata::from_sequences(&records, None);
        coll_meta.compute_ancillary_digests(&records);

        let lvl1 = coll_meta.to_level1();
        assert_eq!(lvl1.names, coll_meta.names_digest);
        assert_eq!(lvl1.lengths, coll_meta.lengths_digest);
        assert_eq!(lvl1.sequences, coll_meta.sequences_digest);
        assert!(lvl1.name_length_pairs.is_some());
        assert!(lvl1.sorted_name_length_pairs.is_some());
        assert!(lvl1.sorted_sequences.is_some());
    }

    #[test]
    fn test_to_level2() {
        let metadata = test_metadata();
        let records: Vec<SequenceRecord> = metadata
            .iter()
            .map(|m| SequenceRecord::Stub(m.clone()))
            .collect();
        let collection = SequenceCollection::from_records(records);

        let lvl2 = collection.to_level2();
        assert_eq!(lvl2.names, vec!["chrX", "chr1"]);
        assert_eq!(lvl2.lengths, vec![8, 4]);
        assert_eq!(lvl2.sequences.len(), 2);
        assert!(lvl2.sequences[0].starts_with("SQ."));
    }

    #[test]
    fn test_compare_same() {
        let records: Vec<SequenceRecord> = test_metadata().into_iter().map(SequenceRecord::Stub).collect();
        let collection = SequenceCollection::from_records(records);
        let result = collection.compare(&collection);

        assert_eq!(result.digests.a, result.digests.b);
        assert_eq!(result.attributes.a_and_b.len(), 3); // core only, no ancillary
        for attr in &result.attributes.a_and_b {
            assert_eq!(result.array_elements.a_and_b_same_order[attr], Some(true));
        }
    }

    #[test]
    fn test_compare_reversed_order() {
        let metadata = test_metadata();
        let reversed: Vec<SequenceMetadata> = metadata.iter().rev().cloned().collect();
        let coll_a = SequenceCollection::from_records(metadata.into_iter().map(SequenceRecord::Stub).collect());
        let coll_b = SequenceCollection::from_records(reversed.into_iter().map(SequenceRecord::Stub).collect());

        let result = coll_a.compare(&coll_b);
        for attr in &result.attributes.a_and_b {
            assert_eq!(result.array_elements.a_and_b_count[attr], 2);
            assert_eq!(result.array_elements.a_and_b_same_order[attr], Some(false));
        }
    }

    #[test]
    fn test_compare_single_element() {
        let meta = SequenceMetadata {
            name: "chr1".to_string(),
            description: None,
            length: 4,
            sha512t24u: "abc".to_string(),
            md5: "md5".to_string(),
            alphabet: AlphabetType::Dna2bit,
            fai: None,
        };
        let coll = SequenceCollection::from_records(vec![SequenceRecord::Stub(meta)]);
        let result = coll.compare(&coll);

        // With only 1 element overlap, same_order should be None
        for attr in &result.attributes.a_and_b {
            assert_eq!(result.array_elements.a_and_b_same_order[attr], None);
        }
    }
}
