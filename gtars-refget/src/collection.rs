use crate::alphabet::{guess_alphabet, AlphabetType};
use crate::digest::{canonicalize_json, md5, sha512t24u};
use crate::fasta::digest_fasta;
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::utils::PathExtension;

/// Shared implementation for writing RGSI (Refget Sequence Index) files.
///
/// This helper function writes the RGSI format:
/// - Collection digest headers (##seqcol_digest, ##names_digest, etc.)
/// - Column header line
/// - Per-sequence metadata lines (if sequences are provided)
///
/// # Arguments
/// * `file_path` - Path to the output RGSI file
/// * `metadata` - Collection-level metadata containing digests
/// * `sequences` - Optional slice of sequence records to write
fn write_rgsi_impl<P: AsRef<Path>>(
    file_path: P,
    metadata: &SequenceCollectionMetadata,
    sequences: Option<&[SequenceRecord]>,
) -> Result<()> {
    let file_path = file_path.as_ref();
    let mut file = std::fs::File::create(file_path)?;

    // Write collection digest headers
    writeln!(file, "##seqcol_digest={}", metadata.digest)?;
    writeln!(file, "##names_digest={}", metadata.names_digest)?;
    writeln!(file, "##sequences_digest={}", metadata.sequences_digest)?;
    writeln!(file, "##lengths_digest={}", metadata.lengths_digest)?;
    writeln!(file, "#name\tlength\talphabet\tsha512t24u\tmd5\tdescription")?;

    // Write sequence metadata if available
    if let Some(seqs) = sequences {
        for seq_record in seqs {
            let seq_meta = seq_record.metadata();
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}",
                seq_meta.name,
                seq_meta.length,
                seq_meta.alphabet,
                seq_meta.sha512t24u,
                seq_meta.md5,
                seq_meta.description.as_deref().unwrap_or("")
            )?;
        }
    }
    Ok(())
}

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
/// use gtars_refget::collection::digest_sequence;
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
    if let SequenceRecord::Full { ref mut metadata, .. } = seq {
        metadata.description = description.map(String::from);
    }
    seq
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
        let metadata_refs: Vec<&SequenceMetadata> = sequences.iter()
            .map(|r| r.metadata())
            .collect();

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
        SequenceCollectionRecord::Full { metadata, sequences }
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
            SequenceCollectionRecord::Full { metadata, sequences } => {
                SequenceCollection {
                    metadata: metadata.clone(),
                    sequences: sequences.clone(),
                }
            }
        }
    }

    /// Write the collection to an RGSI file
    pub fn write_collection_rgsi<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        write_rgsi_impl(file_path, self.metadata(), self.sequences())
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

/// A single Sequence Collection, which may or may not hold data.
#[derive(Clone, Debug)]
pub struct SequenceCollection {
    /// Collection metadata (digest, level 1 digests, n_sequences, file_path)
    pub metadata: SequenceCollectionMetadata,

    /// Vector of SequenceRecords, which contain metadata (name, length, digests, alphabet)
    /// and optionally the actual sequence data.
    pub sequences: Vec<SequenceRecord>,
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

use std::fs::{self, File};
/// Metadata for a single sequence, including its name, length, digests, and alphabet type.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SequenceMetadata {
    pub name: String,
    /// Description from FASTA header (text after first whitespace).
    /// Only populated when strict_seqnames=true during FASTA parsing.
    #[serde(default)]
    pub description: Option<String>,
    pub length: usize,
    pub sha512t24u: String,
    pub md5: String,
    pub alphabet: AlphabetType,
    pub fai: Option<FaiMetadata>,
}

/// FASTA index (FAI) metadata for a sequence.
/// This data is only present when a sequence was loaded from a FASTA file.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FaiMetadata {
    pub offset: u64,       // byte offset to first base of sequence data
    pub line_bases: u32,   // number of bases per line
    pub line_bytes: u32,   // number of bytes per line (including newline chars)
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

impl SequenceMetadata {
    /// Calculate the disk size in bytes for this sequence
    ///
    /// # Arguments
    /// * `mode` - The storage mode (Raw or Encoded)
    ///
    /// # Returns
    /// The number of bytes this sequence occupies on disk
    ///
    /// # Examples
    /// ```ignore
    /// // For a 1000bp DNA sequence in Encoded mode with Dna2bit alphabet:
    /// // disk_size = (1000 * 2 bits).div_ceil(8) = 250 bytes
    /// let size = metadata.disk_size(&StorageMode::Encoded);
    /// ```
    pub fn disk_size(&self, mode: &crate::store::StorageMode) -> usize {
        match mode {
            crate::store::StorageMode::Raw => self.length,
            crate::store::StorageMode::Encoded => {
                let bits_per_symbol = self.alphabet.bits_per_symbol();
                let total_bits = self.length * bits_per_symbol;
                total_bits.div_ceil(8)
            }
        }
    }
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
            SequenceRecord::Full { sequence: existing, .. } => {
                // Just replace the sequence data
                *existing = sequence;
            }
        }
    }

    /// Utility function to write a single sequence to a file
    pub fn to_file<P: AsRef<Path>>(&self, path: P) -> anyhow::Result<()> {
        let data = match self {
            SequenceRecord::Stub(_) => {
                return Err(anyhow::anyhow!("Cannot write file: sequence data not loaded"));
            }
            SequenceRecord::Full { sequence, .. } => sequence,
        };

        // Create parent directories if they don't exist
        if let Some(parent) = path.as_ref().parent() {
            fs::create_dir_all(parent)?;
        }

        let mut file = File::create(path)?;
        file.write_all(data)?;
        Ok(())
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
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use gtars_refget::collection::SequenceRecord;
    /// # let record: SequenceRecord = todo!();
    /// let sequence = record.decode();
    /// if let Some(seq) = sequence {
    ///     println!("Sequence: {}", seq);
    /// }
    /// ```
    pub fn decode(&self) -> Option<String> {
        use crate::alphabet::lookup_alphabet;
        use crate::encoder::decode_substring_from_bytes;

        let (metadata, data) = match self {
            SequenceRecord::Stub(_) => return None,
            SequenceRecord::Full { metadata, sequence } => (metadata, sequence),
        };

        // For ASCII alphabet (8 bits per symbol), the data is always stored raw
        if metadata.alphabet == crate::alphabet::AlphabetType::Ascii {
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
        let decoded_bytes = decode_substring_from_bytes(
            data,
            0,
            metadata.length,
            alphabet
        );

        // Convert to string
        String::from_utf8(decoded_bytes).ok()
    }
}

impl SequenceCollection {
    /// Default behavior: read and write cache
    pub fn from_fasta<P: AsRef<Path>>(file_path: P) -> Result<Self> {
        Self::from_path_with_cache(file_path, true, true)
    }

    pub fn from_rgsi<P: AsRef<Path>>(file_path: P) -> Result<Self> {
        let rgsi_file_path = file_path.as_ref().replace_exts_with("rgsi");

        if rgsi_file_path.exists() {
            read_rgsi_file(&rgsi_file_path)
        } else {
            Err(anyhow::anyhow!(
                "RGSI file does not exist at {:?}",
                rgsi_file_path
            ))
        }
    }

    /// Create a SequenceCollection from a vector of SequenceRecords.
    pub fn from_records(records: Vec<SequenceRecord>) -> Self {
        // Compute metadata from the sequence records
        let metadata = SequenceCollectionMetadata::from_sequences(&records, None);

        SequenceCollection {
            metadata,
            sequences: records,
        }
    }

    /// No caching at all
    pub fn from_path_no_cache<P: AsRef<Path>>(file_path: P) -> Result<Self> {
        Self::from_path_with_cache(file_path, false, false)
    }

    pub fn from_path_with_cache<P: AsRef<Path>>(
        file_path: P,
        read_cache: bool,
        write_cache: bool,
    ) -> Result<Self> {
        // If the rgsi file exists, just use that.
        let fa_file_path = file_path.as_ref();
        let rgsi_file_path = fa_file_path.replace_exts_with("rgsi");

        // Check if the cache file exists and is valid
        if read_cache && rgsi_file_path.exists() {
            // Read the existing rgsi file
            let seqcol = read_rgsi_file(&rgsi_file_path)?;

            // Validate that the cache has sequences - empty cache is invalid
            if !seqcol.sequences.is_empty() {
                return Ok(seqcol);
            }
            // Cache is empty/stale - fall through to re-digest
            // Optionally delete the stale cache file
            let _ = std::fs::remove_file(&rgsi_file_path);
        }

        // If the rgsi file does not exist or was invalid, compute the digests
        // Digest the fasta file (your function)
        let seqcol: SequenceCollection = digest_fasta(file_path.as_ref())?;

        // Write the SequenceCollection to the RGSI file
        if write_cache && !rgsi_file_path.exists() {
            seqcol.write_rgsi()?;
        }
        Ok(seqcol)
    }

    /// Write the SequenceCollection to a collection RGSI file.
    ///
    /// Creates an RGSI file with collection-level digest headers followed by
    /// sequence metadata for all sequences in this collection.
    ///
    /// # Arguments
    /// * `file_path` - The path to the RGSI file to be written
    ///
    /// # Returns
    /// Result indicating success or error
    ///
    /// # Format
    /// The file includes:
    /// - Collection digest headers (##seqcol_digest, ##names_digest, etc.)
    /// - Column header (#name, description, length, alphabet, sha512t24u, md5)
    /// - One line per sequence with metadata
    pub fn write_collection_rgsi<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        write_rgsi_impl(file_path, &self.metadata, Some(&self.sequences))
    }

    /// Write the SequenceCollection to an RGSI file, using the file path stored in the struct.
    pub fn write_rgsi(&self) -> Result<()> {
        if let Some(file_path) = &self.metadata.file_path {
            let rgsi_file_path = file_path.replace_exts_with("rgsi");
            self.write_collection_rgsi(rgsi_file_path)
        } else {
            Err(anyhow::anyhow!(
                "No file path specified for RGSI output. Use `write_collection_rgsi` to specify a file path."
            ))
        }
    }

    /// Convert to a SequenceCollectionRecord for storage in RefgetStore
    pub fn to_record(&self) -> SequenceCollectionRecord {
        SequenceCollectionRecord::from(self.clone())
    }

    /// Write the SequenceCollection to a FASTA file.
    ///
    /// Writes all sequences in this collection to a FASTA file with the specified line width.
    /// Only sequences with loaded data (Full variant) will be written.
    ///
    /// # Arguments
    /// * `file_path` - The path to the output FASTA file
    /// * `line_width` - Number of bases per line (default: 70)
    ///
    /// # Returns
    /// Result indicating success or error
    ///
    /// # Errors
    /// Returns an error if any sequence doesn't have data loaded (Stub variant)
    pub fn write_fasta<P: AsRef<Path>>(&self, file_path: P, line_width: Option<usize>) -> Result<()> {
        let line_width = line_width.unwrap_or(70);
        let mut output_file = File::create(file_path)?;

        for record in &self.sequences {
            // Check if record has data
            if !record.is_loaded() {
                return Err(anyhow::anyhow!(
                    "Cannot write FASTA: sequence '{}' does not have data loaded",
                    record.metadata().name
                ));
            }

            // Decode the sequence
            let decoded_sequence = record.decode().ok_or_else(|| {
                anyhow::anyhow!(
                    "Failed to decode sequence '{}'",
                    record.metadata().name
                )
            })?;

            // Write FASTA header (include description if present)
            let metadata = record.metadata();
            let header = match &metadata.description {
                Some(desc) => format!(">{} {}", metadata.name, desc),
                None => format!(">{}", metadata.name),
            };
            writeln!(output_file, "{}", header)?;

            // Write sequence with line wrapping
            for chunk in decoded_sequence.as_bytes().chunks(line_width) {
                output_file.write_all(chunk)?;
                output_file.write_all(b"\n")?;
            }
        }

        Ok(())
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

// ============================================================================
// RGSI File I/O
// ============================================================================

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
            description: if parts[5].is_empty() { None } else { Some(parts[5].to_string()) },
            length: parts[1].parse().ok()?,
            alphabet: parts[2].parse().unwrap_or(AlphabetType::Unknown),
            sha512t24u: parts[3].to_string(),
            md5: parts[4].to_string(),
            fai: None,
        }),
        _ => None,
    }
}

/// Read an RGSI file and return a SequenceCollection.
///
/// RGSI (Refget Sequence Index) files contain sequence metadata in a tab-separated format.
/// This function reads the metadata without loading actual sequence data.
///
/// # Arguments
/// * `file_path` - The path to the RGSI file to be read.
///
/// # Returns
/// A SequenceCollection with sequence metadata (Stub records, no sequence data).
pub fn read_rgsi_file<T: AsRef<Path>>(file_path: T) -> Result<SequenceCollection> {
    use std::io::BufRead;

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

        // Parse sequence data line
        if let Some(metadata) = parse_rgsi_line(&line) {
            results.push(SequenceRecord::Stub(metadata));
        }
    }

    // If the digests were not found in the file, compute them from sequence records
    // Otherwise, use the digests from the file header
    let collection_metadata = if sequences_digest.is_empty() || names_digest.is_empty() || lengths_digest.is_empty() {
        // Compute from sequence records
        SequenceCollectionMetadata::from_sequences(
            &results,
            Some(file_path.as_ref().to_path_buf())
        )
    } else {
        // Use digests from file header
        let lvl1 = SeqColDigestLvl1 {
            sequences_digest,
            names_digest,
            lengths_digest,
        };

        let digest = if seqcol_digest.is_empty() {
            lvl1.to_digest()
        } else {
            seqcol_digest
        };

        SequenceCollectionMetadata {
            digest,
            n_sequences: results.len(),
            names_digest: lvl1.names_digest,
            sequences_digest: lvl1.sequences_digest,
            lengths_digest: lvl1.lengths_digest,
            file_path: Some(file_path.as_ref().to_path_buf()),
        }
    };

    Ok(SequenceCollection {
        metadata: collection_metadata,
        sequences: results,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::{digest_fasta, load_fasta};
    use crate::encoder::encode_sequence;

    #[test]
    fn test_decode_returns_none_when_no_data() {
        // Test that decode() returns None when data is None
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        for seq_record in &seqcol.sequences {
            assert!(!seq_record.is_loaded(), "digest_fasta should not load sequence data");
            assert!(matches!(seq_record, SequenceRecord::Stub(_)), "digest_fasta should create Stub records");
            assert_eq!(seq_record.decode(), None, "decode() should return None for Stub records");
        }
    }

    #[test]
    fn test_decode_with_loaded_data() {
        // Test that decode() returns the correct sequence when data is loaded
        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        // Expected sequences from base.fa
        let expected_sequences = vec![
            ("chrX", "TTGGGGAA"),
            ("chr1", "GGAA"),
            ("chr2", "GCGC"),
        ];

        for (seq_record, (expected_name, expected_seq)) in seqcol.sequences.iter().zip(expected_sequences.iter()) {
            assert_eq!(seq_record.metadata().name, *expected_name);
            assert!(seq_record.is_loaded(), "load_fasta should load sequence data");
            assert!(matches!(seq_record, SequenceRecord::Full { .. }), "load_fasta should create Full records");

            let decoded = seq_record.decode().expect("decode() should return Some when data is present");
            assert_eq!(decoded, *expected_seq,
                "Decoded sequence for {} should match expected sequence", expected_name);
        }
    }

    #[test]
    fn test_decode_handles_encoded_data() {
        // Test that decode() correctly handles bit-packed encoded data
        use crate::alphabet::{lookup_alphabet, AlphabetType};

        let sequence = b"ACGT";
        let alphabet = lookup_alphabet(&AlphabetType::Dna2bit);
        let encoded_data = encode_sequence(sequence, alphabet);

        let record = SequenceRecord::Full {
            metadata: SequenceMetadata {
                name: "test_seq".to_string(),
                description: None,
                length: 4,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Dna2bit,
                fai: None,
            },
            sequence: encoded_data,
        };

        let decoded = record.decode().expect("Should decode encoded data");
        assert_eq!(decoded, "ACGT", "Should correctly decode bit-packed DNA sequence");
    }

    #[test]
    fn test_decode_handles_raw_utf8_data() {
        // Test that decode() handles raw UTF-8 data (not encoded)
        let raw_sequence = b"ACGTACGT".to_vec();

        let record = SequenceRecord::Full {
            metadata: SequenceMetadata {
                name: "test_seq".to_string(),
                description: None,
                length: 8,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Ascii,
                fai: None,
            },
            sequence: raw_sequence,
        };

        let decoded = record.decode().expect("Should decode raw UTF-8 data");
        assert_eq!(decoded, "ACGTACGT", "Should correctly decode raw sequence data");
    }

    #[test]
    fn test_decode_with_iupac_alphabet() {
        // Test decode with IUPAC DNA alphabet (4-bit encoding)
        use crate::alphabet::{lookup_alphabet, AlphabetType};

        let sequence = b"ACGTRYMK";
        let alphabet = lookup_alphabet(&AlphabetType::DnaIupac);
        let encoded_data = encode_sequence(sequence, alphabet);

        let record = SequenceRecord::Full {
            metadata: SequenceMetadata {
                name: "iupac_test".to_string(),
                description: None,
                length: 8,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::DnaIupac,
                fai: None,
            },
            sequence: encoded_data,
        };

        let decoded = record.decode().expect("Should decode IUPAC encoded data");
        assert_eq!(decoded, "ACGTRYMK", "Should correctly decode IUPAC DNA sequence");
    }

    #[test]
    fn test_decode_with_protein_alphabet() {
        // Test decode with protein alphabet (5-bit encoding)
        use crate::alphabet::{lookup_alphabet, AlphabetType};

        let sequence = b"ACDEFGHIKLMNPQRSTVWY";
        let alphabet = lookup_alphabet(&AlphabetType::Protein);
        let encoded_data = encode_sequence(sequence, alphabet);

        let record = SequenceRecord::Full {
            metadata: SequenceMetadata {
                name: "protein_test".to_string(),
                description: None,
                length: 20,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Protein,
                fai: None,
            },
            sequence: encoded_data,
        };

        let decoded = record.decode().expect("Should decode protein encoded data");
        assert_eq!(decoded, "ACDEFGHIKLMNPQRSTVWY", "Should correctly decode protein sequence");
    }

    #[test]
    fn test_decode_empty_sequence() {
        // Test decode with empty data
        let record = SequenceRecord::Full {
            metadata: SequenceMetadata {
                name: "empty_seq".to_string(),
                description: None,
                length: 0,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Dna2bit,
                fai: None,
            },
            sequence: Vec::new(),
        };

        let decoded = record.decode().expect("Should handle empty sequence");
        assert_eq!(decoded, "", "Empty sequence should decode to empty string");
    }

    #[test]
    fn test_decode_detects_raw_vs_encoded() {
        // Test that decode() correctly distinguishes between raw and encoded data
        use crate::alphabet::{lookup_alphabet, AlphabetType};

        // Create raw UTF-8 data that matches the sequence length (heuristic for raw data)
        let raw_sequence = b"GGGGGGGG";  // 8 bytes for 8-symbol sequence
        let record_raw = SequenceRecord::Full {
            metadata: SequenceMetadata {
                name: "raw_test".to_string(),
                description: None,
                length: 8,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Dna2bit,
                fai: None,
            },
            sequence: raw_sequence.to_vec(),
        };

        let decoded_raw = record_raw.decode().expect("Should decode raw data");
        assert_eq!(decoded_raw, "GGGGGGGG");

        // Create encoded data (2 bits per symbol, so 8 symbols = 2 bytes)
        let alphabet = lookup_alphabet(&AlphabetType::Dna2bit);
        let encoded_data = encode_sequence(b"GGGGGGGG", alphabet);
        assert_eq!(encoded_data.len(), 2, "Encoded 2-bit DNA should be 2 bytes for 8 symbols");

        let record_encoded = SequenceRecord::Full {
            metadata: SequenceMetadata {
                name: "encoded_test".to_string(),
                description: None,
                length: 8,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Dna2bit,
                fai: None,
            },
            sequence: encoded_data,
        };

        let decoded_encoded = record_encoded.decode().expect("Should decode encoded data");
        assert_eq!(decoded_encoded, "GGGGGGGG");
    }

    #[test]
    fn test_sequence_collection_iterator() {
        // Test that SequenceCollection can be iterated over
        let collection = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        // Test borrowing iterator (&collection)
        let mut count = 0;
        for seq in &collection {
            assert!(seq.metadata().length > 0, "Sequence should have length");
            count += 1;
        }
        assert_eq!(count, 3, "base.fa should have 3 sequences");

        // Collection should still be usable after borrowing iteration
        assert_eq!(collection.sequences.len(), 3);

        // Test consuming iterator (collection)
        let names: Vec<String> = collection
            .into_iter()
            .map(|seq| seq.metadata().name.clone())
            .collect();

        assert_eq!(names.len(), 3);
        assert!(names.contains(&"chrX".to_string()));
        assert!(names.contains(&"chr1".to_string()));
        assert!(names.contains(&"chr2".to_string()));
    }

    // ===== SeqColDigestLvl1 Tests =====

    #[test]
    fn test_seqcol_digest_lvl1_from_metadata() {
        // Test computing lvl1 digests from sequence metadata
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let metadata_refs: Vec<&SequenceMetadata> = seqcol.sequences.iter()
            .map(|r| r.metadata())
            .collect();

        let lvl1 = SeqColDigestLvl1::from_metadata(&metadata_refs);

        // Verify digests are non-empty and valid base64url
        assert!(!lvl1.names_digest.is_empty());
        assert!(!lvl1.sequences_digest.is_empty());
        assert!(!lvl1.lengths_digest.is_empty());

        // Verify they match what the collection has
        assert_eq!(lvl1.names_digest, seqcol.metadata.names_digest);
        assert_eq!(lvl1.sequences_digest, seqcol.metadata.sequences_digest);
        assert_eq!(lvl1.lengths_digest, seqcol.metadata.lengths_digest);
    }

    #[test]
    fn test_seqcol_digest_lvl1_to_digest() {
        // Test that to_digest() produces consistent collection digest
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let computed_digest = seqcol.metadata.to_lvl1().to_digest();
        assert_eq!(computed_digest, seqcol.metadata.digest);
    }

    // ===== SequenceCollectionMetadata Tests =====

    #[test]
    fn test_sequence_collection_metadata_from_sequences() {
        use std::path::PathBuf;

        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let meta = SequenceCollectionMetadata::from_sequences(
            &seqcol.sequences,
            Some(PathBuf::from("../tests/data/fasta/base.fa"))
        );

        assert_eq!(meta.digest, seqcol.metadata.digest);
        assert_eq!(meta.n_sequences, 3);
        assert_eq!(meta.names_digest, seqcol.metadata.names_digest);
        assert!(meta.file_path.is_some());
    }

    #[test]
    fn test_sequence_collection_metadata_from_collection() {
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let meta = SequenceCollectionMetadata::from_collection(&seqcol);

        assert_eq!(meta.digest, seqcol.metadata.digest);
        assert_eq!(meta.n_sequences, seqcol.sequences.len());
        assert_eq!(meta.names_digest, seqcol.metadata.names_digest);
    }

    #[test]
    fn test_sequence_collection_metadata_to_lvl1() {
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let meta = SequenceCollectionMetadata::from_collection(&seqcol);
        let lvl1 = meta.to_lvl1();

        assert_eq!(lvl1.names_digest, seqcol.metadata.names_digest);
        assert_eq!(lvl1.sequences_digest, seqcol.metadata.sequences_digest);
        assert_eq!(lvl1.lengths_digest, seqcol.metadata.lengths_digest);
    }

    // ===== SequenceCollectionRecord Tests =====

    #[test]
    fn test_sequence_collection_record_stub() {
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let meta = SequenceCollectionMetadata::from_collection(&seqcol);
        let record = SequenceCollectionRecord::Stub(meta.clone());

        assert_eq!(record.metadata().digest, seqcol.metadata.digest);
        assert!(record.sequences().is_none());
        assert!(!record.has_sequences());
    }

    #[test]
    fn test_sequence_collection_record_full() {
        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let record: SequenceCollectionRecord = seqcol.clone().into();

        assert!(record.has_sequences());
        assert!(record.sequences().is_some());
        assert_eq!(record.sequences().unwrap().len(), 3);
    }

    #[test]
    fn test_sequence_collection_record_with_sequences() {
        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let meta = SequenceCollectionMetadata::from_collection(&seqcol);
        let stub = SequenceCollectionRecord::Stub(meta);

        // Convert stub to full by adding sequences
        let full = stub.with_sequences(seqcol.sequences.clone());

        assert!(full.has_sequences());
        assert_eq!(full.sequences().unwrap().len(), 3);
    }

    #[test]
    fn test_sequence_collection_record_to_collection() {
        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let record: SequenceCollectionRecord = seqcol.clone().into();
        let converted = record.to_collection();

        assert_eq!(converted.metadata.digest, seqcol.metadata.digest);
        assert_eq!(converted.sequences.len(), seqcol.sequences.len());
    }

    // ===== SequenceRecord Methods Tests =====

    #[test]
    fn test_sequence_record_metadata() {
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let record = &seqcol.sequences[0];
        let meta = record.metadata();

        assert_eq!(meta.name, "chrX");
        assert_eq!(meta.length, 8);
        assert!(!meta.sha512t24u.is_empty());
        assert!(!meta.md5.is_empty());
    }

    #[test]
    fn test_sequence_record_sequence() {
        // Stub should return None
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");
        assert!(seqcol.sequences[0].sequence().is_none());

        // Full should return Some
        let seqcol_with_data = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");
        assert!(seqcol_with_data.sequences[0].sequence().is_some());
    }

    #[test]
    fn test_sequence_record_with_data() {
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let stub = seqcol.sequences[0].clone();
        assert!(!stub.is_loaded());

        let full = stub.with_data(b"TTGGGGAA".to_vec());
        assert!(full.is_loaded());
        assert_eq!(full.sequence().unwrap(), b"TTGGGGAA");
    }

    #[test]
    fn test_sequence_record_load_data() {
        // Test Stub -> Full conversion (in-place, no clone)
        let metadata = SequenceMetadata {
            name: "chr1".to_string(),
            description: None,
            length: 8,
            sha512t24u: "test_digest".to_string(),
            md5: "test_md5".to_string(),
            alphabet: AlphabetType::Dna2bit,
            fai: None,
        };

        let mut record = SequenceRecord::Stub(metadata);
        assert!(!record.is_loaded());

        record.load_data(b"ATGCATGC".to_vec());
        assert!(record.is_loaded());
        assert_eq!(record.sequence().unwrap(), b"ATGCATGC");
        assert_eq!(record.metadata().name, "chr1");

        // Test Full -> Full (data replacement)
        record.load_data(b"GGGGAAAA".to_vec());
        assert_eq!(record.sequence().unwrap(), b"GGGGAAAA");
        // Metadata should be preserved
        assert_eq!(record.metadata().name, "chr1");
    }

    #[test]
    fn test_sequence_record_to_file() {
        use std::fs;
        use tempfile::tempdir;

        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let dir = tempdir().expect("Failed to create temp dir");
        let file_path = dir.path().join("test_seq.txt");

        seqcol.sequences[0].to_file(&file_path).expect("Failed to write file");

        let content = fs::read(&file_path).expect("Failed to read file");
        assert!(!content.is_empty());
    }

    // ===== SequenceMetadata disk_size Tests =====

    #[test]
    fn test_sequence_metadata_disk_size() {
        use crate::store::StorageMode;
        use crate::alphabet::AlphabetType;

        let metadata = SequenceMetadata {
            name: "test".to_string(),
            description: None,
            length: 1000,
            sha512t24u: "test".to_string(),
            md5: "test".to_string(),
            alphabet: AlphabetType::Dna2bit,
            fai: None,
        };

        // Raw mode: 1 byte per base
        assert_eq!(metadata.disk_size(&StorageMode::Raw), 1000);

        // Encoded mode for Dna2bit: 2 bits per base = 250 bytes for 1000 bases
        assert_eq!(metadata.disk_size(&StorageMode::Encoded), 250);
    }

    #[test]
    fn test_sequence_metadata_disk_size_protein() {
        use crate::store::StorageMode;
        use crate::alphabet::AlphabetType;

        let metadata = SequenceMetadata {
            name: "protein_test".to_string(),
            description: None,
            length: 100,
            sha512t24u: "test".to_string(),
            md5: "test".to_string(),
            alphabet: AlphabetType::Protein,  // 5 bits per symbol
            fai: None,
        };

        // Raw mode: 1 byte per symbol
        assert_eq!(metadata.disk_size(&StorageMode::Raw), 100);

        // Encoded mode for Protein: 5 bits per symbol = 500 bits = 63 bytes (ceil(500/8))
        assert_eq!(metadata.disk_size(&StorageMode::Encoded), 63);
    }

    // ===== SequenceCollection I/O Tests =====

    #[test]
    fn test_sequence_collection_from_fasta() {
        let seqcol = SequenceCollection::from_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load FASTA");

        assert_eq!(seqcol.sequences.len(), 3);
        assert!(!seqcol.metadata.digest.is_empty());
        assert!(seqcol.metadata.file_path.is_some());
    }

    #[test]
    fn test_sequence_collection_from_records() {
        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let records = seqcol.sequences.clone();
        let reconstructed = SequenceCollection::from_records(records);

        assert_eq!(reconstructed.metadata.digest, seqcol.metadata.digest);
        assert_eq!(reconstructed.sequences.len(), 3);
    }

    #[test]
    fn test_sequence_collection_write_and_read_rgsi() {
        use tempfile::tempdir;

        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let dir = tempdir().expect("Failed to create temp dir");
        let rgsi_path = dir.path().join("test.rgsi");

        seqcol.write_collection_rgsi(&rgsi_path).expect("Failed to write RGSI");

        // Read it back using read_rgsi_file
        let loaded = read_rgsi_file(&rgsi_path)
            .expect("Failed to read RGSI");

        assert_eq!(loaded.metadata.digest, seqcol.metadata.digest);
        assert_eq!(loaded.sequences.len(), seqcol.sequences.len());
    }

    #[test]
    fn test_sequence_collection_write_fasta() {
        use tempfile::tempdir;
        use std::fs;

        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let dir = tempdir().expect("Failed to create temp dir");
        let fasta_path = dir.path().join("output.fa");

        seqcol.write_fasta(&fasta_path, None).expect("Failed to write FASTA");

        let content = fs::read_to_string(&fasta_path).expect("Failed to read file");
        assert!(content.contains(">chrX"));
        assert!(content.contains("TTGGGGAA"));
    }

    #[test]
    fn test_sequence_collection_write_fasta_with_line_width() {
        use tempfile::tempdir;
        use std::fs;

        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let dir = tempdir().expect("Failed to create temp dir");
        let fasta_path = dir.path().join("output.fa");

        // Write with 4 chars per line
        seqcol.write_fasta(&fasta_path, Some(4)).expect("Failed to write FASTA");

        let content = fs::read_to_string(&fasta_path).expect("Failed to read file");
        // TTGGGGAA (8 chars) should be split into 2 lines
        assert!(content.contains("TTGG\n") || content.contains("GGAA\n"));
    }

    // ===== Display Implementation Tests =====

    #[test]
    fn test_sequence_collection_display() {
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let display = format!("{}", seqcol);
        assert!(display.contains("3 sequences"));
        assert!(display.contains(&seqcol.metadata.digest));
    }

    #[test]
    fn test_sequence_record_display() {
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let display = format!("{}", seqcol.sequences[0]);
        assert!(display.contains("chrX"));
        assert!(display.contains("length: 8"));
    }

    // ===== SequenceCollectionRecord write_collection_rgsi Test =====

    #[test]
    fn test_sequence_collection_record_write_rgsi() {
        use tempfile::tempdir;
        use std::fs;

        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        let record: SequenceCollectionRecord = seqcol.clone().into();

        let dir = tempdir().expect("Failed to create temp dir");
        let rgsi_path = dir.path().join("test.rgsi");

        record.write_collection_rgsi(&rgsi_path).expect("Failed to write RGSI");

        let content = fs::read_to_string(&rgsi_path).expect("Failed to read file");
        assert!(content.contains("##seqcol_digest="));
        assert!(content.contains("##names_digest="));
        assert!(content.contains("chrX"));
        assert!(content.contains("chr1"));
        assert!(content.contains("chr2"));
    }

    // ===== digest_sequence() Tests =====

    #[test]
    fn test_digest_sequence_basic() {
        use super::digest_sequence;

        let seq = digest_sequence("test_seq", b"ACGTACGT");
        assert_eq!(seq.metadata().name, "test_seq");
        assert_eq!(seq.metadata().length, 8);
        assert!(!seq.metadata().sha512t24u.is_empty());
        assert!(!seq.metadata().md5.is_empty());
        assert!(seq.is_loaded());

        // Verify data is uppercased
        let data = seq.sequence().unwrap();
        assert_eq!(data, b"ACGTACGT");
    }

    #[test]
    fn test_digest_sequence_empty() {
        use super::digest_sequence;

        let seq = digest_sequence("empty", b"");
        assert_eq!(seq.metadata().name, "empty");
        assert_eq!(seq.metadata().length, 0);
        assert!(seq.is_loaded());
        assert_eq!(seq.sequence().unwrap().len(), 0);
    }

    #[test]
    fn test_digest_sequence_protein() {
        use super::digest_sequence;
        use crate::alphabet::AlphabetType;

        // Sequence with protein-only characters (E, F, I, L, P, Q are protein-only)
        let seq = digest_sequence("protein", b"EFILPQ");
        assert_eq!(seq.metadata().alphabet, AlphabetType::Protein);
        assert_eq!(seq.metadata().length, 6);
    }

    #[test]
    fn test_digest_sequence_digest_matches() {
        use super::digest_sequence;
        use crate::digest::{sha512t24u, md5};

        // Verify digest matches what sha512t24u() returns directly
        let data = b"ACGTACGT";
        let seq = digest_sequence("test", data);
        assert_eq!(seq.metadata().sha512t24u, sha512t24u(data));
        assert_eq!(seq.metadata().md5, md5(data));
    }

    #[test]
    fn test_digest_sequence_lowercase_uppercased() {
        use super::digest_sequence;

        // Verify lowercase input is uppercased for digests
        let seq_lower = digest_sequence("test_lower", b"acgtacgt");
        let seq_upper = digest_sequence("test_upper", b"ACGTACGT");

        // Digests should match because both are uppercased
        assert_eq!(seq_lower.metadata().sha512t24u, seq_upper.metadata().sha512t24u);
        assert_eq!(seq_lower.metadata().md5, seq_upper.metadata().md5);

        // Stored data should be uppercased
        assert_eq!(seq_lower.sequence().unwrap(), b"ACGTACGT");
    }

    #[test]
    fn test_digest_sequence_no_fai() {
        use super::digest_sequence;

        // Programmatically created sequences should not have FAI data
        let seq = digest_sequence("test", b"ACGT");
        assert!(seq.metadata().fai.is_none());
    }

    #[test]
    fn test_digest_sequence_with_description() {
        use super::digest_sequence_with_description;

        let seq = digest_sequence_with_description("test_seq", Some("my description"), b"ACGT");
        assert_eq!(seq.metadata().name, "test_seq");
        assert_eq!(seq.metadata().description, Some("my description".to_string()));
        assert_eq!(seq.metadata().length, 4);

        // Without description
        let seq_no_desc = digest_sequence_with_description("test_seq2", None, b"ACGT");
        assert!(seq_no_desc.metadata().description.is_none());
    }

    #[test]
    fn test_digest_sequence_matches_fasta_digest() {
        use super::digest_sequence;

        // Verify that digest_sequence produces the same digests as digest_fasta
        // for the same sequence data

        // Load a FASTA file and get the first sequence
        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");
        let fasta_seq = &seqcol.sequences[0]; // chrX: TTGGGGAA

        // Create the same sequence with digest_sequence
        let prog_seq = digest_sequence("chrX", b"TTGGGGAA");

        // Digests should match
        assert_eq!(prog_seq.metadata().sha512t24u, fasta_seq.metadata().sha512t24u);
        assert_eq!(prog_seq.metadata().md5, fasta_seq.metadata().md5);
        assert_eq!(prog_seq.metadata().length, fasta_seq.metadata().length);
        assert_eq!(prog_seq.metadata().alphabet, fasta_seq.metadata().alphabet);
    }
}
