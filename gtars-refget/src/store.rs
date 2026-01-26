//! # RefgetStore
//!
//! A store for managing reference genome sequences with support for both
//! in-memory and disk-backed storage.
//!
//! ## Store Creation Patterns
//!
//! ### New stores (empty)
//! - `in_memory()` - All data in RAM, fast but lost on drop
//! - `on_disk(path)` - Sequences written to disk immediately, only metadata in RAM
//!
//! ### Loading existing stores
//! - `load_local(path)` - Load from local directory (lazy-loads sequences)
//! - `load_remote(path, url)` - Load from URL, caches to local directory
//!
//! ## Runtime Configuration
//!
//! ### Persistence control
//! - `enable_persistence(path)` - Start writing to disk, flush in-memory data
//! - `disable_persistence()` - Stop writing to disk (can still read)
//!
//! ### Encoding control
//! - `set_encoding_mode(mode)` - Switch between Raw and Encoded storage
//! - `enable_encoding()` - Use 2-bit encoding (space efficient)
//! - `disable_encoding()` - Use raw bytes

use super::alphabet::{AlphabetType, lookup_alphabet};
use seq_io::fasta::{Reader, Record};
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fmt::{Display, Formatter};
use std::path::{Path, PathBuf};
use std::time::Instant;

use super::encoder::SequenceEncoder;
use super::encoder::{decode_substring_from_bytes, decode_string_from_bytes, encode_sequence};
use crate::collection::{parse_rgsi_line, read_rgsi_file};
use crate::hashkeyable::HashKeyable;
use anyhow::anyhow;
use anyhow::{Context, Result};
use chrono::Utc;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use gtars_core::utils::{get_dynamic_reader, get_file_info, parse_bedlike_file};
use serde::{Deserialize, Serialize};
use std::fs::{self, File, create_dir_all};
use std::io::{BufRead, BufReader, Read, Write};
use std::str;
// Import the HashKeyable trait for converting types to a 32-byte key

// Import collection types
use super::collection::{SequenceCollection, SequenceCollectionMetadata, SequenceCollectionRecord, SequenceMetadata, SequenceRecord};

// const DEFAULT_COLLECTION_ID: [u8; 32] = [0u8; 32]; // Default collection ID for the name lookup table

const DEFAULT_COLLECTION_ID: &str = "DEFAULT_REFGET_SEQUENCE_COLLECTION"; // Default collection ID for the name lookup table
const DEFAULT_SEQDATA_PATH_TEMPLATE: &str = "sequences/%s2/%s.seq"; // Default template for sequence file paths

/// Parse a single line from an RGCI (collection index) file.
///
/// RGCI format is tab-separated with 5 columns:
/// digest, n_sequences, names_digest, sequences_digest, lengths_digest
///
/// Lines starting with '#' are treated as comments and return None.
/// Lines with fewer than 5 columns return None.
fn parse_rgci_line(line: &str) -> Option<SequenceCollectionMetadata> {
    if line.starts_with('#') {
        return None;
    }
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 5 {
        return None;
    }
    Some(SequenceCollectionMetadata {
        digest: parts[0].to_string(),
        n_sequences: parts[1].parse().ok()?,
        names_digest: parts[2].to_string(),
        sequences_digest: parts[3].to_string(),
        lengths_digest: parts[4].to_string(),
        file_path: None,
    })
}

/// Enum storing whether sequences will be stored in Raw or Encoded form
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq)]
pub enum StorageMode {
    Raw,
    Encoded,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RetrievedSequence {
    pub sequence: String,
    pub chrom_name: String,
    pub start: u32,
    pub end: u32,
}

/// Global store handling cross-collection sequence management
/// Holds a global sequence_store, which holds all sequences (across collections) so that
/// sequences are deduplicated.
/// This allows lookup by sequence digest directly (bypassing collection information).
/// The RefgetStore also holds a collections hashmap, to provide lookup by collection+name
#[derive(Debug)]
pub struct RefgetStore {
    /// SHA512t24u digest -> SequenceRecord (metadata + optional data)
    sequence_store: HashMap<[u8; 32], SequenceRecord>,
    /// MD5 digest -> SHA512t24u digest lookup
    md5_lookup: HashMap<[u8; 32], [u8; 32]>,

    /// Collection digest -> {name -> SHA512t24u digest}
    name_lookup: HashMap<[u8; 32], HashMap<String, [u8; 32]>>,
    /// Active sequence collections (now using SequenceCollectionRecord for Stub/Full pattern)
    collections: HashMap<[u8; 32], SequenceCollectionRecord>,
    /// Storage strategy for sequences
    mode: StorageMode,
    /// Where the store lives on disk (local store or cache directory)
    local_path: Option<PathBuf>,
    /// Where to pull sequences from (if remote-backed)
    remote_source: Option<String>,
    /// Template for sequence file paths (e.g., "sequences/%s2/%s.seq")
    seqdata_path_template: Option<String>,
    /// Whether to persist sequences to disk (write-through caching)
    persist_to_disk: bool,
    /// Whether to suppress progress output
    quiet: bool,
}

/// Metadata for the entire store.
/// This is used to serialize metadata to `rgstore.json`, which can be loaded by the application.
#[derive(Serialize, Deserialize, Debug)]
struct StoreMetadata {
    /// Version of the metadata format
    version: u32,
    /// Template for sequence file paths
    seqdata_path_template: String,
    /// Template for collection file paths
    collections_path_template: String,
    /// Path to the sequence metadata index file
    sequence_index: String,
    /// Path to the collection metadata index file (NEW)
    #[serde(default)]
    collection_index: Option<String>,
    /// Storage mode (Raw or Encoded)
    mode: StorageMode,
    /// Creation timestamp
    created_at: String,
}

pub struct SubstringsFromRegions<'a, K>
where
    K: AsRef<[u8]>,
{
    store: &'a mut RefgetStore,
    reader: BufReader<Box<dyn Read>>,
    collection_digest: K,
    previous_parsed_chr: String,
    current_seq_digest: String,
    line_num: usize,
}

impl<K> Iterator for SubstringsFromRegions<'_, K>
where
    K: AsRef<[u8]>,
{
    type Item = Result<RetrievedSequence, Box<dyn std::error::Error>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line_string = String::new();

        let num_bytes = self.reader.read_line(&mut line_string);
        match num_bytes {
            Ok(bytes) => {
                if bytes == 0 {
                    return None;
                }
            }
            Err(err) => return Some(Err(err.into())),
        };

        self.line_num += 1;

        let (parsed_chr, parsed_start, parsed_end) = match parse_bedlike_file(line_string.trim()) {
            Some(coords) => coords,
            None => {
                let err_str = format!(
                    "Error reading line {} because it could not be parsed as a BED-like entry: '{}'",
                    self.line_num + 1,
                    line_string
                );
                return Some(Err(err_str.into()));
            }
        };

        if parsed_start == -1 || parsed_end == -1 {
            let err_str = format!(
                "Error reading line {} due to invalid start or end coordinates: '{}'",
                self.line_num + 1,
                line_string
            );
            return Some(Err(err_str.into()));
        }

        if self.previous_parsed_chr != parsed_chr {
            self.previous_parsed_chr = parsed_chr.clone();

            let result = match self
                .store
                .get_sequence_by_name(&self.collection_digest, &parsed_chr)
            {
                Ok(seq_record) => seq_record,
                Err(e) => {
                    let err_str = format!(
                        "Line {}: sequence '{}' not found in collection '{}': {}",
                        self.line_num + 1,
                        parsed_chr,
                        String::from_utf8_lossy(self.collection_digest.as_ref()),
                        e
                    );
                    return Some(Err(err_str.into()));
                }
            };

            self.current_seq_digest = result.metadata().sha512t24u.clone();
        }

        let retrieved_substring = match self.store.get_substring(
            &self.current_seq_digest,
            parsed_start as usize,
            parsed_end as usize,
        ) {
            Ok(substring) => substring,
            Err(e) => {
                let err_str = format!(
                    "Line {}: failed to get substring for digest '{}' from {} to {}: {}",
                    self.line_num + 1,
                    self.current_seq_digest,
                    parsed_start,
                    parsed_end,
                    e
                );
                return Some(Err(err_str.into()));
            }
        };

        Some(Ok(RetrievedSequence {
            sequence: retrieved_substring,
            chrom_name: parsed_chr,
            start: parsed_start as u32, // Convert i32 to u32
            end: parsed_end as u32,     // Convert i32 to u32
        }))
    }
}

impl RefgetStore {
    /// Generic constructor. Creates a new, empty `RefgetStore`.
    /// This is a private helper - use `on_disk()` or `in_memory()` instead.
    fn new(mode: StorageMode) -> Self {
        // Initialize the name lookup with a default collection
        let mut name_lookup = HashMap::new();
        name_lookup.insert(DEFAULT_COLLECTION_ID.to_key(), HashMap::new());

        RefgetStore {
            sequence_store: HashMap::new(),
            md5_lookup: HashMap::new(),
            name_lookup,
            collections: HashMap::new(),
            mode,
            local_path: None,
            remote_source: None,
            seqdata_path_template: None,
            persist_to_disk: false,  // on_disk() overrides to true
            quiet: false,
        }
    }

    /// Set whether to suppress progress output.
    ///
    /// When quiet is true, operations like add_sequence_collection_from_fasta
    /// will not print progress messages.
    ///
    /// # Arguments
    /// * `quiet` - Whether to suppress progress output
    pub fn set_quiet(&mut self, quiet: bool) {
        self.quiet = quiet;
    }

    /// Returns whether the store is in quiet mode.
    pub fn is_quiet(&self) -> bool {
        self.quiet
    }

    /// Create a disk-backed RefgetStore
    ///
    /// Sequences are written to disk immediately and loaded on-demand (lazy loading).
    /// Only metadata is kept in memory.
    ///
    /// # Arguments
    /// * `cache_path` - Directory for storing sequences and metadata
    /// * `mode` - Storage mode (Raw or Encoded)
    ///
    /// # Returns
    /// Result with a configured disk-backed store
    ///
    /// # Example
    /// ```ignore
    /// let store = RefgetStore::on_disk("/data/store")?;
    /// store.add_sequence_collection_from_fasta("genome.fa")?;
    /// ```
    pub fn on_disk<P: AsRef<Path>>(cache_path: P) -> Result<Self> {
        let cache_path = cache_path.as_ref();
        let index_path = cache_path.join("rgstore.json");

        if index_path.exists() {
            // Load existing store
            Self::open_local(cache_path)
        } else {
            // Create new store with default Encoded mode
            let mode = StorageMode::Encoded;
            create_dir_all(cache_path)?;

            // Use private new() helper
            let mut store = Self::new(mode);
            store.local_path = Some(cache_path.to_path_buf());
            store.seqdata_path_template = Some(DEFAULT_SEQDATA_PATH_TEMPLATE.to_string());
            store.persist_to_disk = true;  // Always true for on_disk

            // Create directory structure
            create_dir_all(cache_path.join("sequences"))?;
            create_dir_all(cache_path.join("collections"))?;

            Ok(store)
        }
    }

    /// Create an in-memory RefgetStore
    ///
    /// All sequences kept in RAM for fast access.
    /// Defaults to Encoded storage mode (2-bit packing for space efficiency).
    /// Use set_encoding_mode() to change storage mode after creation.
    ///
    /// # Example
    /// ```ignore
    /// let store = RefgetStore::in_memory();
    /// store.add_sequence_collection_from_fasta("genome.fa")?;
    /// ```
    pub fn in_memory() -> Self {
        Self::new(StorageMode::Encoded)
    }

    /// Change the storage mode, re-encoding/decoding existing sequences as needed.
    ///
    /// When switching from Raw to Encoded:
    /// - All Full sequences in memory are encoded (2-bit packed)
    ///
    /// When switching from Encoded to Raw:
    /// - All Full sequences in memory are decoded back to raw bytes
    ///
    /// Note: Stub sequences (lazy-loaded from disk) are not affected.
    /// They will be loaded in the NEW mode when accessed.
    ///
    /// # Arguments
    /// * `new_mode` - The storage mode to switch to
    pub fn set_encoding_mode(&mut self, new_mode: StorageMode) {
        if self.mode == new_mode {
            return;  // No change needed
        }

        // Re-encode/decode all Full sequences in memory
        for record in self.sequence_store.values_mut() {
            match record {
                SequenceRecord::Full { metadata, sequence } => {
                    match (self.mode, new_mode) {
                        (StorageMode::Raw, StorageMode::Encoded) => {
                            // Encode: raw bytes -> 2-bit packed
                            let alphabet = lookup_alphabet(&metadata.alphabet);
                            *sequence = encode_sequence(&*sequence, alphabet);
                        }
                        (StorageMode::Encoded, StorageMode::Raw) => {
                            // Decode: 2-bit packed -> raw bytes
                            let alphabet = lookup_alphabet(&metadata.alphabet);
                            *sequence = decode_string_from_bytes(
                                &*sequence,
                                metadata.length,
                                alphabet
                            );
                        }
                        _ => {}  // Same mode, no conversion needed
                    }
                }
                SequenceRecord::Stub(_) => {
                    // Stubs don't hold sequence data, nothing to convert
                }
            }
        }

        self.mode = new_mode;
    }

    /// Enable 2-bit encoding for space efficiency.
    /// Re-encodes any existing Raw sequences in memory.
    pub fn enable_encoding(&mut self) {
        self.set_encoding_mode(StorageMode::Encoded);
    }

    /// Disable encoding, use raw byte storage.
    /// Decodes any existing Encoded sequences in memory.
    pub fn disable_encoding(&mut self) {
        self.set_encoding_mode(StorageMode::Raw);
    }

    /// Enable disk persistence for this store.
    ///
    /// Sets up the store to write sequences to disk. Any in-memory Full sequences
    /// are flushed to disk and converted to Stubs.
    ///
    /// # Arguments
    /// * `path` - Directory for storing sequences and metadata
    ///
    /// # Returns
    /// Result indicating success or error
    pub fn enable_persistence<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let path = path.as_ref();

        // Set up persistence configuration
        self.local_path = Some(path.to_path_buf());
        self.persist_to_disk = true;
        self.seqdata_path_template
            .get_or_insert_with(|| DEFAULT_SEQDATA_PATH_TEMPLATE.to_string());

        // Create directory structure
        create_dir_all(path.join("sequences"))?;
        create_dir_all(path.join("collections"))?;

        // Flush any in-memory Full sequences to disk
        let keys: Vec<[u8; 32]> = self.sequence_store.keys().cloned().collect();
        for key in keys {
            if let Some(SequenceRecord::Full { metadata, sequence }) = self.sequence_store.get(&key) {
                // Write to disk
                self.write_sequence_to_disk_single(metadata, sequence)?;
                // Convert to stub
                let stub = SequenceRecord::Stub(metadata.clone());
                self.sequence_store.insert(key, stub);
            }
        }

        // Write all collections to disk
        for record in self.collections.values() {
            self.write_collection_to_disk_single(record)?;
        }

        // Write index files
        self.write_index_files()?;

        Ok(())
    }

    /// Disable disk persistence for this store.
    ///
    /// New sequences will be kept in memory only. Existing Stub sequences
    /// can still be loaded from disk if local_path is set.
    pub fn disable_persistence(&mut self) {
        self.persist_to_disk = false;
    }

    /// Check if persistence to disk is enabled.
    pub fn is_persisting(&self) -> bool {
        self.persist_to_disk
    }

    /// Adds a sequence to the Store
    /// Ensure that it is added to the appropriate collection.
    /// If no collection is specified, it will be added to the default collection.
    ///
    /// # Arguments
    /// * `sequence_record` - The sequence to add
    /// * `collection_digest` - Collection to add to (or None for default)
    /// * `force` - If true, overwrite existing sequences. If false, skip duplicates.
    // Using Into here  instead of the Option direction allows us to accept
    // either None or [u8; 32], without having to wrap it in Some().
    pub fn add_sequence<T: Into<Option<[u8; 32]>>>(
        &mut self,
        sequence_record: SequenceRecord,
        collection_digest: T,
        force: bool,
    ) -> Result<()> {
        // Ensure collection exists; otherwise use the default collection
        let collection_digest = collection_digest
            .into()
            .unwrap_or(DEFAULT_COLLECTION_ID.to_key());
        self.collections.get(&collection_digest).ok_or_else(|| {
            anyhow::anyhow!("Collection not found for digest: {:?}", collection_digest)
        })?;

        // Get metadata from the record (works for both Stub and Full variants)
        let metadata = sequence_record.metadata();

        // Add to name lookup for the collection
        self.name_lookup
            .entry(collection_digest)
            .or_default()
            .insert(
                metadata.name.clone(),
                metadata.sha512t24u.to_key(),
            );

        // Finally, add SequenceRecord to store (consuming the object)
        self.add_sequence_record(sequence_record, force)?;

        Ok(())
    }

    /// Adds a collection, and all sequences in it, to the store.
    ///
    /// Skips collections and sequences that already exist.
    /// Use `add_sequence_collection_force()` to overwrite existing data.
    ///
    /// # Arguments
    /// * `collection` - The sequence collection to add
    pub fn add_sequence_collection(&mut self, collection: SequenceCollection) -> Result<()> {
        self.add_sequence_collection_internal(collection, false)
    }

    /// Adds a collection, and all sequences in it, to the store, overwriting existing data.
    ///
    /// Forces overwrite of collections and sequences that already exist.
    /// Use `add_sequence_collection()` to skip duplicates (safer default).
    ///
    /// # Arguments
    /// * `collection` - The sequence collection to add
    pub fn add_sequence_collection_force(&mut self, collection: SequenceCollection) -> Result<()> {
        self.add_sequence_collection_internal(collection, true)
    }

    /// Internal implementation for adding a sequence collection.
    fn add_sequence_collection_internal(&mut self, collection: SequenceCollection, force: bool) -> Result<()> {
        let coll_digest = collection.metadata.digest.to_key();

        // Check if collection already exists
        if !force && self.collections.contains_key(&coll_digest) {
            // Skip - collection already exists and force=false
            return Ok(());
        }

        // Convert to SequenceCollectionRecord
        let record = SequenceCollectionRecord::from(collection.clone());

        // Write collection to disk if persist_to_disk is enabled (before moving sequences)
        if self.persist_to_disk && self.local_path.is_some() {
            self.write_collection_to_disk_single(&record)?;
        }

        // Register the collection record
        self.collections.insert(coll_digest, record);

        // Add all sequences in the collection to the store
        for sequence_record in collection.sequences {
            self.add_sequence(sequence_record, coll_digest, force)?;
        }

        // Write index files so store is immediately loadable
        if self.persist_to_disk && self.local_path.is_some() {
            self.write_index_files()?;
        }

        Ok(())
    }

    // Adds SequenceRecord to the store.
    // Should only be used internally, via `add_sequence`, which ensures sequences are added to collections.
    // If the store is disk-backed (persist_to_disk=true), Full records are written to disk and replaced with Stubs.
    fn add_sequence_record(&mut self, sr: SequenceRecord, force: bool) -> Result<()> {
        let metadata = sr.metadata();
        let key = metadata.sha512t24u.to_key();

        // Check if sequence already exists
        if !force && self.sequence_store.contains_key(&key) {
            // Skip - sequence already exists and force=false
            return Ok(());
        }

        self.md5_lookup
            .insert(metadata.md5.to_key(), metadata.sha512t24u.to_key());

        // Check if we should write Full records to disk
        if self.persist_to_disk && self.local_path.is_some() {
            match &sr {
                SequenceRecord::Full { metadata, sequence } => {
                    // Write to disk
                    self.write_sequence_to_disk_single(metadata, sequence)?;
                    // Store as stub instead
                    let stub = SequenceRecord::Stub(metadata.clone());
                    self.sequence_store.insert(key, stub);
                    return Ok(());
                }
                SequenceRecord::Stub(_) => {
                    // Already a stub, just add it normally below
                }
            }
        }

        // Add as-is (either memory-only mode, or already a Stub)
        self.sequence_store.insert(key, sr);
        Ok(())
    }

    /// Add a sequence collection from a FASTA file.
    ///
    /// Skips sequences and collections that already exist in the store.
    /// Use `add_sequence_collection_from_fasta_force()` to overwrite existing data.
    ///
    /// # Arguments
    /// * `file_path` - Path to the FASTA file
    ///
    /// # Returns
    /// A tuple of (SequenceCollectionMetadata, was_new) where was_new indicates
    /// whether the collection was newly added (true) or already existed (false).
    ///
    /// # Notes
    /// Loading sequence data requires 2 passes through the FASTA file:
    /// 1. First pass digests and guesses the alphabet to produce SequenceMetadata
    /// 2. Second pass encodes the sequences based on the detected alphabet
    pub fn add_sequence_collection_from_fasta<P: AsRef<Path>>(&mut self, file_path: P) -> Result<(SequenceCollectionMetadata, bool)> {
        self.add_sequence_collection_from_fasta_internal(file_path, false)
    }

    /// Add a sequence collection from a FASTA file, overwriting existing data.
    ///
    /// Forces overwrite of collections and sequences that already exist in the store.
    /// Use `add_sequence_collection_from_fasta()` to skip duplicates (safer default).
    ///
    /// # Arguments
    /// * `file_path` - Path to the FASTA file
    ///
    /// # Returns
    /// A tuple of (SequenceCollectionMetadata, was_new) where was_new is always true
    /// since force mode always overwrites.
    pub fn add_sequence_collection_from_fasta_force<P: AsRef<Path>>(&mut self, file_path: P) -> Result<(SequenceCollectionMetadata, bool)> {
        self.add_sequence_collection_from_fasta_internal(file_path, true)
    }

    /// Internal implementation for adding a sequence collection from FASTA.
    /// Returns (SequenceCollectionMetadata, was_new) where was_new indicates if the collection was added.
    fn add_sequence_collection_from_fasta_internal<P: AsRef<Path>>(&mut self, file_path: P, force: bool) -> Result<(SequenceCollectionMetadata, bool)> {
        // Print start message
        if !self.quiet {
            println!("Processing {}...", file_path.as_ref().display());
        }

        // Phase 1: Digest computation
        let digest_start = Instant::now();
        let seqcol = SequenceCollection::from_fasta(&file_path)?;
        let digest_elapsed = digest_start.elapsed();

        // Get metadata directly from the collection
        let metadata = seqcol.metadata.clone();

        // Check if collection already exists and skip if not forcing
        if !force && self.collections.contains_key(&seqcol.metadata.digest.to_key()) {
            if !self.quiet {
                println!("Skipped {} (already exists)", seqcol.metadata.digest);
            }
            return Ok((metadata, false));
        }

        // Register the collection
        self.add_sequence_collection_internal(seqcol.clone(), force)?;

        // Local hashmap to store SequenceMetadata (digests)
        let mut seqmeta_hashmap: HashMap<String, SequenceMetadata> = HashMap::new();
        let seqcol_sequences = seqcol.sequences.clone(); // Clone to avoid partial move
        for record in seqcol_sequences {
            let seqmeta = record.metadata().clone();
            seqmeta_hashmap.insert(seqmeta.name.clone(), seqmeta);
        }

        let file_reader = get_dynamic_reader(file_path.as_ref())?;
        let mut fasta_reader = Reader::new(file_reader);

        // Phase 2: Load/encode sequences
        let encode_start = Instant::now();

        let mut seq_count = 0;
        while let Some(record) = fasta_reader.next() {
            let record = record?;
            let header = std::str::from_utf8(record.head())?;
            // Parse header to get name (first word) - same logic as digest_fasta
            let (name, _description) = crate::fasta::parse_fasta_header(header);
            let dr = seqmeta_hashmap.get(&name)
                .ok_or_else(|| {
                    let available_keys: Vec<_> = seqmeta_hashmap.keys().collect();
                    let total = available_keys.len();
                    let sample: Vec<_> = available_keys.iter().take(3).collect();
                    anyhow::anyhow!(
                        "Sequence '{}' not found in metadata. Available ({} total): {:?}{}",
                        name,
                        total,
                        sample,
                        if total > 3 { " ..." } else { "" }
                    )
                })?
                .clone();

            seq_count += 1;

            match self.mode {
                StorageMode::Raw => {
                    let mut raw_sequence = Vec::with_capacity(dr.length);
                    // For raw, just extend with the line content.
                    for seq_line in record.seq_lines() {
                        raw_sequence.extend(seq_line);
                    }

                    // Always replace Stubs with Full sequences from FASTA
                    self.add_sequence(
                        SequenceRecord::Full {
                            metadata: dr,
                            sequence: raw_sequence,
                        },
                        seqcol.metadata.digest.to_key(),
                        true,  // Always replace Stubs with Full
                    )?;
                }
                StorageMode::Encoded => {
                    // Create a SequenceEncoder to handle the encoding of the sequence.
                    let mut encoder = SequenceEncoder::new(dr.alphabet, dr.length);
                    for seq_line in record.seq_lines() {
                        encoder.update(seq_line);
                    }
                    let encoded_sequence = encoder.finalize();

                    // Always replace Stubs with Full sequences from FASTA
                    self.add_sequence(
                        SequenceRecord::Full {
                            metadata: dr,
                            sequence: encoded_sequence,
                        },
                        seqcol.metadata.digest.to_key(),
                        true,  // Always replace Stubs with Full
                    )?;
                }
            }
        }

        let encode_elapsed = encode_start.elapsed();

        // Print summary with timing breakdown
        if !self.quiet {
            println!(
                "Added {} ({} seqs) in {:.1}s [{:.1}s digest + {:.1}s encode]",
                seqcol.metadata.digest,
                seq_count,
                digest_elapsed.as_secs_f64() + encode_elapsed.as_secs_f64(),
                digest_elapsed.as_secs_f64(),
                encode_elapsed.as_secs_f64()
            );
        }

        // Note: If persist_to_disk=true, sequences were already written to disk
        // and replaced with stubs by add_sequence_record()

        Ok((metadata, true))
    }

    /// Returns an iterator over all sequence digests in the store
    pub fn sequence_digests(&self) -> impl Iterator<Item = [u8; 32]> + '_ {
        self.sequence_store.keys().cloned()
    }

    /// Returns an iterator over sequence metadata for all sequences in the store.
    ///
    /// This is a lightweight operation that returns only metadata (name, length, digests)
    /// without loading sequence data.
    ///
    /// # Returns
    /// An iterator over `SequenceMetadata` references.
    ///
    /// # Example
    /// ```ignore
    /// for metadata in store.sequence_metadata() {
    ///     println!("{}: {} bp", metadata.name, metadata.length);
    /// }
    /// ```
    pub fn sequence_metadata(&self) -> impl Iterator<Item = &SequenceMetadata> + '_ {
        self.sequence_store.values().map(|rec| rec.metadata())
    }

    /// Calculate the total disk size of all sequences in the store
    ///
    /// This computes the disk space used by sequence data based on:
    /// - Sequence length
    /// - Alphabet type (bits per symbol)
    /// - Storage mode (Raw or Encoded)
    ///
    /// # Returns
    /// Total bytes used for sequence data on disk
    ///
    /// # Note
    /// This only accounts for sequence data files (.seq), not metadata files
    /// like RGSI files, rgstore.json, or directory overhead.
    ///
    /// # Examples
    /// ```ignore
    /// let store = RefgetStore::on_disk("store");
    /// store.add_sequence_collection_from_fasta("genome.fa")?;
    /// let disk_size = store.total_disk_size();
    /// println!("Sequences use {} bytes on disk", disk_size);
    /// ```
    pub fn total_disk_size(&self) -> usize {
        self.sequence_store
            .values()
            .map(|rec| rec.metadata().disk_size(&self.mode))
            .sum()
    }

    /// Returns the actual disk usage of the store directory.
    ///
    /// Walks the local_path directory (if set) and sums all file sizes.
    /// For in-memory stores without a local_path, returns 0.
    ///
    /// This is useful for stats reporting to show actual disk consumption
    /// regardless of whether sequences are loaded in memory.
    pub fn actual_disk_usage(&self) -> usize {
        let Some(path) = &self.local_path else {
            return 0;
        };

        fn dir_size(path: &std::path::Path) -> usize {
            let mut total = 0;
            if let Ok(entries) = std::fs::read_dir(path) {
                for entry in entries.flatten() {
                    let path = entry.path();
                    if path.is_file() {
                        total += entry.metadata().map(|m| m.len() as usize).unwrap_or(0);
                    } else if path.is_dir() {
                        total += dir_size(&path);
                    }
                }
            }
            total
        }

        dir_size(path)
    }

    // =========================================================================
    // Collection API
    // =========================================================================

    /// List all collections in the store (metadata only, no sequence data).
    ///
    /// Returns metadata for all collections without loading sequence data.
    /// Use this for browsing/inventory operations.
    ///
    /// # Example
    /// ```ignore
    /// for meta in store.list_collections() {
    ///     println!("{}: {} sequences", meta.digest, meta.n_sequences);
    /// }
    /// ```
    pub fn list_collections(&self) -> Vec<SequenceCollectionMetadata> {
        let mut result: Vec<_> = self.collections.values().map(|record| record.metadata().clone()).collect();
        result.sort_by(|a, b| a.digest.cmp(&b.digest));
        result
    }

    /// Get metadata for a single collection by digest (no sequence data).
    ///
    /// Use this for lightweight lookups when you don't need sequence data.
    pub fn get_collection_metadata<K: AsRef<[u8]>>(&self, collection_digest: K) -> Option<&SequenceCollectionMetadata> {
        let key = collection_digest.to_key();
        self.collections.get(&key).map(|record| record.metadata())
    }

    /// Get a collection with all its sequences loaded.
    ///
    /// This loads the collection metadata and all sequence data, returning
    /// a complete `SequenceCollection` ready for use.
    ///
    /// # Example
    /// ```ignore
    /// let collection = store.get_collection("abc123")?;
    /// for seq in &collection.sequences {
    ///     println!("{}: {}", seq.metadata().name, seq.decode()?);
    /// }
    /// ```
    pub fn get_collection(&mut self, collection_digest: &str) -> Result<SequenceCollection> {
        let key = collection_digest.to_key();
        self.ensure_collection_loaded(&key)?;

        // Get all sequence digests for this collection
        let seq_digests: Vec<[u8; 32]> = self.name_lookup
            .get(&key)
            .map(|name_map| name_map.values().cloned().collect())
            .unwrap_or_default();

        // NOTE: We do NOT load sequence data here - that would be too slow for remote stores
        // with hundreds of sequences. Sequences are returned as Stubs with metadata.
        // Call decode() on individual sequences to load their data on demand.

        // Get collection metadata
        let metadata = self.collections.get(&key)
            .ok_or_else(|| anyhow!("Collection not found: {}", collection_digest))?
            .metadata()
            .clone();

        // Build sequences list from sequence_store (as Stubs with metadata only)
        let sequences: Vec<SequenceRecord> = seq_digests
            .iter()
            .filter_map(|seq_key| self.sequence_store.get(seq_key).cloned())
            .collect();

        Ok(SequenceCollection { metadata, sequences })
    }

    // =========================================================================
    // Sequence API
    // =========================================================================

    /// List all sequences in the store (metadata only, no sequence data).
    ///
    /// Returns metadata for all sequences without loading sequence data.
    /// Use this for browsing/inventory operations.
    ///
    /// # Example
    /// ```ignore
    /// for meta in store.list_sequences() {
    ///     println!("{}: {} bp", meta.name, meta.length);
    /// }
    /// ```
    pub fn list_sequences(&self) -> Vec<SequenceMetadata> {
        let mut result: Vec<_> = self.sequence_store.values().map(|rec| rec.metadata().clone()).collect();
        result.sort_by(|a, b| a.sha512t24u.cmp(&b.sha512t24u));
        result
    }

    /// Get metadata for a single sequence by digest (no sequence data).
    ///
    /// Use this for lightweight lookups when you don't need the actual sequence.
    pub fn get_sequence_metadata<K: AsRef<[u8]>>(&self, seq_digest: K) -> Option<&SequenceMetadata> {
        let key = seq_digest.to_key();
        self.sequence_store.get(&key).map(|rec| rec.metadata())
    }

    /// Get a sequence by its SHA512t24u digest, loading data if needed.
    ///
    /// # Example
    /// ```ignore
    /// let seq = store.get_sequence("abc123")?;
    /// println!("{}: {}", seq.metadata().name, seq.decode()?);
    /// ```
    pub fn get_sequence<K: AsRef<[u8]>>(&mut self, seq_digest: K) -> Result<&SequenceRecord> {
        let digest_key = seq_digest.to_key();
        // Try MD5 lookup first, fallback to using digest directly (SHA512t24u)
        let actual_key = self.md5_lookup.get(&digest_key).copied().unwrap_or(digest_key);
        self.ensure_sequence_loaded(&actual_key)?;
        self.sequence_store.get(&actual_key)
            .ok_or_else(|| anyhow!("Sequence not found: {}", String::from_utf8_lossy(seq_digest.as_ref())))
    }

    /// Get a sequence by collection digest and name, loading data if needed.
    ///
    /// # Example
    /// ```ignore
    /// let seq = store.get_sequence_by_name("collection123", "chr1")?;
    /// println!("{}", seq.decode()?);
    /// ```
    pub fn get_sequence_by_name<K: AsRef<[u8]>>(
        &mut self,
        collection_digest: K,
        sequence_name: &str,
    ) -> Result<&SequenceRecord> {
        let collection_key = collection_digest.to_key();
        self.ensure_collection_loaded(&collection_key)?;

        let digest_key = if let Some(name_map) = self.name_lookup.get(&collection_key) {
            name_map.get(sequence_name).cloned()
                .ok_or_else(|| anyhow!("Sequence '{}' not found in collection", sequence_name))?
        } else {
            return Err(anyhow!("Collection not found: {}", String::from_utf8_lossy(collection_digest.as_ref())));
        };

        self.ensure_sequence_loaded(&digest_key)?;
        self.sequence_store.get(&digest_key)
            .ok_or_else(|| anyhow!("Sequence record not found for '{}' after loading", sequence_name))
    }

    /// Iterate over all collections with their sequences loaded.
    ///
    /// This loads all collection data upfront and returns an iterator over
    /// `SequenceCollection` objects with full sequence data.
    ///
    /// # Example
    /// ```ignore
    /// for collection in store.iter_collections() {
    ///     println!("{}: {} sequences", collection.metadata.digest, collection.sequences.len());
    /// }
    /// ```
    ///
    /// Note: For browsing without loading data, use `list_collections()` instead.
    pub fn iter_collections(&mut self) -> impl Iterator<Item = SequenceCollection> {
        // Collect digests first to avoid borrow issues
        let mut digests: Vec<String> = self.collections.values()
            .map(|rec| rec.metadata().digest.clone())
            .collect();
        digests.sort();

        // Load each collection in sorted order
        let mut collections = Vec::new();
        for digest in digests {
            if let Ok(collection) = self.get_collection(&digest) {
                collections.push(collection);
            }
        }
        collections.into_iter()
    }

    /// Iterate over all sequences with their data loaded.
    ///
    /// This ensures all sequence data is loaded and returns an iterator over
    /// `SequenceRecord` objects with full sequence data.
    ///
    /// # Example
    /// ```ignore
    /// for seq in store.iter_sequences() {
    ///     println!("{}: {}", seq.metadata().name, seq.decode().unwrap_or_default());
    /// }
    /// ```
    ///
    /// Note: For browsing without loading data, use `list_sequences()` instead.
    pub fn iter_sequences(&mut self) -> impl Iterator<Item = SequenceRecord> {
        // Collect keys first to avoid borrow issues
        let keys: Vec<[u8; 32]> = self.sequence_store.keys().cloned().collect();

        // Load each sequence
        for key in &keys {
            let _ = self.ensure_sequence_loaded(key);
        }

        // Return cloned records sorted by digest
        let mut records: Vec<_> = self.sequence_store.values().cloned().collect();
        records.sort_by(|a, b| a.metadata().sha512t24u.cmp(&b.metadata().sha512t24u));
        records.into_iter()
    }

    /// Check if a collection is fully loaded (Full) or just metadata (Stub)
    pub fn is_collection_loaded<K: AsRef<[u8]>>(&self, collection_digest: K) -> bool {
        let key = collection_digest.to_key();
        self.collections.get(&key).map_or(false, |record| record.has_sequences())
    }

    /// Returns the local path where the store is located (if any)
    pub fn local_path(&self) -> Option<&PathBuf> {
        self.local_path.as_ref()
    }

    /// Returns the remote source URL (if any)
    pub fn remote_source(&self) -> Option<&str> {
        self.remote_source.as_deref()
    }

    /// Returns the storage mode used by this store
    pub fn storage_mode(&self) -> StorageMode {
        self.mode
    }

    /// Get an iterator over substrings defined by BED file regions.
    ///
    /// Reads a BED file line-by-line and yields substrings for each region.
    /// This is memory-efficient for large BED files as it streams results.
    ///
    /// # Arguments
    /// * `collection_digest` - The collection digest containing the sequences
    /// * `bed_file_path` - Path to the BED file defining regions
    ///
    /// # Returns
    /// Iterator yielding `Result<RetrievedSequence>` for each BED region
    ///
    /// # Example
    /// ```ignore
    /// let iter = store.substrings_from_regions(digest, "regions.bed")?;
    /// for result in iter {
    ///     let seq = result?;
    ///     println!("{}:{}-{}: {}", seq.chrom_name, seq.start, seq.end, seq.sequence);
    /// }
    /// ```
    pub fn substrings_from_regions<'a, K: AsRef<[u8]>>(
        &'a mut self,
        collection_digest: K,
        bed_file_path: &str,
    ) -> Result<SubstringsFromRegions<'a, K>, Box<dyn std::error::Error>> {
        let path = Path::new(bed_file_path);
        let file_info = get_file_info(path);
        let is_gzipped = file_info.is_gzipped;

        let opened_bed_file = File::open(path)?;

        let reader: Box<dyn Read> = match is_gzipped {
            true => Box::new(GzDecoder::new(BufReader::new(opened_bed_file))),
            false => Box::new(opened_bed_file),
        };
        let reader = BufReader::new(reader);

        Ok(SubstringsFromRegions {
            store: self,
            reader,
            collection_digest,
            previous_parsed_chr: String::new(),
            current_seq_digest: String::new(),
            line_num: 0,
        })
    }

    /// Export sequences from BED file regions to a FASTA file.
    ///
    /// Reads a BED file defining genomic regions and exports the sequences
    /// for those regions to a FASTA file. This is useful for extracting
    /// specific regions of interest from a genome.
    ///
    /// # Arguments
    /// * `collection_digest` - The collection digest containing the sequences
    /// * `bed_file_path` - Path to the BED file defining regions
    /// * `output_file_path` - Path to write the output FASTA file
    ///
    /// # Returns
    /// Result indicating success or error
    ///
    /// # Example
    /// ```ignore
    /// store.export_fasta_from_regions(
    ///     digest,
    ///     "regions.bed",
    ///     "output.fa"
    /// )?;
    /// ```
    pub fn export_fasta_from_regions<K: AsRef<[u8]>>(
        &mut self,
        collection_digest: K,
        bed_file_path: &str,
        output_file_path: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Set up the output path and create directories if they don't exist
        let output_path_obj = Path::new(output_file_path);
        if let Some(parent) = output_path_obj.parent() {
            create_dir_all(parent)?;
        }

        // Create output file with optional gzip compression
        let file = File::create(output_file_path)?;

        let mut writer: Box<dyn Write> = if output_path_obj.extension() == Some(OsStr::new("gz")) {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };

        // Pre-fetch all sequence metadata from the collection to avoid borrowing issues
        let collection_key = collection_digest.as_ref().to_key();

        // Ensure collection is loaded (populates name_lookup for lazy-loaded stores)
        self.ensure_collection_loaded(&collection_key)?;

        let name_to_metadata: HashMap<String, (String, usize, AlphabetType, String, String)> = self
            .name_lookup
            .get(&collection_key)
            .map(|name_map| {
                name_map
                    .iter()
                    .filter_map(|(name, seq_digest)| {
                        self.sequence_store.get(seq_digest).map(|record| {
                            let metadata = record.metadata();
                            (
                                name.clone(),
                                (
                                    metadata.name.clone(),
                                    metadata.length,
                                    metadata.alphabet,
                                    metadata.sha512t24u.clone(),
                                    metadata.md5.clone(),
                                ),
                            )
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();

        let seq_iter = self.substrings_from_regions(&collection_digest, bed_file_path)?;

        let mut previous_parsed_chr = String::new();
        let mut current_header: String = String::new();
        let mut previous_header: String = String::new();

        for rs in seq_iter.into_iter() {
            let rs = rs?;

            if previous_parsed_chr != rs.chrom_name {
                previous_parsed_chr = rs.chrom_name.clone();

                // Look up metadata from our pre-fetched map
                if let Some((name, length, alphabet, sha512, md5)) =
                    name_to_metadata.get(&rs.chrom_name)
                {
                    current_header = format!(
                        ">{} {} {} {} {}",
                        name, length, alphabet, sha512, md5
                    );
                }
            }

            let retrieved_substring = rs.sequence;

            if previous_header != current_header {
                let prefix = if previous_header.is_empty() { "" } else { "\n" };

                previous_header = current_header.clone();

                // Combine the prefix, current_header, and a trailing newline
                let header_to_be_written = format!("{}{}\n", prefix, current_header);
                writer.write_all(header_to_be_written.as_bytes())?;
            }

            writer.write_all(retrieved_substring.as_ref())?;
        }

        // Ensure all data is flushed (important for gzip)
        writer.flush()?;

        Ok(())
    }

    /// Retrieves a substring from an encoded sequence by its SHA512t24u digest.
    ///
    /// # Arguments
    ///
    /// * `sha512_digest` - The SHA512t24u digest of the sequence
    /// * `start` - The start index of the substring (inclusive)
    /// * `end` - The end index of the substring (exclusive)
    ///
    /// # Returns
    ///
    /// The substring if the sequence is found, or an error if not found or invalid range
    pub fn get_substring<K: AsRef<[u8]>>(
        &mut self,
        sha512_digest: K,
        start: usize,
        end: usize,
    ) -> Result<String> {
        let digest_key = sha512_digest.to_key();

        // Ensure the sequence data is loaded
        self.ensure_sequence_loaded(&digest_key)?;

        let record = self.sequence_store.get(&digest_key)
            .ok_or_else(|| anyhow!("Sequence not found: {}", String::from_utf8_lossy(sha512_digest.as_ref())))?;
        let (metadata, sequence) = match record {
            SequenceRecord::Stub(_) => return Err(anyhow!("Sequence data not loaded (stub only)")),
            SequenceRecord::Full { metadata, sequence } => (metadata, sequence),
        };

        if start >= metadata.length || end > metadata.length || start >= end {
            return Err(anyhow!(
                "Invalid substring range: start={}, end={}, sequence length={}",
                start, end, metadata.length
            ));
        }

        match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                let decoded_sequence = decode_substring_from_bytes(sequence, start, end, alphabet);
                String::from_utf8(decoded_sequence)
                    .map_err(|e| anyhow!("Failed to decode UTF-8 sequence: {}", e))
            }
            StorageMode::Raw => {
                let raw_slice: &[u8] = &sequence[start..end];
                String::from_utf8(raw_slice.to_vec())
                    .map_err(|e| anyhow!("Failed to decode UTF-8 sequence: {}", e))
            }
        }
    }

    /// Export sequences from a collection to a FASTA file
    ///
    /// # Arguments
    /// * `collection_digest` - The digest of the collection to export from
    /// * `output_path` - Path to write the FASTA file
    /// * `sequence_names` - Optional list of sequence names to export.
    ///                      If None, exports all sequences in the collection.
    /// * `line_width` - Optional line width for wrapping sequences (default: 80)
    ///
    /// # Returns
    /// Result indicating success or error
    pub fn export_fasta<K: AsRef<[u8]>, P: AsRef<Path>>(
        &mut self,
        collection_digest: K,
        output_path: P,
        sequence_names: Option<Vec<&str>>,
        line_width: Option<usize>,
    ) -> Result<()> {
        let line_width = line_width.unwrap_or(80);
        let output_path = output_path.as_ref();
        let collection_key = collection_digest.as_ref().to_key();

        // Ensure collection is loaded (populates name_lookup for lazy-loaded stores)
        self.ensure_collection_loaded(&collection_key)?;

        // Get the name map for this collection and build a map of name -> digest
        let name_to_digest: HashMap<String, [u8; 32]> = self
            .name_lookup
            .get(&collection_key)
            .ok_or_else(|| {
                anyhow!("Collection not found: {:?}", String::from_utf8_lossy(collection_digest.as_ref()))
            })?
            .clone();

        // Determine which sequences to export
        let names_to_export: Vec<String> = if let Some(names) = sequence_names {
            // Filter to only requested names
            names.iter().map(|s| s.to_string()).collect()
        } else {
            // Export all sequences in the collection
            name_to_digest.keys().cloned().collect()
        };

        // Create output file with optional gzip compression
        let file = File::create(output_path)
            .context(format!("Failed to create output file: {}", output_path.display()))?;

        let mut writer: Box<dyn Write> = if output_path.extension() == Some(OsStr::new("gz")) {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };

        // Export each sequence
        for seq_name in names_to_export {
            // Get the sequence digest from the name map
            let seq_digest = name_to_digest.get(&seq_name).ok_or_else(|| {
                anyhow!("Sequence '{}' not found in collection", seq_name)
            })?;

            // Ensure sequence is loaded
            self.ensure_sequence_loaded(seq_digest)?;

            // Get the sequence record
            let record = self.sequence_store.get(seq_digest).ok_or_else(|| {
                anyhow!("Sequence record not found for digest: {:?}", seq_digest)
            })?;

            // Get the sequence data
            let (metadata, sequence_data) = match record {
                SequenceRecord::Stub(_) => {
                    return Err(anyhow!("Sequence data not loaded for '{}'", seq_name));
                }
                SequenceRecord::Full { metadata, sequence } => (metadata, sequence),
            };

            // Decode the sequence based on storage mode
            let decoded_sequence = match self.mode {
                StorageMode::Encoded => {
                    let alphabet = lookup_alphabet(&metadata.alphabet);
                    let decoded = decode_substring_from_bytes(
                        sequence_data,
                        0,
                        metadata.length,
                        alphabet,
                    );
                    String::from_utf8(decoded).context("Failed to decode sequence as UTF-8")?
                }
                StorageMode::Raw => {
                    String::from_utf8(sequence_data.clone())
                        .context("Failed to decode raw sequence as UTF-8")?
                }
            };

            // Write FASTA header (include description if present)
            let header = match &metadata.description {
                Some(desc) => format!(">{} {}", metadata.name, desc),
                None => format!(">{}", metadata.name),
            };
            writeln!(writer, "{}", header)?;

            // Write sequence with line wrapping
            for chunk in decoded_sequence.as_bytes().chunks(line_width) {
                writer.write_all(chunk)?;
                writer.write_all(b"\n")?;
            }
        }

        // Ensure all data is flushed (important for gzip)
        writer.flush()?;

        Ok(())
    }

    /// Export sequences by their sequence digests to a FASTA file
    ///
    /// Bypasses collection information and exports sequences directly via sequence digests.
    /// # Arguments
    /// * `seq_digests` - List of SHA512t24u sequence digests (not collection digests) to export
    /// * `output_path` - Path to write the FASTA file
    /// * `line_width` - Optional line width for wrapping sequences (default: 80)
    ///
    /// # Returns
    /// Result indicating success or error
    pub fn export_fasta_by_digests<P: AsRef<Path>>(
        &mut self,
        seq_digests: Vec<&str>,
        output_path: P,
        line_width: Option<usize>,
    ) -> Result<()> {
        let line_width = line_width.unwrap_or(80);
        let output_path = output_path.as_ref();

        // Create output file with optional gzip compression
        let file = File::create(output_path)
            .context(format!("Failed to create output file: {}", output_path.display()))?;

        let mut writer: Box<dyn Write> = if output_path.extension() == Some(OsStr::new("gz")) {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };

        // Export each sequence
        for digest_str in seq_digests {
            let digest_key = digest_str.as_bytes().to_key();

            // Ensure sequence is loaded
            self.ensure_sequence_loaded(&digest_key)?;

            // Get the sequence record
            let record = self.sequence_store.get(&digest_key).ok_or_else(|| {
                anyhow!("Sequence record not found for digest: {}", digest_str)
            })?;

            // Get the sequence data
            let (metadata, sequence_data) = match record {
                SequenceRecord::Stub(_) => {
                    return Err(anyhow!("Sequence data not loaded for digest: {}", digest_str));
                }
                SequenceRecord::Full { metadata, sequence } => (metadata, sequence),
            };

            // Decode the sequence based on storage mode
            let decoded_sequence = match self.mode {
                StorageMode::Encoded => {
                    let alphabet = lookup_alphabet(&metadata.alphabet);
                    let decoded = decode_substring_from_bytes(
                        sequence_data,
                        0,
                        metadata.length,
                        alphabet,
                    );
                    String::from_utf8(decoded).context("Failed to decode sequence as UTF-8")?
                }
                StorageMode::Raw => {
                    String::from_utf8(sequence_data.clone())
                        .context("Failed to decode raw sequence as UTF-8")?
                }
            };

            // Write FASTA header (include description if present)
            let header = match &metadata.description {
                Some(desc) => format!(">{} {}", metadata.name, desc),
                None => format!(">{}", metadata.name),
            };
            writeln!(writer, "{}", header)?;

            // Write sequence with line wrapping
            for chunk in decoded_sequence.as_bytes().chunks(line_width) {
                writer.write_all(chunk)?;
                writer.write_all(b"\n")?;
            }
        }

        // Ensure all data is flushed (important for gzip)
        writer.flush()?;

        Ok(())
    }

    /// Helper function to get the relative path for a sequence based on its SHA512t24u digest string
    fn get_sequence_path(digest_str: &str, template: &str) -> PathBuf {
        let path_str = template
            .replace("%s2", &digest_str[0..2])
            .replace("%s", digest_str);

        PathBuf::from(path_str)
    }

    /// Write a single sequence to disk using the configured path template
    fn write_sequence_to_disk_single(
        &self,
        metadata: &SequenceMetadata,
        sequence: &[u8],
    ) -> Result<()> {
        let template = self.seqdata_path_template.as_ref()
            .context("seqdata_path_template not set")?;
        let local_path = self.local_path.as_ref()
            .context("local_path not set")?;

        // Build path using template
        let seq_file_path = Self::get_sequence_path(&metadata.sha512t24u, template);
        let full_path = local_path.join(&seq_file_path);

        // Create parent directory
        if let Some(parent) = full_path.parent() {
            create_dir_all(parent)?;
        }

        // Write sequence data
        let mut file = File::create(&full_path)?;
        file.write_all(sequence)?;

        Ok(())
    }

    /// Write a single collection RGSI file to disk
    /// Used when persist_to_disk=true to persist collections incrementally
    fn write_collection_to_disk_single(&self, record: &SequenceCollectionRecord) -> Result<()> {
        let local_path = self.local_path.as_ref()
            .context("local_path not set")?;

        // Build path: collections/{digest}.rgsi
        let coll_file_path = format!("collections/{}.rgsi", record.metadata().digest);
        let full_path = local_path.join(&coll_file_path);

        // Create parent directory
        if let Some(parent) = full_path.parent() {
            create_dir_all(parent)?;
        }

        // Write collection RGSI file
        record.write_collection_rgsi(&full_path)?;

        Ok(())
    }

    /// Write index files (sequences.rgsi, collections.rgci, and rgstore.json) to disk
    ///
    /// This allows the store to be loaded later via load_local().
    /// Called automatically when adding collections in disk-backed mode.
    fn write_index_files(&self) -> Result<()> {
        let local_path = self.local_path.as_ref()
            .context("local_path not set")?;
        let template = self.seqdata_path_template.as_ref()
            .context("seqdata_path_template not set")?;

        // Write the sequence metadata index file (sequences.rgsi)
        let sequence_index_path = local_path.join("sequences.rgsi");
        self.write_sequences_rgsi(&sequence_index_path)?;

        // Write the collection metadata index file (NEW)
        let collection_index_path = local_path.join("collections.rgci");
        self.write_collections_rgci(&collection_index_path)?;

        // Create the metadata structure
        let metadata = StoreMetadata {
            version: 1,
            seqdata_path_template: template.clone(),
            collections_path_template: "collections/%s.rgsi".to_string(),
            sequence_index: "sequences.rgsi".to_string(),
            collection_index: Some("collections.rgci".to_string()),
            mode: self.mode,
            created_at: Utc::now().to_rfc3339(),
        };

        // Write metadata to rgstore.json
        let json = serde_json::to_string_pretty(&metadata)
            .context("Failed to serialize metadata to JSON")?;
        fs::write(local_path.join("rgstore.json"), json)
            .context("Failed to write rgstore.json")?;

        Ok(())
    }

    /// Write collection metadata index (collections.rgci) to disk
    ///
    /// Creates a master index of all collections with their metadata.
    /// Format: TSV with columns: digest, n_sequences, names_digest, sequences_digest, lengths_digest
    fn write_collections_rgci<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        let file_path = file_path.as_ref();
        let mut file = File::create(file_path)?;

        // Write header
        writeln!(file, "#digest\tn_sequences\tnames_digest\tsequences_digest\tlengths_digest")?;

        // Write collection metadata for all collections
        for record in self.collections.values() {
            let meta = record.metadata();
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}",
                meta.digest,
                meta.n_sequences,
                meta.names_digest,
                meta.sequences_digest,
                meta.lengths_digest,
            )?;
        }
        Ok(())
    }

    /// Write all sequence metadata to an RGSI file.
    ///
    /// Creates a global sequence index file containing metadata for all sequences
    /// in the store across all collections.
    pub fn write_sequences_rgsi<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        let file_path = file_path.as_ref();
        let mut file = std::fs::File::create(file_path)?;

        // Write header with column names (6-column format, description at end)
        writeln!(file, "#name\tlength\talphabet\tsha512t24u\tmd5\tdescription")?;

        // Write sequence metadata for all sequences
        for result_sr in self.sequence_store.values() {
            let result = result_sr.metadata().clone();
            let description = result.description.as_deref().unwrap_or("");
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}",
                result.name, result.length, result.alphabet, result.sha512t24u, result.md5, description
            )?;
        }
        Ok(())
    }

    /// Validate a relative path to prevent directory traversal attacks.
    /// Rejects absolute paths, paths with "..", and paths with null bytes.
    fn sanitize_relative_path(path: &str) -> Result<()> {
        if path.starts_with('/') || path.starts_with('\\') {
            return Err(anyhow!("Absolute paths not allowed: {}", path));
        }
        if path.contains("..") {
            return Err(anyhow!("Directory traversal not allowed: {}", path));
        }
        if path.contains('\0') {
            return Err(anyhow!("Null bytes not allowed in path"));
        }
        Ok(())
    }

    /// Helper function to fetch a file from local path or remote source
    /// Returns the file contents as Vec<u8>
    fn fetch_file(
        local_path: &Option<PathBuf>,
        remote_source: &Option<String>,
        relative_path: &str,
        persist_to_disk: bool,
    ) -> Result<Vec<u8>> {
        // Validate the relative path to prevent directory traversal
        Self::sanitize_relative_path(relative_path)?;

        // Check if file exists locally (only if caching is enabled and path exists)
        if persist_to_disk {
            if let Some(local_path) = local_path {
                let full_local_path = local_path.join(relative_path);
                if full_local_path.exists() {
                    return fs::read(&full_local_path)
                        .context(format!("Failed to read local file: {}", full_local_path.display()));
                }
            }
        }

        // If not local and we have a remote source, fetch from remote
        if let Some(remote_url) = remote_source {
            let full_remote_url = if remote_url.ends_with('/') {
                format!("{}{}", remote_url, relative_path)
            } else {
                format!("{}/{}", remote_url, relative_path)
            };

            let response = ureq::get(&full_remote_url)
                .call()
                .map_err(|e| anyhow!("Failed to fetch from remote: {}", e))?;

            let mut data = Vec::new();
            response
                .into_reader()
                .read_to_end(&mut data)
                .context("Failed to read response body")?;

            // Save to local cache only if caching is enabled
            if persist_to_disk {
                if let Some(local_path) = local_path {
                    let full_local_path = local_path.join(relative_path);

                    // Create parent directory if needed
                    if let Some(parent) = full_local_path.parent() {
                        create_dir_all(parent)?;
                    }

                    // Save to disk
                    fs::write(&full_local_path, &data)
                        .context(format!("Failed to cache file to: {}", full_local_path.display()))?;
                }
            }

            Ok(data)
        } else {
            Err(anyhow!(
                "File not found locally and no remote source configured: {}",
                relative_path
            ))
        }
    }

    /// Open a local RefgetStore from a directory.
    ///
    /// This loads only lightweight metadata and stubs. Collections and sequences
    /// remain as stubs until explicitly loaded with load_collection()/load_sequence().
    ///
    /// # Arguments
    /// * `path` - Path to the store directory
    ///
    /// Expects: rgstore.json, sequences.rgsi, collections.rgci, collections/*.rgsi
    pub fn open_local<P: AsRef<Path>>(path: P) -> Result<Self> {
        let root_path = path.as_ref();

        // Read rgstore.json index file
        let index_path = root_path.join("rgstore.json");
        let json = fs::read_to_string(&index_path).context(format!(
            "Failed to read rgstore.json from {}",
            index_path.display()
        ))?;

        let metadata: StoreMetadata =
            serde_json::from_str(&json).context("Failed to parse store metadata")?;

        // Validate paths from metadata to prevent directory traversal
        Self::sanitize_relative_path(&metadata.seqdata_path_template)?;
        Self::sanitize_relative_path(&metadata.sequence_index)?;
        if let Some(ref ci) = metadata.collection_index {
            Self::sanitize_relative_path(ci)?;
        }

        // Create a new empty store with the correct mode
        let mut store = RefgetStore::new(metadata.mode);
        store.local_path = Some(root_path.to_path_buf());
        store.seqdata_path_template = Some(metadata.seqdata_path_template.clone());
        store.persist_to_disk = true;  // Local stores always use disk

        // Load sequence metadata from the sequence index file (metadata only, no data)
        let sequence_index_path = root_path.join(&metadata.sequence_index);
        if sequence_index_path.exists() {
            Self::load_sequences_from_index(&mut store, &sequence_index_path)?;
        }

        // Try to load collection stubs from collections.rgci (new format)
        if let Some(ref collection_index) = metadata.collection_index {
            let collection_index_path = root_path.join(collection_index);
            if collection_index_path.exists() {
                Self::load_collection_stubs_from_rgci(&mut store, &collection_index_path)?;
            }
        }

        // If no collection stubs loaded (missing rgci), load full collections from directory
        if store.collections.is_empty() {
            let collections_dir = root_path.join("collections");
            Self::load_collections_from_directory(&mut store, &collections_dir)?;
        }

        Ok(store)
    }

    /// Load sequence metadata from a sequence index file (sequences.rgsi)
    fn load_sequences_from_index(store: &mut RefgetStore, index_path: &Path) -> Result<()> {
        let file = std::fs::File::open(index_path)?;
        let reader = std::io::BufReader::new(file);

        for line in reader.lines() {
            let line = line?;

            // Skip comment lines
            if line.starts_with('#') {
                continue;
            }

            // Parse sequence metadata line
            if let Some(seq_metadata) = parse_rgsi_line(&line) {
                // Create a SequenceRecord with no data (lazy loading)
                let record = SequenceRecord::Stub(seq_metadata.clone());

                // Add to store
                let sha512_key = seq_metadata.sha512t24u.to_key();
                store.sequence_store.insert(sha512_key, record);

                // Add to MD5 lookup
                let md5_key = seq_metadata.md5.to_key();
                store.md5_lookup.insert(md5_key, sha512_key);
            }
        }

        Ok(())
    }

    /// Load collection stubs from collections.rgci index file (new format)
    fn load_collection_stubs_from_rgci(store: &mut RefgetStore, index_path: &Path) -> Result<()> {
        let file = std::fs::File::open(index_path)?;
        let reader = std::io::BufReader::new(file);

        for line in reader.lines() {
            let line = line?;

            if let Some(metadata) = parse_rgci_line(&line) {
                let key = metadata.digest.to_key();
                // Create a SequenceCollectionRecord::Stub (sequences not loaded)
                // Note: name_lookup is NOT populated for stubs - it will be populated
                // when the collection is loaded via ensure_collection_loaded()
                store.collections.insert(key, SequenceCollectionRecord::Stub(metadata));
            }
        }

        Ok(())
    }

    /// Load full collections from a collections directory (fallback when no RGCI exists).
    ///
    /// Reads all .rgsi files in the directory and loads them as Full collections.
    fn load_collections_from_directory(store: &mut RefgetStore, collections_dir: &Path) -> Result<()> {
        if !collections_dir.exists() {
            return Ok(());
        }

        for entry in fs::read_dir(collections_dir)? {
            let entry = entry?;
            let path = entry.path();

            if path.is_file() && path.extension() == Some(OsStr::new("rgsi")) {
                // Load the collection from the file
                let collection = read_rgsi_file(&path)?;
                let collection_digest = collection.metadata.digest.to_key();

                // Convert to SequenceCollectionRecord::Full
                let record = SequenceCollectionRecord::from(collection.clone());

                // Add collection record to store
                store.collections.insert(collection_digest, record);

                // Build name lookup for this collection
                let mut name_map = HashMap::new();
                for sequence_record in &collection.sequences {
                    let metadata = sequence_record.metadata();
                    let sha512_key = metadata.sha512t24u.to_key();
                    name_map.insert(metadata.name.clone(), sha512_key);
                }
                store.name_lookup.insert(collection_digest, name_map);
            }
        }

        Ok(())
    }

    /// Open a remote RefgetStore with local caching.
    ///
    /// This loads only lightweight metadata and stubs from the remote URL.
    /// Data is fetched on-demand when load_collection()/load_sequence() is called.
    ///
    /// # Arguments
    /// * `cache_path` - Local directory for caching fetched data
    /// * `remote_url` - URL of the remote store
    ///
    /// # Notes
    /// By default, persistence is enabled (sequences are cached to disk).
    /// Call `disable_persistence()` after loading to keep only in memory.
    pub fn open_remote<P: AsRef<Path>, S: AsRef<str>>(
        cache_path: P,
        remote_url: S,
    ) -> Result<Self> {
        let cache_path = cache_path.as_ref();
        let remote_url = remote_url.as_ref().to_string();

        // Create cache directory
        create_dir_all(cache_path)?;

        // Fetch rgstore.json from remote
        let index_data = Self::fetch_file(&Some(cache_path.to_path_buf()), &Some(remote_url.clone()), "rgstore.json", true)?;

        let json = String::from_utf8(index_data)
            .context("Store metadata contains invalid UTF-8")?;

        let metadata: StoreMetadata =
            serde_json::from_str(&json).context("Failed to parse store metadata")?;

        // Validate paths from metadata to prevent directory traversal
        Self::sanitize_relative_path(&metadata.seqdata_path_template)?;
        Self::sanitize_relative_path(&metadata.sequence_index)?;
        if let Some(ref ci) = metadata.collection_index {
            Self::sanitize_relative_path(ci)?;
        }

        // Create a new empty store with the correct mode
        let mut store = RefgetStore::new(metadata.mode);
        store.local_path = Some(cache_path.to_path_buf());
        store.remote_source = Some(remote_url.clone());
        store.seqdata_path_template = Some(metadata.seqdata_path_template.clone());
        store.persist_to_disk = true;  // Default to true; user can call disable_persistence() after

        // Fetch sequence index from remote (always cache metadata - it's small)
        let sequence_index_data = Self::fetch_file(
            &Some(cache_path.to_path_buf()),
            &Some(remote_url.clone()),
            &metadata.sequence_index,
            true,  // Always cache metadata
        )?;
        let sequence_index_str = String::from_utf8(sequence_index_data)
            .context("sequence index contains invalid UTF-8")?;

        // Parse sequence metadata (metadata only, no data)
        for line in sequence_index_str.lines() {
            // Skip comment lines
            if line.starts_with('#') {
                continue;
            }

            // Parse sequence metadata line
            if let Some(seq_metadata) = parse_rgsi_line(line) {
                // Create a SequenceRecord with no data (lazy loading)
                let record = SequenceRecord::Stub(seq_metadata.clone());

                // Add to store
                let sha512_key = seq_metadata.sha512t24u.to_key();
                store.sequence_store.insert(sha512_key, record);

                // Add to MD5 lookup
                let md5_key = seq_metadata.md5.to_key();
                store.md5_lookup.insert(md5_key, sha512_key);
            }
        }

        // Try to fetch and load collection stubs from collections.rgci (new format)
        if let Some(ref collection_index) = metadata.collection_index {
            if let Ok(collection_index_data) = Self::fetch_file(
                &Some(cache_path.to_path_buf()),
                &Some(remote_url.clone()),
                collection_index,
                true,
            ) {
                let collection_index_str = String::from_utf8(collection_index_data)
                    .context("collection index contains invalid UTF-8")?;

                // Parse collection stubs
                for line in collection_index_str.lines() {
                    if let Some(coll_metadata) = parse_rgci_line(line) {
                        let key = coll_metadata.digest.to_key();
                        store.collections.insert(key, SequenceCollectionRecord::Stub(coll_metadata));
                    }
                }
            }
        }

        // If no collection stubs loaded, check for cached collections in local directory
        if store.collections.is_empty() {
            let local_collections_dir = cache_path.join("collections");
            create_dir_all(&local_collections_dir)?;  // Ensure cache dir exists
            Self::load_collections_from_directory(&mut store, &local_collections_dir)?;
        }

        Ok(store)
    }

    /// Ensure a collection is loaded into the store
    /// If the collection is a Stub, try to fetch full data from local or remote
    /// and upgrade it to Full. Also builds name_lookup for the collection.
    fn ensure_collection_loaded(&mut self, collection_digest: &[u8; 32]) -> Result<()> {
        // Check if name_lookup is already populated for this collection
        if self.name_lookup.contains_key(collection_digest) {
            return Ok(());
        }

        // Check if we have a Stub that needs to be loaded
        let needs_fetch = match self.collections.get(collection_digest) {
            Some(SequenceCollectionRecord::Stub(_)) => true,
            Some(SequenceCollectionRecord::Full { .. }) => false,
            None => true,  // Not in collections at all, need to fetch
        };

        if needs_fetch {
            // Get the digest string (either from Stub or from the key)
            let digest_str = if let Some(SequenceCollectionRecord::Stub(meta)) = self.collections.get(collection_digest) {
                meta.digest.clone()
            } else {
                String::from_utf8_lossy(collection_digest).to_string()
            };

            let relative_path = format!("collections/{}.rgsi", digest_str);

            // Fetch the collection file
            // Always cache metadata files (they're small), even when persist_to_disk is false
            if !self.quiet {
                let cached = self.local_path.as_ref()
                    .map(|p| p.join(&relative_path).exists())
                    .unwrap_or(false);
                let verb = if cached { "Loading" } else { "Downloading" };
                eprintln!("{} collection {}...", verb, digest_str);
            }
            let _collection_data = Self::fetch_file(&self.local_path, &self.remote_source, &relative_path, true)?;

            // Read the collection from the cached file
            let local_path = self.local_path
                .as_ref()
                .ok_or_else(|| anyhow!("No local path configured"))?;

            let collection_file_path = local_path.join(&relative_path);

            let collection = read_rgsi_file(&collection_file_path)?;

            // Verify the collection digest matches what we requested
            let loaded_digest = collection.metadata.digest.to_key();
            if loaded_digest != *collection_digest {
                return Err(anyhow!(
                    "Collection digest mismatch: expected {}, got {}",
                    String::from_utf8_lossy(collection_digest),
                    String::from_utf8_lossy(&loaded_digest)
                ));
            }

            // Convert to SequenceCollectionRecord::Full and replace Stub if present
            let record = SequenceCollectionRecord::from(collection.clone());

            // Add collection to store (replacing Stub if present)
            self.collections.insert(*collection_digest, record);

            // Build name lookup and add sequences to sequence_store as Stubs
            let mut name_map = HashMap::new();
            for sequence_record in &collection.sequences {
                let metadata = sequence_record.metadata();
                let sha512_key = metadata.sha512t24u.to_key();
                name_map.insert(metadata.name.clone(), sha512_key);

                // Add to sequence_store if not already present (as Stub for lazy loading)
                if !self.sequence_store.contains_key(&sha512_key) {
                    self.sequence_store.insert(sha512_key, SequenceRecord::Stub(metadata.clone()));
                    // Also add MD5 lookup
                    let md5_key = metadata.md5.to_key();
                    self.md5_lookup.insert(md5_key, sha512_key);
                }
            }
            self.name_lookup.insert(*collection_digest, name_map);
        } else {
            // Collection is Full but name_lookup not built yet - build it now
            // First, collect the data we need to avoid borrow conflicts
            let sequences_data: Vec<(SequenceMetadata, [u8; 32], [u8; 32])> =
                if let Some(SequenceCollectionRecord::Full { sequences, .. }) = self.collections.get(collection_digest) {
                    sequences.iter().map(|seq| {
                        let metadata = seq.metadata().clone();
                        let sha512_key = metadata.sha512t24u.to_key();
                        let md5_key = metadata.md5.to_key();
                        (metadata, sha512_key, md5_key)
                    }).collect()
                } else {
                    Vec::new()
                };

            // Now build name_lookup and add sequences to sequence_store
            let mut name_map = HashMap::new();
            for (metadata, sha512_key, md5_key) in sequences_data {
                name_map.insert(metadata.name.clone(), sha512_key);

                // Add to sequence_store if not already present
                if !self.sequence_store.contains_key(&sha512_key) {
                    self.sequence_store.insert(sha512_key, SequenceRecord::Stub(metadata));
                    self.md5_lookup.insert(md5_key, sha512_key);
                }
            }
            self.name_lookup.insert(*collection_digest, name_map);
        }

        Ok(())
    }

    /// Ensure a sequence is loaded into memory
    /// If the sequence data is not already loaded, fetch it from local or remote
    fn ensure_sequence_loaded(&mut self, digest: &[u8; 32]) -> Result<()> {
        // Check if sequence exists
        let record = self
            .sequence_store
            .get(digest)
            .ok_or_else(|| anyhow!("Sequence not found in store"))?;

        // If data is already loaded, return early
        if matches!(record, SequenceRecord::Full { .. }) {
            return Ok(());
        }

        // Get the necessary information before borrowing mutably
        let digest_str = &record.metadata().sha512t24u;
        let template = self
            .seqdata_path_template
            .as_ref()
            .ok_or_else(|| anyhow!("No sequence data path template configured"))?;

        // Build the relative path using the template
        let relative_path = template
            .replace("%s2", &digest_str[0..2])
            .replace("%s4", &digest_str[0..4])
            .replace("%s", digest_str);

        // Fetch the sequence data
        // Use persist_to_disk flag - this is where memory-only mode saves disk I/O
        if !self.quiet {
            let cached = self.local_path.as_ref().map(|p| p.join(&relative_path).exists()).unwrap_or(false);
            let verb = if cached { "Loading" } else { "Downloading" };
            eprintln!("{} sequence {}...", verb, digest_str);
        }
        let data = Self::fetch_file(&self.local_path, &self.remote_source, &relative_path, self.persist_to_disk)?;

        // Update the record with the loaded data (in-place, no clone needed)
        self.sequence_store.entry(*digest).and_modify(|r| {
            r.load_data(data);
        });

        Ok(())
    }

    /// Write all sequence metadata to an RGSI file (without collection headers).
    ///
    /// Creates a global sequence index file containing metadata for all sequences
    /// in the store across all collections. Does not include collection-level digest headers.
    ///

    /// Write the store using its configured paths
    ///
    /// For disk-backed stores (on_disk), this updates index files only since
    /// sequences/collections are already written incrementally.
    /// For in-memory stores, this is not supported (use write_store_to_dir instead).
    ///
    /// # Returns
    /// Result indicating success or error
    ///
    /// # Errors
    /// Returns an error if `local_path` is not set.
    ///
    /// # Example
    /// ```ignore
    /// let store = RefgetStore::on_disk("/data/store")?;
    /// store.add_sequence_collection_from_fasta("genome.fa")?;
    /// store.write()?;  // Updates index files
    /// ```
    pub fn write(&self) -> Result<()> {
        if !self.persist_to_disk {
            return Err(anyhow!("write() only works with disk-backed stores - use write_store_to_dir() instead"));
        }

        // For disk-backed stores, just update indexes (sequences/collections already written)
        self.write_index_files()
    }

    /// Write a RefgetStore object to a directory
    pub fn write_store_to_dir<P: AsRef<Path>>(
        &self,
        root_path: P,
        seqdata_path_template: Option<&str>,
    ) -> Result<()> {
        let root_path = root_path.as_ref();

        // Use provided template, or store's template, or default
        let template = seqdata_path_template
            .or(self.seqdata_path_template.as_deref())
            .unwrap_or(DEFAULT_SEQDATA_PATH_TEMPLATE);

        println!(
            "Writing store to directory: {}; Using seqdata path template: {}",
            root_path.display(),
            template
        );

        // Create the root directory if it doesn't exist
        fs::create_dir_all(root_path)?;

        // Create sequences directory
        let sequences_dir = root_path.join("sequences");
        fs::create_dir_all(&sequences_dir)?;

        // Create collections directory
        let collections_dir = root_path.join("collections");
        fs::create_dir_all(&collections_dir)?;

        // Write each sequence to its own file
        for record in self.sequence_store.values() {
            match record {
                SequenceRecord::Full { metadata, .. } => {
                    // Get the path for this sequence using the template and base64url-encoded digest
                    let rel_path =
                        Self::get_sequence_path(&metadata.sha512t24u, template);
                    let full_path = root_path.join(&rel_path);

                    // Write the sequence data to file
                    record.to_file(full_path)?;
                }
                SequenceRecord::Stub(_metadata) => {
                    // Stub means sequence already on disk - skip writing
                    continue;
                }
            }
        }

        // Write each collection to its own .rgsi file
        for record in self.collections.values() {
            let collection_file_path =
                root_path.join(format!("collections/{}.rgsi", record.metadata().digest));
            record.write_collection_rgsi(&collection_file_path)?;
        }

        // Write the sequence metadata index file
        let sequence_index_path = root_path.join("sequences.rgsi");
        self.write_sequences_rgsi(&sequence_index_path)?;

        // Write the collection metadata index file
        let collection_index_path = root_path.join("collections.rgci");
        self.write_collections_rgci(&collection_index_path)?;

        // Create the metadata structure
        let metadata = StoreMetadata {
            version: 1,
            seqdata_path_template: template.to_string(),
            collections_path_template: "collections/%s.rgsi".to_string(),
            sequence_index: "sequences.rgsi".to_string(),
            collection_index: Some("collections.rgci".to_string()),
            mode: self.mode,
            created_at: Utc::now().to_rfc3339(),
        };

        // Write metadata to rgstore.json
        let json = serde_json::to_string_pretty(&metadata)
            .context("Failed to serialize metadata to JSON")?;
        fs::write(root_path.join("rgstore.json"), json).context("Failed to write rgstore.json")?;

        Ok(())
    }

    /// Returns statistics about the store
    ///
    /// # Returns
    /// A tuple of (n_sequences, n_collections_loaded, storage_mode_str)
    ///
    /// Note: n_collections_loaded only reflects collections currently loaded in memory.
    /// For remote stores, collections are loaded on-demand when accessed.
    pub fn stats(&self) -> (usize, usize, &'static str) {
        let n_sequences = self.sequence_store.len();
        let n_collections_loaded = self.collections.values()
            .filter(|record| record.has_sequences())
            .count();
        let mode_str = match self.mode {
            StorageMode::Raw => "Raw",
            StorageMode::Encoded => "Encoded",
        };
        (n_sequences, n_collections_loaded, mode_str)
    }

    /// Extended statistics including stub/loaded breakdown for collections
    pub fn stats_extended(&self) -> StoreStats {
        let n_sequences = self.sequence_store.len();
        let n_sequences_loaded = self.sequence_store.values()
            .filter(|record| record.is_loaded())
            .count();
        let n_collections = self.collections.len();
        let n_collections_loaded = self.collections.values()
            .filter(|record| record.has_sequences())
            .count();
        let mode_str = match self.mode {
            StorageMode::Raw => "Raw",
            StorageMode::Encoded => "Encoded",
        };
        let total_disk_size = self.actual_disk_usage();
        StoreStats {
            n_sequences,
            n_sequences_loaded,
            n_collections,
            n_collections_loaded,
            storage_mode: mode_str.to_string(),
            total_disk_size,
        }
    }
}

/// Extended statistics for a RefgetStore
#[derive(Debug, Clone)]
pub struct StoreStats {
    /// Total number of sequences (Stub + Full)
    pub n_sequences: usize,
    /// Number of sequences with data loaded (Full)
    pub n_sequences_loaded: usize,
    /// Total number of collections (Stub + Full)
    pub n_collections: usize,
    /// Number of collections with sequences loaded (Full)
    pub n_collections_loaded: usize,
    /// Storage mode (Raw or Encoded)
    pub storage_mode: String,
    /// Actual disk usage in bytes (all files in store directory)
    pub total_disk_size: usize,
}

/// Format bytes into human-readable size (KB, MB, GB, etc.)
fn format_bytes(bytes: usize) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    let mut size = bytes as f64;
    let mut unit_idx = 0;

    while size >= 1024.0 && unit_idx < UNITS.len() - 1 {
        size /= 1024.0;
        unit_idx += 1;
    }

    if unit_idx == 0 {
        format!("{} {}", bytes, UNITS[0])
    } else {
        format!("{:.2} {}", size, UNITS[unit_idx])
    }
}

impl Display for RefgetStore {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let total_size = self.total_disk_size();
        let size_str = format_bytes(total_size);
        writeln!(f, "SeqColStore object:")?;
        writeln!(f, "  Mode: {:?}", self.mode)?;
        writeln!(f, "  Disk size: {} ({} bytes)", size_str, total_size)?;
        writeln!(f, ">Sequences (n={}):", self.sequence_store.len())?;
        // Print out the sequences in the store
        for (i, (sha512_digest, sequence_record)) in self.sequence_store.iter().take(10).enumerate()
        {
            let metadata = sequence_record.metadata();
            let first_8_chars = match sequence_record {
                SequenceRecord::Stub(_) => "<stub>".to_string(),
                SequenceRecord::Full { metadata, sequence: seq } => {
                    // Extract the first 8 characters of the sequence (or fewer if the sequence is shorter)
                    match self.mode {
                        StorageMode::Encoded => {
                            let alphabet = lookup_alphabet(&metadata.alphabet);
                            let decoded = decode_substring_from_bytes(
                                seq,
                                0,
                                8.min(metadata.length),
                                alphabet,
                            );
                            String::from_utf8(decoded).unwrap_or_else(|_| "???".to_string())
                        }
                        StorageMode::Raw => String::from_utf8(seq[0..8.min(seq.len())].to_vec())
                            .unwrap_or_else(|_| "???".to_string()),
                    }
                }
            };

            writeln!(
                f,
                "   - {}. {:02x?}, MD5: {:02x?}, Length: {}, Alphabet: {:?}, Start: {}",
                i + 1,
                std::str::from_utf8(sha512_digest).unwrap(),
                &metadata.md5,
                &metadata.length,
                &metadata.alphabet,
                first_8_chars
            )?;
        }
        writeln!(f, ">Collections (n={:?}):", self.name_lookup.len())?;
        // Print out the collections in the store
        for (i, (digest, name_map)) in self.name_lookup.iter().enumerate() {
            // Convert the digest to a hex string
            let seqcol_digest_str = String::from_utf8_lossy(digest);
            writeln!(
                f,
                "  {}. Collection Digest: {:02x?} ({} sequences)",
                i + 1,
                seqcol_digest_str,
                name_map.len()
            )?;
            // Only show first 5 sequences in each collection
            for (name, sha512_digest) in name_map.iter().take(5) {
                // Convert the sha512_digest to a hex string
                let sha512_str = String::from_utf8_lossy(sha512_digest);
                writeln!(f, "   - Name: {}, SHA512: {:02x?}", name, sha512_str)?;
            }
            if name_map.len() > 5 {
                writeln!(f, "   - ... and {} more", name_map.len() - 5)?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use std::time::Instant;
    use crate::collection::{
        SequenceCollection, SequenceCollectionMetadata, SequenceMetadata, SequenceRecord,
    };
    use crate::digest::{md5, sha512t24u};
    use tempfile::tempdir;

    // Note: FASTA→RGSI roundtrip testing is in fasta.rs::digests_fa_to_rgsi

    // Helper function to calculate actual digests for testing
    fn calculate_test_digests(sequence: &[u8]) -> (String, String) {
        (sha512t24u(sequence), md5(sequence))
    }

    /// Creates a test store with 3 sequences for export testing
    fn setup_export_test_store(temp_path: &std::path::Path) -> (RefgetStore, [u8; 32]) {
        let fasta_content = ">chr1\nATGCATGCATGC\n>chr2\nGGGGAAAA\n>chr3\nTTTTCCCC\n";
        let temp_fasta_path = temp_path.join("test.fa");
        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        let mut store = RefgetStore::in_memory();
        store.add_sequence_collection_from_fasta(&temp_fasta_path).unwrap();

        let collections: Vec<_> = store.collections.keys().cloned().collect();
        let collection_digest = collections[0];

        (store, collection_digest)
    }

    #[test]
    fn test_mode_basics() {
        // Test default mode and convenience methods (no sequences needed)
        let mut store = RefgetStore::in_memory();

        // Default is Encoded
        assert_eq!(store.mode, StorageMode::Encoded);

        // Convenience methods
        store.disable_encoding();
        assert_eq!(store.mode, StorageMode::Raw);
        store.enable_encoding();
        assert_eq!(store.mode, StorageMode::Encoded);

        // set_encoding_mode() also works
        store.set_encoding_mode(StorageMode::Raw);
        assert_eq!(store.mode, StorageMode::Raw);
        store.set_encoding_mode(StorageMode::Encoded);
        assert_eq!(store.mode, StorageMode::Encoded);
    }

    #[test]
    fn test_mode_switching() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();
        let fasta_content = ">chr1\nATGCATGCATGC\n>chr2\nGGGGAAAA\n";
        let temp_fasta_path = temp_path.join("test.fa");
        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        let (chr1_sha, _) = calculate_test_digests(b"ATGCATGCATGC");
        let chr1_key = chr1_sha.as_bytes().to_key();

        // Test Raw -> Encoded
        {
            let mut store = RefgetStore::in_memory();
            store.disable_encoding();
            store.add_sequence_collection_from_fasta(&temp_fasta_path).unwrap();

            // Verify raw (ASCII)
            if let Some(SequenceRecord::Full { sequence, .. }) = store.sequence_store.get(&chr1_key) {
                assert_eq!(sequence, b"ATGCATGCATGC");
            }
            let seq_before = store.get_sequence(&chr1_sha).unwrap().decode().unwrap();

            // Switch to encoded
            store.set_encoding_mode(StorageMode::Encoded);

            // Verify encoded (smaller)
            if let Some(SequenceRecord::Full { sequence, .. }) = store.sequence_store.get(&chr1_key) {
                assert_eq!(sequence.len(), 3); // 12 bases * 2 bits = 3 bytes
            }
            let seq_after = store.get_sequence(&chr1_sha).unwrap().decode().unwrap();
            assert_eq!(seq_before, seq_after);
        }

        // Test Encoded -> Raw
        {
            let mut store = RefgetStore::in_memory();
            store.add_sequence_collection_from_fasta(&temp_fasta_path).unwrap();

            // Verify encoded
            if let Some(SequenceRecord::Full { sequence, .. }) = store.sequence_store.get(&chr1_key) {
                assert_eq!(sequence.len(), 3);
            }
            let seq_before = store.get_sequence(&chr1_sha).unwrap().decode().unwrap();

            // Switch to raw
            store.disable_encoding();

            // Verify raw (full size)
            if let Some(SequenceRecord::Full { sequence, .. }) = store.sequence_store.get(&chr1_key) {
                assert_eq!(sequence, b"ATGCATGCATGC");
            }
            let seq_after = store.get_sequence(&chr1_sha).unwrap().decode().unwrap();
            assert_eq!(seq_before, seq_after);
        }
    }

    #[test]
    fn test_refget_store_retrieve_seq_and_vec() {
        // Create temporary directory for all test files
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // --- 1. Prepare Test FASTA Data ---
        let fasta_content = "\
>chr1
ATGCATGCATGC
>chr2
GGGGAAAA
";
        let temp_fasta_path = temp_path.join("test.fa");

        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        // --- 2. Initialize RefgetStore and import FASTA ---
        let mut store = RefgetStore::in_memory();
        store.add_sequence_collection_from_fasta(&temp_fasta_path).unwrap();

        let sequence_keys: Vec<[u8; 32]> = store.sequence_store.keys().cloned().collect();

        let _ = sequence_keys[0]; //ww1QMyfFm1f4qa3fRLqqJGafIeEuZR1V
        let _ = sequence_keys[1]; //OyXJErGtjgcIVSdobGkHE3sBdQ5faDTf
        let collection_digest_ref: &str = "uC_UorBNf3YUu1YIDainBhI94CedlNeH";

        // Calculate expected SHA512t24u and MD5 for test sequences
        let (chr1_sha, chr1_md5) = calculate_test_digests(b"ATGCATGCATGC");
        let (chr2_sha, chr2_md5) = calculate_test_digests(b"GGGGAAAA");
        println!("chr1 values: {}  {}", chr1_sha, chr1_md5);
        println!("chr2 values: {}  {}", chr2_sha, chr2_md5);

        // --- 3. Prepare Test BED Data ---
        // Use only valid entries for the success test
        let bed_content = "\
chr1\t0\t5
chr1\t8\t12
chr2\t0\t4
";
        let temp_bed_path = temp_path.join("test.bed");

        fs::write(&temp_bed_path, bed_content).expect("Failed to write test BED file");

        let temp_output_fa_path = temp_path.join("output.fa");

        store
            .export_fasta_from_regions(
                collection_digest_ref,
                temp_bed_path.to_str().unwrap(),
                temp_output_fa_path.to_str().unwrap(),
            )
            .expect("export_fasta_from_regions failed");

        // Read the output FASTA file and verify its content
        let output_fa_content =
            fs::read_to_string(&temp_output_fa_path).expect("Failed to read output FASTA file");

        // Expected output content (headers and sequences should match the logic of the function)
        let expected_fa_content = format!(
            ">chr1 12 dna2bit {} {}\nATGCAATGC\n>chr2 8 dna2bit {} {}\nGGGG\n",
            chr1_sha, chr1_md5, chr2_sha, chr2_md5
        );
        assert_eq!(
            output_fa_content.trim(),
            expected_fa_content.trim(),
            "Output FASTA file content mismatch"
        );
        println!("✓ export_fasta_from_regions test passed.");

        // --- Test substrings_from_regions iterator (returns iterator of RetrievedSequence) ---
        let vec_result: Vec<_> = store
            .substrings_from_regions(collection_digest_ref, temp_bed_path.to_str().unwrap())
            .expect("substrings_from_regions failed")
            .collect::<Result<Vec<_>, _>>()
            .expect("substrings_from_regions had errors");

        // Define the expected vector of RetrievedSequence structs
        let expected_vec = vec![
            RetrievedSequence {
                sequence: "ATGCA".to_string(),
                chrom_name: "chr1".to_string(),
                start: 0,
                end: 5,
            },
            RetrievedSequence {
                sequence: "ATGC".to_string(),
                chrom_name: "chr1".to_string(),
                start: 8,
                end: 12,
            },
            RetrievedSequence {
                sequence: "GGGG".to_string(),
                chrom_name: "chr2".to_string(),
                start: 0,
                end: 4,
            },
        ];

        // Assert that the returned vector matches the expected vector
        assert_eq!(
            vec_result, expected_vec,
            "Retrieved sequence vector mismatch"
        );
        println!("✓ substrings_from_regions test passed.");
    }

    #[test]
    fn test_global_refget_store() {
        let sequence = b"ACGT";
        let name = "test_seq";
        println!("Testing RefgetStore with sequence: {}", name);

        // Create a sequence collection
        let mut collection = SequenceCollection {
            metadata: SequenceCollectionMetadata {
                digest: "test_collection".to_string(),
                n_sequences: 0,
                names_digest: "test".to_string(),
                sequences_digest: "test".to_string(),
                lengths_digest: "test".to_string(),
                file_path: None,
            },
            sequences: Vec::new(),
        };

        // Create a sequence record
        let seq_metadata = SequenceMetadata {
            name: name.to_string(),
            description: None,
            length: sequence.len(),
            sha512t24u: sha512t24u(sequence),
            md5: md5(sequence),
            alphabet: AlphabetType::Dna2bit,
            fai: None,
        };

        let record = SequenceRecord::Full {
            metadata: seq_metadata.clone(),
            sequence: sequence.to_vec(),
        };

        collection.sequences.push(record);

        // Add the sequence to the store
        let mut store = RefgetStore::in_memory();
        store.add_sequence_collection(collection.clone()).unwrap();

        // Verify the store has the sequence
        assert!(!store.sequence_store.is_empty());

        // Test sequence lookup by collection+name (using string digest)
        let retrieved_by_name_str =
            store.get_sequence_by_name(&collection.metadata.digest, name);
        assert!(retrieved_by_name_str.is_ok());
        let retrieved_record = retrieved_by_name_str.unwrap();
        assert_eq!(retrieved_record.metadata().name, name);
        assert_eq!(retrieved_record.sequence().unwrap(), sequence);

        // Test sequence lookup by collection+name (using [u8; 32] digest)
        let retrieved_by_name_key =
            store.get_sequence_by_name(collection.metadata.digest.to_key(), name);
        assert!(retrieved_by_name_key.is_ok());
        let retrieved_record = retrieved_by_name_key.unwrap();
        assert_eq!(retrieved_record.metadata().name, name);
        assert_eq!(retrieved_record.sequence().unwrap(), sequence);

        // Test sequence lookup by SHA512 digest (using string)
        let retrieved_by_sha512_str = store.get_sequence(&seq_metadata.sha512t24u);
        assert!(retrieved_by_sha512_str.is_ok());
        let retrieved_record = retrieved_by_sha512_str.unwrap();
        assert_eq!(retrieved_record.metadata().name, name);
        assert_eq!(retrieved_record.sequence().unwrap(), sequence);

        // Test sequence lookup by SHA512 digest (using [u8; 32])
        let retrieved_by_sha512_key = store.get_sequence(seq_metadata.sha512t24u.to_key());
        assert!(retrieved_by_sha512_key.is_ok());
        let retrieved_record = retrieved_by_sha512_key.unwrap();
        assert_eq!(retrieved_record.metadata().name, name);
        assert_eq!(retrieved_record.sequence().unwrap(), sequence);
    }

    #[test]
    fn test_import_fasta() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Copy test FASTA file to temp directory
        let test_fa = "../tests/data/fasta/base.fa";
        let temp_fa = temp_path.join("base.fa");

        std::fs::copy(test_fa, &temp_fa).expect("Failed to copy test FASTA file");

        let mut store = RefgetStore::in_memory();

        // Import the FASTA file
        store.add_sequence_collection_from_fasta(temp_fa).unwrap();

        // Check that the store has sequences
        assert!(!store.sequence_store.is_empty());

        // Try writing to a file
        let seq_template = "sequences/%s2/%s.seq";
        store
            .write_store_to_dir(temp_path.to_str().unwrap(), Some(seq_template))
            .unwrap();
    }

    #[test]
    fn test_disk_persistence() {
        // Create a temporary directory for the test
        let temp_dir = tempdir().unwrap();
        let temp_path = temp_dir.path();
        let temp_fasta = temp_path.join("base.fa.gz");
        std::fs::copy("../tests/data/fasta/base.fa.gz", &temp_fasta)
            .expect("Failed to copy base.fa.gz to tempdir");

        // Create a new sequence store
        let mut store = RefgetStore::in_memory();

        // Import a FASTA file into the store
        // store.add_sequence_collection_from_fasta("../tests/data/subset.fa.gz").unwrap();
        store.add_sequence_collection_from_fasta(&temp_fasta).unwrap();

        // Get the sequence keys for verification (assuming we know the test file contains 3 sequences)
        let sequence_keys: Vec<[u8; 32]> = store.sequence_store.keys().cloned().collect();
        assert_eq!(
            sequence_keys.len(),
            3,
            "Test file should contain exactly 3 sequences"
        );

        let sha512_key1 = sequence_keys[0];
        let sha512_key2 = sequence_keys[1];

        // Store original sequences for comparison
        let original_seq1 = store.sequence_store.get(&sha512_key1).unwrap().clone();
        let original_seq2 = store.sequence_store.get(&sha512_key2).unwrap().clone();

        // Write the store to the temporary directory
        let seq_template = "sequences/%s2/%s.seq";
        store
            .write_store_to_dir(temp_path, Some(seq_template))
            .unwrap();

        // Verify that the files were created (using new names)
        assert!(temp_path.join("sequences").exists());
        assert!(temp_path.join("sequences").read_dir().unwrap().count() > 0);
        assert!(temp_path.join("rgstore.json").exists());
        assert!(temp_path.join("sequences.rgsi").exists());
        assert!(temp_path.join("collections.rgci").exists());
        assert!(temp_path.join("collections").exists());

        // Load the store from disk
        let mut loaded_store = RefgetStore::open_local(temp_path).unwrap();

        // Verify that the loaded store has the same sequences
        assert_eq!(loaded_store.sequence_store.len(), 3);

        // Verify that we can retrieve sequences by their keys
        assert!(loaded_store.sequence_store.contains_key(&sha512_key1));
        assert!(loaded_store.sequence_store.contains_key(&sha512_key2));

        // Verify the content of the sequences
        let loaded_seq1 = loaded_store.sequence_store.get(&sha512_key1).unwrap();
        let loaded_seq2 = loaded_store.sequence_store.get(&sha512_key2).unwrap();

        // Check metadata equality
        assert_eq!(original_seq1.metadata().name, loaded_seq1.metadata().name);
        assert_eq!(original_seq1.metadata().length, loaded_seq1.metadata().length);
        assert_eq!(
            original_seq1.metadata().sha512t24u,
            loaded_seq1.metadata().sha512t24u
        );
        assert_eq!(original_seq1.metadata().md5, loaded_seq1.metadata().md5);

        assert_eq!(original_seq2.metadata().name, loaded_seq2.metadata().name);
        assert_eq!(original_seq2.metadata().length, loaded_seq2.metadata().length);
        assert_eq!(
            original_seq2.metadata().sha512t24u,
            loaded_seq2.metadata().sha512t24u
        );
        assert_eq!(original_seq2.metadata().md5, loaded_seq2.metadata().md5);

        // Check data is not loaded initially (lazy loading)
        assert_eq!(loaded_seq1.is_loaded(), false, "Data should not be loaded initially with lazy loading");
        assert_eq!(loaded_seq2.is_loaded(), false, "Data should not be loaded initially with lazy loading");

        // Verify MD5 lookup is preserved
        assert_eq!(loaded_store.md5_lookup.len(), 3);

        // Verify collections are preserved
        assert_eq!(loaded_store.collections.len(), store.collections.len());

        // Test sequence retrieval functionality
        for (digest, original_record) in &store.sequence_store {
            let loaded_record = loaded_store.get_sequence(*digest).unwrap();
            assert_eq!(original_record.metadata().name, loaded_record.metadata().name);
            assert_eq!(
                original_record.metadata().length,
                loaded_record.metadata().length
            );

            // Test substring retrieval works on loaded store
            if original_record.metadata().length > 0 {
                let substring_len = std::cmp::min(5, original_record.metadata().length);
                let substring = loaded_store.get_substring(digest, 0, substring_len);
                assert!(
                    substring.is_ok(),
                    "Should be able to retrieve substring from loaded sequence"
                );
            }
        }

        println!("✓ Disk persistence test passed - all data preserved correctly");
    }

    #[test]
    fn test_export_fasta_all_sequences() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let (mut store, collection_digest) = setup_export_test_store(temp_dir.path());

        let output_path = temp_dir.path().join("exported_all.fa");
        store.export_fasta(&collection_digest, &output_path, None, Some(80)).unwrap();

        let exported = fs::read_to_string(&output_path).unwrap();
        assert!(exported.contains(">chr1") && exported.contains(">chr2") && exported.contains(">chr3"));
        assert!(exported.contains("ATGCATGCATGC") && exported.contains("GGGGAAAA") && exported.contains("TTTTCCCC"));
    }

    #[test]
    fn test_export_fasta_subset_sequences() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let (mut store, collection_digest) = setup_export_test_store(temp_dir.path());

        let output_path = temp_dir.path().join("exported_subset.fa");
        store.export_fasta(&collection_digest, &output_path, Some(vec!["chr1", "chr3"]), Some(80)).unwrap();

        let exported = fs::read_to_string(&output_path).unwrap();
        assert!(exported.contains(">chr1") && exported.contains(">chr3"));
        assert!(!exported.contains(">chr2") && !exported.contains("GGGGAAAA"));
    }

    #[test]
    fn test_export_fasta_roundtrip() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Create test FASTA with longer sequences
        let fasta_content = "\
>seq1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>seq2
GGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCC
";
        let temp_fasta_path = temp_path.join("original.fa");
        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        // Import into store
        let mut store1 = RefgetStore::in_memory();
        store1.add_sequence_collection_from_fasta(&temp_fasta_path).unwrap();

        // Get original digests
        let original_digests: Vec<String> = store1
            .sequence_store
            .values()
            .map(|r| r.metadata().sha512t24u.clone())
            .collect();

        // Export to new FASTA
        let collections: Vec<_> = store1.collections.keys().cloned().collect();
        let collection_digest = collections[0];
        let exported_path = temp_path.join("exported.fa");
        store1
            .export_fasta(&collection_digest, &exported_path, None, Some(60))
            .expect("Failed to export FASTA");

        // Re-import the exported FASTA
        let mut store2 = RefgetStore::in_memory();
        store2.add_sequence_collection_from_fasta(&exported_path).unwrap();

        // Verify digests match (same sequences)
        let new_digests: Vec<String> = store2
            .sequence_store
            .values()
            .map(|r| r.metadata().sha512t24u.clone())
            .collect();

        assert_eq!(original_digests.len(), new_digests.len(), "Should have same number of sequences");
        for digest in original_digests {
            assert!(
                new_digests.contains(&digest),
                "Digest {} should be present after roundtrip",
                digest
            );
        }

        println!("✓ Export/import roundtrip test passed");
    }

    #[test]
    fn test_export_fasta_by_digests() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let (mut store, _) = setup_export_test_store(temp_dir.path());

        let digests: Vec<String> = store.sequence_store.values()
            .map(|r| r.metadata().sha512t24u.clone()).collect();
        let digest_refs: Vec<&str> = digests.iter().map(|s| s.as_str()).collect();

        let output_path = temp_dir.path().join("exported_by_digests.fa");
        store.export_fasta_by_digests(digest_refs, &output_path, Some(80)).unwrap();

        let exported = fs::read_to_string(&output_path).unwrap();
        assert!(exported.contains(">chr1") && exported.contains(">chr2") && exported.contains(">chr3"));
    }

    #[test]
    fn test_export_fasta_error_handling() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let (mut store, collection_digest) = setup_export_test_store(temp_dir.path());

        let output_path = temp_dir.path().join("should_fail.fa");

        // Test with non-existent collection
        let fake_collection = b"fake_collection_digest_12345678";
        assert!(store.export_fasta(fake_collection, &output_path, None, Some(80)).is_err());

        // Test with non-existent sequence name
        assert!(store.export_fasta(&collection_digest, &output_path, Some(vec!["nonexistent_chr"]), Some(80)).is_err());
    }

    #[test]
    fn test_export_fasta_after_load_local() {
        // Test that export_fasta works on disk-loaded stores
        // This verifies the lazy loading fix (ensure_collection_loaded is called)
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();
        let store_path = temp_path.join("store");

        // Create test FASTA
        let fasta_content = ">chr1\nACGTACGT\n>chr2\nGGGGAAAA\n";
        let fasta_path = temp_path.join("test.fa");
        fs::write(&fasta_path, fasta_content).unwrap();

        // Create and populate store on disk, save digest before closing
        let collection_digest: [u8; 32];
        {
            let mut store = RefgetStore::on_disk(&store_path).unwrap();
            store.add_sequence_collection_from_fasta(&fasta_path).unwrap();
            let collections: Vec<_> = store.collections.keys().cloned().collect();
            assert_eq!(collections.len(), 1, "Should have exactly one collection");
            collection_digest = collections[0];
        } // store dropped here

        // Load the store fresh from disk
        let mut loaded_store = RefgetStore::open_local(&store_path).unwrap();

        // Verify collection is initially a stub (lazy loaded)
        assert!(
            !loaded_store.is_collection_loaded(&collection_digest),
            "Collection should be Stub after loading from disk"
        );

        // This should work (was failing before fix due to missing ensure_collection_loaded call)
        let output_path = temp_path.join("exported.fa");
        loaded_store
            .export_fasta(&collection_digest, &output_path, None, Some(80))
            .expect("export_fasta should work on disk-loaded stores");

        // Verify output
        let exported = fs::read_to_string(&output_path).unwrap();
        assert!(exported.contains(">chr1"));
        assert!(exported.contains("ACGTACGT"));
        assert!(exported.contains(">chr2"));
        assert!(exported.contains("GGGGAAAA"));

        println!("✓ export_fasta after load_local test passed");
    }

    #[test]
    fn test_sequence_names_with_spaces() {
        // Test FASTA header parsing: name is first word, rest is description
        // Following FASTA standard, we now split on whitespace
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Create test FASTA with headers containing descriptions after the ID
        // This mimics the structure from HPRC pangenome files
        let fasta_content = "\
>JAHKSE010000016.1 unmasked:primary_assembly HG002.alt.pat.f1_v2:JAHKSE010000016.1:1:100:1
ATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGC
>JAHKSE010000012.1 unmasked:primary_assembly HG002.alt.pat.f1_v2:JAHKSE010000012.1:1:100:1
GGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGG
GGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGG
";
        let temp_fasta_path = temp_path.join("spaces_in_names.fa");
        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        // Import FASTA - headers will be split into name (first word) and description (rest)
        let mut store = RefgetStore::in_memory();
        store.add_sequence_collection_from_fasta(&temp_fasta_path)
            .expect("Should parse FASTA headers correctly");

        // Verify the sequences were loaded
        assert_eq!(store.sequence_store.len(), 2);

        // Names are now just the first word (before whitespace)
        let name1 = "JAHKSE010000016.1";
        let name2 = "JAHKSE010000012.1";

        // Get the collection
        let collections: Vec<_> = store.collections.keys().cloned().collect();
        assert_eq!(collections.len(), 1);
        let collection_digest = collections[0];

        // Verify we can retrieve sequences by their short names (first word only)
        // and check the description was captured
        {
            let seq1 = store.get_sequence_by_name(&collection_digest, name1);
            assert!(seq1.is_ok(), "Should retrieve sequence by name (first word)");

            let seq1_meta = seq1.unwrap().metadata();
            assert_eq!(seq1_meta.name, "JAHKSE010000016.1");
            assert_eq!(
                seq1_meta.description,
                Some("unmasked:primary_assembly HG002.alt.pat.f1_v2:JAHKSE010000016.1:1:100:1".to_string())
            );
        }

        {
            let seq2 = store.get_sequence_by_name(&collection_digest, name2);
            assert!(seq2.is_ok(), "Should retrieve sequence by name (first word)");
        }

        println!("✓ FASTA header parsing test passed");
    }

    #[test]
    fn test_rgsi_filename_with_dots() {
        // Test that RGSI filenames preserve dots in the base name
        // Real HPRC files like "HG002.alt.pat.f1_v2.unmasked.fa.gz"
        // should create "HG002.alt.pat.f1_v2.unmasked.rgsi", NOT "HG002.rgsi"

        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Copy test file to temp (so .rgsi file gets created there, not in test data)
        let test_file = "../tests/data/fasta/HG002.alt.pat.f1_v2.unmasked.fa";
        let temp_fasta = temp_path.join("HG002.alt.pat.f1_v2.unmasked.fa");
        fs::copy(test_file, &temp_fasta).expect("Failed to copy test file");

        // Load the FASTA - this creates a .rgsi file
        let mut store = RefgetStore::in_memory();
        store.add_sequence_collection_from_fasta(&temp_fasta)
            .expect("Should load FASTA");

        // Check which .rgsi file was created
        let correct_rgsi = temp_path.join("HG002.alt.pat.f1_v2.unmasked.rgsi");
        let wrong_rgsi = temp_path.join("HG002.rgsi");

        let files: Vec<_> = std::fs::read_dir(temp_path).unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
            .collect();

        assert!(
            correct_rgsi.exists(),
            "Expected 'HG002.alt.pat.f1_v2.unmasked.rgsi' but found: {:?}",
            files
        );

        assert!(
            !wrong_rgsi.exists(),
            "Should NOT create 'HG002.rgsi' (strips too many dots)"
        );

        println!("✓ RGSI filename with dots test passed");
    }

    #[test]
    fn test_on_disk_collection_written_incrementally() {
        // Test that collection RGSI files are written to disk immediately
        // when using on_disk() store, not just when write_store_to_dir() is called
        let temp_dir = tempdir().unwrap();
        let temp_path = temp_dir.path();
        let temp_fasta = temp_path.join("base.fa.gz");
        std::fs::copy("../tests/data/fasta/base.fa.gz", &temp_fasta)
            .expect("Failed to copy base.fa.gz to tempdir");

        let cache_path = temp_path.join("cache");
        let mut store = RefgetStore::on_disk(&cache_path).unwrap();

        // Load FASTA file into the store
        store.add_sequence_collection_from_fasta(&temp_fasta).unwrap();

        // BEFORE calling write_store_to_dir, verify collection RGSI files exist
        let collections_dir = cache_path.join("collections");
        assert!(collections_dir.exists(), "Collections directory should exist");

        let rgsi_files: Vec<_> = std::fs::read_dir(&collections_dir)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
            .collect();

        assert!(!rgsi_files.is_empty(), "Collection RGSI files should be written incrementally, found: {:?}", rgsi_files);
        assert!(rgsi_files.iter().any(|f| f.ends_with(".rgsi")), "Should have .rgsi files");

        println!("✓ On-disk collection incremental write test passed");
    }

    #[test]
    fn test_disk_size_calculation() {
        let mut store = RefgetStore::in_memory();
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz").unwrap();

        let disk_size = store.total_disk_size();
        assert!(disk_size > 0, "Disk size should be greater than 0");

        // Verify against manual calculation
        let manual: usize = store.list_sequences()
            .iter()
            .map(|m| (m.length * m.alphabet.bits_per_symbol()).div_ceil(8))
            .sum();
        assert_eq!(disk_size, manual);
    }

    #[test]
    fn test_incremental_index_writing() {
        let temp_dir = tempdir().unwrap();
        let cache_path = temp_dir.path().join("store");
        let mut store = RefgetStore::on_disk(&cache_path).unwrap();

        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz").unwrap();

        // Index files should exist immediately (using new names)
        assert!(cache_path.join("rgstore.json").exists(), "rgstore.json should exist");
        assert!(cache_path.join("sequences.rgsi").exists(), "sequences.rgsi should exist");
        assert!(cache_path.join("collections.rgci").exists(), "collections.rgci should exist");

        // Store should be loadable (mode ignored for existing store)
        let _loaded = RefgetStore::on_disk(&cache_path).unwrap();
    }

    #[test]
    fn test_write_method() {
        let temp_dir = tempdir().unwrap();
        let cache_path = temp_dir.path().join("store");
        let mut store = RefgetStore::on_disk(&cache_path).unwrap();

        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz").unwrap();
        store.write().unwrap();  // Should succeed

        assert!(cache_path.join("rgstore.json").exists());
    }

    #[test]
    fn test_on_disk_smart_constructor() {
        let temp_dir = tempdir().unwrap();
        let cache_path = temp_dir.path().join("store");

        // Create new store (defaults to Encoded mode)
        let mut store1 = RefgetStore::on_disk(&cache_path).unwrap();
        assert_eq!(store1.mode, StorageMode::Encoded);
        store1.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz").unwrap();

        // Load existing store - should preserve Encoded mode
        let store2 = RefgetStore::on_disk(&cache_path).unwrap();
        assert_eq!(store2.sequence_store.len(), store1.sequence_store.len());
        assert_eq!(store2.mode, StorageMode::Encoded, "Loaded store should preserve Encoded mode");

        // Test with Raw mode
        let cache_path_raw = temp_dir.path().join("store_raw");
        let mut store3 = RefgetStore::on_disk(&cache_path_raw).unwrap();
        store3.disable_encoding(); // Switch to Raw
        assert_eq!(store3.mode, StorageMode::Raw);
        store3.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz").unwrap();

        // Load and verify Raw mode is persisted
        let store4 = RefgetStore::on_disk(&cache_path_raw).unwrap();
        assert_eq!(store4.mode, StorageMode::Raw, "Loaded store should preserve Raw mode");

        // Verify rgstore.json contains the mode
        let index_path = cache_path_raw.join("rgstore.json");
        let json = fs::read_to_string(&index_path).unwrap();
        assert!(json.contains("\"mode\":\"Raw\"") || json.contains("\"mode\": \"Raw\""),
                "rgstore.json should contain mode: Raw");
    }

    #[test]
    fn test_collection_metadata_methods() {
        // Test list_collections, get_collection_metadata, is_collection_loaded
        let temp_dir = tempdir().unwrap();
        let cache_path = temp_dir.path().join("store");
        let mut store = RefgetStore::on_disk(&cache_path).unwrap();

        // Add a FASTA file
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz").unwrap();

        // Test list_collections
        let collections = store.list_collections();
        assert_eq!(collections.len(), 1, "Should have 1 collection");
        let digest = collections[0].digest.clone();

        // Test get_collection_metadata
        let meta = store.get_collection_metadata(&digest);
        assert!(meta.is_some(), "Should get collection metadata");
        let meta = meta.unwrap();
        assert_eq!(meta.n_sequences, 3, "Collection should have 3 sequences");

        // Test is_collection_loaded - should be true since we just added it
        assert!(store.is_collection_loaded(&digest), "Collection should be loaded (Full)");

        // Test stats_extended returns collection counts
        let stats = store.stats_extended();
        assert_eq!(stats.n_collections, 1, "Should have 1 collection total");
        assert_eq!(stats.n_collections_loaded, 1, "Should have 1 collection loaded");
        assert_eq!(stats.n_sequences, 3, "Should have 3 sequences");

        println!("✓ Collection metadata methods test passed");
    }

    #[test]
    fn test_collection_stub_lazy_loading() {
        // Test that collections load as Stubs and upgrade to Full on-demand
        let temp_dir = tempdir().unwrap();
        let cache_path = temp_dir.path().join("store");

        // Create and populate the store
        let mut store = RefgetStore::on_disk(&cache_path).unwrap();
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz").unwrap();
        let digest = store.list_collections()[0].digest.clone();

        // Drop the store and reload from disk
        drop(store);
        let mut loaded_store = RefgetStore::open_local(&cache_path).unwrap();

        // VERIFY: Metadata is available (from collections.rgci)
        let meta = loaded_store.get_collection_metadata(&digest);
        assert!(meta.is_some(), "Metadata should be available for Stub");
        assert_eq!(meta.unwrap().n_sequences, 3, "Stub should know sequence count");

        // VERIFY: Collection is a Stub (not loaded into memory)
        assert!(!loaded_store.is_collection_loaded(&digest),
            "Collection should be Stub after loading from disk");

        // VERIFY: stats shows 0 collections loaded
        let stats_before = loaded_store.stats_extended();
        assert_eq!(stats_before.n_collections, 1, "Should have 1 collection total");
        assert_eq!(stats_before.n_collections_loaded, 0, "Should have 0 collections loaded initially");

        // TRIGGER: Access a sequence by name - this should trigger lazy loading
        let seq = loaded_store.get_sequence_by_name(&digest, "chr1");
        assert!(seq.is_ok(), "Should be able to retrieve sequence after lazy load");
        assert_eq!(seq.unwrap().metadata().name, "chr1");

        // VERIFY: Collection is now Full (loaded into memory)
        assert!(loaded_store.is_collection_loaded(&digest),
            "Collection should be Full after accessing a sequence");

        // VERIFY: stats now shows 1 collection loaded
        let stats_after = loaded_store.stats_extended();
        assert_eq!(stats_after.n_collections_loaded, 1, "Should have 1 collection loaded after access");

        println!("✓ Collection stub lazy loading test passed");
    }

    // Note: open_local is tested in test_disk_persistence which is more comprehensive

    #[test]
    fn test_get_collection() {
        // Test the get_collection method (returns collection with sequence metadata, lazy loading)
        let temp_dir = tempdir().unwrap();
        let cache_path = temp_dir.path().join("store");

        // Create and populate the store
        let mut store = RefgetStore::on_disk(&cache_path).unwrap();
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz").unwrap();
        let digest = store.list_collections()[0].digest.clone();
        drop(store);

        // Reload and test get_collection
        let mut loaded_store = RefgetStore::open_local(&cache_path).unwrap();

        // Before loading - collection is a Stub
        assert!(!loaded_store.is_collection_loaded(&digest));

        // Before loading - no sequences loaded
        let stats_before = loaded_store.stats_extended();
        assert_eq!(stats_before.n_sequences_loaded, 0, "No sequences should be loaded initially");

        // Get the collection (loads metadata only, sequences are lazy)
        let collection = loaded_store.get_collection(&digest).unwrap();
        assert!(!collection.sequences.is_empty(), "Collection should have sequences");
        assert_eq!(collection.sequences.len(), 3);

        // After get_collection - collection is loaded but sequences are still stubs (lazy loading)
        let stats_after = loaded_store.stats_extended();
        assert_eq!(stats_after.n_sequences_loaded, 0, "Sequences not loaded until explicitly fetched");
        assert_eq!(stats_after.n_collections_loaded, 1, "Collection should be loaded");

        // Verify sequences are stubs (not loaded)
        for record in loaded_store.sequence_store.values() {
            assert!(!record.is_loaded(), "Sequences should be stubs after get_collection");
        }

        // Now explicitly load a sequence
        let seq_digest = collection.sequences[0].metadata().sha512t24u.clone();
        let loaded_seq = loaded_store.get_sequence(&seq_digest).unwrap();
        assert!(loaded_seq.is_loaded(), "Sequence should be loaded after get_sequence");

        println!("✓ get_collection test passed");
    }

    #[test]
    fn test_get_sequence() {
        // Test the get_sequence method (loads sequence on demand)
        let temp_dir = tempdir().unwrap();
        let cache_path = temp_dir.path().join("store");

        // Create and populate the store
        let mut store = RefgetStore::on_disk(&cache_path).unwrap();
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz").unwrap();

        // Get a sequence digest from the sequence_store
        let seq_digest = store.sequence_store.values().next().unwrap().metadata().sha512t24u.clone();
        drop(store);

        // Reload and test get_sequence
        let mut loaded_store = RefgetStore::open_local(&cache_path).unwrap();

        // Before loading - sequence is a Stub (no data)
        let seq_before = loaded_store.sequence_store.get(&seq_digest.to_key()).unwrap();
        assert!(!seq_before.is_loaded(), "Sequence should not have data before get_sequence");

        // Get the sequence (loads data on demand)
        let loaded_seq = loaded_store.get_sequence(&seq_digest).unwrap();
        assert!(loaded_seq.is_loaded(), "Sequence should have data after get_sequence");
        assert!(loaded_seq.sequence().is_some(), "Sequence data should be available");

        println!("✓ get_sequence test passed");
    }

    #[test]
    fn test_get_collection_idempotent() {
        // Test that calling get_collection twice is safe (idempotent)
        let temp_dir = tempdir().unwrap();
        let cache_path = temp_dir.path().join("store");

        // Create and populate the store
        let mut store = RefgetStore::on_disk(&cache_path).unwrap();
        store.add_sequence_collection_from_fasta("../tests/data/fasta/base.fa.gz").unwrap();
        let digest = store.list_collections()[0].digest.clone();
        drop(store);

        // Reload and test idempotent loading
        let mut loaded_store = RefgetStore::open_local(&cache_path).unwrap();

        // Get twice - both should succeed
        let result1 = loaded_store.get_collection(&digest);
        assert!(result1.is_ok(), "First get should succeed");

        let result2 = loaded_store.get_collection(&digest);
        assert!(result2.is_ok(), "Second get should also succeed");

        // Store state should be unchanged after second get
        assert_eq!(loaded_store.stats_extended().n_collections_loaded, 1);

        println!("✓ get_collection idempotent test passed");
    }

    #[test]
    fn test_sanitize_relative_path_rejects_traversal() {
        assert!(RefgetStore::sanitize_relative_path("../etc/passwd").is_err());
        assert!(RefgetStore::sanitize_relative_path("foo/../bar").is_err());
        assert!(RefgetStore::sanitize_relative_path("foo/../../bar").is_err());
        assert!(RefgetStore::sanitize_relative_path("..").is_err());
    }

    #[test]
    fn test_sanitize_relative_path_rejects_absolute() {
        assert!(RefgetStore::sanitize_relative_path("/etc/passwd").is_err());
        assert!(RefgetStore::sanitize_relative_path("\\windows\\system32").is_err());
    }

    #[test]
    fn test_sanitize_relative_path_accepts_valid() {
        assert!(RefgetStore::sanitize_relative_path("sequences/ab/abc123.seq").is_ok());
        assert!(RefgetStore::sanitize_relative_path("collections/xyz.rgsi").is_ok());
        assert!(RefgetStore::sanitize_relative_path("rgstore.json").is_ok());
        assert!(RefgetStore::sanitize_relative_path("sequences/%s2/%s.seq").is_ok());
    }

    #[test]
    fn test_stale_rgsi_cache_is_ignored() {
        // Reproduces issue where empty/stale .rgsi cache causes
        // "Sequence not found in metadata. Available (0 total): []"
        use std::io::Write;

        let temp_dir = tempdir().unwrap();

        // Create a test FASTA file
        let fasta_path = temp_dir.path().join("test.fa");
        let mut fasta_file = fs::File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">chr1\nATGCATGC\n>chr2\nGGGGAAAA").unwrap();

        // Create an EMPTY .rgsi cache file (simulating stale/corrupt cache)
        let rgsi_path = temp_dir.path().join("test.rgsi");
        let mut rgsi_file = fs::File::create(&rgsi_path).unwrap();
        writeln!(rgsi_file, "#name\tlength\talphabet\tsha512t24u\tmd5\tdescription").unwrap();

        // Create on-disk store
        let store_path = temp_dir.path().join("store");
        let mut store = RefgetStore::on_disk(&store_path).unwrap();

        // Before fix: Failed with "Sequence 'chr1' not found in metadata. Available (0 total): []"
        // After fix: Detects empty cache, deletes it, re-digests FASTA
        let result = store.add_sequence_collection_from_fasta(&fasta_path);
        assert!(result.is_ok(), "Should handle stale cache: {:?}", result.err());

        // Verify sequences were loaded
        assert_eq!(store.sequence_store.len(), 2, "Should have 2 sequences");

        println!("✓ Stale RGSI cache test passed");
    }
}
