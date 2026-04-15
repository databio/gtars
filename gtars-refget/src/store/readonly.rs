//! ReadonlyRefgetStore struct definition and core methods.

use super::*;
use super::alias::AliasManager;

use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::path::{Path, PathBuf};

use indexmap::IndexMap;

use anyhow::{anyhow, Context, Result};

use crate::collection::{read_rgsi_file, SequenceMetadataExt, SequenceRecordExt};
use crate::digest::lookup_alphabet;
use crate::digest::{
    SequenceCollectionMetadata, SequenceCollectionRecord, SequenceMetadata,
    SequenceRecord,
};
use crate::digest::{decode_string_from_bytes, decode_substring_from_bytes, encode_sequence, StreamingDecoder};
use crate::hashkeyable::{DigestKey, HashKeyable, key_to_digest_string};
use crate::seqcol::metadata_matches_attribute;

use std::fs::{self, create_dir_all, File};
use std::io::{BufReader, Read, Seek, SeekFrom};

/// Open a bounded HTTP byte range against a remote refget store and return
/// a `Read` yielding exactly `byte_end - byte_start` bytes.
#[cfg(feature = "http")]
fn open_remote_range(
    remote: &str,
    relpath: &Path,
    byte_start: u64,
    byte_end: u64,
) -> Result<Box<dyn Read + Send>> {
    let relpath_str = relpath
        .to_str()
        .ok_or_else(|| anyhow!("Non-UTF8 sequence path"))?;
    let url = if remote.ends_with('/') {
        format!("{}{}", remote, relpath_str)
    } else {
        format!("{}/{}", remote, relpath_str)
    };

    let byte_len = byte_end - byte_start;
    let last_byte = byte_end.saturating_sub(1);
    let response = ureq::get(&url)
        .set("Range", &format!("bytes={}-{}", byte_start, last_byte))
        .call()
        .map_err(|e| anyhow!("Failed to open remote byte range: {}", e))?;

    Ok(Box::new(response.into_reader().take(byte_len)))
}

#[cfg(not(feature = "http"))]
fn open_remote_range(
    _remote: &str,
    _relpath: &Path,
    _byte_start: u64,
    _byte_end: u64,
) -> Result<Box<dyn Read + Send>> {
    Err(anyhow!(
        "Remote streaming requires the 'http' feature to be enabled"
    ))
}

/// Core refget store with `&self` read methods, suitable for `Arc` sharing in servers.
///
/// Mutating methods are used during the setup/loading phase; once wrapped in `Arc`,
/// only `&self` reads are accessible, making concurrent access thread-safe.
///
/// Holds a global sequence_store with all sequences (across collections) deduplicated.
/// This allows lookup by sequence digest directly (bypassing collection information).
/// Also holds a collections hashmap, to provide lookup by collection+name.
#[derive(Debug)]
pub struct ReadonlyRefgetStore {
    /// SHA512t24u digest -> SequenceRecord (metadata + optional data)
    pub(crate) sequence_store: HashMap<DigestKey, SequenceRecord>,
    /// MD5 digest -> SHA512t24u digest lookup
    pub(crate) md5_lookup: HashMap<DigestKey, DigestKey>,

    /// Collection digest -> {name -> SHA512t24u digest} (IndexMap preserves FASTA insertion order)
    pub(crate) name_lookup: HashMap<DigestKey, IndexMap<String, DigestKey>>,
    /// Active sequence collections (now using SequenceCollectionRecord for Stub/Full pattern)
    pub(crate) collections: HashMap<DigestKey, SequenceCollectionRecord>,
    /// Storage strategy for sequences
    pub(crate) mode: StorageMode,
    /// Where the store lives on disk (local store or cache directory)
    pub(crate) local_path: Option<PathBuf>,
    /// Where to pull sequences from (if remote-backed)
    pub(crate) remote_source: Option<String>,
    /// Template for sequence file paths (e.g., "sequences/%s2/%s.seq")
    pub(crate) seqdata_path_template: Option<String>,
    /// Whether to persist sequences to disk (write-through caching)
    pub(crate) persist_to_disk: bool,
    /// Whether to suppress progress output
    pub(crate) quiet: bool,
    /// Whether to compute ancillary digests (nlp, snlp, sorted_sequences).
    /// Default: true for new stores.
    pub(crate) ancillary_digests: bool,
    /// Whether on-disk attribute reverse index is enabled.
    /// Default: false. Part 2 implements the indexed path.
    pub(crate) attribute_index: bool,
    /// Human-readable aliases for sequences and collections.
    pub(crate) aliases: AliasManager,
    /// FHR metadata for collections, keyed by collection digest.
    pub(crate) fhr_metadata: HashMap<DigestKey, super::fhr_metadata::FhrMetadata>,
    /// Available sequence alias namespaces (from manifest, for remote discovery).
    pub(crate) available_sequence_alias_namespaces: Vec<String>,
    /// Available collection alias namespaces (from manifest, for remote discovery).
    pub(crate) available_collection_alias_namespaces: Vec<String>,
    /// Cache of decoded sequence bytes, keyed by SHA512t24u digest.
    /// Populated by ensure_decoded(), read by sequence_bytes().
    pub(crate) decoded_cache: HashMap<DigestKey, Vec<u8>>,
}

impl ReadonlyRefgetStore {
    /// Generic constructor. Creates a new, empty `ReadonlyRefgetStore`.
    /// Internal only - users should go through RefgetStore.
    pub(crate) fn new(mode: StorageMode) -> Self {
        ReadonlyRefgetStore {
            sequence_store: HashMap::new(),
            md5_lookup: HashMap::new(),
            name_lookup: HashMap::new(),
            collections: HashMap::new(),
            mode,
            local_path: None,
            remote_source: None,
            seqdata_path_template: None,
            persist_to_disk: false,
            quiet: false,
            ancillary_digests: true,
            attribute_index: false,
            aliases: AliasManager::default(),
            fhr_metadata: HashMap::new(),
            decoded_cache: HashMap::new(),
            available_sequence_alias_namespaces: Vec::new(),
            available_collection_alias_namespaces: Vec::new(),
        }
    }

    /// Set whether to suppress progress output.
    pub fn set_quiet(&mut self, quiet: bool) {
        self.quiet = quiet;
    }

    /// Returns whether the store is in quiet mode.
    pub fn is_quiet(&self) -> bool {
        self.quiet
    }

    /// Check whether a valid RefgetStore exists at the given path.
    pub fn store_exists<P: AsRef<Path>>(path: P) -> bool {
        path.as_ref().join("rgstore.json").exists()
    }

    /// Change the storage mode, re-encoding/decoding existing sequences as needed.
    pub fn set_encoding_mode(&mut self, new_mode: StorageMode) {
        if self.mode == new_mode {
            return;
        }

        for record in self.sequence_store.values_mut() {
            match record {
                SequenceRecord::Full { metadata, sequence } => {
                    match (self.mode, new_mode) {
                        (StorageMode::Raw, StorageMode::Encoded) => {
                            let alphabet = lookup_alphabet(&metadata.alphabet);
                            *sequence = encode_sequence(&*sequence, alphabet);
                        }
                        (StorageMode::Encoded, StorageMode::Raw) => {
                            let alphabet = lookup_alphabet(&metadata.alphabet);
                            *sequence =
                                decode_string_from_bytes(&*sequence, metadata.length, alphabet);
                        }
                        _ => {}
                    }
                }
                SequenceRecord::Stub(_) => {}
            }
        }

        self.mode = new_mode;
    }

    /// Enable 2-bit encoding for space efficiency.
    pub fn enable_encoding(&mut self) {
        self.set_encoding_mode(StorageMode::Encoded);
    }

    /// Disable encoding, use raw byte storage.
    pub fn disable_encoding(&mut self) {
        self.set_encoding_mode(StorageMode::Raw);
    }

    /// Enable disk persistence for this store.
    pub fn enable_persistence<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let path = path.as_ref();

        self.local_path = Some(path.to_path_buf());
        self.persist_to_disk = true;
        self.seqdata_path_template
            .get_or_insert_with(|| DEFAULT_SEQDATA_PATH_TEMPLATE.to_string());

        create_dir_all(path.join("sequences"))?;
        create_dir_all(path.join("collections"))?;

        let keys: Vec<DigestKey> = self.sequence_store.keys().cloned().collect();
        for key in keys {
            if let Some(SequenceRecord::Full { metadata, sequence }) = self.sequence_store.get(&key)
            {
                self.write_sequence_to_disk_single(metadata, sequence)?;
                let stub = SequenceRecord::Stub(metadata.clone());
                self.sequence_store.insert(key, stub);
            }
        }

        for record in self.collections.values() {
            self.write_collection_to_disk_single(record)?;
        }

        self.write_index_files()?;

        Ok(())
    }

    /// Disable disk persistence for this store.
    pub fn disable_persistence(&mut self) {
        self.persist_to_disk = false;
    }

    /// Check if persistence to disk is enabled.
    pub fn is_persisting(&self) -> bool {
        self.persist_to_disk
    }

    /// Adds a sequence to the Store
    pub fn add_sequence<T: Into<Option<DigestKey>>>(
        &mut self,
        sequence_record: SequenceRecord,
        collection_digest: T,
        force: bool,
    ) -> Result<()> {
        let collection_digest = collection_digest
            .into()
            .ok_or_else(|| anyhow::anyhow!("Collection digest is required"))?;
        self.collections.get(&collection_digest).ok_or_else(|| {
            anyhow::anyhow!("Collection not found for digest: {:?}", collection_digest)
        })?;

        let metadata = sequence_record.metadata();

        self.name_lookup
            .entry(collection_digest)
            .or_default()
            .insert(metadata.name.clone(), metadata.sha512t24u.to_key());

        self.add_sequence_record(sequence_record, force)?;

        Ok(())
    }

    /// Adds a collection, and all sequences in it, to the store.
    pub fn add_sequence_collection(
        &mut self,
        collection: crate::digest::SequenceCollection,
    ) -> Result<()> {
        self.add_sequence_collection_internal(collection, false)
    }

    /// Adds a collection, overwriting existing data.
    pub fn add_sequence_collection_force(
        &mut self,
        collection: crate::digest::SequenceCollection,
    ) -> Result<()> {
        self.add_sequence_collection_internal(collection, true)
    }

    /// Internal implementation for adding a sequence collection.
    pub(crate) fn add_sequence_collection_internal(
        &mut self,
        collection: crate::digest::SequenceCollection,
        force: bool,
    ) -> Result<()> {
        let coll_digest = collection.metadata.digest.to_key();

        if !force && self.collections.contains_key(&coll_digest) {
            return Ok(());
        }

        let crate::digest::SequenceCollection { metadata, sequences } = collection;

        let record = SequenceCollectionRecord::Full {
            metadata: metadata.clone(),
            sequences: sequences.iter().map(|s| SequenceRecord::Stub(s.metadata().clone())).collect(),
        };

        if self.persist_to_disk && self.local_path.is_some() {
            self.write_collection_to_disk_single(&record)?;
        }

        self.collections.insert(coll_digest, record);

        for sequence_record in sequences {
            self.add_sequence(sequence_record, coll_digest, force)?;
        }

        if self.persist_to_disk && self.local_path.is_some() {
            self.write_index_files()?;
        }

        Ok(())
    }

    /// Adds a SequenceRecord directly to the store without collection association.
    pub fn add_sequence_record(&mut self, sr: SequenceRecord, force: bool) -> Result<()> {
        let metadata = sr.metadata();
        let key = metadata.sha512t24u.to_key();

        if !force && self.sequence_store.contains_key(&key) {
            return Ok(());
        }

        self.md5_lookup
            .insert(metadata.md5.to_key(), metadata.sha512t24u.to_key());

        if self.persist_to_disk && self.local_path.is_some() {
            match &sr {
                SequenceRecord::Full { metadata, sequence } => {
                    self.write_sequence_to_disk_single(metadata, sequence)?;
                    let stub = SequenceRecord::Stub(metadata.clone());
                    self.sequence_store.insert(key, stub);
                    return Ok(());
                }
                SequenceRecord::Stub(_) => {}
            }
        }

        self.sequence_store.insert(key, sr);
        Ok(())
    }

    // =========================================================================
    // Sequence query methods
    // =========================================================================

    /// Returns an iterator over all sequence digests in the store
    pub fn sequence_digests(&self) -> impl Iterator<Item = DigestKey> + '_ {
        self.sequence_store.keys().cloned()
    }

    /// Returns an iterator over sequence metadata for all sequences in the store.
    pub fn sequence_metadata(&self) -> impl Iterator<Item = &SequenceMetadata> + '_ {
        self.sequence_store.values().map(|rec| rec.metadata())
    }

    /// Calculate the total disk size of all sequences in the store
    pub fn total_disk_size(&self) -> usize {
        self.sequence_store
            .values()
            .map(|rec| rec.metadata().disk_size(&self.mode))
            .sum()
    }

    /// Returns the actual disk usage of the store directory.
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

    /// List collections with pagination and optional attribute filtering.
    pub fn list_collections(
        &self,
        page: usize,
        page_size: usize,
        filters: &[(&str, &str)],
    ) -> Result<PagedResult<SequenceCollectionMetadata>> {
        let mut filtered: Vec<SequenceCollectionMetadata> = Vec::new();
        for record in self.collections.values() {
            let meta = record.metadata();
            let mut passes = true;
            for &(attr_name, attr_digest) in filters {
                if !metadata_matches_attribute(meta, attr_name, attr_digest)? {
                    passes = false;
                    break;
                }
            }
            if passes {
                filtered.push(meta.clone());
            }
        }

        filtered.sort_by(|a, b| a.digest.cmp(&b.digest));

        let total = filtered.len();
        let start = page * page_size;
        let results = if start < total {
            filtered.into_iter().skip(start).take(page_size).collect()
        } else {
            Vec::new()
        };

        Ok(PagedResult {
            results,
            pagination: Pagination {
                page,
                page_size,
                total,
            },
        })
    }

    /// Get metadata for a single collection by digest (no sequence data).
    pub fn get_collection_metadata<K: AsRef<[u8]>>(
        &self,
        collection_digest: K,
    ) -> Option<&SequenceCollectionMetadata> {
        let key = collection_digest.to_key();
        self.collections.get(&key).map(|record| record.metadata())
    }

    /// Get a collection with all its sequences loaded.
    pub fn get_collection(&self, collection_digest: &str) -> Result<crate::digest::SequenceCollection> {
        let key = collection_digest.to_key();

        if !self.name_lookup.contains_key(&key) {
            return Err(anyhow!(
                "Collection not loaded: {}. Call load_collection() or load_all_collections() first.",
                collection_digest
            ));
        }

        let metadata = self
            .collections
            .get(&key)
            .ok_or_else(|| anyhow!("Collection not found: {}", collection_digest))?
            .metadata()
            .clone();

        // Iterate name_lookup for (name, digest) pairs so each record gets the
        // correct per-collection name, not the last-written global name.
        let sequences: Vec<SequenceRecord> = self
            .name_lookup
            .get(&key)
            .map(|name_map| {
                name_map
                    .iter()
                    .map(|(name, seq_key)| {
                        let record = self.sequence_store.get(seq_key).ok_or_else(|| {
                            anyhow!(
                                "Sequence {} not found in store for collection {}",
                                key_to_digest_string(seq_key),
                                collection_digest,
                            )
                        })?;
                        let mut meta = record.metadata().clone();
                        meta.name = name.clone();
                        Ok(match record.sequence() {
                            Some(seq) => SequenceRecord::Full {
                                metadata: meta,
                                sequence: seq.to_vec(),
                            },
                            None => SequenceRecord::Stub(meta),
                        })
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .transpose()?
            .unwrap_or_default();

        Ok(crate::digest::SequenceCollection {
            metadata,
            sequences,
        })
    }

    /// Remove a collection from the store.
    pub fn remove_collection(
        &mut self,
        digest: &str,
        remove_orphan_sequences: bool,
    ) -> Result<bool> {
        let key = digest.to_key();

        if self.collections.remove(&key).is_none() {
            return Ok(false);
        }

        let orphan_candidates: Vec<DigestKey> = self
            .name_lookup
            .get(&key)
            .map(|name_map| name_map.values().cloned().collect())
            .unwrap_or_default();

        self.name_lookup.remove(&key);
        self.fhr_metadata.remove(&key);

        // Remove collection aliases pointing to this digest
        let alias_pairs = self.aliases.reverse_lookup_collection(digest);
        let affected_namespaces: std::collections::HashSet<String> = alias_pairs
            .iter()
            .map(|(ns, _)| ns.clone())
            .collect();
        for (ns, alias) in &alias_pairs {
            self.aliases.remove_collection(ns, alias);
        }
        for ns in &affected_namespaces {
            self.persist_alias_namespace(AliasKind::Collection, ns)?;
        }

        if remove_orphan_sequences && !orphan_candidates.is_empty() {
            let mut still_referenced: std::collections::HashSet<DigestKey> =
                std::collections::HashSet::new();
            for name_map in self.name_lookup.values() {
                for seq_key in name_map.values() {
                    still_referenced.insert(*seq_key);
                }
            }

            let orphans: Vec<DigestKey> = orphan_candidates
                .into_iter()
                .filter(|k| !still_referenced.contains(k))
                .collect();

            for orphan_key in &orphans {
                self.sequence_store.remove(orphan_key);
                self.md5_lookup.retain(|_, v| v != orphan_key);
                self.decoded_cache.remove(orphan_key);
            }

            if self.persist_to_disk {
                if let (Some(local_path), Some(template)) =
                    (&self.local_path, &self.seqdata_path_template)
                {
                    for orphan_key in &orphans {
                        let orphan_digest = key_to_digest_string(orphan_key);
                        let seq_file_path = Self::expand_template(&orphan_digest, template);
                        let full_path = local_path.join(&seq_file_path);
                        let _ = fs::remove_file(&full_path);
                        if let Some(parent) = full_path.parent() {
                            let _ = fs::remove_dir(parent);
                        }
                    }
                }
            }
        }

        if self.persist_to_disk {
            if let Some(local_path) = &self.local_path {
                let rgsi_path = local_path.join(format!("collections/{}.rgsi", digest));
                let _ = fs::remove_file(&rgsi_path);
                let fhr_path = local_path.join(format!("fhr/{}.fhr.json", digest));
                let _ = fs::remove_file(&fhr_path);
            }
            self.write_index_files()?;
        }

        Ok(true)
    }

    // =========================================================================
    // Import from another store
    // =========================================================================

    /// Import a single collection (with all its sequences, aliases, and FHR
    /// metadata) from another store into this store.
    ///
    /// Both stores must be disk-backed with matching storage modes.
    /// The source store must have the collection loaded (call
    /// `load_collection()` or `load_all_collections()` first).
    pub fn import_collection(&mut self, source: &ReadonlyRefgetStore, digest: &str) -> Result<()> {
        // Both stores must be disk-backed with same mode
        if source.local_path.is_none() || self.local_path.is_none() || !self.persist_to_disk {
            return Err(anyhow!("import_collection requires both stores to be disk-backed"));
        }
        if source.mode != self.mode {
            return Err(anyhow!(
                "import_collection requires matching storage modes (source={:?}, dest={:?})",
                source.mode,
                self.mode,
            ));
        }

        self.import_collection_file_copy(source, digest)?;

        // Copy sequence aliases for every sequence in the imported collection
        let coll_key = digest.to_key();
        if let Some(name_map) = source.name_lookup.get(&coll_key) {
            for seq_key in name_map.values() {
                let seq_digest = key_to_digest_string(seq_key);
                for (ns, alias) in source.aliases.reverse_lookup_sequence(&seq_digest) {
                    self.add_sequence_alias(&ns, &alias, &seq_digest)?;
                }
            }
        }

        // Copy collection aliases
        for (ns, alias) in source.aliases.reverse_lookup_collection(digest) {
            self.add_collection_alias(&ns, &alias, digest)?;
        }

        // Copy FHR metadata
        if let Some(fhr) = source.get_fhr_metadata(digest) {
            self.set_fhr_metadata(digest, fhr.clone())?;
        }

        Ok(())
    }

    /// File-copy based import: copies RGSI and .seq files directly from
    /// source to destination, then registers the collection in memory.
    fn import_collection_file_copy(
        &mut self,
        source: &ReadonlyRefgetStore,
        digest: &str,
    ) -> Result<()> {
        use crate::collection::read_rgsi_file;

        // 1. Read the source RGSI file to get collection metadata
        let src_rgsi_path = source
            .collection_file_path(digest)
            .ok_or_else(|| anyhow!("Source store has no local path for collection {}", digest))?;
        let collection = read_rgsi_file(&src_rgsi_path)
            .with_context(|| format!("Failed to read source RGSI file: {}", src_rgsi_path.display()))?;

        let mut metadata = collection.metadata.clone();

        // 2. Handle ancillary digests and copy/write the RGSI file
        let dst_rgsi_path = self
            .collection_file_path(digest)
            .ok_or_else(|| anyhow!("Dest store has no local path for collection {}", digest))?;
        if let Some(parent) = dst_rgsi_path.parent() {
            create_dir_all(parent)?;
        }

        let needs_ancillary_rewrite = self.ancillary_digests
            && metadata.name_length_pairs_digest.is_none();

        if needs_ancillary_rewrite {
            // Source lacks ancillary digests but destination wants them --
            // compute them and write a new RGSI file.
            metadata.compute_ancillary_digests(&collection.sequences);
            use crate::collection::SequenceCollectionRecordExt;
            let record = SequenceCollectionRecord::Full {
                metadata: metadata.clone(),
                sequences: collection.sequences.iter()
                    .map(|s| SequenceRecord::Stub(s.metadata().clone()))
                    .collect(),
            };
            record.write_collection_rgsi(&dst_rgsi_path)?;
        } else {
            // Byte-for-byte copy of the RGSI file
            fs::copy(&src_rgsi_path, &dst_rgsi_path)
                .with_context(|| format!(
                    "Failed to copy RGSI {} -> {}",
                    src_rgsi_path.display(),
                    dst_rgsi_path.display(),
                ))?;
        }

        // 3. Copy sequence data files
        for seq_record in &collection.sequences {
            let seq_meta = seq_record.metadata();
            let seq_digest = &seq_meta.sha512t24u;

            let src_seq_path = source
                .sequence_file_path(seq_digest)
                .ok_or_else(|| anyhow!("Source has no path for sequence {}", seq_digest))?;
            let dst_seq_path = self
                .sequence_file_path(seq_digest)
                .ok_or_else(|| anyhow!("Dest has no path for sequence {}", seq_digest))?;

            // Skip if destination already has this sequence (dedup across collections)
            if dst_seq_path.exists() {
                // Still need to register in memory below
            } else {
                if let Some(parent) = dst_seq_path.parent() {
                    create_dir_all(parent)?;
                }
                fs::copy(&src_seq_path, &dst_seq_path).with_context(|| {
                    format!(
                        "Failed to copy sequence {} -> {}",
                        src_seq_path.display(),
                        dst_seq_path.display(),
                    )
                })?;
            }
        }

        // 4. Register in memory
        let coll_key = digest.to_key();

        // Build the collection record with stub sequences
        let stub_sequences: Vec<SequenceRecord> = collection
            .sequences
            .iter()
            .map(|s| SequenceRecord::Stub(s.metadata().clone()))
            .collect();

        let record = SequenceCollectionRecord::Full {
            metadata: metadata.clone(),
            sequences: stub_sequences,
        };
        self.collections.insert(coll_key, record);

        // Register sequences and populate name_lookup
        let mut name_map = IndexMap::new();
        for seq_record in &collection.sequences {
            let seq_meta = seq_record.metadata();
            let seq_key = seq_meta.sha512t24u.to_key();

            name_map.insert(seq_meta.name.clone(), seq_key);

            // Insert stub into sequence_store (skip if already present -- dedup)
            if !self.sequence_store.contains_key(&seq_key) {
                self.sequence_store
                    .insert(seq_key, SequenceRecord::Stub(seq_meta.clone()));
                self.md5_lookup
                    .insert(seq_meta.md5.to_key(), seq_key);
            }
        }
        self.name_lookup.insert(coll_key, name_map);

        // 5. Update index files
        self.write_index_files()?;

        Ok(())
    }

    // =========================================================================
    // Sequence API
    // =========================================================================

    /// List all sequences in the store (metadata only, no sequence data).
    pub fn list_sequences(&self) -> Vec<SequenceMetadata> {
        let mut result: Vec<_> = self
            .sequence_store
            .values()
            .map(|rec| rec.metadata().clone())
            .collect();
        result.sort_by(|a, b| a.sha512t24u.cmp(&b.sha512t24u));
        result
    }

    /// Get metadata for a single sequence by digest (no sequence data).
    pub fn get_sequence_metadata<K: AsRef<[u8]>>(
        &self,
        seq_digest: K,
    ) -> Option<&SequenceMetadata> {
        let key = seq_digest.to_key();
        self.sequence_store.get(&key).map(|rec| rec.metadata())
    }

    /// Get a sequence by its SHA512t24u digest.
    pub fn get_sequence<K: AsRef<[u8]>>(&self, seq_digest: K) -> Result<&SequenceRecord> {
        let digest_key = seq_digest.to_key();
        let actual_key = self
            .md5_lookup
            .get(&digest_key)
            .copied()
            .unwrap_or(digest_key);
        self.sequence_store.get(&actual_key).ok_or_else(|| {
            anyhow!(
                "Sequence not found: {}",
                String::from_utf8_lossy(seq_digest.as_ref())
            )
        })
    }

    /// Ensure a sequence is loaded and decoded into the decoded cache.
    pub fn ensure_decoded<K: AsRef<[u8]>>(&mut self, seq_digest: K) -> Result<()> {
        let digest_key = seq_digest.to_key();
        let actual_key = self
            .md5_lookup
            .get(&digest_key)
            .copied()
            .unwrap_or(digest_key);

        if self.decoded_cache.contains_key(&actual_key) {
            return Ok(());
        }

        let record = self
            .sequence_store
            .get(&actual_key)
            .ok_or_else(|| anyhow!("Sequence not found"))?;
        let decoded = record
            .decode()
            .ok_or_else(|| anyhow!("Sequence not loaded (stub). Call load_sequence() first."))?;

        self.decoded_cache.insert(actual_key, decoded.into_bytes());
        Ok(())
    }

    /// Clear the decoded sequence cache to reclaim memory.
    pub fn clear_decoded_cache(&mut self) {
        self.decoded_cache.clear();
    }

    /// Clear sequence data from the store to free memory.
    pub fn clear(&mut self) {
        self.sequence_store.clear();
        self.decoded_cache.clear();
    }

    /// Get decoded sequence bytes from the cache.
    pub fn sequence_bytes<K: AsRef<[u8]>>(&self, seq_digest: K) -> Option<&[u8]> {
        let digest_key = seq_digest.to_key();
        let actual_key = self
            .md5_lookup
            .get(&digest_key)
            .copied()
            .unwrap_or(digest_key);
        self.decoded_cache.get(&actual_key).map(|v| v.as_slice())
    }

    /// Get a sequence by collection digest and name.
    pub fn get_sequence_by_name<K: AsRef<[u8]>>(
        &self,
        collection_digest: K,
        sequence_name: &str,
    ) -> Result<&SequenceRecord> {
        let collection_key = collection_digest.to_key();

        if !self.name_lookup.contains_key(&collection_key) {
            return Err(anyhow!(
                "Collection not loaded. Call load_collection() or load_all_collections() first."
            ));
        }

        let digest_key = self.name_lookup.get(&collection_key)
            .and_then(|name_map| name_map.get(sequence_name).cloned())
            .ok_or_else(|| anyhow!("Sequence '{}' not found in collection", sequence_name))?;

        let record = self.sequence_store.get(&digest_key).ok_or_else(|| {
            anyhow!("Sequence record not found for '{}'. Call load_sequence() first.", sequence_name)
        })?;

        Ok(record)
    }

    // =========================================================================
    // Loading methods
    // =========================================================================

    /// Eagerly load all Stub collections to Full.
    pub fn load_all_collections(&mut self) -> Result<()> {
        let keys: Vec<DigestKey> = self.collections.keys().cloned().collect();
        for key in keys {
            self.ensure_collection_loaded(&key)?;
        }
        Ok(())
    }

    /// Eagerly load all Stub sequences to Full.
    pub fn load_all_sequences(&mut self) -> Result<()> {
        let keys: Vec<DigestKey> = self.sequence_store.keys().cloned().collect();
        for key in keys {
            self.ensure_sequence_loaded(&key)?;
        }
        Ok(())
    }

    /// Load a single collection by digest.
    pub fn load_collection(&mut self, digest: &str) -> Result<()> {
        let key = digest.to_key();
        self.ensure_collection_loaded(&key)
    }

    /// Load a single sequence by digest.
    pub fn load_sequence(&mut self, digest: &str) -> Result<()> {
        let key = digest.to_key();
        self.ensure_sequence_loaded(&key)
    }

    /// Iterate over all collections with their sequences loaded.
    pub fn iter_collections(&self) -> impl Iterator<Item = crate::digest::SequenceCollection> + '_ {
        let mut digests: Vec<String> = self
            .collections
            .values()
            .map(|rec| rec.metadata().digest.clone())
            .collect();
        digests.sort();

        digests.into_iter().filter_map(move |digest| {
            self.get_collection(&digest).ok()
        })
    }

    /// Iterate over all sequences with their data loaded.
    pub fn iter_sequences(&self) -> impl Iterator<Item = SequenceRecord> + '_ {
        let mut records: Vec<_> = self.sequence_store.values().cloned().collect();
        records.sort_by(|a, b| a.metadata().sha512t24u.cmp(&b.metadata().sha512t24u));
        records.into_iter()
    }

    /// Check if a collection is fully loaded.
    pub fn is_collection_loaded<K: AsRef<[u8]>>(&self, collection_digest: K) -> bool {
        let key = collection_digest.to_key();
        self.collections
            .get(&key)
            .map_or(false, |record| record.has_sequences())
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

    // =========================================================================
    // Substring retrieval
    // =========================================================================

    /// Retrieves a substring from an encoded sequence by its SHA512t24u digest.
    pub fn get_substring<K: AsRef<[u8]>>(
        &self,
        sha512_digest: K,
        start: usize,
        end: usize,
    ) -> Result<String> {
        let digest_key = sha512_digest.to_key();

        let record = self.sequence_store.get(&digest_key).ok_or_else(|| {
            anyhow!(
                "Sequence not found: {}",
                String::from_utf8_lossy(sha512_digest.as_ref())
            )
        })?;
        let (metadata, sequence) = match record {
            SequenceRecord::Stub(_) => return Err(anyhow!("Sequence data not loaded (stub only)")),
            SequenceRecord::Full { metadata, sequence } => (metadata, sequence),
        };

        if start >= metadata.length || end > metadata.length || start >= end {
            return Err(anyhow!(
                "Invalid substring range: start={}, end={}, sequence length={}",
                start,
                end,
                metadata.length
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

    // =========================================================================
    // Streaming retrieval
    // =========================================================================

    /// Resolve the relative on-disk/remote path for a sequence's encoded bytes.
    fn resolve_seq_file_relpath(&self, digest: &str) -> Result<PathBuf> {
        let template = self
            .seqdata_path_template
            .as_deref()
            .unwrap_or(DEFAULT_SEQDATA_PATH_TEMPLATE);
        let relpath = Self::expand_template(digest, template);
        let relpath_str = relpath
            .to_str()
            .ok_or_else(|| anyhow!("Non-UTF8 sequence path"))?;
        Self::sanitize_relative_path(relpath_str)?;
        Ok(relpath)
    }

    /// Stream a (sub)sequence as decoded ASCII bases without loading the full
    /// sequence into memory.
    ///
    /// Opens a bounded byte window on the backing store (local file via
    /// `seek`+`take`, or remote via HTTP `Range: bytes=`), and wraps it in a
    /// [`StreamingDecoder`] when the store is in `Encoded` mode. In `Raw`
    /// mode the bounded byte source is returned directly, since bytes and
    /// bases are 1:1.
    ///
    /// Memory use is O(1) in the sequence length.
    pub fn stream_sequence<K: AsRef<[u8]>>(
        &self,
        sha512_digest: K,
        start: Option<u32>,
        end: Option<u32>,
    ) -> Result<Box<dyn Read + Send>> {
        // 1. Resolve digest (MD5 -> SHA512t24u fallback).
        let digest_key = sha512_digest.to_key();
        let actual_key = self
            .md5_lookup
            .get(&digest_key)
            .copied()
            .unwrap_or(digest_key);

        // 2. Fetch metadata only (do not load bytes).
        let record = self.sequence_store.get(&actual_key).ok_or_else(|| {
            anyhow!(
                "Sequence not found: {}",
                String::from_utf8_lossy(sha512_digest.as_ref())
            )
        })?;
        let metadata = record.metadata();

        // 3. Resolve bounds.
        let length = metadata.length as u64;
        let start = start.map(|s| s as u64).unwrap_or(0);
        let end = end.map(|e| e as u64).unwrap_or(length);
        if start > end || end > length {
            return Err(anyhow!(
                "Invalid stream range: start={}, end={}, length={}",
                start,
                end,
                length
            ));
        }

        // 4. Look up alphabet and compute bits-per-base.
        let alphabet = crate::digest::lookup_alphabet(&metadata.alphabet);
        let bps = match self.mode {
            StorageMode::Encoded => alphabet.bits_per_symbol as u64,
            StorageMode::Raw => 8,
        };

        // 5. Compute byte range + leading-skip bits.
        let start_bit = start * bps;
        let end_bit = end * bps;
        let byte_start = start_bit / 8;
        let byte_end = end_bit.div_ceil(8);
        let leading_skip_bits = (start_bit - byte_start * 8) as u8;
        let bases_to_emit = end - start;
        let byte_len = byte_end - byte_start;

        // Fast path for empty range.
        if bases_to_emit == 0 {
            return Ok(Box::new(std::io::empty()));
        }

        let sha = metadata.sha512t24u.clone();
        let relpath = self.resolve_seq_file_relpath(&sha)?;

        // 6. Build the source reader.
        let source: Box<dyn Read + Send> = if let Some(local) = self.local_path.as_ref() {
            let full = local.join(&relpath);
            if full.exists() {
                let mut file = File::open(&full)
                    .with_context(|| format!("Failed to open local seq file: {}", full.display()))?;
                file.seek(SeekFrom::Start(byte_start))
                    .context("Failed to seek in local seq file")?;
                Box::new(BufReader::new(file).take(byte_len))
            } else if let Some(remote) = self.remote_source.as_ref() {
                open_remote_range(remote, &relpath, byte_start, byte_end)?
            } else {
                return Err(anyhow!(
                    "Sequence file missing locally and no remote source: {}",
                    full.display()
                ));
            }
        } else if let Some(remote) = self.remote_source.as_ref() {
            open_remote_range(remote, &relpath, byte_start, byte_end)?
        } else {
            return Err(anyhow!("No backing source configured for sequence {}", sha));
        };

        // 7. Dispatch on mode.
        match self.mode {
            StorageMode::Encoded => Ok(Box::new(StreamingDecoder::new(
                source,
                alphabet,
                leading_skip_bits,
                bases_to_emit,
            ))),
            StorageMode::Raw => {
                debug_assert_eq!(leading_skip_bits, 0);
                Ok(source)
            }
        }
    }

    // =========================================================================
    // Internal helpers
    // =========================================================================

    /// Expand a path template by substituting digest-based placeholders.
    pub(crate) fn expand_template(digest_str: &str, template: &str) -> PathBuf {
        debug_assert!(
            digest_str.len() >= 4,
            "Digest string must be at least 4 characters for template expansion, got {} chars",
            digest_str.len()
        );
        let path_str = template
            .replace("%s2", digest_str.get(0..2).unwrap_or(digest_str))
            .replace("%s4", digest_str.get(0..4).unwrap_or(digest_str))
            .replace("%s", digest_str);
        PathBuf::from(path_str)
    }

    /// Return the full filesystem path to a sequence `.seq` file for the given digest.
    ///
    /// Returns `None` if the store has no local path or no seqdata path template.
    pub fn sequence_file_path(&self, seq_digest: &str) -> Option<PathBuf> {
        let local_path = self.local_path.as_ref()?;
        let template = self.seqdata_path_template.as_ref()?;
        Some(local_path.join(Self::expand_template(seq_digest, template)))
    }

    /// Return the full filesystem path to a collection RGSI file for the given digest.
    ///
    /// Returns `None` if the store has no local path.
    pub fn collection_file_path(&self, coll_digest: &str) -> Option<PathBuf> {
        let local_path = self.local_path.as_ref()?;
        Some(local_path.join(format!("collections/{}.rgsi", coll_digest)))
    }

    /// Validate a relative path to prevent directory traversal attacks.
    pub(crate) fn sanitize_relative_path(path: &str) -> Result<()> {
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
    pub(crate) fn fetch_file(
        local_path: &Option<PathBuf>,
        remote_source: &Option<String>,
        relative_path: &str,
        persist_to_disk: bool,
        force_refresh: bool,
    ) -> Result<Vec<u8>> {
        Self::sanitize_relative_path(relative_path)?;

        if persist_to_disk && !force_refresh {
            if let Some(local_path) = local_path {
                let full_local_path = local_path.join(relative_path);
                if full_local_path.exists() {
                    return fs::read(&full_local_path).context(format!(
                        "Failed to read local file: {}",
                        full_local_path.display()
                    ));
                }
            }
        }

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

            if persist_to_disk {
                if let Some(local_path) = local_path {
                    let full_local_path = local_path.join(relative_path);

                    if let Some(parent) = full_local_path.parent() {
                        create_dir_all(parent)?;
                    }

                    fs::write(&full_local_path, &data).context(format!(
                        "Failed to cache file to: {}",
                        full_local_path.display()
                    ))?;
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

    /// Ensure a collection is loaded into the store
    pub(crate) fn ensure_collection_loaded(&mut self, collection_digest: &DigestKey) -> Result<()> {
        if self.name_lookup.contains_key(collection_digest) {
            return Ok(());
        }

        let needs_fetch = match self.collections.get(collection_digest) {
            Some(SequenceCollectionRecord::Stub(_)) => true,
            Some(SequenceCollectionRecord::Full { .. }) => false,
            None => true,
        };

        if needs_fetch {
            let digest_str = if let Some(SequenceCollectionRecord::Stub(meta)) =
                self.collections.get(collection_digest)
            {
                meta.digest.clone()
            } else {
                key_to_digest_string(collection_digest)
            };

            let relative_path = format!("collections/{}.rgsi", digest_str);

            if !self.quiet {
                let cached = self
                    .local_path
                    .as_ref()
                    .map(|p| p.join(&relative_path).exists())
                    .unwrap_or(false);
                let verb = if cached { "Loading" } else { "Downloading" };
                eprintln!("{} collection metadata {}...", verb, digest_str);
            }
            let _collection_data =
                Self::fetch_file(&self.local_path, &self.remote_source, &relative_path, true, false)?;

            let local_path = self
                .local_path
                .as_ref()
                .ok_or_else(|| anyhow!("No local path configured"))?;

            let collection_file_path = local_path.join(&relative_path);

            let collection = read_rgsi_file(&collection_file_path)?;

            let loaded_digest = collection.metadata.digest.to_key();
            if loaded_digest != *collection_digest {
                return Err(anyhow!(
                    "Collection digest mismatch: expected {}, got {}",
                    key_to_digest_string(collection_digest),
                    key_to_digest_string(&loaded_digest)
                ));
            }

            let mut name_map = IndexMap::new();
            for sequence_record in &collection.sequences {
                let metadata = sequence_record.metadata();
                let sha512_key = metadata.sha512t24u.to_key();
                name_map.insert(metadata.name.clone(), sha512_key);

                if !self.sequence_store.contains_key(&sha512_key) {
                    self.sequence_store
                        .insert(sha512_key, SequenceRecord::Stub(metadata.clone()));
                    let md5_key = metadata.md5.to_key();
                    self.md5_lookup.insert(md5_key, sha512_key);
                }
            }
            self.name_lookup.insert(*collection_digest, name_map);

            let record = SequenceCollectionRecord::from(collection);
            self.collections.insert(*collection_digest, record);
        } else {
            let sequences_data: Vec<(SequenceMetadata, DigestKey, DigestKey)> =
                if let Some(SequenceCollectionRecord::Full { sequences, .. }) =
                    self.collections.get(collection_digest)
                {
                    sequences
                        .iter()
                        .map(|seq| {
                            let metadata = seq.metadata().clone();
                            let sha512_key = metadata.sha512t24u.to_key();
                            let md5_key = metadata.md5.to_key();
                            (metadata, sha512_key, md5_key)
                        })
                        .collect()
                } else {
                    Vec::new()
                };

            let mut name_map = IndexMap::new();
            for (metadata, sha512_key, md5_key) in sequences_data {
                name_map.insert(metadata.name.clone(), sha512_key);

                if !self.sequence_store.contains_key(&sha512_key) {
                    self.sequence_store
                        .insert(sha512_key, SequenceRecord::Stub(metadata));
                    self.md5_lookup.insert(md5_key, sha512_key);
                }
            }
            self.name_lookup.insert(*collection_digest, name_map);
        }

        Ok(())
    }

    /// Ensure a sequence is loaded into memory
    pub(crate) fn ensure_sequence_loaded(&mut self, digest: &DigestKey) -> Result<()> {
        let record = self
            .sequence_store
            .get(digest)
            .ok_or_else(|| anyhow!("Sequence not found in store"))?;

        if matches!(record, SequenceRecord::Full { .. }) {
            return Ok(());
        }

        let digest_str = &record.metadata().sha512t24u;
        let template = self
            .seqdata_path_template
            .as_ref()
            .ok_or_else(|| anyhow!("No sequence data path template configured"))?;

        let relative_path = Self::expand_template(digest_str, template)
            .to_string_lossy()
            .into_owned();

        if !self.quiet {
            let cached = self
                .local_path
                .as_ref()
                .map(|p| p.join(&relative_path).exists())
                .unwrap_or(false);
            let verb = if cached { "Loading" } else { "Downloading" };
            eprintln!("{} sequence {}...", verb, digest_str);
        }
        let data = Self::fetch_file(
            &self.local_path,
            &self.remote_source,
            &relative_path,
            self.persist_to_disk,
            false,
        )?;

        self.sequence_store.entry(*digest).and_modify(|r| {
            r.load_data(data);
        });

        Ok(())
    }

    // =========================================================================
    // Write methods
    // =========================================================================

    /// Write the store using its configured paths.
    pub fn write(&self) -> Result<()> {
        if !self.persist_to_disk {
            return Err(anyhow!(
                "write() only works with disk-backed stores - use write_store_to_dir() instead"
            ));
        }
        self.write_index_files()
    }

    /// Write a RefgetStore object to a directory
    pub fn write_store_to_dir<P: AsRef<Path>>(
        &self,
        root_path: P,
        seqdata_path_template: Option<&str>,
    ) -> Result<()> {
        let root_path = root_path.as_ref();

        let template = seqdata_path_template
            .or(self.seqdata_path_template.as_deref())
            .unwrap_or(DEFAULT_SEQDATA_PATH_TEMPLATE);

        if !self.quiet {
            eprintln!(
                "Writing store to directory: {}; Using seqdata path template: {}",
                root_path.display(),
                template
            );
        }

        fs::create_dir_all(root_path)?;

        let sequences_dir = root_path.join("sequences");
        fs::create_dir_all(&sequences_dir)?;

        let collections_dir = root_path.join("collections");
        fs::create_dir_all(&collections_dir)?;

        for record in self.sequence_store.values() {
            match record {
                SequenceRecord::Full { metadata, .. } => {
                    let rel_path = Self::expand_template(&metadata.sha512t24u, template);
                    let full_path = root_path.join(&rel_path);
                    record.to_file(full_path)?;
                }
                SequenceRecord::Stub(_) => {
                    continue;
                }
            }
        }

        for record in self.collections.values() {
            let collection_file_path =
                root_path.join(format!("collections/{}.rgsi", record.metadata().digest));
            record.write_collection_rgsi(&collection_file_path)?;
        }

        let sequence_index_path = root_path.join("sequences.rgsi");
        self.write_sequences_rgsi(&sequence_index_path)?;

        let collection_index_path = root_path.join("collections.rgci");
        self.write_collections_rgci(&collection_index_path)?;

        let aliases_dir = root_path.join("aliases");
        self.aliases.write_to_dir(&aliases_dir)?;

        super::fhr_metadata::write_sidecars(&root_path.join("fhr"), &self.fhr_metadata)?;

        self.write_rgstore_json(root_path, template)?;

        Ok(())
    }

    /// Returns statistics about the store
    pub fn stats(&self) -> StoreStats {
        let n_sequences = self.sequence_store.len();
        let n_sequences_loaded = self
            .sequence_store
            .values()
            .filter(|record| record.is_loaded())
            .count();
        let n_collections = self.collections.len();
        let n_collections_loaded = self
            .collections
            .values()
            .filter(|record| record.has_sequences())
            .count();
        let mode_str = match self.mode {
            StorageMode::Raw => "Raw",
            StorageMode::Encoded => "Encoded",
        };
        StoreStats {
            n_sequences,
            n_sequences_loaded,
            n_collections,
            n_collections_loaded,
            storage_mode: mode_str.to_string(),
        }
    }

    /// List alias namespaces available on this store (from manifest).
    pub fn available_alias_namespaces(&self) -> AvailableAliases<'_> {
        AvailableAliases {
            sequences: &self.available_sequence_alias_namespaces,
            collections: &self.available_collection_alias_namespaces,
        }
    }
}

impl Display for ReadonlyRefgetStore {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let total_size = self.total_disk_size();
        let size_str = format_bytes(total_size);
        writeln!(f, "ReadonlyRefgetStore object:")?;
        writeln!(f, "  Mode: {:?}", self.mode)?;
        writeln!(f, "  Disk size: {} ({} bytes)", size_str, total_size)?;
        writeln!(f, ">Sequences (n={}):", self.sequence_store.len())?;
        for (i, (sha512_digest, sequence_record)) in self.sequence_store.iter().take(10).enumerate()
        {
            let metadata = sequence_record.metadata();
            let first_8_chars = match sequence_record {
                SequenceRecord::Stub(_) => "<stub>".to_string(),
                SequenceRecord::Full {
                    metadata,
                    sequence: seq,
                } => {
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
                key_to_digest_string(sha512_digest),
                &metadata.md5,
                &metadata.length,
                &metadata.alphabet,
                first_8_chars
            )?;
        }
        writeln!(f, ">Collections (n={:?}):", self.name_lookup.len())?;
        for (i, (digest, name_map)) in self.name_lookup.iter().enumerate() {
            let seqcol_digest_str = key_to_digest_string(digest);
            writeln!(
                f,
                "  {}. Collection Digest: {:02x?} ({} sequences)",
                i + 1,
                seqcol_digest_str,
                name_map.len()
            )?;
            for (name, sha512_digest) in name_map.iter().take(5) {
                let sha512_str = key_to_digest_string(sha512_digest);
                writeln!(f, "   - Name: {}, SHA512: {:02x?}", name, sha512_str)?;
            }
            if name_map.len() > 5 {
                writeln!(f, "   - ... and {} more", name_map.len() - 5)?;
            }
        }

        Ok(())
    }
}

// Extension traits used by collection.rs
use crate::collection::SequenceCollectionRecordExt;
