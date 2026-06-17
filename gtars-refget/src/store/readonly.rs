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
use crate::digest::{byte_range_for_bases, decode_substring_from_bytes_at_offset};
use crate::hashkeyable::{DigestKey, HashKeyable, key_to_digest_string};
use crate::seqcol::metadata_matches_attribute;

use std::fs::{self, create_dir_all, File};
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::sync::{Arc, Mutex};

/// Default capacity of the per-store open-file-descriptor cache used by the
/// partial-read path. MUST stay well under typical OS fd limits (often 1024):
/// stores can have huge sequence counts (e.g. ~122k for a transcriptome), so an
/// unbounded cache would exhaust file descriptors. Mirrors seqrepo's
/// `fd_cache_size`.
const SEQ_FD_CACHE_CAP: usize = 256;

/// Batch size at which `get_substrings` against a *remote-only* sequence stops
/// issuing one HTTP byte-range request per range and instead downloads the whole
/// `.seq` once (cached locally, like `load_sequence`), then serves every range
/// from the local file. Each remote range is a separate round-trip, so for bulk
/// extraction (e.g. BED region sets) the per-request latency dominates and a
/// single cached download wins; below the threshold, fetching only the touched
/// bytes is cheaper than pulling a whole chromosome. Heuristic, count-based:
/// small/ad-hoc reads stay lean, bulk reads promote.
const REMOTE_BULK_FETCH_THRESHOLD: usize = 16;

/// A tiny bounded LRU cache mapping a sequence digest to an open, shared
/// `.seq` file handle. `.seq` content for a digest is immutable, so a cached
/// handle never goes stale and needs no invalidation. When the cache is full,
/// inserting a new handle evicts (closes) the least-recently-used one.
///
/// Recency is tracked with a monotonically increasing counter stamped on each
/// access; eviction scans for the minimum stamp. With a small cap (default 256)
/// this O(n) eviction is negligible compared to the avoided `open`/`seek`
/// syscalls on the hot path.
#[derive(Debug)]
struct FdCache {
    cap: usize,
    tick: u64,
    /// digest key -> (shared file handle, last-use stamp)
    entries: HashMap<DigestKey, (Arc<File>, u64)>,
}

impl FdCache {
    fn new(cap: usize) -> Self {
        FdCache {
            cap: cap.max(1),
            tick: 0,
            entries: HashMap::with_capacity(cap.max(1)),
        }
    }

    /// Return a cloned handle on hit, bumping recency. `None` on miss.
    fn get(&mut self, key: &DigestKey) -> Option<Arc<File>> {
        self.tick += 1;
        let tick = self.tick;
        if let Some(slot) = self.entries.get_mut(key) {
            slot.1 = tick;
            Some(Arc::clone(&slot.0))
        } else {
            None
        }
    }

    /// Insert a freshly opened handle, evicting (closing) the LRU entry if full.
    fn insert(&mut self, key: DigestKey, file: Arc<File>) {
        self.tick += 1;
        let tick = self.tick;
        if !self.entries.contains_key(&key) && self.entries.len() >= self.cap {
            // Evict the least-recently-used entry. Dropping the Arc closes the
            // fd once no in-flight read still holds a clone.
            if let Some(lru_key) = self
                .entries
                .iter()
                .min_by_key(|(_, (_, stamp))| *stamp)
                .map(|(k, _)| *k)
            {
                self.entries.remove(&lru_key);
            }
        }
        self.entries.insert(key, (file, tick));
    }
}


/// A `Read` adapter that yields a byte range from a shared `Arc<Vec<u8>>`
/// buffer without cloning the backing data.
///
/// Used by `stream_sequence` when the target sequence is already resident
/// in memory as `SequenceRecord::Full`. The alternative of slicing+cloning
/// would allocate a second copy of the (possibly chromosome-sized) buffer
/// and defeat the O(1) memory guarantee of streaming.
pub(crate) struct ArcSliceReader {
    buf: Arc<Vec<u8>>,
    pos: usize,
    end: usize,
}

impl ArcSliceReader {
    pub(crate) fn new(buf: Arc<Vec<u8>>, start: usize, end: usize) -> Self {
        debug_assert!(start <= end);
        debug_assert!(end <= buf.len());
        Self { buf, pos: start, end }
    }
}

impl Read for ArcSliceReader {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        let remaining = self.end - self.pos;
        if remaining == 0 {
            return Ok(0);
        }
        let n = remaining.min(out.len());
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

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
    let agent: ureq::Agent = ureq::Agent::config_builder()
        .timeout_connect(Some(std::time::Duration::from_secs(30)))
        .timeout_recv_response(Some(std::time::Duration::from_secs(60)))
        .build()
        .into();
    let response = agent
        .get(&url)
        .header("Range", &format!("bytes={}-{}", byte_start, last_byte))
        .call()
        .map_err(|e| anyhow!("Failed to open remote byte range: {}", e))?;

    // Guard against servers that ignore the Range header and return the full
    // body with HTTP 200 instead of a 206 Partial Content slice. Wrapping the
    // full body here would silently produce wrong bases, so refuse outright.
    let status = response.status();
    if status.as_u16() != 206 {
        return Err(anyhow!(
            "Remote server did not honor Range header (got HTTP {}); \
             refusing to stream to avoid silent data corruption. URL: {}",
            status,
            url
        ));
    }

    Ok(Box::new(BufReader::with_capacity(
        64 * 1024,
        response.into_body().into_reader().take(byte_len),
    )))
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
    /// Bounded LRU cache of open `.seq` file handles for the partial-read path,
    /// behind a `Mutex` because `get_substring` takes `&self`. Avoids re-opening
    /// the same chromosome file on every per-region `get_substring` call.
    seq_fd_cache: Mutex<FdCache>,
    /// Whether the sequence index (sequences.rgsi) has been loaded.
    /// For remote stores, this starts as `false` and is lazily loaded on first
    /// sequence access, avoiding the costly download when only browsing collections.
    pub(crate) sequence_index_loaded: bool,
    /// The relative path to the sequence index file (from rgstore.json),
    /// stored for deferred loading in remote stores.
    pub(crate) sequence_index_path: Option<String>,
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
            available_sequence_alias_namespaces: Vec::new(),
            available_collection_alias_namespaces: Vec::new(),
            seq_fd_cache: Mutex::new(FdCache::new(SEQ_FD_CACHE_CAP)),
            sequence_index_loaded: true,
            sequence_index_path: None,
        }
    }

    /// Test-only: shrink the partial-read fd-cache to a tiny capacity so the
    /// eviction/reopen path is exercised deterministically. Only used by the
    /// filesystem-gated partial-read tests.
    #[cfg(all(test, feature = "filesystem"))]
    pub(crate) fn set_seq_fd_cache_cap_for_test(&self, cap: usize) {
        let mut cache = self.seq_fd_cache.lock().unwrap();
        *cache = FdCache::new(cap);
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
                            *sequence = std::sync::Arc::new(encode_sequence(&sequence[..], alphabet));
                        }
                        (StorageMode::Encoded, StorageMode::Raw) => {
                            let alphabet = lookup_alphabet(&metadata.alphabet);
                            *sequence = std::sync::Arc::new(decode_string_from_bytes(
                                &sequence[..],
                                metadata.length,
                                alphabet,
                            ));
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

        // CONTRACT: ingested sequence bytes are ASCII (the refget digest and the
        // on-the-fly decode paths assume one byte == one residue). Enforced in
        // debug builds only to keep the ingestion hot path allocation/branch
        // free in release.
        debug_assert!(
            match &sr {
                SequenceRecord::Full { sequence, .. } => sequence.is_ascii(),
                SequenceRecord::Stub(_) => true,
            },
            "add_sequence_record: sequence bytes must be ASCII"
        );

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

    /// Inserter-side half of [`add_sequence_record`](Self::add_sequence_record)
    /// for the parallel streaming import: do the FAST in-memory dedup decision
    /// and map/`md5_lookup` insert here (single inserter, so no race on the
    /// dedup truth), but DO NOT write the `.seq` to disk. Instead return the
    /// expanded full path + bytes so the caller can dispatch the (independent,
    /// per-digest) disk write to a bounded writer pool.
    ///
    /// Returns `Some((full_path, bytes))` when a `.seq` write must be performed
    /// (disk-backed store, `Full` record, and either a first occurrence or a
    /// forced overwrite), or `None` when nothing needs to be written to disk
    /// (in-memory mode, a `Stub`, or a dedup hit). The in-memory state is fully
    /// updated before returning, so dedup stays the single source of truth on
    /// the inserter thread.
    #[cfg_attr(not(feature = "filesystem"), allow(dead_code))]
    pub(crate) fn add_sequence_record_deferred_write(
        &mut self,
        sr: SequenceRecord,
        force: bool,
    ) -> Result<Option<(PathBuf, Vec<u8>)>> {
        let metadata = sr.metadata();
        let key = metadata.sha512t24u.to_key();

        if !force && self.sequence_store.contains_key(&key) {
            return Ok(None);
        }

        self.md5_lookup
            .insert(metadata.md5.to_key(), metadata.sha512t24u.to_key());

        if self.persist_to_disk && self.local_path.is_some() {
            match sr {
                SequenceRecord::Full { metadata, sequence } => {
                    // Insert the stub into the map NOW (reserve), so the dedup
                    // decision is committed on the inserter thread. The actual
                    // byte write is dispatched to the writer pool by the caller.
                    let full_path = self
                        .sequence_file_path(&metadata.sha512t24u)
                        .context("seqdata_path_template/local_path not set")?;
                    let stub = SequenceRecord::Stub(metadata);
                    self.sequence_store.insert(key, stub);
                    // INVARIANT: strong_count == 1 here — the Arc was just
                    // constructed in the import pipeline (Arc::new in import.rs)
                    // and no other owner exists at the point of this call, so
                    // try_unwrap always succeeds and reclaims the buffer with zero
                    // copy. The clone fallback is dead code / future-proofing only.
                    debug_assert_eq!(
                        std::sync::Arc::strong_count(&sequence),
                        1,
                        "sequence Arc must be uniquely owned for zero-copy reclaim"
                    );
                    let bytes = std::sync::Arc::try_unwrap(sequence)
                        .unwrap_or_else(|arc| (*arc).clone());
                    return Ok(Some((full_path, bytes)));
                }
                SequenceRecord::Stub(s) => {
                    self.sequence_store.insert(key, SequenceRecord::Stub(s));
                    return Ok(None);
                }
            }
        }

        self.sequence_store.insert(key, sr);
        Ok(None)
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
                        Ok(match record.sequence_arc() {
                            Some(seq) => SequenceRecord::Full {
                                metadata: meta,
                                sequence: seq,
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

    /// Remove any sequence from `sequence_store` (and its `.seq` file on disk)
    /// that is not referenced by any entry in `self.collections`.
    ///
    /// Used as a best-effort cleanup after a failed import:
    /// [`add_sequence_collections_from_fastas`](Self::add_sequence_collections_from_fastas)
    /// may have inserted sequences for a partially-built collection that never
    /// received an `End` message (and was therefore never added to
    /// `self.collections`).  Those sequences are orphans — they have no
    /// owning collection — and this method removes them so the store stays
    /// internally consistent. The operation is O(sequences + collections).
    ///
    /// Sequences shared with a successfully-finalized collection are preserved.
    #[cfg_attr(not(feature = "filesystem"), allow(dead_code))]
    pub(crate) fn remove_orphan_seq_files(&mut self) {
        // Build the set of all sequence digest keys that are referenced by at
        // least one surviving (finalized) collection.
        let mut referenced: std::collections::HashSet<DigestKey> =
            std::collections::HashSet::new();
        for name_map in self.name_lookup.values() {
            for &seq_key in name_map.values() {
                referenced.insert(seq_key);
            }
        }

        // Collect the orphan keys (present in sequence_store but not referenced).
        let orphans: Vec<DigestKey> = self
            .sequence_store
            .keys()
            .filter(|k| !referenced.contains(*k))
            .copied()
            .collect();

        if orphans.is_empty() {
            return;
        }

        // Remove from in-memory structures.
        for key in &orphans {
            self.sequence_store.remove(key);
            self.md5_lookup.retain(|_, v| v != key);
        }

        // Best-effort remove on-disk `.seq` files.
        if self.persist_to_disk {
            if let (Some(local_path), Some(template)) =
                (&self.local_path, &self.seqdata_path_template)
            {
                for key in &orphans {
                    let digest_str = key_to_digest_string(key);
                    let rel = Self::expand_template(&digest_str, template);
                    let full = local_path.join(&rel);
                    let _ = fs::remove_file(&full);
                    if let Some(parent) = full.parent() {
                        let _ = fs::remove_dir(parent); // ignore if non-empty
                    }
                }
            }
        }
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

    /// Clear sequence data from the store to free memory.
    ///
    /// Removes all sequence records. Metadata is preserved:
    /// collections, name lookups, MD5 lookups, aliases, and FHR
    /// metadata remain intact.
    ///
    /// For on-disk stores, data on disk is unaffected.
    pub fn clear(&mut self) {
        self.sequence_store.clear();
        if let Ok(mut cache) = self.seq_fd_cache.lock() {
            cache.entries.clear();
        }
    }

    /// Check whether a sequence's record is loaded (Full) rather
    /// than a metadata-only Stub.
    pub fn is_sequence_loaded<K: AsRef<[u8]>>(&self, seq_digest: K) -> bool {
        let digest_key = seq_digest.to_key();
        let actual_key = self
            .md5_lookup
            .get(&digest_key)
            .copied()
            .unwrap_or(digest_key);
        match self.sequence_store.get(&actual_key) {
            Some(SequenceRecord::Stub(_)) | None => false,
            Some(SequenceRecord::Full { .. }) => true,
        }
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

    /// Ensure the sequence index (sequences.rgsi) is loaded.
    /// For remote stores this is deferred until first sequence access.
    /// Returns Ok(()) immediately if already loaded.
    pub(crate) fn ensure_sequence_index_loaded(&mut self) -> Result<()> {
        if self.sequence_index_loaded {
            return Ok(());
        }

        let seq_index_path = self.sequence_index_path.clone()
            .ok_or_else(|| anyhow!("No sequence_index_path set for deferred loading"))?;

        if !self.quiet {
            eprintln!("Downloading sequence index {}...", seq_index_path);
        }

        let data = Self::fetch_file(
            &self.local_path,
            &self.remote_source,
            &seq_index_path,
            true,
            false,
        )?;
        let data_str = String::from_utf8(data)
            .context("sequence index contains invalid UTF-8")?;

        Self::load_sequences_from_reader(self, data_str.as_bytes())?;
        self.sequence_index_loaded = true;
        Ok(())
    }

    /// Explicitly load the sequence index. For servers that want to preload
    /// all data during startup before converting to `ReadonlyRefgetStore`.
    pub fn load_sequence_index(&mut self) -> Result<()> {
        self.ensure_sequence_index_loaded()
    }

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
        let actual_key = self
            .md5_lookup
            .get(&digest_key)
            .copied()
            .unwrap_or(digest_key);

        let record = self.sequence_store.get(&actual_key).ok_or_else(|| {
            anyhow!(
                "Sequence not found: {}",
                String::from_utf8_lossy(sha512_digest.as_ref())
            )
        })?;
        let (metadata, sequence): (&SequenceMetadata, &[u8]) = match record {
            // Not resident: for a disk-backed store, read only the bytes covering
            // [start, end) directly from the `.seq` file instead of loading the
            // whole sequence into RAM. This is the random-access extract path:
            // sparse region extraction touches a tiny fraction of each sequence,
            // so a per-query partial read avoids the cost of `load_sequence`
            // pulling entire chromosomes into memory.
            SequenceRecord::Stub(meta) => {
                // Partial-read resolution for a non-resident sequence:
                //   local `.seq` (if present) -> remote byte-range (if configured).
                // A missing local file falls through to the remote fallback
                // rather than erroring.
                if self.seqdata_path_template.is_some() {
                    if let Some(path) = self.sequence_file_path(&meta.sha512t24u) {
                        if path.exists() {
                            return self.get_substring_from_disk(meta, start, end);
                        }
                    }
                }
                if self.remote_source.is_some() {
                    return self.get_substring_from_remote(meta, start, end);
                }
                return Err(anyhow!("Sequence data not loaded (stub only)"));
            }
            SequenceRecord::Full { metadata, sequence } => (metadata, sequence.as_slice()),
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
                // SAFETY: the Encoded decoder emits only bytes drawn from
                // `alphabet.decoding_array`, which are ASCII sequence symbols
                // (A/C/G/T/N/IUPAC/protein/ASCII), hence valid UTF-8. Skipping
                // the O(n) `from_utf8` validation matters on multi-GB extracts.
                Ok(unsafe { String::from_utf8_unchecked(decoded_sequence) })
            }
            StorageMode::Raw => {
                let raw_slice: &[u8] = &sequence[start..end];
                // SAFETY: Raw mode stores one ASCII sequence byte per base per the
                // refget data model, so the slice is valid UTF-8.
                Ok(unsafe { String::from_utf8_unchecked(raw_slice.to_vec()) })
            }
        }
    }

    /// Batch substring retrieval for many ranges of ONE sequence.
    ///
    /// Resolves the record once and validates all ranges up front, then:
    /// - **Resident `Full`**: decode/slice each range from the resident bytes.
    /// - **Disk-backed `Stub`**: partial-read each range independently via the
    ///   fd-cached [`get_substring_from_disk`](Self::get_substring_from_disk)
    ///   path (reads only the bytes covering each range; reuses the open `.seq`
    ///   handle across ranges).
    ///
    /// The main benefit over a loop of [`get_substring`](Self::get_substring) is
    /// for FFI callers (e.g. the Python binding): resolving the record once and
    /// crossing the boundary once per sequence instead of once per range.
    ///
    /// A whole-sequence "decode once, slice all" strategy was tried for
    /// wide/overlapping range sets but measured *slower* on real chromosome-scale
    /// sequences (decoding a 250 Mbp chromosome into one buffer, plus the large
    /// allocation, outweighs the avoided re-decode of scattered spans), so the
    /// partial path is used uniformly. Callers wanting resident behavior should
    /// load the sequence (`load_sequence`) first; the `Full` arm then serves it.
    ///
    /// Output strings are byte-identical to a loop of
    /// [`get_substring`](Self::get_substring).
    pub fn get_substrings<K: AsRef<[u8]>>(
        &self,
        sha512_digest: K,
        ranges: &[(usize, usize)],
    ) -> Result<Vec<String>> {
        let digest_key = sha512_digest.to_key();
        let actual_key = self
            .md5_lookup
            .get(&digest_key)
            .copied()
            .unwrap_or(digest_key);

        let record = self.sequence_store.get(&actual_key).ok_or_else(|| {
            anyhow!(
                "Sequence not found: {}",
                String::from_utf8_lossy(sha512_digest.as_ref())
            )
        })?;

        let metadata = record.metadata();
        let seq_length = metadata.length;

        // Validate every range up front (same rule as get_substring).
        for &(start, end) in ranges {
            if start >= seq_length || end > seq_length || start >= end {
                return Err(anyhow!(
                    "Invalid substring range: start={}, end={}, sequence length={}",
                    start,
                    end,
                    seq_length
                ));
            }
        }

        match record {
            SequenceRecord::Full { metadata, sequence } => {
                let seq: &[u8] = sequence.as_slice();
                let mut out = Vec::with_capacity(ranges.len());
                match self.mode {
                    StorageMode::Encoded => {
                        let alphabet = lookup_alphabet(&metadata.alphabet);
                        for &(start, end) in ranges {
                            let decoded =
                                decode_substring_from_bytes(seq, start, end, alphabet);
                            // SAFETY: see get_substring (ASCII alphabet bytes).
                            out.push(unsafe { String::from_utf8_unchecked(decoded) });
                        }
                    }
                    StorageMode::Raw => {
                        for &(start, end) in ranges {
                            // SAFETY: see get_substring (Raw stores ASCII bases).
                            out.push(unsafe {
                                String::from_utf8_unchecked(seq[start..end].to_vec())
                            });
                        }
                    }
                }
                Ok(out)
            }
            SequenceRecord::Stub(meta) => {
                // Resolve once for the whole batch: local `.seq` (if present) ->
                // remote byte-range (if configured). A missing local file falls
                // through to the remote fallback rather than erroring.
                let local_seq_exists = self.seqdata_path_template.is_some()
                    && self
                        .sequence_file_path(&meta.sha512t24u)
                        .map(|p| p.exists())
                        .unwrap_or(false);
                let has_remote = self.remote_source.is_some();
                if !local_seq_exists && !has_remote {
                    return Err(anyhow!("Sequence data not loaded (stub only)"));
                }

                // Bulk remote extraction: one HTTP round-trip per range is slow
                // for large batches. Past REMOTE_BULK_FETCH_THRESHOLD, download and
                // cache the whole `.seq` once (needs a local cache dir + path
                // template to write into), then read every range from the local
                // file. Single / small-batch reads fall through to pure byte-range.
                if !local_seq_exists
                    && has_remote
                    && self.local_path.is_some()
                    && self.seqdata_path_template.is_some()
                    && ranges.len() >= REMOTE_BULK_FETCH_THRESHOLD
                {
                    let relpath = self.resolve_seq_file_relpath(&meta.sha512t24u)?;
                    let relpath_str = relpath
                        .to_str()
                        .ok_or_else(|| anyhow!("Non-UTF8 sequence path"))?;
                    // persist_to_disk=true caches the file for this and future
                    // reads; force_refresh=false reuses it if already present.
                    Self::fetch_file(&self.local_path, &self.remote_source, relpath_str, true, false)?;
                    let mut out = Vec::with_capacity(ranges.len());
                    for &(start, end) in ranges {
                        out.push(self.get_substring_from_disk(meta, start, end)?);
                    }
                    return Ok(out);
                }

                let mut out = Vec::with_capacity(ranges.len());
                for &(start, end) in ranges {
                    if local_seq_exists {
                        out.push(self.get_substring_from_disk(meta, start, end)?);
                    } else {
                        out.push(self.get_substring_from_remote(meta, start, end)?);
                    }
                }
                Ok(out)
            }
        }
    }

    /// Read a substring directly from the on-disk `.seq` file, fetching only the
    /// bytes that cover `[start, end)`.
    ///
    /// `.seq` files are headerless raw byte arrays (bit-packed in Encoded mode,
    /// one ASCII byte per base in Raw mode), so byte offsets map directly to base
    /// positions. We compute the covering byte range with [`byte_range_for_bases`],
    /// read just that window with a positioned read, and decode it with
    /// [`decode_substring_from_bytes_at_offset`]. This is the partial-read path
    /// that makes random-access extraction cheap: it never materializes the whole
    /// sequence in memory.
    fn get_substring_from_disk(
        &self,
        metadata: &SequenceMetadata,
        start: usize,
        end: usize,
    ) -> Result<String> {
        if start >= metadata.length || end > metadata.length || start >= end {
            return Err(anyhow!(
                "Invalid substring range: start={}, end={}, sequence length={}",
                start,
                end,
                metadata.length
            ));
        }

        // Fetch (or open + cache) the shared `.seq` file handle. The critical
        // section is kept tiny: lock -> get-or-open -> clone Arc -> unlock. The
        // actual positioned read happens OUTSIDE the lock, so concurrent reads
        // are never serialized by the cache.
        let file = self.get_cached_seq_file(&metadata.sha512t24u)?;

        // Compute the covering byte window, then read it with a positioned read.
        let (byte_start, byte_end) = match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                byte_range_for_bases(start, end, alphabet.bits_per_symbol)
            }
            StorageMode::Raw => (start, end),
        };
        let buf = crate::posread::read_exact_window(&file, byte_start, byte_end - byte_start)?;

        let decoded: Vec<u8> = match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                decode_substring_from_bytes_at_offset(&buf, byte_start, start, end, alphabet)
            }
            StorageMode::Raw => buf,
        };

        // SAFETY: Encoded decode emits only ASCII alphabet bytes; Raw stores
        // ASCII sequence bytes. Either way `decoded` is valid UTF-8, so we skip
        // the O(n) validation. See `get_substring` for the full justification.
        Ok(unsafe { String::from_utf8_unchecked(decoded) })
    }

    /// Return a shared handle to the `.seq` file for `digest`, opening and
    /// caching it on a miss. Bounded LRU; the LRU handle is evicted (closed)
    /// when the cache is full.
    fn get_cached_seq_file(&self, digest: &str) -> Result<Arc<File>> {
        let key = digest.to_key();

        // Fast path: hit. Hold the lock only for the tiny map op.
        {
            let mut cache = self
                .seq_fd_cache
                .lock()
                .map_err(|_| anyhow!("seq fd cache mutex poisoned"))?;
            if let Some(f) = cache.get(&key) {
                return Ok(f);
            }
        }

        // Miss: build the path and open OUTSIDE the lock (the open syscall must
        // not serialize other readers).
        let local_path = self
            .local_path
            .as_ref()
            .ok_or_else(|| anyhow!("Partial read requires a local disk-backed store"))?;
        let template = self
            .seqdata_path_template
            .as_ref()
            .ok_or_else(|| anyhow!("No sequence data path template configured"))?;
        let full_path = local_path.join(Self::expand_template(digest, template));
        let file = Arc::new(
            File::open(&full_path)
                .with_context(|| format!("Failed to open seq file {}", full_path.display()))?,
        );

        // Insert under the lock. A concurrent miss may have inserted the same
        // digest meanwhile; that is harmless (both handles read the same
        // immutable file) -- prefer the already-cached one to avoid duplicating.
        let mut cache = self
            .seq_fd_cache
            .lock()
            .map_err(|_| anyhow!("seq fd cache mutex poisoned"))?;
        if let Some(existing) = cache.get(&key) {
            return Ok(existing);
        }
        cache.insert(key, Arc::clone(&file));
        Ok(file)
    }

    /// Read a substring from a remote `.seq` file via an HTTP byte-range
    /// request, fetching only the bytes that cover `[start, end)` and persisting
    /// nothing locally.
    ///
    /// This is the remote fallback for the partial-read path (flow 1): when a
    /// sequence is a remote-only `Stub` — no resident bytes and no cached local
    /// `.seq` — [`get_substring`](Self::get_substring) /
    /// [`get_substrings`](Self::get_substrings) fall through to here instead of
    /// erroring. It mirrors
    /// [`get_substring_from_disk`](Self::get_substring_from_disk) but fetches the
    /// covering byte window over HTTP `Range:` (reusing [`open_remote_range`], the
    /// same machinery as [`stream_sequence`](Self::stream_sequence)) rather than
    /// reading from a local file. Decoded with
    /// [`decode_substring_from_bytes_at_offset`].
    ///
    /// Kept `&self` so the `Arc`-shared readonly store stays usable; this path
    /// never mutates the store and never writes to disk (use `load_sequence`,
    /// flow 3, for caching).
    fn get_substring_from_remote(
        &self,
        metadata: &SequenceMetadata,
        start: usize,
        end: usize,
    ) -> Result<String> {
        if start >= metadata.length || end > metadata.length || start >= end {
            return Err(anyhow!(
                "Invalid substring range: start={}, end={}, sequence length={}",
                start,
                end,
                metadata.length
            ));
        }

        let remote = self
            .remote_source
            .as_ref()
            .ok_or_else(|| anyhow!("Remote partial read requires a remote source"))?;

        let relpath = self.resolve_seq_file_relpath(&metadata.sha512t24u)?;

        // Compute the covering byte window. `.seq` files are headerless raw byte
        // arrays, so byte offsets map directly to base positions.
        let (byte_start, byte_end) = match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                byte_range_for_bases(start, end, alphabet.bits_per_symbol)
            }
            StorageMode::Raw => (start, end),
        };

        let mut reader =
            open_remote_range(remote, &relpath, byte_start as u64, byte_end as u64)?;
        let mut buf = Vec::with_capacity(byte_end - byte_start);
        reader
            .read_to_end(&mut buf)
            .context("Failed to read remote byte range")?;

        let decoded: Vec<u8> = match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                decode_substring_from_bytes_at_offset(&buf, byte_start, start, end, alphabet)
            }
            StorageMode::Raw => buf,
        };

        // SAFETY: see `get_substring` — Encoded decode emits only ASCII alphabet
        // bytes; Raw stores ASCII sequence bytes. Either way valid UTF-8.
        Ok(unsafe { String::from_utf8_unchecked(decoded) })
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
        start: Option<u64>,
        end: Option<u64>,
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
        let start = start.unwrap_or(0);
        let end = end.unwrap_or(length);
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

        // 6. Build the source reader. In-memory Full records are streamed
        // directly from their encoded byte buffer. Otherwise try the
        // on-disk path (local_path / remote_source).
        if let SequenceRecord::Full { sequence, .. } = record {
            // Share ownership of the already-resident buffer instead of
            // cloning the byte range. This keeps peak allocation O(1) in
            // the sequence size.
            let shared = Arc::clone(sequence);
            let source: Box<dyn Read + Send> = Box::new(ArcSliceReader::new(
                shared,
                byte_start as usize,
                byte_end as usize,
            ));
            return match self.mode {
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
            };
        }

        let relpath = self.resolve_seq_file_relpath(&sha)?;

        let source: Box<dyn Read + Send> = if let Some(local) = self.local_path.as_ref() {
            let full = local.join(&relpath);
            if full.exists() {
                let mut file = File::open(&full)
                    .with_context(|| format!("Failed to open local seq file: {}", full.display()))?;
                file.seek(SeekFrom::Start(byte_start))
                    .context("Failed to seek in local seq file")?;
                let file_len = file.metadata()
                    .with_context(|| format!("Failed to stat local seq file: {}", full.display()))?.len();
                if file_len < byte_start + byte_len {
                    return Err(anyhow!(
                        "Local seq file is too short to satisfy the requested byte range \
                         [{}..{}) (file has {} bytes): {}",
                        byte_start, byte_start + byte_len,
                        file_len,
                        full.display()
                    ));
                }
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

        #[cfg(feature = "http")]
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
                .into_body()
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

            return Ok(data);
        }

        let _ = remote_source;
        Err(anyhow!(
            "File not found locally and no remote source configured: {}",
            relative_path
        ))
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

        // Already loaded (Full or Mmap)
        if record.is_loaded() {
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

#[cfg(all(test, feature = "http"))]
mod open_remote_range_tests {
    use super::*;
    use std::net::TcpListener;
    use std::sync::{Arc, atomic::{AtomicBool, Ordering}};

    /// Spin up an HTTP server that DELIBERATELY IGNORES the Range header,
    /// responding 200 OK with the full body. Models a misbehaving CDN/origin.
    /// Returns (base_url, shutdown_fn).
    fn start_range_ignoring_server(body: Vec<u8>) -> (String, impl FnOnce()) {
        use std::io::{Read as _, Write as _};

        let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
        let port = listener.local_addr().unwrap().port();
        let base_url = format!("http://127.0.0.1:{}", port);
        let stop = Arc::new(AtomicBool::new(false));
        let stop_clone = Arc::clone(&stop);

        std::thread::spawn(move || {
            listener.set_nonblocking(false).ok();
            while !stop_clone.load(Ordering::Relaxed) {
                match listener.accept() {
                    Ok((mut stream, _)) => {
                        let mut buf = [0u8; 4096];
                        let _ = stream.read(&mut buf).unwrap_or(0);
                        // Respond 200 OK with the full body, ignoring the Range header.
                        let header = format!(
                            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                            body.len()
                        );
                        let _ = stream.write_all(header.as_bytes());
                        let _ = stream.write_all(&body);
                    }
                    Err(_) => break,
                }
            }
        });

        let shutdown = move || {
            stop.store(true, Ordering::Relaxed);
            let _ = std::net::TcpStream::connect(format!("127.0.0.1:{}", port));
        };

        (base_url, shutdown)
    }

    /// When a remote server ignores the Range header and returns 200 + full body,
    /// `open_remote_range` MUST return an error rather than silently wrapping the
    /// wrong bytes in a StreamingDecoder. Otherwise callers get silent corruption.
    #[test]
    fn test_open_remote_range_rejects_non_206_response() {
        let full_body: Vec<u8> = (0u8..=255).collect();
        let (base_url, shutdown) = start_range_ignoring_server(full_body.clone());

        // Request only a middle slice: bytes 10..20.
        let result = open_remote_range(
            &base_url,
            std::path::Path::new("seq.dat"),
            10,
            20,
        );

        shutdown();

        match result {
            Err(e) => {
                let msg = format!("{}", e);
                assert!(
                    msg.contains("Range") || msg.contains("206") || msg.contains("status"),
                    "expected a Range/206/status error, got: {}",
                    msg
                );
            }
            Ok(mut reader) => {
                use std::io::Read as _;
                let mut got = Vec::new();
                reader.read_to_end(&mut got).unwrap();
                // If we got here, the bug is live: the reader yielded bytes
                // that are NOT the requested 10..20 slice. The `.take(byte_len)`
                // in the buggy code returns the FIRST `byte_len` bytes of the
                // full body, i.e. 0..10, instead of the requested 10..20.
                assert_eq!(
                    got,
                    full_body[10..20],
                    "BUG: open_remote_range silently returned wrong bytes \
                     when server ignored the Range header. Got {:?}, expected {:?}. \
                     The function should have returned an error instead.",
                    got,
                    &full_body[10..20]
                );
            }
        }
    }
}

/// Tests for the flow-1 remote byte-range fallback: a cold `get_substring` /
/// `get_substrings` against a remote-only `Stub` must fetch just the covering
/// bytes over HTTP `Range:`, decode correctly, and persist nothing.
#[cfg(all(test, feature = "http"))]
mod remote_partial_read_tests {
    use super::*;
    use crate::digest::digest_sequence;
    use std::io::{Read as _, Write as _};
    use std::net::TcpListener;
    use std::sync::{
        atomic::{AtomicBool, Ordering},
        Arc,
    };

    /// Spin up an HTTP server that HONORS the `Range` header, replying 206 with
    /// the requested byte slice (and 200 + full body if no Range is sent). The
    /// same `body` is served for every path. Returns (base_url, shutdown_fn).
    fn start_range_honoring_server(body: Vec<u8>) -> (String, impl FnOnce()) {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
        let port = listener.local_addr().unwrap().port();
        let base_url = format!("http://127.0.0.1:{}", port);
        let stop = Arc::new(AtomicBool::new(false));
        let stop_clone = Arc::clone(&stop);

        std::thread::spawn(move || {
            while !stop_clone.load(Ordering::Relaxed) {
                match listener.accept() {
                    Ok((mut stream, _)) => {
                        let mut buf = [0u8; 4096];
                        let n = stream.read(&mut buf).unwrap_or(0);
                        let req = String::from_utf8_lossy(&buf[..n]);
                        // Parse "Range: bytes=START-END" (inclusive END).
                        let range = req.lines().find_map(|l| {
                            let l = l.trim();
                            if !l.to_ascii_lowercase().starts_with("range:") {
                                return None;
                            }
                            let v = l.split('=').nth(1)?.trim();
                            let mut parts = v.split('-');
                            let s: usize = parts.next()?.trim().parse().ok()?;
                            let e: usize = parts.next()?.trim().parse().ok()?;
                            Some((s, e))
                        });
                        if let Some((s, e)) = range {
                            let s = s.min(body.len());
                            let end_excl = (e + 1).min(body.len());
                            let slice = &body[s..end_excl];
                            let header = format!(
                                "HTTP/1.1 206 Partial Content\r\nContent-Range: bytes {}-{}/{}\r\n\
                                 Content-Length: {}\r\nConnection: close\r\n\r\n",
                                s,
                                end_excl.saturating_sub(1),
                                body.len(),
                                slice.len()
                            );
                            let _ = stream.write_all(header.as_bytes());
                            let _ = stream.write_all(slice);
                        } else {
                            let header = format!(
                                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                                body.len()
                            );
                            let _ = stream.write_all(header.as_bytes());
                            let _ = stream.write_all(&body);
                        }
                    }
                    Err(_) => break,
                }
            }
        });

        let shutdown = move || {
            stop.store(true, Ordering::Relaxed);
            let _ = std::net::TcpStream::connect(format!("127.0.0.1:{}", port));
        };

        (base_url, shutdown)
    }

    /// A remote-only store: a `Stub` whose encoded bytes live solely on the
    /// remote server. No `local_path`, so nothing can be (or is) persisted.
    fn remote_only_store(remote: &str, record: &SequenceRecord) -> ReadonlyRefgetStore {
        let mut store = ReadonlyRefgetStore::new(StorageMode::Encoded);
        store.remote_source = Some(remote.to_string());
        store.seqdata_path_template = Some(DEFAULT_SEQDATA_PATH_TEMPLATE.to_string());
        let meta = record.metadata().clone();
        store
            .sequence_store
            .insert(meta.sha512t24u.to_key(), SequenceRecord::Stub(meta));
        store
    }

    #[test]
    fn test_get_substring_cold_remote_byte_range_and_no_persist() {
        // Mixed ACGT + IUPAC (N) bases exercise the variable-bit decode path.
        let bases = b"ACGTACGTACGTTTGGCCAANNNNACGTACGTACGTACGTACGT";
        let record = digest_sequence("chrTest", bases);
        let meta = record.metadata().clone();
        let alphabet = lookup_alphabet(&meta.alphabet);
        let encoded = encode_sequence(&bases[..], alphabet);

        let (base_url, shutdown) = start_range_honoring_server(encoded.clone());
        let store = remote_only_store(&base_url, &record);

        // Cold read of an interior sub-range: must fetch via byte-range + decode.
        let got = store
            .get_substring(meta.sha512t24u.as_bytes(), 4, 12)
            .expect("cold remote get_substring should succeed");
        let expected = String::from_utf8(bases[4..12].to_vec()).unwrap();
        assert_eq!(got, expected, "remote byte-range decode mismatch");

        // Batch read of several disjoint ranges via the same remote fallback.
        let batch = store
            .get_substrings(meta.sha512t24u.as_bytes(), &[(0, 4), (20, 24), (40, 44)])
            .expect("cold remote get_substrings should succeed");
        shutdown();
        assert_eq!(batch[0], String::from_utf8(bases[0..4].to_vec()).unwrap());
        assert_eq!(batch[1], String::from_utf8(bases[20..24].to_vec()).unwrap());
        assert_eq!(batch[2], String::from_utf8(bases[40..44].to_vec()).unwrap());

        // Persists nothing: no local path, and the record stays a Stub (never
        // promoted to a resident Full by the partial-read path).
        assert!(store.local_path.is_none());
        assert!(
            matches!(
                store.sequence_store.get(&meta.sha512t24u.to_key()).unwrap(),
                SequenceRecord::Stub(_)
            ),
            "remote partial read must not promote the record to Full (no caching)"
        );
    }

    #[test]
    fn test_get_substrings_bulk_remote_promotes_to_whole_download() {
        // 64-base sequence so we can request many disjoint ranges.
        let bases: &[u8] =
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let record = digest_sequence("chrBulk", bases);
        let meta = record.metadata().clone();
        let alphabet = lookup_alphabet(&meta.alphabet);
        let encoded = encode_sequence(bases, alphabet);

        let (base_url, shutdown) = start_range_honoring_server(encoded.clone());

        // Remote store WITH a local cache dir, so the bulk path can promote into it.
        let cache_dir =
            std::env::temp_dir().join(format!("gtars_bulk_fetch_test_{}", meta.sha512t24u));
        let _ = std::fs::remove_dir_all(&cache_dir);
        std::fs::create_dir_all(&cache_dir).expect("create cache dir");

        let mut store = ReadonlyRefgetStore::new(StorageMode::Encoded);
        store.remote_source = Some(base_url);
        store.seqdata_path_template = Some(DEFAULT_SEQDATA_PATH_TEMPLATE.to_string());
        store.local_path = Some(cache_dir.clone());
        store
            .sequence_store
            .insert(meta.sha512t24u.to_key(), SequenceRecord::Stub(meta.clone()));

        let cached_path = cache_dir.join(ReadonlyRefgetStore::expand_template(
            &meta.sha512t24u,
            DEFAULT_SEQDATA_PATH_TEMPLATE,
        ));

        // Sub-threshold batch: served per-range by byte-range; nothing cached.
        let small: Vec<(usize, usize)> = vec![(0, 3), (8, 11), (60, 63)];
        let small_out = store
            .get_substrings(meta.sha512t24u.as_bytes(), &small)
            .expect("small remote get_substrings should succeed");
        for (i, &(s, e)) in small.iter().enumerate() {
            assert_eq!(small_out[i], String::from_utf8(bases[s..e].to_vec()).unwrap());
        }
        assert!(
            !cached_path.exists(),
            "a sub-threshold batch must not download/cache the whole .seq"
        );

        // Bulk batch (>= REMOTE_BULK_FETCH_THRESHOLD): promote to one whole
        // download, cache it, and serve every range from the local file.
        let bulk: Vec<(usize, usize)> = (0..REMOTE_BULK_FETCH_THRESHOLD)
            .map(|i| (i * 4, i * 4 + 3))
            .collect();
        assert_eq!(bulk.len(), REMOTE_BULK_FETCH_THRESHOLD);
        let bulk_out = store
            .get_substrings(meta.sha512t24u.as_bytes(), &bulk)
            .expect("bulk remote get_substrings should succeed");
        shutdown();
        for (i, &(s, e)) in bulk.iter().enumerate() {
            assert_eq!(
                bulk_out[i],
                String::from_utf8(bases[s..e].to_vec()).unwrap(),
                "bulk range {} decode mismatch",
                i
            );
        }

        // The promotion cached the WHOLE encoded `.seq` (not a partial window).
        assert!(
            cached_path.exists(),
            "bulk batch must cache the whole .seq locally"
        );
        let on_disk = std::fs::read(&cached_path).expect("read cached .seq");
        assert_eq!(on_disk, encoded, "cached file must be the full encoded sequence");

        let _ = std::fs::remove_dir_all(&cache_dir);
    }
}
