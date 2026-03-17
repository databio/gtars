//! RefgetStore wrapper struct with lazy-loading read methods.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::fmt::{Display, Formatter};
use std::path::Path;

use anyhow::Result;
use std::fs::create_dir_all;

use crate::digest::{SequenceCollection, SequenceCollectionMetadata, SequenceRecord};

/// User-facing store with lazy-loading read methods.
///
/// Wraps `ReadonlyRefgetStore` and provides `&mut self` read methods that
/// automatically load data on first access. Use `into_readonly()` to convert
/// to an immutable store suitable for `Arc<ReadonlyRefgetStore>` in servers.
#[derive(Debug)]
pub struct RefgetStore {
    pub(crate) inner: ReadonlyRefgetStore,
}

impl std::ops::Deref for RefgetStore {
    type Target = ReadonlyRefgetStore;
    fn deref(&self) -> &ReadonlyRefgetStore {
        &self.inner
    }
}

impl Display for RefgetStore {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.inner)
    }
}

impl RefgetStore {
    // =====================================================================
    // Constructors
    // =====================================================================

    /// Check whether a valid RefgetStore exists at the given path.
    pub fn store_exists<P: AsRef<Path>>(path: P) -> bool {
        ReadonlyRefgetStore::store_exists(path)
    }

    /// Create a disk-backed RefgetStore.
    pub fn on_disk<P: AsRef<Path>>(cache_path: P) -> Result<Self> {
        let cache_path = cache_path.as_ref();
        let index_path = cache_path.join("rgstore.json");

        if index_path.exists() {
            Self::open_local(cache_path)
        } else {
            let mode = StorageMode::Encoded;
            create_dir_all(cache_path)?;
            let mut inner = ReadonlyRefgetStore::new(mode);
            inner.local_path = Some(cache_path.to_path_buf());
            inner.seqdata_path_template = Some(DEFAULT_SEQDATA_PATH_TEMPLATE.to_string());
            inner.persist_to_disk = true;
            create_dir_all(cache_path.join("sequences"))?;
            create_dir_all(cache_path.join("collections"))?;
            create_dir_all(cache_path.join("fhr"))?;
            Ok(Self { inner })
        }
    }

    /// Create an in-memory RefgetStore.
    pub fn in_memory() -> Self {
        Self {
            inner: ReadonlyRefgetStore::new(StorageMode::Encoded),
        }
    }

    /// Open a local RefgetStore from a directory.
    pub fn open_local<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self {
            inner: ReadonlyRefgetStore::open_local(path)?,
        })
    }

    /// Open a remote RefgetStore with local caching.
    pub fn open_remote<P: AsRef<Path>, S: AsRef<str>>(
        cache_path: P,
        remote_url: S,
    ) -> Result<Self> {
        Ok(Self {
            inner: ReadonlyRefgetStore::open_remote(cache_path, remote_url)?,
        })
    }

    /// Convert to a ReadonlyRefgetStore for concurrent read access.
    pub fn into_readonly(self) -> ReadonlyRefgetStore {
        self.inner
    }

    // =====================================================================
    // Preload methods (delegate to inner)
    // =====================================================================

    /// Preload all collections (delegates to inner).
    pub fn load_all_collections(&mut self) -> Result<()> {
        self.inner.load_all_collections()
    }

    /// Preload all sequences (delegates to inner).
    pub fn load_all_sequences(&mut self) -> Result<()> {
        self.inner.load_all_sequences()
    }

    /// Load a single collection by digest.
    pub fn load_collection(&mut self, digest: &str) -> Result<()> {
        self.inner.load_collection(digest)
    }

    /// Load a single sequence by digest.
    pub fn load_sequence(&mut self, digest: &str) -> Result<()> {
        self.inner.load_sequence(digest)
    }

    // =====================================================================
    // Lazy-loading read methods
    // =====================================================================

    /// Lazy-loading get_collection. Loads on first access.
    pub fn get_collection(&mut self, digest: &str) -> Result<SequenceCollection> {
        if let Ok(coll) = self.inner.get_collection(digest) {
            return Ok(coll);
        }
        self.inner.load_collection(digest)?;
        self.inner.get_collection(digest)
    }

    /// Lazy-loading get_collection_level2.
    pub fn get_collection_level2(&mut self, digest: &str) -> Result<crate::digest::CollectionLevel2> {
        if let Ok(lvl2) = self.inner.get_collection_level2(digest) {
            return Ok(lvl2);
        }
        self.inner.load_collection(digest)?;
        self.inner.get_collection_level2(digest)
    }

    /// Lazy-loading compare.
    pub fn compare(&mut self, digest_a: &str, digest_b: &str) -> Result<crate::digest::SeqColComparison> {
        if !self.inner.is_collection_loaded(digest_a) {
            self.inner.load_collection(digest_a)?;
        }
        if !self.inner.is_collection_loaded(digest_b) {
            self.inner.load_collection(digest_b)?;
        }
        self.inner.compare(digest_a, digest_b)
    }

    /// Lazy-loading get_attribute.
    pub fn get_attribute(
        &mut self,
        attr_name: &str,
        attr_digest: &str,
    ) -> Result<Option<serde_json::Value>> {
        let collections = self.inner.find_collections_by_attribute(attr_name, attr_digest)?;
        if collections.is_empty() {
            return Ok(None);
        }
        if !self.inner.is_collection_loaded(&collections[0]) {
            self.inner.load_collection(&collections[0])?;
        }
        self.inner.get_attribute(attr_name, attr_digest)
    }

    // =====================================================================
    // Write/mutation methods (delegate to inner)
    // =====================================================================

    /// Set whether to suppress progress output.
    pub fn set_quiet(&mut self, quiet: bool) {
        self.inner.set_quiet(quiet);
    }

    /// Change the storage mode.
    pub fn set_encoding_mode(&mut self, new_mode: StorageMode) {
        self.inner.set_encoding_mode(new_mode);
    }

    /// Enable 2-bit encoding for space efficiency.
    pub fn enable_encoding(&mut self) {
        self.inner.enable_encoding();
    }

    /// Disable encoding, use raw byte storage.
    pub fn disable_encoding(&mut self) {
        self.inner.disable_encoding();
    }

    /// Enable disk persistence for this store.
    pub fn enable_persistence<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        self.inner.enable_persistence(path)
    }

    /// Disable disk persistence for this store.
    pub fn disable_persistence(&mut self) {
        self.inner.disable_persistence();
    }

    /// Add a sequence to the store.
    pub fn add_sequence<T: Into<Option<[u8; 32]>>>(
        &mut self,
        sequence_record: SequenceRecord,
        collection_digest: T,
        force: bool,
    ) -> Result<()> {
        self.inner.add_sequence(sequence_record, collection_digest, force)
    }

    /// Add a collection and all sequences in it to the store.
    pub fn add_sequence_collection(&mut self, collection: SequenceCollection) -> Result<()> {
        self.inner.add_sequence_collection(collection)
    }

    /// Add a collection, overwriting existing data.
    pub fn add_sequence_collection_force(&mut self, collection: SequenceCollection) -> Result<()> {
        self.inner.add_sequence_collection_force(collection)
    }

    /// Add a sequence collection from a FASTA file.
    pub fn add_sequence_collection_from_fasta<P: AsRef<Path>>(
        &mut self,
        file_path: P,
        opts: FastaImportOptions<'_>,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
        self.inner.add_sequence_collection_from_fasta(file_path, opts)
    }

    /// Add a SequenceRecord directly to the store.
    pub fn add_sequence_record(&mut self, sr: SequenceRecord, force: bool) -> Result<()> {
        self.inner.add_sequence_record(sr, force)
    }

    /// Remove a collection from the store.
    pub fn remove_collection(&mut self, digest: &str, remove_orphan_sequences: bool) -> Result<bool> {
        self.inner.remove_collection(digest, remove_orphan_sequences)
    }

    /// Import a collection (with sequences, aliases, FHR) from another store.
    /// The source must already have the collection loaded.
    pub fn import_collection_from_readonly(&mut self, source: &ReadonlyRefgetStore, digest: &str) -> Result<()> {
        self.inner.import_collection(source, digest)
    }

    /// Import a collection from another RefgetStore, lazy-loading the source if needed.
    pub fn import_collection(&mut self, source: &mut RefgetStore, digest: &str) -> Result<()> {
        if !source.inner.is_collection_loaded(digest) {
            source.inner.load_collection(digest)?;
        }
        self.inner.import_collection(&source.inner, digest)
    }

    /// Ensure a sequence is decoded into the decoded cache.
    pub fn ensure_decoded<K: AsRef<[u8]>>(&mut self, seq_digest: K) -> Result<()> {
        self.inner.ensure_decoded(seq_digest)
    }

    /// Clear the decoded sequence cache.
    pub fn clear_decoded_cache(&mut self) {
        self.inner.clear_decoded_cache();
    }

    /// Clear in-memory sequence data and decoded cache, preserving metadata.
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    // --- Seqcol config methods ---

    /// Enable computation of ancillary digests.
    pub fn enable_ancillary_digests(&mut self) {
        self.inner.enable_ancillary_digests();
    }

    /// Disable computation of ancillary digests.
    pub fn disable_ancillary_digests(&mut self) {
        self.inner.disable_ancillary_digests();
    }

    /// Enable indexed attribute lookup.
    pub fn enable_attribute_index(&mut self) {
        self.inner.enable_attribute_index();
    }

    /// Disable indexed attribute lookup.
    pub fn disable_attribute_index(&mut self) {
        self.inner.disable_attribute_index();
    }

    // --- Write methods ---

    /// Write the store using its configured paths.
    pub fn write(&self) -> Result<()> {
        self.inner.write()
    }

    /// Write the store to a directory.
    pub fn write_store_to_dir<P: AsRef<Path>>(
        &self,
        root_path: P,
        seqdata_path_template: Option<&str>,
    ) -> Result<()> {
        self.inner.write_store_to_dir(root_path, seqdata_path_template)
    }

    /// Write all sequence metadata to an RGSI file.
    pub fn write_sequences_rgsi<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        self.inner.write_sequences_rgsi(file_path)
    }
}
