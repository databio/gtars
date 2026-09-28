//! RefgetStore wrapper struct with lazy-loading read methods.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::fmt::{Display, Formatter};
use std::path::Path;

use anyhow::Result;
use std::fs::create_dir_all;

use crate::digest::{SequenceCollection, SequenceRecord};
// Only the filesystem FASTA-import wrappers use this return type.
#[cfg(feature = "filesystem")]
use crate::digest::SequenceCollectionMetadata;

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

    /// Load the sequence index if not already loaded.
    /// Servers should call this during startup to preload the index.
    pub fn load_sequence_index(&mut self) -> Result<()> {
        self.inner.load_sequence_index()
    }

    /// Preload all collections (delegates to inner).
    pub fn load_all_collections(&mut self) -> Result<()> {
        self.inner.load_all_collections()
    }

    /// Preload all sequences (delegates to inner).
    pub fn load_all_sequences(&mut self) -> Result<()> {
        self.inner.ensure_sequence_index_loaded()?;
        self.inner.load_all_sequences()
    }

    /// Load a single collection by digest.
    pub fn load_collection(&mut self, digest: &str) -> Result<()> {
        self.inner.load_collection(digest)
    }

    /// Load a single sequence by digest.
    pub fn load_sequence(&mut self, digest: &str) -> Result<()> {
        self.inner.ensure_sequence_index_loaded()?;
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

    /// List sequence metadata, loading the deferred sequence index if needed.
    ///
    /// This only populates metadata-only records; sequence bodies remain lazy.
    pub fn list_sequences(&mut self) -> Result<Vec<crate::digest::SequenceMetadata>> {
        self.inner.ensure_sequence_index_loaded()?;
        Ok(self.inner.list_sequences())
    }

    /// Retrieve one sequence, loading only its body when necessary.
    ///
    /// The readonly lookup accepts MD5 aliases, but the disk loader needs the
    /// canonical SHA-512/24u key, so resolve before loading.
    pub fn get_sequence<K: AsRef<[u8]>>(&mut self, digest: K) -> Result<SequenceRecord> {
        self.inner.ensure_sequence_index_loaded()?;
        let canonical_digest = self
            .inner
            .get_sequence(digest.as_ref())
            .map(|record| record.metadata().sha512t24u.clone())?;
        if !self.inner.is_sequence_loaded(&canonical_digest) {
            self.inner.load_sequence(&canonical_digest)?;
        }
        Ok(self.inner.get_sequence(&canonical_digest)?.clone())
    }

    /// Retrieve a named sequence from a collection, loading only that
    /// collection and sequence body when necessary.
    pub fn get_sequence_by_name<K: AsRef<[u8]>>(
        &mut self,
        collection_digest: K,
        sequence_name: &str,
    ) -> Result<SequenceRecord> {
        let collection_digest = String::from_utf8_lossy(collection_digest.as_ref());
        if !self.inner.is_collection_loaded(collection_digest.as_bytes()) {
            self.inner.load_collection(collection_digest.as_ref())?;
        }
        let canonical_digest = self
            .inner
            .get_sequence_by_name(collection_digest.as_bytes(), sequence_name)
            .map(|record| record.metadata().sha512t24u.clone())?;
        if !self.inner.is_sequence_loaded(&canonical_digest) {
            self.inner.load_sequence(&canonical_digest)?;
        }
        self.inner
            .get_sequence_by_name(collection_digest.as_bytes(), sequence_name)
    }

    /// Lazily load two collection records and match their local names by content.
    pub fn match_sequence_names(
        &mut self,
        digest_a: &str,
        digest_b: &str,
    ) -> Result<crate::digest::CollectionNameMatch> {
        if !self.inner.is_collection_loaded(digest_a) {
            self.inner.load_collection(digest_a)?;
        }
        if !self.inner.is_collection_loaded(digest_b) {
            self.inner.load_collection(digest_b)?;
        }
        self.inner.match_sequence_names(digest_a, digest_b)
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

    /// Verify sequences against their stored digests. `digests: None` checks
    /// every sequence in the store.
    ///
    /// See [`ReadonlyRefgetStore::verify_sequences`] for what this actually
    /// checks and how parallelism is controlled by `opts.jobs`.
    pub fn verify_sequences(
        &mut self,
        digests: Option<&[String]>,
        opts: &VerifyOptions,
    ) -> Result<VerifyReport> {
        self.inner.ensure_sequence_index_loaded()?;
        self.inner.verify_sequences(digests, opts)
    }

    /// Verify only the sequences belonging to one collection.
    ///
    /// Loads the collection's metadata (not sequence bytes) via the usual
    /// lazy-load path, then verifies each of its (deduped) sequence digests.
    pub fn verify_collection(
        &mut self,
        collection_digest: &str,
        opts: &VerifyOptions,
    ) -> Result<VerifyReport> {
        use crate::hashkeyable::HashKeyable;

        self.inner.ensure_sequence_index_loaded()?;
        let collection_key = collection_digest.to_key();
        self.inner.ensure_collection_loaded(&collection_key)?;
        let mut digests: Vec<String> = self
            .inner
            .collection_sequence_metadata(&collection_key)?
            .iter()
            .map(|record| record.metadata().sha512t24u.clone())
            .collect();
        digests.sort();
        digests.dedup();
        self.inner.verify_sequences(Some(&digests), opts)
    }

    /// Lazy-loading export_fasta. Loads the collection and only the sequence
    /// bodies it will write (every record whose name is requested, or all of
    /// the collection's when `sequence_names` is `None`), never the whole
    /// store. Name validation is left to the readonly implementation.
    #[cfg(feature = "filesystem")]
    pub fn export_fasta<K: AsRef<[u8]>, P: AsRef<Path>>(
        &mut self,
        collection_digest: K,
        output_path: P,
        sequence_names: Option<Vec<&str>>,
        line_width: Option<usize>,
    ) -> Result<()> {
        use crate::hashkeyable::HashKeyable;

        let collection_key = collection_digest.as_ref().to_key();
        self.inner.ensure_collection_loaded(&collection_key)?;
        let wanted: Option<std::collections::HashSet<&str>> =
            sequence_names.as_ref().map(|names| names.iter().copied().collect());
        let to_load: Vec<String> = self
            .inner
            .collection_sequence_metadata(&collection_key)?
            .iter()
            .map(|record| record.metadata())
            .filter(|meta| wanted.as_ref().is_none_or(|w| w.contains(meta.name.as_str())))
            .map(|meta| meta.sha512t24u.clone())
            .collect();
        for digest in &to_load {
            self.inner.load_sequence(digest)?;
        }

        self.inner
            .export_fasta(collection_digest, output_path, sequence_names, line_width)
    }

    /// Lazy-loading export_fasta_by_digests. Loads any requested sequence
    /// digest that is not yet [`Self::is_sequence_loaded`], then delegates to
    /// the readonly implementation. There is no collection context here, so
    /// each digest is loaded directly (no name resolution needed).
    #[cfg(feature = "filesystem")]
    pub fn export_fasta_by_digests<P: AsRef<Path>>(
        &mut self,
        seq_digests: Vec<&str>,
        output_path: P,
        line_width: Option<usize>,
    ) -> Result<()> {
        self.inner.ensure_sequence_index_loaded()?;
        for digest in seq_digests.iter() {
            if !self.inner.is_sequence_loaded(*digest) {
                self.inner.load_sequence(digest)?;
            }
        }

        self.inner
            .export_fasta_by_digests(seq_digests, output_path, line_width)
    }

    // =====================================================================
    // Write/mutation methods (delegate to inner)
    // =====================================================================

    /// Set whether to suppress progress output.
    pub fn set_quiet(&mut self, quiet: bool) {
        self.inner.set_quiet(quiet);
    }

    /// Change the storage mode.
    ///
    /// Switching to Encoded fails, leaving the store unchanged, if a sequence
    /// cannot be encoded with its stored alphabet. The store must then be
    /// re-imported.
    pub fn set_encoding_mode(&mut self, new_mode: StorageMode) -> Result<()> {
        self.inner.set_encoding_mode(new_mode)
    }

    /// Enable 2-bit encoding for space efficiency.
    pub fn enable_encoding(&mut self) -> Result<()> {
        self.inner.enable_encoding()
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
    pub fn add_sequence<T: Into<Option<DigestKey>>>(
        &mut self,
        sequence_record: SequenceRecord,
        collection_digest: T,
        force: bool,
    ) -> Result<()> {
        self.inner.add_sequence(sequence_record, collection_digest, force)
    }

    /// Digest and insert raw sequence bytes without a collection association.
    ///
    /// Unlike `add_sequence_record`, this method can reliably encode a
    /// single-base sequence because its input is known to be raw bytes.
    /// Returns the sequence's sha512t24u digest.
    pub fn ingest_sequence(&mut self, name: &str, bytes: &[u8]) -> Result<String> {
        self.inner.ingest_sequence(name, bytes)
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
    #[cfg(feature = "filesystem")]
    pub fn add_sequence_collection_from_fasta<P: AsRef<Path>>(
        &mut self,
        file_path: P,
        opts: FastaImportOptions<'_>,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
        self.inner.add_sequence_collection_from_fasta(file_path, opts)
    }

    /// Add sequence collections from multiple FASTA files, decoding up to
    /// `opts.file_jobs` files concurrently. Inserts in fixed input order so the
    /// resulting store is byte-identical to a serial build.
    ///
    /// Returns an [`ImportReport`] with per-file results and per-run ingest
    /// counters.
    #[cfg(feature = "filesystem")]
    pub fn add_sequence_collections_from_fastas(
        &mut self,
        files: &[std::path::PathBuf],
        opts: FastaImportOptions<'_>,
    ) -> Result<ImportReport> {
        self.inner.add_sequence_collections_from_fastas(files, opts)
    }

    /// Add a SequenceRecord directly to the store.
    pub fn add_sequence_record(&mut self, sr: SequenceRecord, force: bool) -> Result<()> {
        self.inner.add_sequence_record(sr, force)
    }

    /// Remove a collection from the store.
    pub fn remove_collection(&mut self, digest: &str, remove_orphan_sequences: bool) -> Result<bool> {
        self.inner.remove_collection(digest, remove_orphan_sequences)
    }

    /// Dry-run of the orphan cleanup in [`Self::remove_collection`].
    pub fn plan_orphan_removal(&self, digest: &str) -> Result<Vec<String>> {
        self.inner.plan_orphan_removal(digest)
    }

    // --- Write locking (delegates; `Deref` only yields `&`) ---

    /// Hold the store's exclusive writer lock across several mutations.
    pub fn lock_for_batch(&mut self, operation: &str) -> Result<()> {
        self.inner.lock_for_batch(operation)
    }

    /// Release a lock taken by [`Self::lock_for_batch`].
    pub fn release_batch_lock(&mut self) {
        self.inner.release_batch_lock();
    }

    /// Override the timeout/staleness settings used when acquiring the write lock.
    pub fn set_lock_options(&mut self, options: super::LockOptions) {
        self.inner.set_lock_options(options);
    }

    /// Allow a commit to overwrite an alias another writer already published.
    pub fn set_force_alias(&mut self, force: bool) {
        self.inner.set_force_alias(force);
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

    /// Clear in-memory sequence data, preserving metadata.
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    /// Check whether a sequence is loaded (Full).
    pub fn is_sequence_loaded<K: AsRef<[u8]>>(&self, seq_digest: K) -> bool {
        self.inner.is_sequence_loaded(seq_digest)
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
