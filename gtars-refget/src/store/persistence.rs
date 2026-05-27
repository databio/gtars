//! Disk I/O for RefgetStore: reading and writing index files,
//! opening stores from disk, and loading sequences/collections.

use super::*;
use super::readonly::ReadonlyRefgetStore;
use super::fhr_metadata;

use std::collections::HashMap;
use std::ffi::OsStr;

use indexmap::IndexMap;
use std::fs::{self, File, create_dir_all};
use std::io::{BufRead, Write};
use std::path::Path;

use anyhow::{Context, Result};
use sha2::{Sha256, Digest};

use crate::collection::{
    SequenceCollectionRecordExt,
    read_rgsi_file,
};
use crate::digest::{
    SequenceCollectionRecord, SequenceMetadata, SequenceRecord,
    parse_rgci_line, parse_rgsi_line,
};
use crate::hashkeyable::HashKeyable;

use chrono::Utc;

// ============================================================================
// ReadonlyRefgetStore disk I/O methods
// ============================================================================

impl ReadonlyRefgetStore {
    /// Write a single sequence to disk using the configured path template
    pub(crate) fn write_sequence_to_disk_single(
        &self,
        metadata: &SequenceMetadata,
        sequence: &[u8],
    ) -> Result<()> {
        let template = self
            .seqdata_path_template
            .as_ref()
            .context("seqdata_path_template not set")?;
        let local_path = self.local_path.as_ref().context("local_path not set")?;

        let seq_file_path = Self::expand_template(&metadata.sha512t24u, template);
        let full_path = local_path.join(&seq_file_path);

        if let Some(parent) = full_path.parent() {
            create_dir_all(parent)?;
        }

        let mut file = File::create(&full_path)?;
        file.write_all(sequence)?;

        Ok(())
    }

    /// Write a single collection RGSI file to disk.
    /// Used when persist_to_disk=true to persist collections incrementally.
    pub(crate) fn write_collection_to_disk_single(&self, record: &SequenceCollectionRecord) -> Result<()> {
        let local_path = self.local_path.as_ref().context("local_path not set")?;

        let coll_file_path = format!("collections/{}.rgsi", record.metadata().digest);
        let full_path = local_path.join(&coll_file_path);

        if let Some(parent) = full_path.parent() {
            create_dir_all(parent)?;
        }

        record.write_collection_rgsi(&full_path)?;

        Ok(())
    }

    /// Write index files (sequences.rgsi, collections.rgci, and rgstore.json) to disk.
    ///
    /// This allows the store to be loaded later via open_local().
    /// Called automatically when adding collections in disk-backed mode.
    pub(crate) fn write_index_files(&self) -> Result<()> {
        let local_path = self.local_path.as_ref().context("local_path not set")?;
        let template = self
            .seqdata_path_template
            .as_ref()
            .context("seqdata_path_template not set")?;

        let sequence_index_path = local_path.join("sequences.rgsi");
        self.write_sequences_rgsi(&sequence_index_path)?;

        let collection_index_path = local_path.join("collections.rgci");
        self.write_collections_rgci(&collection_index_path)?;

        // Compute digests of the files we just wrote
        let sequences_digest = Self::sha256_file(&sequence_index_path).ok();
        let collections_digest = Self::sha256_file(&collection_index_path).ok();
        let aliases_digest = self.compute_aliases_digest();
        let fhr_digest = self.compute_fhr_digest();

        self.write_rgstore_json_with_digests(
            local_path, template,
            collections_digest, sequences_digest,
            aliases_digest, fhr_digest,
        )?;

        Ok(())
    }

    /// Write the rgstore.json metadata file to the given directory.
    pub(crate) fn write_rgstore_json(&self, dir: &Path, seqdata_template: &str) -> Result<()> {
        self.write_rgstore_json_with_digests(dir, seqdata_template, None, None, None, None)
    }

    /// Write the rgstore.json metadata file with state digests and modified timestamp.
    pub(crate) fn write_rgstore_json_with_digests(
        &self,
        dir: &Path,
        seqdata_template: &str,
        collections_digest: Option<String>,
        sequences_digest: Option<String>,
        aliases_digest: Option<String>,
        fhr_digest: Option<String>,
    ) -> Result<()> {
        let metadata = StoreMetadata {
            version: 1,
            seqdata_path_template: seqdata_template.to_string(),
            collections_path_template: "collections/%s.rgsi".to_string(),
            sequence_index: "sequences.rgsi".to_string(),
            collection_index: Some("collections.rgci".to_string()),
            mode: self.mode,
            created_at: Utc::now().to_rfc3339(),
            ancillary_digests: self.ancillary_digests,
            attribute_index: self.attribute_index,
            sequence_alias_namespaces: self.aliases.sequence_namespaces(),
            collection_alias_namespaces: self.aliases.collection_namespaces(),
            modified: Some(Utc::now().to_rfc3339()),
            collections_digest,
            sequences_digest,
            aliases_digest,
            fhr_digest,
        };

        let json = serde_json::to_string_pretty(&metadata)
            .context("Failed to serialize metadata to JSON")?;
        fs::write(dir.join("rgstore.json"), json).context("Failed to write rgstore.json")?;

        Ok(())
    }

    /// Write collection metadata index (collections.rgci) to disk.
    ///
    /// Creates a master index of all collections with their metadata.
    pub(crate) fn write_collections_rgci<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        let file_path = file_path.as_ref();
        let mut file = File::create(file_path)?;

        writeln!(
            file,
            "#digest\tn_sequences\tnames_digest\tsequences_digest\tlengths_digest\tname_length_pairs_digest\tsorted_name_length_pairs_digest\tsorted_sequences_digest"
        )?;

        for record in self.collections.values() {
            let meta = record.metadata();
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                meta.digest,
                meta.n_sequences,
                meta.names_digest,
                meta.sequences_digest,
                meta.lengths_digest,
                meta.name_length_pairs_digest.as_deref().unwrap_or(""),
                meta.sorted_name_length_pairs_digest.as_deref().unwrap_or(""),
                meta.sorted_sequences_digest.as_deref().unwrap_or(""),
            )?;
        }
        Ok(())
    }

    /// Write all sequence metadata to an RGSI file.
    pub fn write_sequences_rgsi<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        let file_path = file_path.as_ref();
        let mut file = std::fs::File::create(file_path)?;

        writeln!(
            file,
            "#name\tlength\talphabet\tsha512t24u\tmd5\tdescription"
        )?;

        // Sort by sha512t24u digest for deterministic output (sequence_store is a HashMap).
        let mut entries: Vec<&SequenceRecord> = self.sequence_store.values().collect();
        entries.sort_by(|a, b| a.metadata().sha512t24u.cmp(&b.metadata().sha512t24u));

        for result_sr in entries {
            let result = result_sr.metadata();
            let description = result.description.as_deref().unwrap_or("");
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}",
                result.name,
                result.length,
                result.alphabet,
                result.sha512t24u,
                result.md5,
                description
            )?;
        }
        Ok(())
    }

    /// Read the store metadata from rgstore.json, returning state digests and timestamp.
    ///
    /// Returns a HashMap with keys: modified, collections_digest, sequences_digest,
    /// aliases_digest, fhr_digest. Missing values are omitted from the map.
    pub fn store_metadata(&self) -> Result<HashMap<String, String>> {
        let local_path = self.local_path.as_ref().context("local_path not set")?;
        let json = fs::read_to_string(local_path.join("rgstore.json"))
            .context("Failed to read rgstore.json")?;
        let metadata: StoreMetadata =
            serde_json::from_str(&json).context("Failed to parse store metadata")?;

        let mut map = HashMap::new();
        if let Some(v) = metadata.modified { map.insert("modified".to_string(), v); }
        if let Some(v) = metadata.collections_digest { map.insert("collections_digest".to_string(), v); }
        if let Some(v) = metadata.sequences_digest { map.insert("sequences_digest".to_string(), v); }
        if let Some(v) = metadata.aliases_digest { map.insert("aliases_digest".to_string(), v); }
        if let Some(v) = metadata.fhr_digest { map.insert("fhr_digest".to_string(), v); }
        Ok(map)
    }

    /// Compute SHA256 digest of a file's contents.
    fn sha256_file(path: &Path) -> Result<String> {
        let bytes = fs::read(path)?;
        let hash = Sha256::digest(&bytes);
        Ok(format!("{:x}", hash))
    }

    /// Compute a combined SHA256 digest of all alias namespace files (sorted).
    /// Alias files live under aliases/sequences/*.tsv and aliases/collections/*.tsv.
    fn compute_aliases_digest(&self) -> Option<String> {
        let local_path = self.local_path.as_ref()?;
        let aliases_dir = local_path.join("aliases");
        if !aliases_dir.exists() { return None; }

        let mut paths: Vec<_> = Vec::new();
        for subdir in &["sequences", "collections"] {
            let sub = aliases_dir.join(subdir);
            if sub.exists() {
                if let Ok(entries) = fs::read_dir(&sub) {
                    for entry in entries.filter_map(|e| e.ok()) {
                        let p = entry.path();
                        if p.extension().map_or(false, |e| e == "tsv") {
                            paths.push(p);
                        }
                    }
                }
            }
        }
        if paths.is_empty() { return None; }
        paths.sort();

        let mut hasher = Sha256::new();
        for path in &paths {
            hasher.update(fs::read(path).ok()?);
        }
        Some(format!("{:x}", hasher.finalize()))
    }

    /// Compute a combined SHA256 digest of all FHR sidecar files (sorted).
    fn compute_fhr_digest(&self) -> Option<String> {
        let local_path = self.local_path.as_ref()?;
        let fhr_dir = local_path.join("fhr");
        if !fhr_dir.exists() { return None; }

        let mut paths: Vec<_> = fs::read_dir(&fhr_dir).ok()?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .collect();
        if paths.is_empty() { return None; }
        paths.sort();

        let mut hasher = Sha256::new();
        for path in paths {
            hasher.update(fs::read(&path).ok()?);
        }
        Some(format!("{:x}", hasher.finalize()))
    }

    // =========================================================================
    // Open methods
    // =========================================================================

    /// Open a local store (internal). Users should use RefgetStore::open_local().
    pub(crate) fn open_local<P: AsRef<Path>>(path: P) -> Result<Self> {
        let root_path = path.as_ref();

        let index_path = root_path.join("rgstore.json");
        let json = fs::read_to_string(&index_path).context(format!(
            "Failed to read rgstore.json from {}",
            index_path.display()
        ))?;

        let metadata: StoreMetadata =
            serde_json::from_str(&json).context("Failed to parse store metadata")?;

        Self::sanitize_relative_path(&metadata.seqdata_path_template)?;
        Self::sanitize_relative_path(&metadata.sequence_index)?;
        if let Some(ref ci) = metadata.collection_index {
            Self::sanitize_relative_path(ci)?;
        }

        let mut store = ReadonlyRefgetStore::new(metadata.mode);
        store.local_path = Some(root_path.to_path_buf());
        store.seqdata_path_template = Some(metadata.seqdata_path_template.clone());
        store.persist_to_disk = true;
        store.ancillary_digests = metadata.ancillary_digests;
        store.attribute_index = metadata.attribute_index;

        let sequence_index_path = root_path.join(&metadata.sequence_index);
        if sequence_index_path.exists() {
            Self::load_sequences_from_index(&mut store, &sequence_index_path)?;
        }

        if let Some(ref collection_index) = metadata.collection_index {
            let collection_index_path = root_path.join(collection_index);
            if collection_index_path.exists() {
                Self::load_collection_stubs_from_rgci(&mut store, &collection_index_path)?;
            }
        }

        if store.collections.is_empty() {
            let collections_dir = root_path.join("collections");
            Self::load_collections_from_directory(&mut store, &collections_dir)?;
        }

        let aliases_dir = root_path.join("aliases");
        store.aliases.load_from_dir(&aliases_dir)?;

        store.available_sequence_alias_namespaces = if metadata.sequence_alias_namespaces.is_empty() {
            store.aliases.sequence_namespaces()
        } else {
            metadata.sequence_alias_namespaces
        };
        store.available_collection_alias_namespaces = if metadata.collection_alias_namespaces.is_empty() {
            store.aliases.collection_namespaces()
        } else {
            metadata.collection_alias_namespaces
        };

        store.fhr_metadata =
            fhr_metadata::load_sidecars(&root_path.join("fhr"));

        Ok(store)
    }

    /// Open a remote store (internal). Users should use RefgetStore::open_remote().
    pub(crate) fn open_remote<P: AsRef<Path>, S: AsRef<str>>(
        cache_path: P,
        remote_url: S,
    ) -> Result<Self> {
        let cache_path = cache_path.as_ref();
        let remote_url = remote_url.as_ref().to_string();

        create_dir_all(cache_path)?;

        let index_data = Self::fetch_file(
            &Some(cache_path.to_path_buf()),
            &Some(remote_url.clone()),
            "rgstore.json",
            true,
            false,
        )?;

        let json =
            String::from_utf8(index_data).context("Store metadata contains invalid UTF-8")?;

        let metadata: StoreMetadata =
            serde_json::from_str(&json).context("Failed to parse store metadata")?;

        Self::sanitize_relative_path(&metadata.seqdata_path_template)?;
        Self::sanitize_relative_path(&metadata.sequence_index)?;
        if let Some(ref ci) = metadata.collection_index {
            Self::sanitize_relative_path(ci)?;
        }

        let mut store = ReadonlyRefgetStore::new(metadata.mode);
        store.local_path = Some(cache_path.to_path_buf());
        store.remote_source = Some(remote_url.clone());
        store.seqdata_path_template = Some(metadata.seqdata_path_template.clone());
        store.persist_to_disk = true;
        store.ancillary_digests = metadata.ancillary_digests;
        store.attribute_index = metadata.attribute_index;
        store.available_sequence_alias_namespaces = metadata.sequence_alias_namespaces;
        store.available_collection_alias_namespaces = metadata.collection_alias_namespaces;

        // Defer sequence index loading — it can be 66+ MB and is only needed
        // when accessing individual sequences, not for browsing collections.
        store.sequence_index_loaded = false;
        store.sequence_index_path = Some(metadata.sequence_index.clone());

        if let Some(ref collection_index) = metadata.collection_index {
            if let Ok(collection_index_data) = Self::fetch_file(
                &Some(cache_path.to_path_buf()),
                &Some(remote_url.clone()),
                collection_index,
                true,
                false,
            ) {
                let collection_index_str = String::from_utf8(collection_index_data)
                    .context("collection index contains invalid UTF-8")?;

                Self::load_collection_stubs_from_reader(
                    &mut store,
                    collection_index_str.as_bytes(),
                )?;
            }
        }

        if store.collections.is_empty() {
            let local_collections_dir = cache_path.join("collections");
            create_dir_all(&local_collections_dir)?;
            Self::load_collections_from_directory(&mut store, &local_collections_dir)?;
        }

        Ok(store)
    }

    // =========================================================================
    // Loading helpers
    // =========================================================================

    /// Parse RGSI lines from a reader and load as Stub sequence records.
    pub(crate) fn load_sequences_from_reader<R: BufRead>(store: &mut ReadonlyRefgetStore, reader: R) -> Result<()> {
        for line in reader.lines() {
            let line = line?;

            if line.starts_with('#') {
                continue;
            }

            if let Some(seq_metadata) = parse_rgsi_line(&line) {
                let record = SequenceRecord::Stub(seq_metadata.clone());

                let sha512_key = seq_metadata.sha512t24u.to_key();
                store.sequence_store.insert(sha512_key, record);

                let md5_key = seq_metadata.md5.to_key();
                store.md5_lookup.insert(md5_key, sha512_key);
            }
        }

        Ok(())
    }

    /// Load sequence metadata from a sequence index file (sequences.rgsi).
    pub(crate) fn load_sequences_from_index(store: &mut ReadonlyRefgetStore, index_path: &Path) -> Result<()> {
        let file = std::fs::File::open(index_path)?;
        let reader = std::io::BufReader::new(file);
        Self::load_sequences_from_reader(store, reader)
    }

    /// Parse RGCI lines from a reader and load as Stub collection records.
    pub(crate) fn load_collection_stubs_from_reader<R: BufRead>(
        store: &mut ReadonlyRefgetStore,
        reader: R,
    ) -> Result<()> {
        for line in reader.lines() {
            let line = line?;

            if let Some(metadata) = parse_rgci_line(&line) {
                let key = metadata.digest.to_key();
                store
                    .collections
                    .insert(key, SequenceCollectionRecord::Stub(metadata));
            }
        }

        Ok(())
    }

    /// Load collection stubs from collections.rgci index file (new format).
    pub(crate) fn load_collection_stubs_from_rgci(store: &mut ReadonlyRefgetStore, index_path: &Path) -> Result<()> {
        let file = std::fs::File::open(index_path)?;
        let reader = std::io::BufReader::new(file);
        Self::load_collection_stubs_from_reader(store, reader)
    }

    /// Load full collections from a collections directory (fallback when no RGCI exists).
    pub(crate) fn load_collections_from_directory(
        store: &mut ReadonlyRefgetStore,
        collections_dir: &Path,
    ) -> Result<()> {
        if !collections_dir.exists() {
            return Ok(());
        }

        for entry in fs::read_dir(collections_dir)? {
            let entry = entry?;
            let path = entry.path();

            if path.is_file() && path.extension() == Some(OsStr::new("rgsi")) {
                let collection = read_rgsi_file(&path)?;
                let collection_digest = collection.metadata.digest.to_key();

                let mut name_map = IndexMap::new();
                for sequence_record in &collection.sequences {
                    let metadata = sequence_record.metadata();
                    name_map.insert(metadata.name.clone(), metadata.sha512t24u.to_key());
                }
                store.name_lookup.insert(collection_digest, name_map);

                let record = SequenceCollectionRecord::from(collection);
                store.collections.insert(collection_digest, record);
            }
        }

        Ok(())
    }
}
