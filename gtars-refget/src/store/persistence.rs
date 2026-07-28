//! Disk I/O for RefgetStore: reading and writing index files,
//! opening stores from disk, and loading sequences/collections.
//!
//! # Commit model
//!
//! [`ReadonlyRefgetStore::write_index_files`] is the ONE place the shared,
//! whole-file store artifacts are published. It is not a "serialize my memory to
//! disk" operation — it is a MERGE:
//!
//! ```text
//! committed = (rows currently on disk  ∪  rows in memory)  \  tombstones
//! ```
//!
//! This matters because a store handle holds a snapshot taken at `open_local`
//! time and nothing re-reads it afterwards. For a FASTA import that snapshot is
//! hours old. A writer that serialized its own map would silently drop every row
//! another writer committed in the meantime — which is exactly how a genome went
//! missing from a production store on 2026-07-23 (its `collections/<digest>.rgsi`
//! was on disk, but no index referenced it).
//!
//! The merge is well-defined because both index files are digest-keyed sets of
//! content-derived rows: the same key implies the same content, so a union cannot
//! conflict on anything load-bearing. See [`merge_sequence_rows`] and
//! [`merge_collection_rows`] for the two places it is not purely mechanical.
//!
//! Removal has to be explicit, or a blind union would resurrect whatever
//! `remove_collection` just deleted. Hence tombstones (see
//! [`super::readonly::Tombstones`]) — in-memory only, cleared on a successful
//! commit, never written to disk.
//!
//! Everything expensive stays OUTSIDE the lock: FASTA parsing, digesting, and
//! every per-`.seq` and per-`collections/*.rgsi` write. Those are
//! digest-addressed and concurrency-safe by construction (see
//! [`ReadonlyRefgetStore::write_seq_bytes_to_full_path`]). The exclusive section
//! is O(index size): parse two TSVs, union two hashmaps, publish. Seconds, not
//! hours.

use super::*;
use super::readonly::{ReadonlyRefgetStore, Tombstones};
use super::alias::AliasKind;
use super::atomic::atomic_write;
use super::fhr_metadata;
use super::lock::StoreLock;

use std::collections::{BTreeMap, HashMap, HashSet};
use std::ffi::OsStr;

use indexmap::IndexMap;
use std::fs::{self, File, create_dir_all};
use std::io::{BufRead, Write};
use std::path::Path;

use anyhow::{Context, Result};
use sha2::{Sha256, Digest};

use crate::collection::{
    SequenceCollectionRecordExt, SequenceMetadataExt,
    read_rgsi_file,
};
use crate::digest::{
    SequenceCollectionMetadata, SequenceCollectionRecord, SequenceMetadata, SequenceRecord,
    parse_rgci_line, parse_rgsi_line,
};
use crate::hashkeyable::{HashKeyable, key_to_digest_string};

use chrono::Utc;

/// Header line of a `sequences.rgsi` file.
const RGSI_HEADER: &str = "#name\tlength\talphabet\tsha512t24u\tmd5\tdescription";

/// Header line of a `collections.rgci` file.
const RGCI_HEADER: &str = "#digest\tn_sequences\tnames_digest\tsequences_digest\tlengths_digest\tname_length_pairs_digest\tsorted_name_length_pairs_digest\tsorted_sequences_digest";

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

    /// Write a `.seq` file for a single sequence, given an already-expanded full
    /// path, skipping `create_dir_all` when the shard directory is known to exist.
    ///
    /// Per-digest `.seq` files are independent, so this is safe to call
    /// concurrently from a writer pool (see [`import`](super::import)). The
    /// `created_shards` cache holds shard dirs already created in this run so the
    /// common case avoids a `create_dir_all` syscall per (tiny) sequence -- the
    /// throughput bottleneck for transcriptomes with hundreds of thousands of
    /// records.
    #[cfg_attr(not(feature = "filesystem"), allow(dead_code))]
    pub(crate) fn write_seq_bytes_to_full_path(
        full_path: &Path,
        sequence: &[u8],
        created_shards: &std::sync::Mutex<HashSet<std::path::PathBuf>>,
    ) -> Result<()> {
        if let Some(parent) = full_path.parent() {
            // Fast path: shard dir already created this run -> skip the syscall.
            let known = {
                let guard = created_shards.lock().unwrap();
                guard.contains(parent)
            };
            if !known {
                create_dir_all(parent)?;
                created_shards.lock().unwrap().insert(parent.to_path_buf());
            }
        }

        let mut file = File::create(full_path)?;
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

    // =========================================================================
    // The commit path
    // =========================================================================

    /// Acquire the store's exclusive writer lock for a commit, unless this store
    /// already holds one via [`Self::lock_for_batch`].
    ///
    /// Returns `Ok(None)` when the lock is already held (re-entrancy: a batch
    /// guard wrapping several mutations must not deadlock against the
    /// `write_index_files` each mutation triggers) or when the store is not
    /// disk-backed (nothing to serialize against).
    pub(crate) fn acquire_commit_lock(&self, operation: &str) -> Result<Option<StoreLock>> {
        if self.commit_lock.is_some() {
            return Ok(None);
        }
        let Some(local_path) = self.local_path.as_ref() else {
            return Ok(None);
        };
        Ok(Some(StoreLock::acquire(
            local_path,
            operation,
            self.lock_options.clone(),
        )?))
    }

    /// Commit the store's shared index/manifest/alias artifacts.
    ///
    /// This MERGES with what is currently on disk rather than overwriting from
    /// the open-time snapshot — see the module docs for why. The publish order is
    /// load-bearing: indexes and alias files first, `rgstore.json` LAST, because
    /// the manifest advertises digests of the files it describes. A reader that
    /// catches the manifest mid-commit would otherwise be told about bytes that
    /// have not landed.
    ///
    /// Called automatically when adding or removing collections in disk-backed
    /// mode, and by [`Self::write`].
    pub(crate) fn write_index_files(&self) -> Result<()> {
        let local_path = self
            .local_path
            .as_ref()
            .context("local_path not set")?
            .clone();
        let template = self
            .seqdata_path_template
            .as_ref()
            .context("seqdata_path_template not set")?
            .clone();

        let _guard = self.acquire_commit_lock("write_index_files")?;

        // Snapshot the tombstones under the store lock, so the merge and the
        // subsequent clear agree on exactly which rows were suppressed.
        let tombstones = self.tombstones.lock().unwrap().clone();

        let sequence_index_path = local_path.join("sequences.rgsi");
        let collection_index_path = local_path.join("collections.rgci");

        // (1) Re-read the CURRENT on-disk rows -- not the open-time snapshot --
        // and union them with memory, minus tombstones.
        let merged_sequences = self.merged_sequence_rows(&sequence_index_path, &tombstones)?;
        let merged_collections =
            self.merged_collection_rows(&collection_index_path, &tombstones)?;

        // (2) Publish the indexes atomically.
        write_sequences_rgsi_rows(&sequence_index_path, &merged_sequences)?;
        write_collections_rgci_rows(&collection_index_path, &merged_collections)?;

        // (3) Publish the alias TSVs (also merged, per namespace).
        self.commit_all_alias_namespaces(&local_path, &tombstones)?;

        // (4) Recompute manifest fields from what we just wrote -- NOT from
        // memory, which is a subset after a merge.
        let logical_bytes: u64 = merged_sequences
            .values()
            .map(|m| m.disk_size(&self.mode) as u64)
            .sum();

        let mut metadata = self.manifest_base(&local_path, &template)?;
        metadata.sequences_digest = Self::sha256_file(&sequence_index_path).ok();
        metadata.collections_digest = Self::sha256_file(&collection_index_path).ok();
        metadata.aliases_digest = self.compute_aliases_digest();
        metadata.fhr_digest = self.compute_fhr_digest();
        metadata.logical_sequence_bytes = Some(logical_bytes);
        let (seq_ns, coll_ns) = published_alias_namespaces(&local_path);
        metadata.sequence_alias_namespaces = seq_ns;
        metadata.collection_alias_namespaces = coll_ns;

        // (5) Manifest LAST.
        write_manifest(&local_path, &metadata)?;

        // (6) The tombstones have been applied to disk; forget them.
        self.tombstones.lock().unwrap().clear_matching(&tombstones);

        Ok(())
    }

    /// Union the on-disk `sequences.rgsi` rows with this store's in-memory rows,
    /// minus tombstoned digests.
    fn merged_sequence_rows(
        &self,
        index_path: &Path,
        tombstones: &Tombstones,
    ) -> Result<HashMap<DigestKey, SequenceMetadata>> {
        let mut merged = read_sequences_rgsi_rows(index_path)?;
        for (key, record) in &self.sequence_store {
            merge_sequence_rows(&mut merged, *key, record.metadata());
        }
        merged.retain(|k, _| !tombstones.sequences.contains(k));
        Ok(merged)
    }

    /// Union the on-disk `collections.rgci` rows with this store's in-memory
    /// rows, minus tombstoned digests.
    fn merged_collection_rows(
        &self,
        index_path: &Path,
        tombstones: &Tombstones,
    ) -> Result<HashMap<DigestKey, SequenceCollectionMetadata>> {
        let mut merged = read_collections_rgci_rows(index_path)?;
        for (key, record) in &self.collections {
            merge_collection_rows(&mut merged, *key, record.metadata());
        }
        merged.retain(|k, _| !tombstones.collections.contains(k));
        Ok(merged)
    }

    /// Build the manifest to publish, preserving fields that must survive a
    /// commit.
    ///
    /// `created_at` is read back from the existing manifest rather than stamped
    /// with `now()`. It previously meant "last written" (it was overwritten on
    /// every commit), which made it a duplicate of `modified`; since the commit
    /// re-reads the manifest under the lock anyway, it can mean what it says.
    fn manifest_base(&self, local_path: &Path, seqdata_template: &str) -> Result<StoreMetadata> {
        let now = Utc::now().to_rfc3339();
        let existing = read_manifest(local_path)?;
        let created_at = existing
            .as_ref()
            .map(|m| m.created_at.clone())
            .unwrap_or_else(|| now.clone());

        Ok(StoreMetadata {
            version: 1,
            seqdata_path_template: seqdata_template.to_string(),
            collections_path_template: "collections/%s.rgsi".to_string(),
            sequence_index: "sequences.rgsi".to_string(),
            collection_index: Some("collections.rgci".to_string()),
            mode: self.mode,
            created_at,
            ancillary_digests: self.ancillary_digests,
            attribute_index: self.attribute_index,
            sequence_alias_namespaces: self.aliases.sequence_namespaces(),
            collection_alias_namespaces: self.aliases.collection_namespaces(),
            modified: Some(now),
            collections_digest: None,
            sequences_digest: None,
            aliases_digest: None,
            fhr_digest: None,
            logical_sequence_bytes: None,
        })
    }

    /// Write a fresh rgstore.json for an EXPORT to `dir` (a directory this store
    /// is being copied into, not committed to). No merge, no lock: the caller is
    /// producing a new store image.
    pub(crate) fn write_rgstore_json(&self, dir: &Path, seqdata_template: &str) -> Result<()> {
        let now = Utc::now().to_rfc3339();
        let metadata = StoreMetadata {
            version: 1,
            seqdata_path_template: seqdata_template.to_string(),
            collections_path_template: "collections/%s.rgsi".to_string(),
            sequence_index: "sequences.rgsi".to_string(),
            collection_index: Some("collections.rgci".to_string()),
            mode: self.mode,
            created_at: now.clone(),
            ancillary_digests: self.ancillary_digests,
            attribute_index: self.attribute_index,
            sequence_alias_namespaces: self.aliases.sequence_namespaces(),
            collection_alias_namespaces: self.aliases.collection_namespaces(),
            modified: Some(now),
            collections_digest: None,
            sequences_digest: None,
            aliases_digest: None,
            fhr_digest: None,
            logical_sequence_bytes: Some(self.logical_sequence_bytes() as u64),
        };
        write_manifest(dir, &metadata)
    }

    /// Refresh only the alias-namespace lists (and aliases_digest) in an existing
    /// rgstore.json, preserving created_at and the index digests.
    ///
    /// This is the cheap manifest update used after an alias is added/removed: it
    /// re-reads the current manifest and rewrites only the alias-derived fields,
    /// so we do not rehash the (potentially 60+ MB) sequence/collection indexes
    /// on every alias mutation. Keeps `sequence_alias_namespaces` /
    /// `collection_alias_namespaces` in lock-step with what is actually on disk,
    /// so served stores stay self-describing. No-op if the manifest does not yet
    /// exist (the next full index write will include the namespaces).
    ///
    /// The namespace lists come from a DIRECTORY SCAN of what was just
    /// published, not from `self.aliases`. Writing memory's view here is how a
    /// namespace another writer created gets dropped from the manifest: the TSV
    /// stays on disk but stops being advertised, and over HTTP -- where you
    /// cannot list a directory -- that makes it unreachable.
    ///
    /// Callers must already hold the store lock.
    pub(crate) fn refresh_manifest_alias_namespaces(&self) -> Result<()> {
        let local_path = match self.local_path.as_ref() {
            Some(p) => p,
            None => return Ok(()),
        };
        let Some(mut metadata) = read_manifest(local_path)? else {
            return Ok(());
        };

        let (seq_ns, coll_ns) = published_alias_namespaces(local_path);
        metadata.sequence_alias_namespaces = seq_ns;
        metadata.collection_alias_namespaces = coll_ns;
        metadata.aliases_digest = self.compute_aliases_digest();
        metadata.modified = Some(Utc::now().to_rfc3339());

        write_manifest(local_path, &metadata)
    }

    // =========================================================================
    // Whole-store serialization (export paths -- no merge)
    // =========================================================================

    /// Write this store's in-memory collection metadata to a `.rgci` file.
    ///
    /// EXPORT ONLY: serializes exactly what is in memory. The commit path uses
    /// [`write_collections_rgci_rows`] with merged rows instead — do not route a
    /// commit through here, or rows another writer added are lost.
    pub(crate) fn write_collections_rgci<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        let rows: HashMap<DigestKey, SequenceCollectionMetadata> = self
            .collections
            .iter()
            .map(|(k, v)| (*k, v.metadata().clone()))
            .collect();
        write_collections_rgci_rows(file_path.as_ref(), &rows)
    }

    /// Write this store's in-memory sequence metadata to an `.rgsi` file.
    ///
    /// EXPORT ONLY — see [`Self::write_collections_rgci`].
    pub fn write_sequences_rgsi<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        let rows: HashMap<DigestKey, SequenceMetadata> = self
            .sequence_store
            .iter()
            .map(|(k, v)| (*k, v.metadata().clone()))
            .collect();
        write_sequences_rgsi_rows(file_path.as_ref(), &rows)
    }

    // =========================================================================
    // Alias commit
    // =========================================================================

    /// Merge and publish every alias namespace this store or the directory knows
    /// about. Caller must hold the store lock.
    fn commit_all_alias_namespaces(
        &self,
        local_path: &Path,
        tombstones: &Tombstones,
    ) -> Result<()> {
        let aliases_dir = local_path.join("aliases");
        for kind in [AliasKind::Sequence, AliasKind::Collection] {
            let mut namespaces: HashSet<String> = self
                .aliases
                .namespaces_for(kind)
                .into_iter()
                .collect();
            namespaces.extend(scan_alias_namespaces(&aliases_dir.join(kind.subdir())));
            namespaces.extend(
                tombstones
                    .aliases
                    .iter()
                    .filter(|(k, _, _)| *k == kind)
                    .map(|(_, ns, _)| ns.clone()),
            );
            namespaces.extend(
                tombstones
                    .emptied_alias_namespaces
                    .iter()
                    .filter(|(k, _)| *k == kind)
                    .map(|(_, ns)| ns.clone()),
            );

            for ns in namespaces {
                self.commit_alias_namespace(&aliases_dir, kind, &ns, tombstones)?;
            }
        }
        Ok(())
    }

    /// Merge one alias namespace (disk ∪ memory, minus tombstones) and publish
    /// it atomically. Deletes the TSV only when the MERGED result is empty.
    ///
    /// The old behavior — delete the TSV whenever the namespace was absent from
    /// memory — is what orphaned namespaces: `open_local` loads only the
    /// namespaces the manifest advertises, so a writer holding a stale manifest
    /// snapshot simply did not have another writer's brand-new namespace in
    /// memory, and deleted it.
    pub(crate) fn commit_alias_namespace(
        &self,
        aliases_dir: &Path,
        kind: AliasKind,
        namespace: &str,
        tombstones: &Tombstones,
    ) -> Result<()> {
        let tsv_path = aliases_dir.join(kind.subdir()).join(format!("{}.tsv", namespace));

        let disk = read_alias_tsv(&tsv_path)?;
        let memory = self.aliases.namespace_map(kind, namespace);

        let mut merged = disk;
        for (alias, digest) in memory {
            match merged.get(&alias) {
                Some(existing) if *existing != digest => {
                    // A genuine semantic conflict: two writers bound the same
                    // human-readable name to different content. Silently picking
                    // one is the failure class that lost a genome in July.
                    if !self.force_alias {
                        return Err(anyhow::anyhow!(
                            "alias conflict in {} namespace '{}': '{}' is already published as {} \
                             on disk but this writer has it as {}. Another writer bound it first. \
                             Re-run with force-alias enabled (CLI: --force-alias) to overwrite the \
                             published value.",
                            kind.subdir(),
                            namespace,
                            alias,
                            key_to_digest_string(existing),
                            key_to_digest_string(&digest),
                        ));
                    }
                    merged.insert(alias, digest);
                }
                Some(_) => {}
                None => {
                    merged.insert(alias, digest);
                }
            }
        }

        for (k, ns, alias) in &tombstones.aliases {
            if *k == kind && ns == namespace {
                merged.remove(alias);
            }
        }

        write_alias_namespace(&tsv_path, &merged)
    }

    // =========================================================================
    // Orphan GC support
    // =========================================================================

    /// The set of sequence digests still referenced by SOME collection on disk,
    /// excluding `exclude` (the collection being removed).
    ///
    /// Derived from disk, not from `name_lookup`. `name_lookup` is populated only
    /// by `load_collections_from_directory`, `ensure_collection_loaded`, and the
    /// import path — NOT by the normal `open_local` stub load. On a
    /// partially-loaded store it is empty, so a `name_lookup`-derived live set
    /// says "nothing is referenced" and the GC unlinks sequences other
    /// collections still need. That difference is invisible from inside the
    /// function: a partially-loaded store looks exactly like a fully-loaded small
    /// one.
    ///
    /// FAILS CLOSED. If any collection listed in the index has no readable
    /// `.rgsi` and is not fully resident in memory, this returns `Err` and the
    /// caller must delete nothing. A missing input must never be read as "nothing
    /// references this".
    pub(crate) fn live_sequence_digests_from_disk(
        &self,
        exclude: &DigestKey,
    ) -> Result<HashSet<DigestKey>> {
        let local_path = self.local_path.as_ref().context("local_path not set")?;
        let tombstones = self.tombstones.lock().unwrap().clone();
        let index_path = local_path.join("collections.rgci");
        let collections = self.merged_collection_rows(&index_path, &tombstones)?;

        let mut live = HashSet::new();
        for (key, meta) in &collections {
            if key == exclude {
                continue;
            }

            // Prefer a fully-resident in-memory record; it needs no I/O and is
            // authoritative for collections that were never persisted.
            if let Some(record) = self.collections.get(key) {
                if let Some(sequences) = record.sequences() {
                    for seq in sequences {
                        live.insert(seq.metadata().sha512t24u.to_key());
                    }
                    continue;
                }
            }

            let rgsi_path = local_path.join(format!("collections/{}.rgsi", meta.digest));
            let collection = read_rgsi_file(&rgsi_path).with_context(|| {
                format!(
                    "refusing to remove orphan sequences: collection {} is listed in \
                     collections.rgci but {} could not be read. Treating that as \
                     'references nothing' would delete live sequence data.",
                    meta.digest,
                    rgsi_path.display()
                )
            })?;
            for seq in &collection.sequences {
                live.insert(seq.metadata().sha512t24u.to_key());
            }
        }
        Ok(live)
    }

    /// The sequence digests belonging to one collection, read from its
    /// `collections/<digest>.rgsi` (falling back to the in-memory record).
    pub(crate) fn collection_sequence_digests(&self, digest: &str) -> Result<Vec<DigestKey>> {
        let key = digest.to_key();
        if let Some(record) = self.collections.get(&key) {
            if let Some(sequences) = record.sequences() {
                return Ok(sequences
                    .iter()
                    .map(|s| s.metadata().sha512t24u.to_key())
                    .collect());
            }
        }
        if let Some(local_path) = self.local_path.as_ref() {
            let rgsi_path = local_path.join(format!("collections/{}.rgsi", digest));
            if rgsi_path.exists() {
                let collection = read_rgsi_file(&rgsi_path)?;
                return Ok(collection
                    .sequences
                    .iter()
                    .map(|s| s.metadata().sha512t24u.to_key())
                    .collect());
            }
        }
        // Last resort: whatever the name map happens to hold.
        Ok(self
            .name_lookup
            .get(&key)
            .map(|m| m.values().cloned().collect())
            .unwrap_or_default())
    }

    /// Read the store metadata from rgstore.json, returning state digests and timestamp.
    ///
    /// Returns a HashMap with keys: modified, collections_digest, sequences_digest,
    /// aliases_digest, fhr_digest, logical_sequence_bytes. Missing values are omitted
    /// from the map.
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
        if let Some(v) = metadata.logical_sequence_bytes { map.insert("logical_sequence_bytes".to_string(), v.to_string()); }
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

        // Manifest is the single source of truth: load exactly the alias
        // namespaces rgstore.json advertises, rather than scanning the aliases/
        // directory. This unifies local and remote behavior — an alias file the
        // manifest omits is not loaded, so a stale manifest fails the same way on
        // disk as it does over HTTP (where directory listing is impossible),
        // surfacing manifest-write bugs instead of silently masking them.
        let aliases_dir = root_path.join("aliases");
        store.aliases.load_namespaces_from_dir(
            &aliases_dir,
            &metadata.sequence_alias_namespaces,
            &metadata.collection_alias_namespaces,
        )?;
        store.available_sequence_alias_namespaces = metadata.sequence_alias_namespaces;
        store.available_collection_alias_namespaces = metadata.collection_alias_namespaces;

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

// ============================================================================
// Row readers / writers (free functions)
//
// These read and write index ROWS without touching store state, which is what
// makes a merge possible: `load_sequences_from_index` and
// `load_collection_stubs_from_rgci` mutate a store, so they cannot be used to
// look at what is currently on disk during a commit.
// ============================================================================

/// Read `sequences.rgsi` into digest-keyed rows. A missing file is an empty map
/// (a brand-new store), not an error.
pub(crate) fn read_sequences_rgsi_rows(path: &Path) -> Result<HashMap<DigestKey, SequenceMetadata>> {
    let mut rows = HashMap::new();
    let file = match File::open(path) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(rows),
        Err(e) => return Err(e).with_context(|| format!("reading {}", path.display())),
    };
    for line in std::io::BufReader::new(file).lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        if let Some(meta) = parse_rgsi_line(&line) {
            rows.insert(meta.sha512t24u.to_key(), meta);
        }
    }
    Ok(rows)
}

/// Read `collections.rgci` into digest-keyed rows. A missing file is an empty map.
pub(crate) fn read_collections_rgci_rows(
    path: &Path,
) -> Result<HashMap<DigestKey, SequenceCollectionMetadata>> {
    let mut rows = HashMap::new();
    let file = match File::open(path) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(rows),
        Err(e) => return Err(e).with_context(|| format!("reading {}", path.display())),
    };
    for line in std::io::BufReader::new(file).lines() {
        let line = line?;
        if let Some(meta) = parse_rgci_line(&line) {
            rows.insert(meta.digest.to_key(), meta);
        }
    }
    Ok(rows)
}

/// Fold one in-memory sequence row into the merged set.
///
/// `length`, `alphabet`, and `md5` are determined by the sha512t24u key, so a
/// key collision cannot disagree about them. `name` and `description` are NOT:
/// the same sequence is `chr1` in one collection and `1` in another. We keep the
/// row already on disk — never rewrite a published row — which makes the outcome
/// stable but order-dependent across processes: whichever writer commits first
/// sets the name.
///
/// The honest framing is that `sequences.rgsi`'s `name` column is ADVISORY. The
/// authoritative per-collection name mapping lives in `name_lookup` and in each
/// `collections/<digest>.rgsi`. A lexicographic tie-break would be
/// order-independent but would silently rename sequences on rebuild, which is
/// worse.
fn merge_sequence_rows(
    merged: &mut HashMap<DigestKey, SequenceMetadata>,
    key: DigestKey,
    memory: &SequenceMetadata,
) {
    merged.entry(key).or_insert_with(|| memory.clone());
}

/// Fold one in-memory collection row into the merged set.
///
/// The key is the collection digest and every other column is a digest of the
/// collection's CONTENT, so the same key implies the same row. The one asymmetry
/// is vintage: rows written by older gtars may have empty ancillary columns
/// (`name_length_pairs_digest`, `sorted_name_length_pairs_digest`,
/// `sorted_sequences_digest`). Prefer whichever row carries more of them, so a
/// merge upgrades an old row instead of pinning it.
fn merge_collection_rows(
    merged: &mut HashMap<DigestKey, SequenceCollectionMetadata>,
    key: DigestKey,
    memory: &SequenceCollectionMetadata,
) {
    let ancillary_count = |m: &SequenceCollectionMetadata| {
        m.name_length_pairs_digest.is_some() as u8
            + m.sorted_name_length_pairs_digest.is_some() as u8
            + m.sorted_sequences_digest.is_some() as u8
    };
    match merged.get(&key) {
        Some(existing) if ancillary_count(existing) >= ancillary_count(memory) => {}
        _ => {
            merged.insert(key, memory.clone());
        }
    }
}

/// Atomically publish `sequences.rgsi` from a row set.
///
/// Sorted by `sha512t24u` so output is deterministic regardless of build order
/// (the row set is a `HashMap`).
pub(crate) fn write_sequences_rgsi_rows(
    path: &Path,
    rows: &HashMap<DigestKey, SequenceMetadata>,
) -> Result<()> {
    let mut entries: Vec<&SequenceMetadata> = rows.values().collect();
    entries.sort_by(|a, b| a.sha512t24u.cmp(&b.sha512t24u));

    atomic_write(path, |w| {
        writeln!(w, "{}", RGSI_HEADER)?;
        for meta in entries {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}",
                meta.name,
                meta.length,
                meta.alphabet,
                meta.sha512t24u,
                meta.md5,
                meta.description.as_deref().unwrap_or("")
            )?;
        }
        Ok(())
    })
}

/// Atomically publish `collections.rgci` from a row set.
///
/// Sorted by collection digest for deterministic output.
pub(crate) fn write_collections_rgci_rows(
    path: &Path,
    rows: &HashMap<DigestKey, SequenceCollectionMetadata>,
) -> Result<()> {
    let mut entries: Vec<&SequenceCollectionMetadata> = rows.values().collect();
    entries.sort_by(|a, b| a.digest.cmp(&b.digest));

    atomic_write(path, |w| {
        writeln!(w, "{}", RGCI_HEADER)?;
        for meta in entries {
            writeln!(
                w,
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
    })
}

// ============================================================================
// Manifest and alias file helpers
// ============================================================================

/// Read `<dir>/rgstore.json`, or `None` if it does not exist yet.
pub(crate) fn read_manifest(dir: &Path) -> Result<Option<StoreMetadata>> {
    let path = dir.join("rgstore.json");
    match fs::read_to_string(&path) {
        Ok(json) => Ok(Some(serde_json::from_str(&json).with_context(|| {
            format!("Failed to parse {}", path.display())
        })?)),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(e).with_context(|| format!("Failed to read {}", path.display())),
    }
}

/// Atomically publish `<dir>/rgstore.json`.
pub(crate) fn write_manifest(dir: &Path, metadata: &StoreMetadata) -> Result<()> {
    let json = serde_json::to_string_pretty(metadata)
        .context("Failed to serialize metadata to JSON")?;
    super::atomic::atomic_write_bytes(&dir.join("rgstore.json"), json.as_bytes())
        .context("Failed to write rgstore.json")
}

/// Namespace names (file stems) of the `*.tsv` files actually present in `dir`.
fn scan_alias_namespaces(dir: &Path) -> Vec<String> {
    let Ok(entries) = fs::read_dir(dir) else {
        return Vec::new();
    };
    entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|x| x.to_str()) == Some("tsv"))
        .filter_map(|p| p.file_stem().and_then(|s| s.to_str()).map(String::from))
        .collect()
}

/// What the store's alias directory actually advertises, as
/// `(sequence_namespaces, collection_namespaces)`.
///
/// This is the source of truth for the manifest's namespace lists at commit
/// time. `open_local` deliberately trusts the MANIFEST rather than the directory
/// (so a stale manifest fails identically on disk and over HTTP, surfacing bugs
/// instead of masking them) — correct for reading, wrong for writing, where
/// trusting a stale snapshot is precisely how namespaces get dropped.
pub(crate) fn published_alias_namespaces(local_path: &Path) -> (Vec<String>, Vec<String>) {
    let aliases_dir = local_path.join("aliases");
    let mut seq = scan_alias_namespaces(&aliases_dir.join(AliasKind::Sequence.subdir()));
    let mut coll = scan_alias_namespaces(&aliases_dir.join(AliasKind::Collection.subdir()));
    seq.sort();
    coll.sort();
    (seq, coll)
}

/// Read one alias namespace TSV. A missing file is an empty map.
pub(crate) fn read_alias_tsv(path: &Path) -> Result<HashMap<String, DigestKey>> {
    let mut map = HashMap::new();
    let file = match File::open(path) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(map),
        Err(e) => return Err(e).with_context(|| format!("reading {}", path.display())),
    };
    for line in std::io::BufReader::new(file).lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if let Some((alias, digest)) = line.split_once('\t') {
            map.insert(alias.to_string(), digest.to_key());
        }
    }
    Ok(map)
}

/// Atomically publish one alias namespace TSV, or delete it if the merged
/// namespace is empty.
///
/// Sorted by alias so the file is deterministic (the map is a `HashMap`), which
/// keeps `aliases_digest` stable across writers that committed the same content.
pub(crate) fn write_alias_namespace(path: &Path, aliases: &HashMap<String, DigestKey>) -> Result<()> {
    if aliases.is_empty() {
        match fs::remove_file(path) {
            Ok(()) => Ok(()),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
            Err(e) => Err(e).with_context(|| format!("removing empty {}", path.display())),
        }
    } else {
        let sorted: BTreeMap<&String, &DigestKey> = aliases.iter().collect();
        atomic_write(path, |w| {
            for (alias, digest) in sorted {
                writeln!(w, "{}\t{}", alias, key_to_digest_string(digest))?;
            }
            Ok(())
        })
    }
}
