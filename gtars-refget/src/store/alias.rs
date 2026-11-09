//! Alias management for human-readable sequence and collection names.
//!
//! Contains AliasManager types, store bridge methods, RefgetStore delegates,
//! and remote pull methods.

use super::*;
use super::readonly::ReadonlyRefgetStore;
use super::core::RefgetStore;

use std::collections::{HashMap, HashSet};
use std::fs::{self, File, create_dir_all};
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};

use crate::hashkeyable::{DigestKey, HashKeyable, key_to_digest_string};

// =========================================================================
// Types
// =========================================================================

/// Identifies whether an alias targets a sequence or a collection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AliasKind {
    Sequence,
    Collection,
}

impl AliasKind {
    /// Returns the subdirectory name for this alias kind.
    pub fn subdir(&self) -> &'static str {
        match self {
            Self::Sequence => "sequences",
            Self::Collection => "collections",
        }
    }
}

/// Generic alias store: namespace -> { alias -> digest_key }
type AliasStore = HashMap<String, HashMap<String, DigestKey>>;

/// Manages human-readable aliases for sequences and collections.
#[derive(Debug)]
pub struct AliasManager {
    sequence_aliases: AliasStore,
    collection_aliases: AliasStore,
}

// =========================================================================
// Private helpers (operate on any AliasStore)
// =========================================================================

/// Error if `namespace` already binds `alias` to a digest other than `digest`,
/// unless `force`. Read-only: the shared check behind every checked add.
fn alias_check_existing(
    store: &AliasStore,
    kind: AliasKind,
    namespace: &str,
    alias: &str,
    digest: &DigestKey,
    force: bool,
) -> Result<()> {
    match alias_resolve(store, namespace, alias) {
        Some(existing) if existing != *digest && !force => Err(anyhow::anyhow!(
            "alias conflict in {} namespace '{}': '{}' already points to {}, \
             refusing to point it to {}. Aliases must be unique within a namespace; \
             use force-alias (CLI: --force-alias) only to deliberately re-point it.",
            kind.subdir(),
            namespace,
            alias,
            key_to_digest_string(&existing),
            key_to_digest_string(digest),
        )),
        _ => Ok(()),
    }
}

/// Bind `alias` -> `digest` in `namespace`. Same digest: no-op. Different
/// digest: error unless `force`. A refused add leaves `store` untouched (no
/// empty namespace is created).
fn alias_add(
    store: &mut AliasStore,
    kind: AliasKind,
    namespace: &str,
    alias: &str,
    digest: DigestKey,
    force: bool,
) -> Result<()> {
    alias_check_existing(store, kind, namespace, alias, &digest, force)?;
    store
        .entry(namespace.to_string())
        .or_default()
        .insert(alias.to_string(), digest);
    Ok(())
}

fn alias_resolve(store: &AliasStore, namespace: &str, alias: &str) -> Option<DigestKey> {
    store.get(namespace).and_then(|ns| ns.get(alias)).copied()
}

fn alias_reverse_scan(store: &AliasStore, digest: &DigestKey) -> Vec<(String, String)> {
    let mut results = Vec::new();
    for (namespace, aliases) in store {
        for (alias, d) in aliases {
            if d == digest {
                results.push((namespace.clone(), alias.clone()));
            }
        }
    }
    results
}

fn alias_namespaces(store: &AliasStore) -> Vec<String> {
    store.keys().cloned().collect()
}

fn alias_list(store: &AliasStore, namespace: &str) -> Option<Vec<String>> {
    store
        .get(namespace)
        .map(|ns| ns.keys().cloned().collect())
}

fn alias_remove(store: &mut AliasStore, namespace: &str, alias: &str) -> bool {
    if let Some(ns) = store.get_mut(namespace) {
        let removed = ns.remove(alias).is_some();
        if ns.is_empty() {
            store.remove(namespace);
        }
        removed
    } else {
        false
    }
}

/// Load a user-supplied TSV into `store` under `namespace`, all or nothing.
/// Returns the distinct alias names bound, in first-seen order.
///
/// Returns the names rather than a count because the store-level loaders have to
/// record each one as a pending change — a commit only publishes aliases this
/// handle actually touched.
///
/// Refuses the whole file, leaving `store` untouched, when:
/// * the file maps one alias to two different digests. Always an error, even
///   with `force`: one batch can never mean two things by one name, and letting
///   the last row win is how per-genome names like `chr1` silently collapsed
///   onto one genome.
/// * an alias already in `store` points elsewhere, unless `force`.
///
/// Exact duplicate rows are fine and count once.
fn alias_load_tsv(
    store: &mut AliasStore,
    kind: AliasKind,
    namespace: &str,
    path: &Path,
    force: bool,
) -> Result<Vec<String>> {
    // Phase 1: parse every row.
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows: Vec<(usize, String, DigestKey)> = Vec::new();
    for (idx, line) in reader.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.splitn(2, '\t').collect();
        if parts.len() == 2 {
            rows.push((idx + 1, parts[0].to_string(), parts[1].to_key()));
        }
    }

    // Phase 2: within this file, one alias must map to one digest.
    let mut first_seen: HashMap<&str, (DigestKey, usize)> = HashMap::new();
    let mut distinct: Vec<(&str, DigestKey)> = Vec::new();
    let mut clashing: HashSet<&str> = HashSet::new();
    let mut first_clash: Option<(&str, DigestKey, usize, DigestKey, usize)> = None;
    for (line_no, alias, digest) in &rows {
        match first_seen.get(alias.as_str()) {
            None => {
                first_seen.insert(alias, (*digest, *line_no));
                distinct.push((alias, *digest));
            }
            Some(&(first_digest, first_line)) if first_digest != *digest => {
                clashing.insert(alias);
                if first_clash.is_none() {
                    first_clash = Some((alias, first_digest, first_line, *digest, *line_no));
                }
            }
            Some(_) => {}
        }
    }
    if let Some((alias, digest_a, line_a, digest_b, line_b)) = first_clash {
        return Err(anyhow::anyhow!(
            "alias conflict in {} namespace '{}': '{}' maps '{}' to two different digests \
             ({} on line {}, {} on line {}); {}. A namespace maps each alias to one digest; \
             per-genome names such as chromosome names are not aliases.",
            kind.subdir(),
            namespace,
            path.display(),
            alias,
            key_to_digest_string(&digest_a),
            line_a,
            key_to_digest_string(&digest_b),
            line_b,
            clash_scale(clashing.len()),
        ));
    }

    // Phase 3: no alias may silently re-point an existing binding.
    for (alias, digest) in &distinct {
        alias_check_existing(store, kind, namespace, alias, digest, force)?;
    }

    // Phase 4: apply. Validated above, so nothing below can fail.
    if !distinct.is_empty() {
        let ns = store.entry(namespace.to_string()).or_default();
        for (alias, digest) in &distinct {
            ns.insert(alias.to_string(), *digest);
        }
    }
    Ok(distinct.into_iter().map(|(alias, _)| alias.to_string()).collect())
}

/// "N aliases in this file clash", reporting the scale of a clash: thousands
/// means a design mistake, not a typo.
fn clash_scale(n: usize) -> String {
    if n == 1 {
        "1 alias in this file clashes".to_string()
    } else {
        format!("{} aliases in this file clash", n)
    }
}

/// Read an alias TSV without clash checks, last row winning. Only for files
/// gtars wrote itself from a map (so they cannot hold a clash); being strict
/// here would make a store unopenable over one hand-edited file.
fn alias_read_tsv_unchecked(store: &mut AliasStore, namespace: &str, path: &Path) -> Result<()> {
    let map = super::persistence::read_alias_tsv(path)?;
    if !map.is_empty() {
        store.entry(namespace.to_string()).or_default().extend(map);
    }
    Ok(())
}

fn load_aliases_from_dir(store: &mut AliasStore, dir: &Path) -> Result<()> {
    if !dir.exists() {
        return Ok(());
    }
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) == Some("tsv") {
            let namespace = path
                .file_stem()
                .and_then(|s| s.to_str())
                .context("Invalid filename")?
                .to_string();
            alias_read_tsv_unchecked(store, &namespace, &path)?;
        }
    }
    Ok(())
}

fn write_all_aliases(store: &AliasStore, dir: &Path) -> Result<()> {
    if store.is_empty() {
        return Ok(());
    }
    create_dir_all(dir)?;
    for (namespace, aliases) in store {
        let tsv_path = dir.join(format!("{}.tsv", namespace));
        super::persistence::write_alias_namespace(&tsv_path, aliases)?;
    }
    Ok(())
}

// =========================================================================
// AliasManager public API
// =========================================================================

impl AliasManager {
    pub fn new() -> Self {
        AliasManager {
            sequence_aliases: HashMap::new(),
            collection_aliases: HashMap::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.sequence_aliases.is_empty() && self.collection_aliases.is_empty()
    }

    // --- Sequence aliases ---

    /// Bind a sequence alias. Errors if it already points to a different
    /// digest, unless `force`.
    pub fn add_sequence(&mut self, namespace: &str, alias: &str, digest: &str, force: bool) -> Result<()> {
        let key = digest.to_key();
        alias_add(&mut self.sequence_aliases, AliasKind::Sequence, namespace, alias, key, force)
    }

    pub fn resolve_sequence(&self, namespace: &str, alias: &str) -> Option<DigestKey> {
        alias_resolve(&self.sequence_aliases, namespace, alias)
    }

    pub fn reverse_lookup_sequence(&self, digest: &str) -> Vec<(String, String)> {
        let key = digest.to_key();
        alias_reverse_scan(&self.sequence_aliases, &key)
    }

    pub fn sequence_namespaces(&self) -> Vec<String> {
        alias_namespaces(&self.sequence_aliases)
    }

    pub fn sequence_aliases(&self, namespace: &str) -> Option<Vec<String>> {
        alias_list(&self.sequence_aliases, namespace)
    }

    pub fn remove_sequence(&mut self, namespace: &str, alias: &str) -> bool {
        alias_remove(&mut self.sequence_aliases, namespace, alias)
    }

    /// Load a TSV of sequence aliases, all or nothing, returning the distinct
    /// alias names bound. See [`alias_load_tsv`] for what is refused.
    pub fn load_sequence_tsv(&mut self, namespace: &str, path: &Path, force: bool) -> Result<Vec<String>> {
        alias_load_tsv(&mut self.sequence_aliases, AliasKind::Sequence, namespace, path, force)
    }

    // --- Collection aliases ---

    /// Bind a collection alias. Errors if it already points to a different
    /// digest, unless `force`.
    pub fn add_collection(&mut self, namespace: &str, alias: &str, digest: &str, force: bool) -> Result<()> {
        let key = digest.to_key();
        alias_add(&mut self.collection_aliases, AliasKind::Collection, namespace, alias, key, force)
    }

    pub fn resolve_collection(&self, namespace: &str, alias: &str) -> Option<DigestKey> {
        alias_resolve(&self.collection_aliases, namespace, alias)
    }

    pub fn reverse_lookup_collection(&self, digest: &str) -> Vec<(String, String)> {
        let key = digest.to_key();
        alias_reverse_scan(&self.collection_aliases, &key)
    }

    pub fn collection_namespaces(&self) -> Vec<String> {
        alias_namespaces(&self.collection_aliases)
    }

    pub fn collection_aliases(&self, namespace: &str) -> Option<Vec<String>> {
        alias_list(&self.collection_aliases, namespace)
    }

    pub fn remove_collection(&mut self, namespace: &str, alias: &str) -> bool {
        alias_remove(&mut self.collection_aliases, namespace, alias)
    }

    /// Load a TSV of collection aliases, all or nothing, returning the distinct
    /// alias names bound. See [`alias_load_tsv`] for what is refused.
    pub fn load_collection_tsv(&mut self, namespace: &str, path: &Path, force: bool) -> Result<Vec<String>> {
        alias_load_tsv(&mut self.collection_aliases, AliasKind::Collection, namespace, path, force)
    }

    // --- Persistence ---

    pub fn load_from_dir(&mut self, aliases_dir: &Path) -> Result<()> {
        load_aliases_from_dir(&mut self.sequence_aliases, &aliases_dir.join("sequences"))?;
        load_aliases_from_dir(&mut self.collection_aliases, &aliases_dir.join("collections"))?;
        Ok(())
    }

    /// Load only the alias namespaces explicitly listed (the manifest-driven
    /// path). Unlike `load_from_dir`, which scans the directory for any `*.tsv`
    /// present, this loads exactly `<aliases_dir>/{sequences,collections}/<ns>.tsv`
    /// for each listed namespace. This makes `rgstore.json` the single source of
    /// truth for what a store advertises: an alias file on disk that the manifest
    /// does not list is simply not loaded, so a stale manifest fails identically
    /// on a local filesystem and over HTTP (where directory listing is impossible).
    pub fn load_namespaces_from_dir(
        &mut self,
        aliases_dir: &Path,
        sequence_namespaces: &[String],
        collection_namespaces: &[String],
    ) -> Result<()> {
        let seq_dir = aliases_dir.join("sequences");
        for ns in sequence_namespaces {
            let path = seq_dir.join(format!("{}.tsv", ns));
            if path.exists() {
                alias_read_tsv_unchecked(&mut self.sequence_aliases, ns, &path)?;
            }
        }
        let coll_dir = aliases_dir.join("collections");
        for ns in collection_namespaces {
            let path = coll_dir.join(format!("{}.tsv", ns));
            if path.exists() {
                alias_read_tsv_unchecked(&mut self.collection_aliases, ns, &path)?;
            }
        }
        Ok(())
    }

    pub fn write_to_dir(&self, aliases_dir: &Path) -> Result<()> {
        write_all_aliases(&self.sequence_aliases, &aliases_dir.join("sequences"))?;
        write_all_aliases(&self.collection_aliases, &aliases_dir.join("collections"))?;
        Ok(())
    }

    // --- Kind-generic accessors (used by the commit path) ---

    /// This manager's in-memory namespaces for one alias kind.
    pub(crate) fn namespaces_for(&self, kind: AliasKind) -> Vec<String> {
        alias_namespaces(self.store_for(kind))
    }

    /// This manager's in-memory `alias -> digest` map for one namespace.
    /// Empty when the namespace is not in memory — which is NOT the same as
    /// "the namespace has no aliases on disk"; see
    /// [`ReadonlyRefgetStore::commit_alias_namespace`](super::readonly::ReadonlyRefgetStore).
    pub(crate) fn namespace_map(&self, kind: AliasKind, namespace: &str) -> HashMap<String, DigestKey> {
        self.store_for(kind)
            .get(namespace)
            .cloned()
            .unwrap_or_default()
    }

    /// Kind-generic single-alias resolve, used by the commit path to look up the
    /// target of an alias it has recorded as a pending change.
    pub(crate) fn resolve(&self, kind: AliasKind, namespace: &str, alias: &str) -> Option<DigestKey> {
        alias_resolve(self.store_for(kind), namespace, alias)
    }

    fn store_for(&self, kind: AliasKind) -> &AliasStore {
        match kind {
            AliasKind::Sequence => &self.sequence_aliases,
            AliasKind::Collection => &self.collection_aliases,
        }
    }
}

impl Default for AliasManager {
    fn default() -> Self {
        Self::new()
    }
}

// =========================================================================
// ReadonlyRefgetStore alias bridge methods
// =========================================================================

impl ReadonlyRefgetStore {
    /// Add a sequence alias and persist to disk if applicable.
    ///
    /// Re-adding the same binding is a no-op. Pointing an existing alias at a
    /// different digest errors unless `force_alias` is set; a refused add
    /// leaves memory, pending changes and disk untouched.
    pub fn add_sequence_alias(&mut self, namespace: &str, alias: &str, digest: &str) -> Result<()> {
        self.aliases.add_sequence(namespace, alias, digest, self.force_alias)?;
        self.record(|p| p.add_alias(AliasKind::Sequence, namespace, alias));
        self.persist_alias_namespace(AliasKind::Sequence, namespace)?;
        Ok(())
    }

    /// Resolve a sequence alias to sequence metadata (no data loading).
    pub fn get_sequence_metadata_by_alias(&self, namespace: &str, alias: &str) -> Option<&crate::digest::SequenceMetadata> {
        let key = self.aliases.resolve_sequence(namespace, alias)?;
        self.sequence_store.get(&key).map(|rec| rec.metadata())
    }

    /// Resolve a sequence alias and return the loaded sequence record.
    pub fn get_sequence_by_alias(&self, namespace: &str, alias: &str) -> Result<&crate::digest::SequenceRecord> {
        let key = self.aliases.resolve_sequence(namespace, alias)
            .ok_or_else(|| anyhow::anyhow!("Sequence alias not found: {}/{}", namespace, alias))?;
        self.sequence_store.get(&key)
            .ok_or_else(|| anyhow::anyhow!("Sequence not found for alias {}/{}", namespace, alias))
    }

    /// Reverse lookup: find all aliases pointing to this sequence digest.
    pub fn get_aliases_for_sequence(&self, digest: &str) -> Vec<(String, String)> {
        self.aliases.reverse_lookup_sequence(digest)
    }

    /// List all sequence alias namespaces.
    pub fn list_sequence_alias_namespaces(&self) -> Vec<String> {
        self.aliases.sequence_namespaces()
    }

    /// List all aliases in a sequence alias namespace.
    pub fn list_sequence_aliases(&self, namespace: &str) -> Option<Vec<String>> {
        self.aliases.sequence_aliases(namespace)
    }

    /// Remove a single sequence alias.
    pub fn remove_sequence_alias(&mut self, namespace: &str, alias: &str) -> Result<bool> {
        let removed = self.aliases.remove_sequence(namespace, alias);
        if removed {
            // Record BEFORE persisting: the commit starts from the TSV on disk,
            // so without an explicit removal the alias is simply left alone.
            self.record_alias_removal(AliasKind::Sequence, namespace, alias);
            self.persist_alias_namespace(AliasKind::Sequence, namespace)?;
        }
        Ok(removed)
    }

    /// Load sequence aliases from a TSV file into a namespace, all or nothing.
    ///
    /// Returns the number of distinct aliases bound. Errors, loading nothing,
    /// if the file maps one alias to two different digests (always, even with
    /// `force_alias`), or if it would re-point an existing alias (unless
    /// `force_alias`).
    pub fn load_sequence_aliases(&mut self, namespace: &str, path: &str) -> Result<usize> {
        let loaded = self.aliases.load_sequence_tsv(namespace, Path::new(path), self.force_alias)?;
        self.record(|p| {
            for alias in &loaded {
                p.add_alias(AliasKind::Sequence, namespace, alias);
            }
        });
        self.persist_alias_namespace(AliasKind::Sequence, namespace)?;
        Ok(loaded.len())
    }

    /// Add a collection alias and persist to disk if applicable.
    ///
    /// Re-adding the same binding is a no-op. Pointing an existing alias at a
    /// different digest errors unless `force_alias` is set; a refused add
    /// leaves memory, pending changes and disk untouched.
    pub fn add_collection_alias(&mut self, namespace: &str, alias: &str, digest: &str) -> Result<()> {
        self.aliases.add_collection(namespace, alias, digest, self.force_alias)?;
        self.record(|p| p.add_alias(AliasKind::Collection, namespace, alias));
        self.persist_alias_namespace(AliasKind::Collection, namespace)?;
        Ok(())
    }

    /// Check whether an import-time collection alias request would conflict,
    /// WITHOUT mutating anything.
    ///
    /// This is the read-only half of [`Self::register_import_collection_alias`]:
    /// it returns exactly the error registration would return, so callers can
    /// validate up front and bail out before committing any part of an import.
    /// Returns `Ok(())` when `spec` is `None`, when `force` is set, when the
    /// alias is unused, or when it already points at this same collection (the
    /// idempotent re-import case).
    pub(crate) fn check_import_collection_alias(
        &self,
        spec: Option<(&str, &str)>,
        collection_digest: &str,
        force: bool,
    ) -> Result<()> {
        let Some((namespace, alias)) = spec else {
            return Ok(());
        };
        if force {
            return Ok(());
        }

        let new_key = collection_digest.to_key();
        match self.aliases.resolve_collection(namespace, alias) {
            Some(existing) if existing != new_key => Err(anyhow::anyhow!(
                "collection alias {}:{} already resolves to collection {} \
                 but this import would point it at {}; \
                 pass force to overwrite, or choose a different alias",
                namespace,
                alias,
                key_to_digest_string(&existing),
                collection_digest,
            )),
            _ => Ok(()),
        }
    }

    /// Register a collection alias requested at import time.
    ///
    /// No-op when `spec` is `None`. Idempotent when the alias already points at
    /// this digest. Errors on a conflicting existing alias unless `force`,
    /// because silently remapping a name like `ucsc:hg38` from one assembly to
    /// another is exactly the misdirection this option exists to prevent.
    ///
    /// The alias is recorded as pending and published by the import's final
    /// `write_index_files`, together with the collection it names.
    pub(crate) fn register_import_collection_alias(
        &mut self,
        spec: Option<(&str, &str)>,
        collection_digest: &str,
        force: bool,
    ) -> Result<()> {
        let Some((namespace, alias)) = spec else {
            return Ok(());
        };

        // Same error, evaluated the same way, as the up-front validation
        // callers run before mutating anything.
        self.check_import_collection_alias(spec, collection_digest, force)?;

        let new_key = collection_digest.to_key();
        match self.aliases.resolve_collection(namespace, alias) {
            // Same collection: already registered, nothing to do (and no
            // redundant re-persist). This is what makes re-import idempotent.
            Some(existing) if existing == new_key => return Ok(()),
            Some(existing) => {
                let existing_digest = key_to_digest_string(&existing);
                if !self.quiet {
                    println!(
                        "Warning: remapping collection alias {}:{} from {} to {} (force)",
                        namespace, alias, existing_digest, collection_digest,
                    );
                }
            }
            None => {}
        }

        // Deferred, not published here: the collection itself only reaches
        // `collections.rgci` in the import's final `write_index_files`, and an
        // alias visible on disk before its target is a dangling alias to every
        // concurrent reader (permanently so, if that final commit fails).
        // Publishing both in the same commit closes that window. `force` is
        // carried along so the on-disk conflict check honors it too, instead
        // of only the handle-wide `force_alias` flag. Forced in memory because
        // the check above already applied this call's own `force`.
        self.aliases.add_collection(namespace, alias, collection_digest, true)?;
        self.record(|p| {
            if force {
                p.add_forced_alias(AliasKind::Collection, namespace, alias);
            } else {
                p.add_alias(AliasKind::Collection, namespace, alias);
            }
        });
        Ok(())
    }

    /// Check FASTA header-tag sequence aliases `(namespace, alias, digest)`
    /// from one file, WITHOUT mutating anything, so an import can refuse them
    /// before it touches any state.
    ///
    /// The same alias tagged onto two different sequences in one file is
    /// always an error (it ignores `force_alias`): that is a per-genome name,
    /// not an alias. An alias that already points to a different digest
    /// errors unless `force_alias`.
    pub(crate) fn check_import_sequence_aliases(
        &self,
        aliases: &[(String, String, String)],
        source: &Path,
    ) -> Result<()> {
        let mut first_seen: HashMap<(&str, &str), DigestKey> = HashMap::new();
        let mut clashing: HashSet<(&str, &str)> = HashSet::new();
        let mut first_clash: Option<(&str, &str, DigestKey, DigestKey)> = None;
        for (namespace, alias, digest) in aliases {
            let key = digest.to_key();
            match first_seen.get(&(namespace.as_str(), alias.as_str())) {
                None => {
                    first_seen.insert((namespace, alias), key);
                }
                Some(&first) if first != key => {
                    clashing.insert((namespace, alias));
                    if first_clash.is_none() {
                        first_clash = Some((namespace, alias, first, key));
                    }
                }
                Some(_) => {}
            }
        }
        if let Some((namespace, alias, digest_a, digest_b)) = first_clash {
            return Err(anyhow::anyhow!(
                "alias conflict in {} namespace '{}': '{}' tags '{}' onto two different \
                 sequences ({} and {}); {}. A namespace maps each alias to one digest; \
                 per-genome names such as chromosome names are not aliases.",
                AliasKind::Sequence.subdir(),
                namespace,
                source.display(),
                alias,
                key_to_digest_string(&digest_a),
                key_to_digest_string(&digest_b),
                clash_scale(clashing.len()),
            ));
        }
        // In FASTA order, so the error names the first offending header.
        for (namespace, alias, digest) in aliases {
            alias_check_existing(
                &self.aliases.sequence_aliases,
                AliasKind::Sequence,
                namespace,
                alias,
                &digest.to_key(),
                self.force_alias,
            )?;
        }
        Ok(())
    }

    /// Register a sequence alias in memory and mark it pending, WITHOUT
    /// publishing the namespace TSV. Used by the FASTA import path, which
    /// commits all of its changes together at the end; publishing per alias
    /// would rewrite the TSV (and the manifest) once per header.
    pub(crate) fn add_sequence_alias_pending(&mut self, namespace: &str, alias: &str, digest: &str) -> Result<()> {
        self.aliases.add_sequence(namespace, alias, digest, self.force_alias)?;
        self.record(|p| p.add_alias(AliasKind::Sequence, namespace, alias));
        Ok(())
    }

    /// Resolve a collection alias to collection metadata.
    pub fn get_collection_metadata_by_alias(&self, namespace: &str, alias: &str) -> Option<&crate::digest::SequenceCollectionMetadata> {
        let key = self.aliases.resolve_collection(namespace, alias)?;
        self.collections.get(&key).map(|rec| rec.metadata())
    }

    /// Resolve a collection alias and return the loaded collection.
    pub fn get_collection_by_alias(&self, namespace: &str, alias: &str) -> Result<crate::digest::SequenceCollection> {
        let key = self.aliases.resolve_collection(namespace, alias)
            .ok_or_else(|| anyhow::anyhow!("Collection alias not found: {}/{}", namespace, alias))?;
        let digest_str = key_to_digest_string(&key);
        self.get_collection(&digest_str)
    }

    /// Reverse lookup: find all aliases pointing to this collection digest.
    pub fn get_aliases_for_collection(&self, digest: &str) -> Vec<(String, String)> {
        self.aliases.reverse_lookup_collection(digest)
    }

    /// List all collection alias namespaces.
    pub fn list_collection_alias_namespaces(&self) -> Vec<String> {
        self.aliases.collection_namespaces()
    }

    /// List all aliases in a collection alias namespace.
    pub fn list_collection_aliases(&self, namespace: &str) -> Option<Vec<String>> {
        self.aliases.collection_aliases(namespace)
    }

    /// Remove a single collection alias.
    pub fn remove_collection_alias(&mut self, namespace: &str, alias: &str) -> Result<bool> {
        let removed = self.aliases.remove_collection(namespace, alias);
        if removed {
            self.record_alias_removal(AliasKind::Collection, namespace, alias);
            self.persist_alias_namespace(AliasKind::Collection, namespace)?;
        }
        Ok(removed)
    }

    /// Load collection aliases from a TSV file into a namespace, all or nothing.
    ///
    /// Returns the number of distinct aliases bound. Errors, loading nothing,
    /// if the file maps one alias to two different digests (always, even with
    /// `force_alias`), or if it would re-point an existing alias (unless
    /// `force_alias`).
    pub fn load_collection_aliases(&mut self, namespace: &str, path: &str) -> Result<usize> {
        let loaded = self.aliases.load_collection_tsv(namespace, Path::new(path), self.force_alias)?;
        self.record(|p| {
            for alias in &loaded {
                p.add_alias(AliasKind::Collection, namespace, alias);
            }
        });
        self.persist_alias_namespace(AliasKind::Collection, namespace)?;
        Ok(loaded.len())
    }

    /// Commit a single alias namespace to disk (if disk-backed).
    ///
    /// The cheap alias-mutation path: it applies this handle's pending changes
    /// for ONE namespace and patches only the manifest's alias-derived fields,
    /// rather than rehashing the (potentially 60+ MB) sequence and collection
    /// indexes. Takes the store write lock for the duration, like any other
    /// commit.
    pub(crate) fn persist_alias_namespace(&self, kind: AliasKind, namespace: &str) -> Result<()> {
        if !self.persist_to_disk {
            return Ok(());
        }
        let Some(local_path) = self.local_path.clone() else {
            return Ok(());
        };

        let _guard = self.acquire_commit_lock("persist_alias_namespace")?;

        let pending = self.pending.lock().unwrap().clone();
        let aliases_dir = local_path.join("aliases");
        self.commit_alias_namespace(&aliases_dir, kind, namespace, &pending)?;

        // Keep rgstore.json in lock-step with the alias files on disk, so a
        // served store always advertises its true alias namespaces. Without
        // this, aliases added after the index was written (e.g. a post-build
        // alias step) leave the manifest stale and remote/HTTP clients —
        // which cannot list the aliases/ directory — never discover them.
        self.refresh_manifest_alias_namespaces()?;

        // Only this namespace has been published; leave the rest of the pending
        // changes for the next full commit.
        self.pending
            .lock()
            .unwrap()
            .clear_alias_namespace(kind, namespace);

        Ok(())
    }

    /// Record that an alias was deliberately removed, so the commit drops it
    /// from the TSV it reads back from disk instead of leaving it alone.
    pub(crate) fn record_alias_removal(&self, kind: AliasKind, namespace: &str, alias: &str) {
        let emptied = self.aliases.namespace_map(kind, namespace).is_empty();
        self.record(|p| {
            p.remove_alias(kind, namespace, alias);
            // A namespace that just lost its last in-memory alias is explicitly
            // emptied, as distinct from "we never loaded it". This only makes
            // the commit VISIT the namespace; whether the TSV is deleted depends
            // on the merged result being empty, because another writer may have
            // added aliases we never saw.
            if emptied {
                p.emptied_alias_namespaces
                    .insert((kind, namespace.to_string()));
            }
        });
    }
}

// =========================================================================
// RefgetStore alias delegates
// =========================================================================

impl RefgetStore {
    /// Add a sequence alias.
    pub fn add_sequence_alias(&mut self, namespace: &str, alias: &str, digest: &str) -> Result<()> {
        self.inner.add_sequence_alias(namespace, alias, digest)
    }

    /// Remove a single sequence alias.
    pub fn remove_sequence_alias(&mut self, namespace: &str, alias: &str) -> Result<bool> {
        self.inner.remove_sequence_alias(namespace, alias)
    }

    /// Load sequence aliases from a TSV file.
    pub fn load_sequence_aliases(&mut self, namespace: &str, path: &str) -> Result<usize> {
        self.inner.load_sequence_aliases(namespace, path)
    }

    /// Add a collection alias.
    pub fn add_collection_alias(&mut self, namespace: &str, alias: &str, digest: &str) -> Result<()> {
        self.inner.add_collection_alias(namespace, alias, digest)
    }

    /// Remove a single collection alias.
    pub fn remove_collection_alias(&mut self, namespace: &str, alias: &str) -> Result<bool> {
        self.inner.remove_collection_alias(namespace, alias)
    }

    /// Load collection aliases from a TSV file.
    pub fn load_collection_aliases(&mut self, namespace: &str, path: &str) -> Result<usize> {
        self.inner.load_collection_aliases(namespace, path)
    }

    /// Lazy-loading get_collection_by_alias.
    pub fn get_collection_by_alias(
        &mut self,
        namespace: &str,
        alias: &str,
    ) -> Result<crate::digest::SequenceCollection> {
        if let Some(meta) = self.inner.get_collection_metadata_by_alias(namespace, alias) {
            let digest = meta.digest.clone();
            if !self.inner.is_collection_loaded(&digest) {
                self.inner.load_collection(&digest)?;
            }
            return self.inner.get_collection_by_alias(namespace, alias);
        }
        Err(anyhow::anyhow!("Collection alias not found: {}:{}", namespace, alias))
    }

    // --- Sidecar pull: aliases ---

    /// Pull alias sidecars from the remote store.
    pub fn pull_aliases(
        &mut self,
        namespace: Option<&str>,
        strategy: SyncStrategy,
    ) -> Result<PullResult> {
        let seq_namespaces: Vec<String> = match namespace {
            Some(ns) => vec![ns.to_string()],
            None => self.inner.available_sequence_alias_namespaces.clone(),
        };
        let coll_namespaces: Vec<String> = match namespace {
            Some(ns) => vec![ns.to_string()],
            None => self.inner.available_collection_alias_namespaces.clone(),
        };
        self.pull_alias_namespaces(&seq_namespaces, &coll_namespaces, strategy)
    }

    /// Pull the given sequence and collection alias namespaces in one
    /// parallel batch.
    pub(crate) fn pull_alias_namespaces(
        &mut self,
        seq_namespaces: &[String],
        coll_namespaces: &[String],
        strategy: SyncStrategy,
    ) -> Result<PullResult> {
        let rels: Vec<String> = seq_namespaces
            .iter()
            .map(|ns| format!("aliases/sequences/{}.tsv", ns))
            .chain(coll_namespaces.iter().map(|ns| format!("aliases/collections/{}.tsv", ns)))
            .collect();
        let (result, _) = self.inner.sync_sidecars(&rels, &strategy)?;

        if result.pulled > 0 {
            if let Some(ref local_path) = self.inner.local_path {
                let aliases_dir = local_path.join("aliases");
                self.inner.aliases = AliasManager::default();
                self.inner.aliases.load_from_dir(&aliases_dir)?;
            }
        }

        Ok(result)
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    // -----------------------------------------------------------------------
    // Minimal HTTP test server helper
    // -----------------------------------------------------------------------

    /// Spin up a single-threaded HTTP file server on a random port.
    ///
    /// The server serves files from `serve_dir` for as long as the returned
    /// `JoinHandle` is alive (call `drop()` or let it go out of scope, but
    /// the thread will block on the next request; `shutdown_flag` is used to
    /// stop it cleanly after the test).
    ///
    /// Returns `(base_url, shutdown_fn)` where `shutdown_fn()` signals the
    /// server to stop.
    #[cfg(all(feature = "filesystem", feature = "http"))]
    fn start_file_server(serve_dir: std::path::PathBuf) -> (String, impl FnOnce()) {
        use std::io::{Read as _, Write as _};
        use std::net::TcpListener;
        use std::sync::{Arc, atomic::{AtomicBool, Ordering}};

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
                        let n = stream.read(&mut buf).unwrap_or(0);
                        let request = std::str::from_utf8(&buf[..n]).unwrap_or("");
                        // Parse GET /path HTTP/1.x
                        let path = request
                            .lines()
                            .next()
                            .and_then(|l| l.split_whitespace().nth(1))
                            .unwrap_or("/");
                        // Strip leading '/'
                        let rel = path.trim_start_matches('/');
                        let file_path = serve_dir.join(rel);
                        if file_path.exists() && file_path.is_file() {
                            let data = std::fs::read(&file_path).unwrap_or_default();
                            let header = format!(
                                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                                data.len()
                            );
                            let _ = stream.write_all(header.as_bytes());
                            let _ = stream.write_all(&data);
                        } else {
                            let body = b"Not Found";
                            let header = format!(
                                "HTTP/1.1 404 Not Found\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                                body.len()
                            );
                            let _ = stream.write_all(header.as_bytes());
                            let _ = stream.write_all(body);
                        }
                    }
                    Err(_) => break,
                }
            }
        });

        let shutdown = move || {
            stop.store(true, Ordering::Relaxed);
            // Connect once to unblock the accept() call
            let _ = std::net::TcpStream::connect(format!("127.0.0.1:{}", port));
        };

        (base_url, shutdown)
    }

    // --- AliasManager unit tests ---

    #[test]
    fn test_add_and_resolve_sequence() {
        let mut mgr = AliasManager::new();
        assert!(mgr.is_empty());

        mgr.add_sequence("ncbi", "NC_000001.11", "abc123", false).unwrap();
        assert!(!mgr.is_empty());

        let key = mgr.resolve_sequence("ncbi", "NC_000001.11");
        assert!(key.is_some());
        assert_eq!(key.unwrap(), "abc123".to_key());

        assert!(mgr.resolve_sequence("ncbi", "missing").is_none());
        assert!(mgr.resolve_sequence("missing_ns", "NC_000001.11").is_none());
    }

    #[test]
    fn test_reverse_lookup() {
        let mut mgr = AliasManager::new();
        mgr.add_sequence("ncbi", "NC_000001.11", "digest1", false).unwrap();
        mgr.add_sequence("ucsc", "chr1", "digest1", false).unwrap();
        mgr.add_sequence("ncbi", "NC_000002.12", "digest2", false).unwrap();

        let aliases = mgr.reverse_lookup_sequence("digest1");
        assert_eq!(aliases.len(), 2);
        assert!(aliases.contains(&("ncbi".to_string(), "NC_000001.11".to_string())));
        assert!(aliases.contains(&("ucsc".to_string(), "chr1".to_string())));

        let aliases2 = mgr.reverse_lookup_sequence("digest2");
        assert_eq!(aliases2.len(), 1);
    }

    #[test]
    fn test_namespaces_and_list() {
        let mut mgr = AliasManager::new();
        mgr.add_sequence("ncbi", "NC_000001.11", "d1", false).unwrap();
        mgr.add_sequence("ucsc", "chr1", "d1", false).unwrap();

        let ns = mgr.sequence_namespaces();
        assert!(ns.contains(&"ncbi".to_string()));
        assert!(ns.contains(&"ucsc".to_string()));

        let aliases = mgr.sequence_aliases("ncbi").unwrap();
        assert!(aliases.contains(&"NC_000001.11".to_string()));
    }

    #[test]
    fn test_remove() {
        let mut mgr = AliasManager::new();
        mgr.add_sequence("ncbi", "NC_000001.11", "d1", false).unwrap();

        assert!(mgr.remove_sequence("ncbi", "NC_000001.11"));
        assert!(mgr.resolve_sequence("ncbi", "NC_000001.11").is_none());
        assert!(mgr.sequence_namespaces().is_empty());

        assert!(!mgr.remove_sequence("ncbi", "NC_000001.11"));
    }

    #[test]
    fn test_collection_aliases() {
        let mut mgr = AliasManager::new();
        mgr.add_collection("ucsc", "hg38", "coll_digest", false).unwrap();
        mgr.add_collection("gencode", "GRCh38.p14", "coll_digest", false).unwrap();

        assert!(mgr.resolve_collection("ucsc", "hg38").is_some());
        assert_eq!(mgr.collection_namespaces().len(), 2);

        let aliases = mgr.reverse_lookup_collection("coll_digest");
        assert_eq!(aliases.len(), 2);

        assert!(mgr.remove_collection("ucsc", "hg38"));
        assert!(mgr.resolve_collection("ucsc", "hg38").is_none());
    }

    #[test]
    fn test_persistence_roundtrip() {
        let dir = tempdir().unwrap();
        let aliases_dir = dir.path().join("aliases");

        let mut mgr = AliasManager::new();
        mgr.add_sequence("ncbi", "NC_000001.11", "seq_digest", false).unwrap();
        mgr.add_collection("ucsc", "hg38", "coll_digest", false).unwrap();

        mgr.write_to_dir(&aliases_dir).unwrap();

        assert!(aliases_dir.join("sequences/ncbi.tsv").exists());
        assert!(aliases_dir.join("collections/ucsc.tsv").exists());

        let mut mgr2 = AliasManager::new();
        mgr2.load_from_dir(&aliases_dir).unwrap();

        assert!(mgr2.resolve_sequence("ncbi", "NC_000001.11").is_some());
        assert!(mgr2.resolve_collection("ucsc", "hg38").is_some());
    }

    #[test]
    fn test_load_from_missing_dir() {
        let mut mgr = AliasManager::new();
        mgr.load_from_dir(std::path::Path::new("/nonexistent/aliases")).unwrap();
        assert!(mgr.is_empty());
    }

    #[test]
    fn test_write_namespace_single() {
        let dir = tempdir().unwrap();
        let aliases_dir = dir.path().join("aliases");

        let mut mgr = AliasManager::new();
        mgr.add_sequence("ncbi", "NC_000001.11", "d1", false).unwrap();
        mgr.add_sequence("ucsc", "chr1", "d2", false).unwrap();

        let tsv = aliases_dir.join("sequences/ncbi.tsv");
        super::super::persistence::write_alias_namespace(
            &tsv,
            &mgr.namespace_map(AliasKind::Sequence, "ncbi"),
        )
        .unwrap();

        assert!(tsv.exists());
        assert!(!aliases_dir.join("sequences/ucsc.tsv").exists());
    }

    #[test]
    fn test_write_namespace_absent_from_memory_does_not_delete() {
        // The old `write_namespace` deleted the TSV whenever the namespace was
        // missing from memory, which is how another writer's brand-new namespace
        // got orphaned: `open_local` loads only what the manifest advertises, so
        // "not in memory" routinely means "not mine", not "deleted".
        let dir = tempdir().unwrap();
        let tsv = dir.path().join("aliases/sequences/ncbi.tsv");
        std::fs::create_dir_all(tsv.parent().unwrap()).unwrap();
        std::fs::write(&tsv, "NC_000001.11\tother_writers_digest\n").unwrap();

        let mgr = AliasManager::new();
        let mut merged = super::super::persistence::read_alias_tsv(&tsv).unwrap();
        merged.extend(mgr.namespace_map(AliasKind::Sequence, "ncbi"));
        super::super::persistence::write_alias_namespace(&tsv, &merged).unwrap();

        assert!(tsv.exists(), "another writer's namespace was deleted");
        assert_eq!(
            std::fs::read_to_string(&tsv).unwrap(),
            "NC_000001.11\tother_writers_digest\n"
        );
    }

    #[test]
    fn test_empty_write_is_noop() {
        let dir = tempdir().unwrap();
        let aliases_dir = dir.path().join("aliases");

        let mgr = AliasManager::new();
        mgr.write_to_dir(&aliases_dir).unwrap();

        assert!(!aliases_dir.join("sequences").exists());
        assert!(!aliases_dir.join("collections").exists());
    }

    #[test]
    fn test_load_tsv() {
        let dir = tempdir().unwrap();
        let tsv_path = dir.path().join("ncbi.tsv");
        std::fs::write(&tsv_path, "NC_000001.11\tsome_digest\n# comment\n\nNC_000002.12\tanother_digest\n").unwrap();

        let mut mgr = AliasManager::new();
        let loaded = mgr.load_sequence_tsv("ncbi", &tsv_path, false).unwrap();
        assert_eq!(loaded, vec!["NC_000001.11", "NC_000002.12"]);
        assert!(mgr.resolve_sequence("ncbi", "NC_000001.11").is_some());
        assert!(mgr.resolve_sequence("ncbi", "NC_000002.12").is_some());
    }

    // --- Store-level alias integration tests ---

    #[cfg(feature = "filesystem")]
    fn copy_test_fasta(temp_dir: &std::path::Path, name: &str) -> std::path::PathBuf {
        let src = format!("../tests/data/fasta/{}", name);
        let dst = temp_dir.join(name);
        std::fs::copy(&src, &dst)
            .unwrap_or_else(|e| panic!("Failed to copy {} to tempdir: {}", src, e));
        dst
    }

    #[test]
    fn test_store_sequence_aliases() {
        use crate::digest::digest_sequence;

        let mut store = RefgetStore::in_memory();
        let record = digest_sequence("chr1", b"ACGT");
        store.add_sequence_record(record.clone(), false).unwrap();

        let digest = record.metadata().sha512t24u.clone();
        store.add_sequence_alias("ncbi", "NC_000001.11", &digest).unwrap();
        store.add_sequence_alias("ucsc", "chr1", &digest).unwrap();

        let found = store.get_sequence_metadata_by_alias("ncbi", "NC_000001.11").unwrap();
        assert_eq!(found.name, "chr1");

        let aliases = store.get_aliases_for_sequence(&digest);
        assert_eq!(aliases.len(), 2);
        assert!(aliases.contains(&("ncbi".to_string(), "NC_000001.11".to_string())));
        assert!(aliases.contains(&("ucsc".to_string(), "chr1".to_string())));

        let ns = store.list_sequence_alias_namespaces();
        assert!(ns.contains(&"ncbi".to_string()));
        assert!(ns.contains(&"ucsc".to_string()));

        let aliases = store.list_sequence_aliases("ncbi").unwrap();
        assert!(aliases.contains(&"NC_000001.11".to_string()));
    }

    #[test]
    #[cfg(feature = "filesystem")]
    fn test_store_collection_aliases() {
        let temp = tempdir().unwrap();
        let fasta_path = copy_test_fasta(temp.path(), "base.fa");

        let mut store = RefgetStore::in_memory();
        let (meta, _) = store
            .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
            .unwrap();

        store.add_collection_alias("ucsc", "hg38", &meta.digest).unwrap();
        store.add_collection_alias("gencode", "GRCh38.p14", &meta.digest).unwrap();

        let coll = store.get_collection_metadata_by_alias("ucsc", "hg38").unwrap();
        assert_eq!(coll.digest, meta.digest);

        let aliases = store.get_aliases_for_collection(&meta.digest);
        assert_eq!(aliases.len(), 2);
    }

    #[test]
    fn test_store_alias_remove() {
        use crate::digest::digest_sequence;

        let mut store = RefgetStore::in_memory();
        let record = digest_sequence("chr1", b"ACGT");
        store.add_sequence_record(record.clone(), false).unwrap();
        let digest = record.metadata().sha512t24u.clone();

        store.add_sequence_alias("ncbi", "NC_000001.11", &digest).unwrap();
        assert!(store.get_sequence_metadata_by_alias("ncbi", "NC_000001.11").is_some());

        assert!(store.remove_sequence_alias("ncbi", "NC_000001.11").unwrap());
        assert!(store.get_sequence_metadata_by_alias("ncbi", "NC_000001.11").is_none());

        assert!(store.list_sequence_alias_namespaces().is_empty());
    }

    #[test]
    #[cfg(feature = "filesystem")]
    fn test_store_alias_persistence() {
        let dir = tempdir().unwrap();
        let store_path = dir.path().join("store");

        let fasta_temp = tempdir().unwrap();
        let fasta_path = copy_test_fasta(fasta_temp.path(), "base.fa");

        let digest: String;
        let seq_digest: String;
        {
            let mut store = RefgetStore::on_disk(&store_path).unwrap();
            let (meta, _) = store
                .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
                .unwrap();
            digest = meta.digest.clone();
            seq_digest = store.list_sequences().unwrap()[0].sha512t24u.clone();

            store.add_sequence_alias("ncbi", "NC_000001.11", &seq_digest).unwrap();
            store.add_collection_alias("ucsc", "hg38", &digest).unwrap();
        }

        {
            let store = RefgetStore::open_local(&store_path).unwrap();
            assert!(store.get_sequence_metadata_by_alias("ncbi", "NC_000001.11").is_some());
            assert!(store.get_collection_metadata_by_alias("ucsc", "hg38").is_some());

            assert!(store_path.join("aliases/sequences/ncbi.tsv").exists());
            assert!(store_path.join("aliases/collections/ucsc.tsv").exists());
        }
    }

    #[test]
    fn test_store_alias_load_tsv() {
        use crate::digest::digest_sequence;

        let dir = tempdir().unwrap();
        let tsv_path = dir.path().join("ncbi.tsv");

        let mut store = RefgetStore::in_memory();
        let record = digest_sequence("chr1", b"ACGT");
        store.add_sequence_record(record.clone(), false).unwrap();
        let digest = record.metadata().sha512t24u.clone();

        std::fs::write(&tsv_path, format!("NC_000001.11\t{}\n", digest)).unwrap();

        let count = store.load_sequence_aliases("ncbi", tsv_path.to_str().unwrap()).unwrap();
        assert_eq!(count, 1);
        assert!(store.get_sequence_metadata_by_alias("ncbi", "NC_000001.11").is_some());
    }

    #[test]
    fn test_store_alias_reverse_multiple_sequences() {
        use crate::digest::digest_sequence;

        let mut store = RefgetStore::in_memory();
        let r1 = digest_sequence("chr1", b"ACGT");
        let r2 = digest_sequence("chr2", b"TTTT");
        store.add_sequence_record(r1.clone(), false).unwrap();
        store.add_sequence_record(r2.clone(), false).unwrap();

        let d1 = r1.metadata().sha512t24u.clone();
        let d2 = r2.metadata().sha512t24u.clone();

        store.add_sequence_alias("ncbi", "NC_000001.11", &d1).unwrap();
        store.add_sequence_alias("ucsc", "chr1", &d1).unwrap();
        store.add_sequence_alias("ncbi", "NC_000002.12", &d2).unwrap();

        let aliases = store.get_aliases_for_sequence(&d1);
        assert_eq!(aliases.len(), 2);

        let aliases = store.get_aliases_for_sequence(&d2);
        assert_eq!(aliases.len(), 1);
    }

    #[test]
    fn test_store_alias_write_store_to_dir() {
        use crate::digest::digest_sequence;

        let dir = tempdir().unwrap();
        let store_path = dir.path().join("store");

        let mut store = RefgetStore::in_memory();
        let record = digest_sequence("chr1", b"ACGT");
        store.add_sequence_record(record.clone(), false).unwrap();
        let digest = record.metadata().sha512t24u.clone();

        store.add_sequence_alias("ncbi", "NC_000001.11", &digest).unwrap();

        store.write_store_to_dir(&store_path, None).unwrap();

        assert!(store_path.join("aliases/sequences/ncbi.tsv").exists());

        let store2 = RefgetStore::open_local(&store_path).unwrap();
        assert!(store2.get_sequence_metadata_by_alias("ncbi", "NC_000001.11").is_some());
    }

    #[test]
    fn test_get_sequence_metadata_by_alias() {
        use crate::digest::digest_sequence;

        let mut store = RefgetStore::in_memory();
        let record = digest_sequence("chr1", b"ACGT");
        store.add_sequence_record(record.clone(), false).unwrap();
        let digest = record.metadata().sha512t24u.clone();

        store.add_sequence_alias("ncbi", "NC_000001.11", &digest).unwrap();

        let meta = store.get_sequence_metadata_by_alias("ncbi", "NC_000001.11").unwrap();
        assert_eq!(meta.name, "chr1");
        assert_eq!(meta.length, 4);
    }

    #[test]
    fn test_get_sequence_by_alias_loads_data() {
        use crate::digest::digest_sequence;

        let mut store = RefgetStore::in_memory();
        let record = digest_sequence("chr1", b"ACGT");
        store.add_sequence_record(record.clone(), false).unwrap();
        let digest = record.metadata().sha512t24u.clone();

        store.add_sequence_alias("ncbi", "NC_000001.11", &digest).unwrap();

        let rec = store.get_sequence_by_alias("ncbi", "NC_000001.11").unwrap();
        assert_eq!(rec.metadata().name, "chr1");
    }

    #[test]
    #[cfg(feature = "filesystem")]
    fn test_get_collection_metadata_by_alias() {
        let temp = tempdir().unwrap();
        let fasta_path = copy_test_fasta(temp.path(), "base.fa");

        let mut store = RefgetStore::in_memory();
        let (meta, _) = store.add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new()).unwrap();

        store.add_collection_alias("ucsc", "hg38", &meta.digest).unwrap();

        let coll_meta = store.get_collection_metadata_by_alias("ucsc", "hg38").unwrap();
        assert_eq!(coll_meta.digest, meta.digest);
    }

    #[test]
    #[cfg(feature = "filesystem")]
    fn test_get_collection_by_alias_loads() {
        let temp = tempdir().unwrap();
        let fasta_path = copy_test_fasta(temp.path(), "base.fa");

        let mut store = RefgetStore::in_memory();
        let (meta, _) = store.add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new()).unwrap();

        store.add_collection_alias("ucsc", "hg38", &meta.digest).unwrap();

        let coll = store.get_collection_by_alias("ucsc", "hg38").unwrap();
        assert_eq!(coll.metadata.digest, meta.digest);
        assert!(!coll.sequences.is_empty());
    }

    #[test]
    fn test_get_sequence_by_alias_not_found() {
        let store = RefgetStore::in_memory();
        assert!(store.get_sequence_metadata_by_alias("ncbi", "nonexistent").is_none());
    }

    #[test]
    fn test_get_sequence_by_alias_error_not_found() {
        let store = RefgetStore::in_memory();
        assert!(store.get_sequence_by_alias("ncbi", "nonexistent").is_err());
    }

    #[test]
    #[cfg(feature = "filesystem")]
    fn test_fasta_load_with_namespace_aliases() {
        let dir = tempdir().unwrap();
        let fasta = dir.path().join("test.fa");
        fs::write(
            &fasta,
            ">chr1 ncbi:NC_000001.11 refseq:NC_000001.11\nACGT\n>chr2 ncbi:NC_000002.12\nTGCA\n",
        )
        .unwrap();

        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new().namespaces(&["ncbi", "refseq"]))
            .unwrap();

        let result = store.get_sequence_by_alias("ncbi", "NC_000001.11");
        assert!(result.is_ok());
        assert_eq!(result.unwrap().metadata().name, "chr1");

        let result = store.get_sequence_by_alias("refseq", "NC_000001.11");
        assert!(result.is_ok());

        let result = store.get_sequence_by_alias("ncbi", "NC_000002.12");
        assert!(result.is_ok());
        assert_eq!(result.unwrap().metadata().name, "chr2");

        let result = store.get_sequence_by_alias("ncbi", "NC_999999.1");
        assert!(result.is_err());
    }

    #[test]
    #[cfg(feature = "filesystem")]
    fn test_fasta_load_without_namespaces_no_aliases() {
        let dir = tempdir().unwrap();
        let fasta = dir.path().join("test.fa");
        fs::write(&fasta, ">chr1 ncbi:NC_000001.11\nACGT\n").unwrap();

        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
            .unwrap();

        let result = store.get_sequence_by_alias("ncbi", "NC_000001.11");
        assert!(result.is_err());
    }

    #[test]
    #[cfg(feature = "filesystem")]
    fn test_remove_collection_cleans_up_aliases() {
        let dir = tempdir().unwrap();
        let fasta = dir.path().join("test.fa");
        fs::write(&fasta, ">chr1\nACGT\n").unwrap();

        let mut store = RefgetStore::in_memory();
        let (meta, _) = store
            .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
            .unwrap();
        let digest = meta.digest;

        store.add_collection_alias("ucsc", "hg38", &digest).unwrap();
        assert!(store.get_collection_metadata_by_alias("ucsc", "hg38").is_some());

        store.remove_collection(&digest, false).unwrap();
        assert!(store.get_collection_metadata_by_alias("ucsc", "hg38").is_none());
    }

    #[test]
    #[cfg(feature = "filesystem")]
    fn test_manifest_namespace_roundtrip() {
        let dir = tempdir().unwrap();
        let store_path = dir.path().join("store");

        let mut store = RefgetStore::in_memory();
        let fasta_path = dir.path().join("test.fa");
        fs::write(&fasta_path, ">seq1\nACGT\n").unwrap();
        store
            .add_sequence_collection_from_fasta(fasta_path.to_str().unwrap(), FastaImportOptions::new())
            .unwrap();

        let seq_digest = key_to_digest_string(&store.sequence_digests().next().unwrap());
        store.add_sequence_alias("ncbi", "NC_000001.11", &seq_digest).unwrap();
        let coll_digest = {
            let key = *store.collections.keys().next().unwrap();
            key_to_digest_string(&key)
        };
        store.add_collection_alias("ucsc", "hg38", &coll_digest).unwrap();

        store.write_store_to_dir(&store_path, None).unwrap();

        let json_str = fs::read_to_string(store_path.join("rgstore.json")).unwrap();
        let metadata: serde_json::Value = serde_json::from_str(&json_str).unwrap();
        assert!(metadata["sequence_alias_namespaces"].as_array().unwrap().iter().any(|v| v.as_str() == Some("ncbi")));
        assert!(metadata["collection_alias_namespaces"].as_array().unwrap().iter().any(|v| v.as_str() == Some("ucsc")));

        let store2 = RefgetStore::open_local(&store_path).unwrap();
        let available = store2.available_alias_namespaces();
        assert!(available.sequences.contains(&"ncbi".to_string()));
        assert!(available.collections.contains(&"ucsc".to_string()));
    }

    #[test]
    #[cfg(feature = "filesystem")]
    fn test_manifest_empty_namespaces_not_serialized() {
        let dir = tempdir().unwrap();
        let store_path = dir.path().join("store");

        let mut store = RefgetStore::in_memory();
        let fasta_path = dir.path().join("test.fa");
        fs::write(&fasta_path, ">seq1\nACGT\n").unwrap();
        store
            .add_sequence_collection_from_fasta(fasta_path.to_str().unwrap(), FastaImportOptions::new())
            .unwrap();
        store.write_store_to_dir(&store_path, None).unwrap();

        let json_str = fs::read_to_string(store_path.join("rgstore.json")).unwrap();
        assert!(!json_str.contains("sequence_alias_namespaces"));
        assert!(!json_str.contains("collection_alias_namespaces"));
    }

    #[test]
    #[cfg(feature = "filesystem")]
    fn test_old_rgstore_json_without_namespaces() {
        let dir = tempdir().unwrap();
        let store_path = dir.path().join("store");
        fs::create_dir_all(&store_path).unwrap();

        let mut store = RefgetStore::in_memory();
        let fasta_path = dir.path().join("test.fa");
        fs::write(&fasta_path, ">seq1\nACGT\n").unwrap();
        store
            .add_sequence_collection_from_fasta(fasta_path.to_str().unwrap(), FastaImportOptions::new())
            .unwrap();
        store.write_store_to_dir(&store_path, None).unwrap();

        let store2 = RefgetStore::open_local(&store_path).unwrap();
        let available = store2.available_alias_namespaces();
        assert!(available.sequences.is_empty());
        assert!(available.collections.is_empty());
    }

    /// Fix #1: adding an alias to a disk-backed store re-writes rgstore.json so
    /// the served manifest advertises the new namespace (no separate index write
    /// required). Previously the alias TSV was persisted but the manifest was
    /// left stale, so remote/HTTP clients could never discover the namespace.
    #[test]
    #[cfg(feature = "filesystem")]
    fn test_alias_add_refreshes_manifest() {
        let dir = tempdir().unwrap();
        let store_path = dir.path().join("store");

        let mut store = RefgetStore::on_disk(&store_path).unwrap();
        let fasta_path = dir.path().join("test.fa");
        fs::write(&fasta_path, ">chr1\nACGT\n").unwrap();
        store
            .add_sequence_collection_from_fasta(fasta_path.to_str().unwrap(), FastaImportOptions::new())
            .unwrap();

        let seq_digest = key_to_digest_string(&store.sequence_digests().next().unwrap());
        // The manifest exists (written when the collection was added) but lists no
        // sequence namespaces yet.
        let before: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(store_path.join("rgstore.json")).unwrap()).unwrap();
        assert!(before.get("sequence_alias_namespaces").is_none());

        // Adding the alias must refresh the manifest in place.
        store.add_sequence_alias("ucsc", "chr1", &seq_digest).unwrap();

        let after: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(store_path.join("rgstore.json")).unwrap()).unwrap();
        assert!(after["sequence_alias_namespaces"]
            .as_array()
            .unwrap()
            .iter()
            .any(|v| v.as_str() == Some("ucsc")));

        // And a fresh manifest-driven open discovers and resolves it.
        let store2 = RefgetStore::open_local(&store_path).unwrap();
        assert!(store2.available_alias_namespaces().sequences.contains(&"ucsc".to_string()));
        assert!(store2.get_sequence_metadata_by_alias("ucsc", "chr1").is_some());
    }

    /// Fix #2: with the manifest as the single source of truth, a stale manifest
    /// (alias files present on disk but not advertised) yields the SAME result on
    /// a local filesystem as over HTTP — the aliases are not discovered and do not
    /// resolve. This is the failure that used to be silently masked by scanning
    /// the aliases/ directory; making it visible locally is the point.
    #[test]
    #[cfg(feature = "filesystem")]
    fn test_stale_manifest_hides_aliases_locally() {
        let dir = tempdir().unwrap();
        let store_path = dir.path().join("store");

        let mut store = RefgetStore::in_memory();
        let fasta_path = dir.path().join("test.fa");
        fs::write(&fasta_path, ">chr1\nACGT\n").unwrap();
        store
            .add_sequence_collection_from_fasta(fasta_path.to_str().unwrap(), FastaImportOptions::new())
            .unwrap();
        let seq_digest = key_to_digest_string(&store.sequence_digests().next().unwrap());
        store.add_sequence_alias("ucsc", "chr1", &seq_digest).unwrap();
        store.write_store_to_dir(&store_path, None).unwrap();

        // The alias TSV is on disk and the manifest correctly lists it.
        assert!(store_path.join("aliases/sequences/ucsc.tsv").exists());

        // Simulate the stale-manifest bug: drop the advertised namespaces while
        // leaving the alias files in place.
        let mut metadata: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(store_path.join("rgstore.json")).unwrap()).unwrap();
        metadata.as_object_mut().unwrap().remove("sequence_alias_namespaces");
        fs::write(store_path.join("rgstore.json"), serde_json::to_string_pretty(&metadata).unwrap()).unwrap();

        // Manifest-driven open: the un-advertised alias is invisible and unresolvable,
        // exactly as it would be for a remote client that cannot list the directory.
        let store2 = RefgetStore::open_local(&store_path).unwrap();
        assert!(store2.available_alias_namespaces().sequences.is_empty());
        assert!(store2.get_sequence_metadata_by_alias("ucsc", "chr1").is_none());
    }

    // -----------------------------------------------------------------------
    // KeepOurs sync strategy tests (regression test for was_local ordering bug)
    // -----------------------------------------------------------------------

    /// Pull an alias file that does NOT exist locally yet.
    /// KeepOurs: first pull should count as `pulled`, second pull as `skipped`.
    #[test]
    #[cfg(all(feature = "filesystem", feature = "http"))]
    fn test_keep_ours_alias_first_pull_counts_as_pulled() {
        // "Remote" store: a directory with pre-built alias TSV files.
        // pull_aliases(Some("ncbi"), ...) pulls BOTH sequences and collections
        // namespace "ncbi", so we create both files on the remote.
        let remote_dir = tempdir().unwrap();
        let seq_dir = remote_dir.path().join("aliases").join("sequences");
        let coll_dir = remote_dir.path().join("aliases").join("collections");
        fs::create_dir_all(&seq_dir).unwrap();
        fs::create_dir_all(&coll_dir).unwrap();
        fs::write(seq_dir.join("ncbi.tsv"), "NC_000001.11\tsome_digest\n").unwrap();
        fs::write(coll_dir.join("ncbi.tsv"), "hg38\tcoll_digest\n").unwrap();

        // Start a local HTTP server serving the remote_dir.
        let (base_url, shutdown) = start_file_server(remote_dir.path().to_path_buf());

        // "Local" store: empty disk-backed store that will pull from the server.
        let local_dir = tempdir().unwrap();
        let local_store_path = local_dir.path().join("store");
        fs::create_dir_all(&local_store_path).unwrap();

        let mut store = RefgetStore::on_disk(&local_store_path).unwrap();
        store.inner.fetcher = Some(std::sync::Arc::new(
            super::fetch::Fetcher::new(&base_url, super::FetchConfig::default()).unwrap(),
        ));
        store.inner.remote_source = Some(base_url);
        // Advertise "ncbi" as an available namespace for both sequences and collections.
        // pull_aliases(Some("ncbi"), ...) pulls both kinds, so we need both files on remote.
        store.inner.available_sequence_alias_namespaces = vec!["ncbi".to_string()];
        store.inner.available_collection_alias_namespaces = vec!["ncbi".to_string()];

        // First pull: both alias files do not exist locally → both should be pulled.
        // pull_aliases(Some("ncbi"), ...) pulls sequences/ncbi.tsv AND collections/ncbi.tsv.
        let result = store.pull_aliases(Some("ncbi"), SyncStrategy::KeepOurs).unwrap();
        assert_eq!(result.pulled, 2, "first pull should count both files as pulled, not skipped");
        assert_eq!(result.skipped, 0, "first pull should not be skipped");
        assert_eq!(result.not_found, 0);

        // Second pull: both files now exist locally → both should be skipped.
        let result2 = store.pull_aliases(Some("ncbi"), SyncStrategy::KeepOurs).unwrap();
        assert_eq!(result2.skipped, 2, "second pull should skip both files (already local)");
        assert_eq!(result2.pulled, 0, "second pull should not count any files as pulled");

        shutdown();
    }
}
