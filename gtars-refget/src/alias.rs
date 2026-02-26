//! Alias management for human-readable sequence and collection names.
//!
//! This module is self-contained: it defines the `AliasManager` struct and
//! provides all alias CRUD and persistence operations. RefgetStore holds an
//! `AliasManager` field and delegates alias logic here.

use std::collections::HashMap;
use std::fs::{self, File, create_dir_all};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::hashkeyable::{HashKeyable, key_to_digest_string};

/// Identifies whether an alias targets a sequence or a collection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

/// Generic alias store: namespace → { alias → digest_key }
type AliasStore = HashMap<String, HashMap<String, [u8; 32]>>;

/// Manages human-readable aliases for sequences and collections.
///
/// Each alias lives in a namespace (e.g. "ncbi", "ucsc") and maps to a
/// binary digest key. Supports forward lookup, reverse lookup, add/remove,
/// TSV file loading, and directory-based persistence.
#[derive(Debug)]
pub struct AliasManager {
    sequence_aliases: AliasStore,
    collection_aliases: AliasStore,
}

// =========================================================================
// Private helpers (operate on any AliasStore)
// =========================================================================

/// Add an alias to an AliasStore.
fn alias_add(store: &mut AliasStore, namespace: &str, alias: &str, digest: [u8; 32]) {
    store
        .entry(namespace.to_string())
        .or_default()
        .insert(alias.to_string(), digest);
}

/// Forward lookup: resolve alias → digest key.
fn alias_resolve(store: &AliasStore, namespace: &str, alias: &str) -> Option<[u8; 32]> {
    store.get(namespace).and_then(|ns| ns.get(alias)).copied()
}

/// Reverse lookup (scan): find all aliases pointing to a given digest.
/// Returns Vec<(namespace, alias)>.
fn alias_reverse_scan(store: &AliasStore, digest: &[u8; 32]) -> Vec<(String, String)> {
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

/// List all namespaces in an AliasStore.
fn alias_namespaces(store: &AliasStore) -> Vec<String> {
    store.keys().cloned().collect()
}

/// List all aliases in a namespace.
fn alias_list(store: &AliasStore, namespace: &str) -> Option<Vec<String>> {
    store
        .get(namespace)
        .map(|ns| ns.keys().cloned().collect())
}

/// Remove a single alias. Returns true if it existed.
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


/// Load aliases from a TSV file into an AliasStore namespace.
/// Format: alias\tdigest per line. Lines starting with '#' are comments.
fn alias_load_tsv(store: &mut AliasStore, namespace: &str, path: &Path) -> Result<usize> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut count = 0;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.splitn(2, '\t').collect();
        if parts.len() == 2 {
            let key = parts[1].to_key();
            alias_add(store, namespace, parts[0], key);
            count += 1;
        }
    }
    Ok(count)
}

/// Load all alias TSV files from a directory into an AliasStore.
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
            alias_load_tsv(store, &namespace, &path)?;
        }
    }
    Ok(())
}

/// Write all aliases in an AliasStore to TSV files in a directory.
fn write_all_aliases(store: &AliasStore, dir: &Path) -> Result<()> {
    if store.is_empty() {
        return Ok(());
    }
    create_dir_all(dir)?;
    for (namespace, aliases) in store {
        let tsv_path = dir.join(format!("{}.tsv", namespace));
        let mut file = File::create(&tsv_path)?;
        for (alias, digest) in aliases {
            writeln!(file, "{}\t{}", alias, key_to_digest_string(digest))?;
        }
    }
    Ok(())
}

// =========================================================================
// AliasManager public API
// =========================================================================

impl AliasManager {
    /// Create a new, empty AliasManager.
    pub fn new() -> Self {
        AliasManager {
            sequence_aliases: HashMap::new(),
            collection_aliases: HashMap::new(),
        }
    }

    /// Returns true if there are no aliases at all.
    pub fn is_empty(&self) -> bool {
        self.sequence_aliases.is_empty() && self.collection_aliases.is_empty()
    }

    // --- Sequence aliases ---

    /// Add a sequence alias: namespace/alias → digest key.
    pub fn add_sequence(&mut self, namespace: &str, alias: &str, digest: &str) {
        let key = digest.to_key();
        alias_add(&mut self.sequence_aliases, namespace, alias, key);
    }

    /// Forward lookup: resolve sequence alias → digest key.
    pub fn resolve_sequence(&self, namespace: &str, alias: &str) -> Option<[u8; 32]> {
        alias_resolve(&self.sequence_aliases, namespace, alias)
    }

    /// Reverse lookup: find all aliases pointing to this sequence digest.
    pub fn reverse_lookup_sequence(&self, digest: &str) -> Vec<(String, String)> {
        let key = digest.to_key();
        alias_reverse_scan(&self.sequence_aliases, &key)
    }

    /// List all sequence alias namespaces.
    pub fn sequence_namespaces(&self) -> Vec<String> {
        alias_namespaces(&self.sequence_aliases)
    }

    /// List all aliases in a sequence alias namespace.
    pub fn sequence_aliases(&self, namespace: &str) -> Option<Vec<String>> {
        alias_list(&self.sequence_aliases, namespace)
    }

    /// Remove a single sequence alias. Returns true if it existed.
    pub fn remove_sequence(&mut self, namespace: &str, alias: &str) -> bool {
        alias_remove(&mut self.sequence_aliases, namespace, alias)
    }

    /// Load sequence aliases from a TSV file into a namespace.
    pub fn load_sequence_tsv(&mut self, namespace: &str, path: &Path) -> Result<usize> {
        alias_load_tsv(&mut self.sequence_aliases, namespace, path)
    }

    // --- Collection aliases ---

    /// Add a collection alias: namespace/alias → digest key.
    pub fn add_collection(&mut self, namespace: &str, alias: &str, digest: &str) {
        let key = digest.to_key();
        alias_add(&mut self.collection_aliases, namespace, alias, key);
    }

    /// Forward lookup: resolve collection alias → digest key.
    pub fn resolve_collection(&self, namespace: &str, alias: &str) -> Option<[u8; 32]> {
        alias_resolve(&self.collection_aliases, namespace, alias)
    }

    /// Reverse lookup: find all aliases pointing to this collection digest.
    pub fn reverse_lookup_collection(&self, digest: &str) -> Vec<(String, String)> {
        let key = digest.to_key();
        alias_reverse_scan(&self.collection_aliases, &key)
    }

    /// List all collection alias namespaces.
    pub fn collection_namespaces(&self) -> Vec<String> {
        alias_namespaces(&self.collection_aliases)
    }

    /// List all aliases in a collection alias namespace.
    pub fn collection_aliases(&self, namespace: &str) -> Option<Vec<String>> {
        alias_list(&self.collection_aliases, namespace)
    }

    /// Remove a single collection alias. Returns true if it existed.
    pub fn remove_collection(&mut self, namespace: &str, alias: &str) -> bool {
        alias_remove(&mut self.collection_aliases, namespace, alias)
    }

    /// Load collection aliases from a TSV file into a namespace.
    pub fn load_collection_tsv(&mut self, namespace: &str, path: &Path) -> Result<usize> {
        alias_load_tsv(&mut self.collection_aliases, namespace, path)
    }

    // --- Persistence ---

    /// Load aliases from an aliases directory (with sequences/ and collections/ subdirs).
    pub fn load_from_dir(&mut self, aliases_dir: &Path) -> Result<()> {
        load_aliases_from_dir(&mut self.sequence_aliases, &aliases_dir.join("sequences"))?;
        load_aliases_from_dir(&mut self.collection_aliases, &aliases_dir.join("collections"))?;
        Ok(())
    }

    /// Write all aliases to an aliases directory (with sequences/ and collections/ subdirs).
    pub fn write_to_dir(&self, aliases_dir: &Path) -> Result<()> {
        write_all_aliases(&self.sequence_aliases, &aliases_dir.join("sequences"))?;
        write_all_aliases(&self.collection_aliases, &aliases_dir.join("collections"))?;
        Ok(())
    }

    /// Write a single alias namespace to its TSV sidecar file on disk.
    pub fn write_namespace(
        &self,
        aliases_dir: &Path,
        kind: AliasKind,
        namespace: &str,
    ) -> Result<()> {
        let dir = aliases_dir.join(kind.subdir());
        create_dir_all(&dir)?;

        let store = match kind {
            AliasKind::Sequence => &self.sequence_aliases,
            AliasKind::Collection => &self.collection_aliases,
        };

        let tsv_path = dir.join(format!("{}.tsv", namespace));
        if let Some(ns) = store.get(namespace) {
            let mut file = File::create(&tsv_path)?;
            for (alias, digest) in ns {
                writeln!(file, "{}\t{}", alias, key_to_digest_string(digest))?;
            }
        } else {
            // Namespace was cleared — remove the file
            let _ = fs::remove_file(&tsv_path);
        }
        Ok(())
    }
}

impl Default for AliasManager {
    fn default() -> Self {
        Self::new()
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_add_and_resolve_sequence() {
        let mut mgr = AliasManager::new();
        assert!(mgr.is_empty());

        mgr.add_sequence("ncbi", "NC_000001.11", "abc123");
        assert!(!mgr.is_empty());

        let key = mgr.resolve_sequence("ncbi", "NC_000001.11");
        assert!(key.is_some());
        assert_eq!(key.unwrap(), "abc123".to_key());

        // Non-existent
        assert!(mgr.resolve_sequence("ncbi", "missing").is_none());
        assert!(mgr.resolve_sequence("missing_ns", "NC_000001.11").is_none());
    }

    #[test]
    fn test_reverse_lookup() {
        let mut mgr = AliasManager::new();
        mgr.add_sequence("ncbi", "NC_000001.11", "digest1");
        mgr.add_sequence("ucsc", "chr1", "digest1");
        mgr.add_sequence("ncbi", "NC_000002.12", "digest2");

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
        mgr.add_sequence("ncbi", "NC_000001.11", "d1");
        mgr.add_sequence("ucsc", "chr1", "d1");

        let ns = mgr.sequence_namespaces();
        assert!(ns.contains(&"ncbi".to_string()));
        assert!(ns.contains(&"ucsc".to_string()));

        let aliases = mgr.sequence_aliases("ncbi").unwrap();
        assert!(aliases.contains(&"NC_000001.11".to_string()));
    }

    #[test]
    fn test_remove() {
        let mut mgr = AliasManager::new();
        mgr.add_sequence("ncbi", "NC_000001.11", "d1");

        assert!(mgr.remove_sequence("ncbi", "NC_000001.11"));
        assert!(mgr.resolve_sequence("ncbi", "NC_000001.11").is_none());
        // Namespace should be cleaned up
        assert!(mgr.sequence_namespaces().is_empty());

        // Remove non-existent
        assert!(!mgr.remove_sequence("ncbi", "NC_000001.11"));
    }

    #[test]
    fn test_collection_aliases() {
        let mut mgr = AliasManager::new();
        mgr.add_collection("ucsc", "hg38", "coll_digest");
        mgr.add_collection("gencode", "GRCh38.p14", "coll_digest");

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
        mgr.add_sequence("ncbi", "NC_000001.11", "seq_digest");
        mgr.add_collection("ucsc", "hg38", "coll_digest");

        mgr.write_to_dir(&aliases_dir).unwrap();

        // Verify files exist
        assert!(aliases_dir.join("sequences/ncbi.tsv").exists());
        assert!(aliases_dir.join("collections/ucsc.tsv").exists());

        // Load into new manager
        let mut mgr2 = AliasManager::new();
        mgr2.load_from_dir(&aliases_dir).unwrap();

        assert!(mgr2.resolve_sequence("ncbi", "NC_000001.11").is_some());
        assert!(mgr2.resolve_collection("ucsc", "hg38").is_some());
    }

    #[test]
    fn test_load_from_missing_dir() {
        let mut mgr = AliasManager::new();
        // Should succeed silently
        mgr.load_from_dir(std::path::Path::new("/nonexistent/aliases")).unwrap();
        assert!(mgr.is_empty());
    }

    #[test]
    fn test_write_namespace_single() {
        let dir = tempdir().unwrap();
        let aliases_dir = dir.path().join("aliases");

        let mut mgr = AliasManager::new();
        mgr.add_sequence("ncbi", "NC_000001.11", "d1");
        mgr.add_sequence("ucsc", "chr1", "d2");

        // Write only ncbi namespace
        mgr.write_namespace(&aliases_dir, AliasKind::Sequence, "ncbi").unwrap();

        assert!(aliases_dir.join("sequences/ncbi.tsv").exists());
        // ucsc should NOT have been written
        assert!(!aliases_dir.join("sequences/ucsc.tsv").exists());
    }

    #[test]
    fn test_empty_write_is_noop() {
        let dir = tempdir().unwrap();
        let aliases_dir = dir.path().join("aliases");

        let mgr = AliasManager::new();
        mgr.write_to_dir(&aliases_dir).unwrap();

        // Directory should not have been created (empty stores skip writing)
        assert!(!aliases_dir.join("sequences").exists());
        assert!(!aliases_dir.join("collections").exists());
    }

    #[test]
    fn test_load_tsv() {
        let dir = tempdir().unwrap();
        let tsv_path = dir.path().join("ncbi.tsv");
        std::fs::write(&tsv_path, "NC_000001.11\tsome_digest\n# comment\n\nNC_000002.12\tanother_digest\n").unwrap();

        let mut mgr = AliasManager::new();
        let count = mgr.load_sequence_tsv("ncbi", &tsv_path).unwrap();
        assert_eq!(count, 2);
        assert!(mgr.resolve_sequence("ncbi", "NC_000001.11").is_some());
        assert!(mgr.resolve_sequence("ncbi", "NC_000002.12").is_some());
    }
}
