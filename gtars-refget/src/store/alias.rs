//! Alias management for human-readable sequence and collection names.
//!
//! Contains AliasManager types, store bridge methods, RefgetStore delegates,
//! and remote pull methods.

use super::*;
use super::readonly::ReadonlyRefgetStore;
use super::core::RefgetStore;

use std::collections::HashMap;
use std::fs::{self, File, create_dir_all};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::hashkeyable::{HashKeyable, key_to_digest_string};

// =========================================================================
// Types
// =========================================================================

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

/// Generic alias store: namespace -> { alias -> digest_key }
type AliasStore = HashMap<String, HashMap<String, [u8; 32]>>;

/// Manages human-readable aliases for sequences and collections.
#[derive(Debug)]
pub struct AliasManager {
    sequence_aliases: AliasStore,
    collection_aliases: AliasStore,
}

// =========================================================================
// Private helpers (operate on any AliasStore)
// =========================================================================

fn alias_add(store: &mut AliasStore, namespace: &str, alias: &str, digest: [u8; 32]) {
    store
        .entry(namespace.to_string())
        .or_default()
        .insert(alias.to_string(), digest);
}

fn alias_resolve(store: &AliasStore, namespace: &str, alias: &str) -> Option<[u8; 32]> {
    store.get(namespace).and_then(|ns| ns.get(alias)).copied()
}

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

    pub fn add_sequence(&mut self, namespace: &str, alias: &str, digest: &str) {
        let key = digest.to_key();
        alias_add(&mut self.sequence_aliases, namespace, alias, key);
    }

    pub fn resolve_sequence(&self, namespace: &str, alias: &str) -> Option<[u8; 32]> {
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

    pub fn load_sequence_tsv(&mut self, namespace: &str, path: &Path) -> Result<usize> {
        alias_load_tsv(&mut self.sequence_aliases, namespace, path)
    }

    // --- Collection aliases ---

    pub fn add_collection(&mut self, namespace: &str, alias: &str, digest: &str) {
        let key = digest.to_key();
        alias_add(&mut self.collection_aliases, namespace, alias, key);
    }

    pub fn resolve_collection(&self, namespace: &str, alias: &str) -> Option<[u8; 32]> {
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

    pub fn load_collection_tsv(&mut self, namespace: &str, path: &Path) -> Result<usize> {
        alias_load_tsv(&mut self.collection_aliases, namespace, path)
    }

    // --- Persistence ---

    pub fn load_from_dir(&mut self, aliases_dir: &Path) -> Result<()> {
        load_aliases_from_dir(&mut self.sequence_aliases, &aliases_dir.join("sequences"))?;
        load_aliases_from_dir(&mut self.collection_aliases, &aliases_dir.join("collections"))?;
        Ok(())
    }

    pub fn write_to_dir(&self, aliases_dir: &Path) -> Result<()> {
        write_all_aliases(&self.sequence_aliases, &aliases_dir.join("sequences"))?;
        write_all_aliases(&self.collection_aliases, &aliases_dir.join("collections"))?;
        Ok(())
    }

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
// ReadonlyRefgetStore alias bridge methods
// =========================================================================

impl ReadonlyRefgetStore {
    /// Add a sequence alias and persist to disk if applicable.
    pub fn add_sequence_alias(&mut self, namespace: &str, alias: &str, digest: &str) -> Result<()> {
        self.aliases.add_sequence(namespace, alias, digest);
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
            self.persist_alias_namespace(AliasKind::Sequence, namespace)?;
        }
        Ok(removed)
    }

    /// Load sequence aliases from a TSV file into a namespace.
    pub fn load_sequence_aliases(&mut self, namespace: &str, path: &str) -> Result<usize> {
        let count = self.aliases.load_sequence_tsv(namespace, Path::new(path))?;
        self.persist_alias_namespace(AliasKind::Sequence, namespace)?;
        Ok(count)
    }

    /// Add a collection alias and persist to disk if applicable.
    pub fn add_collection_alias(&mut self, namespace: &str, alias: &str, digest: &str) -> Result<()> {
        self.aliases.add_collection(namespace, alias, digest);
        self.persist_alias_namespace(AliasKind::Collection, namespace)?;
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
            self.persist_alias_namespace(AliasKind::Collection, namespace)?;
        }
        Ok(removed)
    }

    /// Load collection aliases from a TSV file into a namespace.
    pub fn load_collection_aliases(&mut self, namespace: &str, path: &str) -> Result<usize> {
        let count = self.aliases.load_collection_tsv(namespace, Path::new(path))?;
        self.persist_alias_namespace(AliasKind::Collection, namespace)?;
        Ok(count)
    }

    /// Write a single alias namespace to disk (if disk-backed).
    pub(crate) fn persist_alias_namespace(&self, kind: AliasKind, namespace: &str) -> Result<()> {
        if self.persist_to_disk {
            if let Some(ref local_path) = self.local_path {
                let aliases_dir = local_path.join("aliases");
                self.aliases.write_namespace(&aliases_dir, kind, namespace)?;
            }
        }
        Ok(())
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
        let mut result = PullResult::default();

        let seq_namespaces: Vec<String> = match namespace {
            Some(ns) => vec![ns.to_string()],
            None => self.inner.available_sequence_alias_namespaces.clone(),
        };
        let coll_namespaces: Vec<String> = match namespace {
            Some(ns) => vec![ns.to_string()],
            None => self.inner.available_collection_alias_namespaces.clone(),
        };

        for ns in &seq_namespaces {
            self.pull_alias_file(ns, "sequences", &strategy, &mut result)?;
        }

        for ns in &coll_namespaces {
            self.pull_alias_file(ns, "collections", &strategy, &mut result)?;
        }

        if result.pulled > 0 {
            if let Some(ref local_path) = self.inner.local_path {
                let aliases_dir = local_path.join("aliases");
                self.inner.aliases = AliasManager::default();
                self.inner.aliases.load_from_dir(&aliases_dir)?;
            }
        }

        Ok(result)
    }

    /// Pull a single alias TSV file.
    fn pull_alias_file(
        &self,
        namespace: &str,
        kind: &str,
        strategy: &SyncStrategy,
        result: &mut PullResult,
    ) -> Result<()> {
        let relative_path = format!("aliases/{}/{}.tsv", kind, namespace);

        match strategy {
            SyncStrategy::KeepOurs => {
                let was_local = self
                    .inner
                    .local_path
                    .as_ref()
                    .map(|p| p.join(&relative_path).exists())
                    .unwrap_or(false);
                match ReadonlyRefgetStore::fetch_file(
                    &self.inner.local_path,
                    &self.inner.remote_source,
                    &relative_path,
                    self.inner.persist_to_disk,
                    false,
                ) {
                    Ok(_) => {
                        if was_local {
                            result.skipped += 1;
                        } else {
                            result.pulled += 1;
                        }
                    }
                    Err(_) => {
                        result.not_found += 1;
                    }
                }
            }
            SyncStrategy::KeepTheirs => {
                match ReadonlyRefgetStore::fetch_file(
                    &self.inner.local_path,
                    &self.inner.remote_source,
                    &relative_path,
                    self.inner.persist_to_disk,
                    true,
                ) {
                    Ok(_) => {
                        result.pulled += 1;
                    }
                    Err(_) => {
                        result.not_found += 1;
                    }
                }
            }
            SyncStrategy::Notify => {
                let local_exists = self
                    .inner
                    .local_path
                    .as_ref()
                    .map(|p| p.join(&relative_path).exists())
                    .unwrap_or(false);

                if local_exists {
                    match ReadonlyRefgetStore::fetch_file(
                        &None,
                        &self.inner.remote_source,
                        &relative_path,
                        false,
                        false,
                    ) {
                        Ok(remote_data) => {
                            let local_path = self
                                .inner
                                .local_path
                                .as_ref()
                                .unwrap()
                                .join(&relative_path);
                            let local_data = fs::read(&local_path)?;
                            if local_data != remote_data {
                                result.conflicts.push(relative_path);
                            } else {
                                result.skipped += 1;
                            }
                        }
                        Err(_) => {
                            result.not_found += 1;
                        }
                    }
                } else {
                    result.conflicts.push(relative_path);
                }
            }
        }

        Ok(())
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

        mgr.add_sequence("ncbi", "NC_000001.11", "abc123");
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
        assert!(mgr.sequence_namespaces().is_empty());

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
        mgr.add_sequence("ncbi", "NC_000001.11", "d1");
        mgr.add_sequence("ucsc", "chr1", "d2");

        mgr.write_namespace(&aliases_dir, AliasKind::Sequence, "ncbi").unwrap();

        assert!(aliases_dir.join("sequences/ncbi.tsv").exists());
        assert!(!aliases_dir.join("sequences/ucsc.tsv").exists());
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
        let count = mgr.load_sequence_tsv("ncbi", &tsv_path).unwrap();
        assert_eq!(count, 2);
        assert!(mgr.resolve_sequence("ncbi", "NC_000001.11").is_some());
        assert!(mgr.resolve_sequence("ncbi", "NC_000002.12").is_some());
    }

    // --- Store-level alias integration tests ---

    fn copy_test_fasta(temp_dir: &std::path::Path, name: &str) -> std::path::PathBuf {
        let src = format!("../tests/data/fasta/{}", name);
        let dst = temp_dir.join(name);
        std::fs::copy(&src, &dst)
            .unwrap_or_else(|e| panic!("Failed to copy {} to tempdir: {}", src, e));
        dst
    }

    #[test]
    fn test_store_sequence_aliases() {
        use crate::collection::digest_sequence;

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
        use crate::collection::digest_sequence;

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
            seq_digest = store.list_sequences()[0].sha512t24u.clone();

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
        use crate::collection::digest_sequence;

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
        use crate::collection::digest_sequence;

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
        use crate::collection::digest_sequence;

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
        use crate::collection::digest_sequence;

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
        use crate::collection::digest_sequence;

        let mut store = RefgetStore::in_memory();
        let record = digest_sequence("chr1", b"ACGT");
        store.add_sequence_record(record.clone(), false).unwrap();
        let digest = record.metadata().sha512t24u.clone();

        store.add_sequence_alias("ncbi", "NC_000001.11", &digest).unwrap();

        let rec = store.get_sequence_by_alias("ncbi", "NC_000001.11").unwrap();
        assert_eq!(rec.metadata().name, "chr1");
    }

    #[test]
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

    // -----------------------------------------------------------------------
    // KeepOurs sync strategy tests (regression test for was_local ordering bug)
    // -----------------------------------------------------------------------

    /// Pull an alias file that does NOT exist locally yet.
    /// KeepOurs: first pull should count as `pulled`, second pull as `skipped`.
    #[test]
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
