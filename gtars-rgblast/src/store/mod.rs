//! Opening a RefgetStore (local directory or http(s) prefix) through gtars-refget,
//! resolving a selection, loading its metadata, and reading sequence bytes.
//!
//! A remote store is never streamed: the sequences a command needs are first pulled into the
//! cache with gtars-refget's own pull steps (`load_collections`, `pull_fhr`,
//! `cache_sequences`, as `refget store pull` does), then read from local files like a local
//! store's.

pub mod access;
pub mod decode;
pub mod fetch;

use crate::error::Error;
use crate::plan::ResolvedSelection;
use crate::select::resolve::{Names, resolve_selector};
pub use fetch::{SeqFetch, StoreFiles, StoredBytes};
use gtars_refget::store::{
    FetchConfig, FetchCounters, FhrMetadata, PrefetchOptions, ReadonlyRefgetStore, RefgetStore,
    RemoteMissing, StorageMode, SyncStrategy,
};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StoreLocation {
    Local(PathBuf),
    /// `cache` is the shared cache root the store downloads into.
    Remote {
        url: String,
        cache: PathBuf,
    },
}

/// One row of a collection, in `.rgsi` order.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CollectionRow {
    pub name: String,
    pub length: u64,
    pub digest: String,
}

/// How a remote store is downloaded (`--fetch-jobs`). Local stores ignore it. Pulled
/// sequences always stay in the cache, so a rerun downloads only what is missing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct RemoteOptions {
    /// Parallel HTTP requests.
    pub fetch_jobs: usize,
}

impl Default for RemoteOptions {
    fn default() -> Self {
        RemoteOptions { fetch_jobs: 32 }
    }
}

impl RemoteOptions {
    /// gtars-refget's fetch configuration (everything else at its defaults).
    pub fn config(&self) -> FetchConfig {
        FetchConfig {
            jobs: self.fetch_jobs.max(1),
            ..FetchConfig::default()
        }
    }
}

/// Download sequences of a remote store into its cache with gtars-refget's
/// `cache_sequences`; sequences already there are not downloaded again.
pub(crate) fn pull(
    store: &ReadonlyRefgetStore,
    digests: &[String],
) -> Result<FetchCounters, Error> {
    store
        .cache_sequences(digests, PrefetchOptions::default())
        .map_err(|e| match e.downcast_ref::<RemoteMissing>() {
            Some(m) => Error::remote_not_found(&m.url, m.status),
            None => Error::Fetch {
                what: "sequence download".into(),
                msg: format!("{e:#}"),
            },
        })
}

pub struct OpenedStore {
    pub store: RefgetStore,
    pub location: StoreLocation,
    cache_root: PathBuf,
}

/// True for a store spec that names a remote store.
pub fn is_remote_spec(spec: &str) -> bool {
    spec.starts_with("http://") || spec.starts_with("https://")
}

impl OpenedStore {
    /// [`open_with`](Self::open_with) and the default fetch settings.
    pub fn open(spec: &str, cache_dir: Option<&Path>) -> Result<Self, Error> {
        Self::open_with(spec, cache_dir, &RemoteOptions::default())
    }

    /// Open a local directory, or an http(s) URL that downloads into the shared cache root
    /// (`cache_dir`, else gtars-refget's default root). A missing remote store is
    /// `Error::NotFound`.
    pub fn open_with(
        spec: &str,
        cache_dir: Option<&Path>,
        remote: &RemoteOptions,
    ) -> Result<Self, Error> {
        let cache_root = crate::cache_root::cache_root(cache_dir);
        if is_remote_spec(spec) {
            let url = gtars_refget::store::normalize_remote_url(spec)
                .unwrap_or_else(|_| spec.trim_end_matches('/').to_string());
            let cache = cache_root.clone();
            let mut store =
                RefgetStore::open_remote(&cache, &url, &remote.config()).map_err(|e| {
                    match e.downcast_ref::<RemoteMissing>() {
                        Some(m) => Error::remote_not_found(&m.url, m.status),
                        None => Error::Store(e),
                    }
                })?;
            store.set_quiet(true);
            Ok(OpenedStore {
                store,
                location: StoreLocation::Remote { url, cache },
                cache_root,
            })
        } else {
            let p = Path::new(spec);
            let not_found =
                |detail: &str| Error::NotFound(format!("no RefgetStore at {spec} ({detail})"));
            if !p.is_dir() {
                return Err(not_found("no such directory"));
            }
            if !p.join("rgstore.json").is_file() {
                return Err(not_found("rgstore.json not found"));
            }
            let mut store = RefgetStore::open_local(p)?;
            store.set_quiet(true);
            Ok(OpenedStore {
                store,
                location: StoreLocation::Local(p.to_path_buf()),
                cache_root,
            })
        }
    }

    /// The cache root this store was opened with (`--cache-dir`, else the shared default),
    /// local stores included: gtars-rgblast's own caches live under it.
    pub fn cache_root(&self) -> &Path {
        &self.cache_root
    }

    pub fn is_remote(&self) -> bool {
        matches!(self.location, StoreLocation::Remote { .. })
    }

    /// Pull sequences of a remote store into the cache (the sequence index must be loaded);
    /// `None` for a local store, which has its files already.
    pub fn pull(&self, digests: &[String]) -> Result<Option<FetchCounters>, Error> {
        if !self.is_remote() {
            return Ok(None);
        }
        pull(&self.store, digests).map(Some)
    }

    /// Directory name, or the URL's last path segment.
    pub fn label(&self) -> String {
        match &self.location {
            StoreLocation::Local(p) => {
                let p = std::fs::canonicalize(p).unwrap_or_else(|_| p.clone());
                p.file_name()
                    .map(|n| n.to_string_lossy().into_owned())
                    .unwrap_or_else(|| p.display().to_string())
            }
            StoreLocation::Remote { url, .. } => url.rsplit('/').next().unwrap_or(url).to_string(),
        }
    }

    /// The path or URL the store was opened from.
    pub fn location_string(&self) -> String {
        match &self.location {
            StoreLocation::Local(p) => p.display().to_string(),
            StoreLocation::Remote { url, .. } => url.clone(),
        }
    }

    pub fn mode(&self) -> StorageMode {
        self.store.storage_mode()
    }

    /// A collection digest resolves to itself; `NS:ALIAS` is looked up in that namespace when
    /// `NS` is a collection alias namespace of the store; otherwise the key is a bare alias that
    /// must be unique across namespaces.
    pub fn resolve_collection(&self, key: &str) -> Result<String, Error> {
        resolve_selector(self, key)
    }

    /// FHR metadata of a collection: the local sidecar, or (remote stores) the sidecar pulled
    /// into the cache (`KeepOurs` reads one already cached). None when the collection has none.
    pub fn fhr(&mut self, digest: &str) -> Result<Option<FhrMetadata>, Error> {
        if self.store.get_fhr_metadata(digest).is_none() && self.is_remote() {
            self.store
                .pull_fhr(Some(digest), SyncStrategy::KeepOurs)
                .map_err(Error::Store)?;
        }
        Ok(self.store.get_fhr_metadata(digest).cloned())
    }

    /// Load what `SubjectPlan::build` reads: the sequence index (md5 and alphabet of every
    /// digest; alias lookups need it too) and the selected collections (a remote store
    /// downloads the missing ones in parallel). FHR is made available for each selected
    /// collection (remote stores pull it into the cache).
    pub fn load_for(&mut self, sel: &ResolvedSelection) -> Result<(), Error> {
        self.store.load_sequence_index()?;
        self.store.load_collections(&sel.collections)?;
        for d in &sel.collections {
            self.fhr(d)?;
        }
        Ok(())
    }

    /// The content digests `rgstore.json` has (any of `sequences_digest`, `collections_digest`,
    /// `aliases_digest`, `fhr_digest`; a store without FHR has no `fhr_digest`).
    pub fn digests(&self) -> std::collections::BTreeMap<String, String> {
        let Ok(m) = self.store.store_metadata() else {
            return Default::default();
        };
        [
            "sequences_digest",
            "collections_digest",
            "aliases_digest",
            "fhr_digest",
        ]
        .iter()
        .filter_map(|k| m.get(*k).map(|v| (k.to_string(), v.clone())))
        .collect()
    }

    /// A collection's rows in `.rgsi` order, from metadata only. Works for any collection in
    /// the store, loading it first when needed.
    pub fn collection_rows(&mut self, digest: &str) -> Result<Vec<CollectionRow>, Error> {
        self.store.load_sequence_index()?;
        if !self.store.is_collection_loaded(digest) {
            self.store.load_collection(digest)?;
        }
        let c = self.store.get_collection(digest)?;
        Ok(c.sequences
            .iter()
            .map(|r| {
                let m = r.metadata();
                CollectionRow {
                    name: m.name.clone(),
                    length: m.length as u64,
                    digest: m.sha512t24u.clone(),
                }
            })
            .collect())
    }
}

impl Names for OpenedStore {
    fn label(&self) -> String {
        self.location_string()
    }

    fn has_collection(&self, digest: &str) -> bool {
        self.store.get_collection_metadata(digest).is_some()
    }

    fn namespaces(&self) -> Vec<String> {
        self.store.available_alias_namespaces().collections.to_vec()
    }

    fn lookup_alias(&self, ns: &str, alias: &str) -> Option<String> {
        self.store
            .get_collection_metadata_by_alias(ns, alias)
            .map(|m| m.digest.clone())
    }
}
