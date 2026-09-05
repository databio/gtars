//! # RefgetStore
//!
//! A store for managing reference genome sequences with support for both
//! in-memory and disk-backed storage.
//!
//! ## Two-Type Design
//!
//! - **`RefgetStore`** (wrapper): User-facing type with `&mut self` read methods
//!   that automatically lazy-load data on first access. Use for CLI and scripts.
//! - **`ReadonlyRefgetStore`** (inner): All reads are `&self`. Suitable for
//!   `Arc<ReadonlyRefgetStore>` in servers. Requires explicit preloading.
//!
//! ## Sequence retrieval flows
//!
//! Sequence *bytes* can be served three ways, depending on where the data lives
//! (resident in RAM, a local `.seq` file, or a remote HTTP store) and how much
//! you want to move over the wire. All three are first-class; choose by access
//! pattern.
//!
//! 1. **Partial read** — `get_substring` / `get_substrings`. Returns the bases
//!    for `[start, end)` as a `String`, reading only the bytes that cover the
//!    region. Source resolution: resident `Full` bytes -> local `.seq`
//!    (positioned read; the whole sequence never enters RAM) -> remote
//!    byte-range (HTTP `Range:` via `open_remote_range`). The whole sequence
//!    never enters RAM and nothing is persisted, so this is the smart default
//!    for sparse, random-access extraction regardless of where the data lives.
//!    A remote-only `Stub` is served by fetching just the covering bytes; it
//!    does NOT trigger a whole-chromosome download (use flow 3 for that).
//!    Repeated remote reads re-fetch, so promote to flow 3 for repeat-heavy
//!    workloads.
//!
//! 2. **Streaming** — `stream_sequence`. Returns a `Read`er over `[start, end)`
//!    and decodes on the fly, so peak memory stays O(1) in the region length.
//!    Source resolution: resident `Full` bytes -> local `.seq` (seek) -> remote
//!    byte-range (HTTP `Range:` via `open_remote_range`). This is the only flow
//!    that pulls a *region* straight from a remote store, fetching just the
//!    covering bytes; nothing is persisted, so repeated remote reads re-fetch.
//!
//! 3. **Load & cache** — `load_sequence` / `load_all_sequences` (via
//!    `ensure_sequence_loaded`). Downloads the *whole* `.seq` once, persists it
//!    under the local cache (when `persist_to_disk` is set) and loads it into
//!    RAM. Subsequent reads are served resident (flow 1's first branch) with no
//!    further I/O. Best when the same sequence is read many times; costs a
//!    full-sequence download + allocation up front (a whole chromosome can be
//!    tens of MB encoded), so avoid it for one-off region pulls.
//!
//! Flows 1 and 2 need only metadata loaded (a `Stub`); flow 3 promotes a `Stub`
//! to a `Full` record. Remote byte-range (flows 1-2) requires the `http`
//! feature; without it, only resident and local sources are available.
//!
//! ## Write discipline (concurrency)
//!
//! A store directory may be written by several processes at once — separate
//! `gtars refget build` invocations, SLURM jobs, Python callers. Three rules make
//! that safe; all three live in Rust, because callers routinely bypass any
//! coordination arranged above it.
//!
//! 1. **Per-item files are digest-addressed and lock-free.** A `.seq` file, a
//!    `collections/<digest>.rgsi`, and an FHR sidecar are named by the digest of
//!    their own content, so two writers producing the same file produce the same
//!    bytes. These are written outside any lock — which is what keeps the
//!    expensive part of an import (parsing, digesting, writing hundreds of
//!    thousands of sequences) fully parallel.
//!
//! 2. **Writers commit a DELTA under an exclusive lock.** The shared artifacts —
//!    `sequences.rgsi`, `collections.rgci`, the alias TSVs, `rgstore.json` — are
//!    rewritten wholesale, so a writer that serialized its open-time snapshot
//!    would drop every row another writer added in the meantime, and resurrect
//!    every row another writer removed. Instead,
//!    [`ReadonlyRefgetStore::write_index_files`] takes [`lock::StoreLock`] on
//!    `<store>/.rgstore.lock`, re-reads the current on-disk rows, applies only
//!    the changes this handle actually made (see [`PendingChanges`]), and
//!    publishes. A row nobody
//!    touched is never rewritten. The lock is held for the commit only (seconds),
//!    never for the lifetime of the handle (hours) — with one deliberate
//!    exception, the orphan scan in
//!    [`ReadonlyRefgetStore::remove_collection`].
//!
//! 3. **Readers never lock.** Every whole-file write goes through
//!    [`atomic::atomic_write`] (temp file, `fsync`, `rename(2)`, `fsync` the
//!    directory), so a reader always sees a complete file. `rgstore.json` is
//!    published LAST in the commit sequence, after the files whose digests it
//!    advertises, so a manifest never describes bytes that have not landed.
//!
//! Transient files (`.rgstore.lock`, `.rgstore.lock.stale.*`, `.rgstore.tmp.*`)
//! may briefly appear in a store directory; anything mirroring or auditing the
//! directory must skip them (see [`atomic::is_transient_store_file`]).

mod readonly;
mod core;
mod alias;
pub(crate) mod atomic;
mod fhr_metadata;
mod lock;
// FASTA import (crossbeam-channel) and export (gtars-core) are filesystem-only.
#[cfg(feature = "filesystem")]
mod import;
mod persistence;
#[cfg(feature = "filesystem")]
mod export;

// The bulk of the store tests import FASTA via the filesystem-only API.
#[cfg(all(test, feature = "filesystem"))]
mod tests;

// Non-filesystem store tests: pure path/template/mode logic that needs no FASTA
// import, so they run under `--no-default-features`.
#[cfg(test)]
mod nofs_tests {
    use super::*;

    #[test]
    fn test_expand_template() {
        let digest = "ABCDEFghijklmnop";

        let result = ReadonlyRefgetStore::expand_template(digest, "sequences/%s2/%s.seq");
        assert_eq!(result, std::path::PathBuf::from("sequences/AB/ABCDEFghijklmnop.seq"));

        let result = ReadonlyRefgetStore::expand_template(digest, "sequences/%s2/%s4/%s.seq");
        assert_eq!(result, std::path::PathBuf::from("sequences/AB/ABCD/ABCDEFghijklmnop.seq"));

        let result = ReadonlyRefgetStore::expand_template(digest, "sequences/%s.seq");
        assert_eq!(result, std::path::PathBuf::from("sequences/ABCDEFghijklmnop.seq"));
    }

    #[test]
    fn test_sanitize_relative_path() {
        // Rejects traversal
        assert!(ReadonlyRefgetStore::sanitize_relative_path("../etc/passwd").is_err());
        assert!(ReadonlyRefgetStore::sanitize_relative_path("foo/../bar").is_err());
        assert!(ReadonlyRefgetStore::sanitize_relative_path("foo/../../bar").is_err());
        assert!(ReadonlyRefgetStore::sanitize_relative_path("..").is_err());

        // Rejects absolute
        assert!(ReadonlyRefgetStore::sanitize_relative_path("/etc/passwd").is_err());
        assert!(ReadonlyRefgetStore::sanitize_relative_path("\\windows\\system32").is_err());

        // Accepts valid
        assert!(ReadonlyRefgetStore::sanitize_relative_path("sequences/ab/abc123.seq").is_ok());
        assert!(ReadonlyRefgetStore::sanitize_relative_path("collections/xyz.rgsi").is_ok());
        assert!(ReadonlyRefgetStore::sanitize_relative_path("rgstore.json").is_ok());
        assert!(ReadonlyRefgetStore::sanitize_relative_path("sequences/%s2/%s.seq").is_ok());
    }

    #[test]
    fn test_mode_basics() {
        let mut store = RefgetStore::in_memory();

        assert_eq!(store.mode, StorageMode::Encoded);

        store.disable_encoding();
        assert_eq!(store.mode, StorageMode::Raw);
        store.enable_encoding();
        assert_eq!(store.mode, StorageMode::Encoded);

        store.set_encoding_mode(StorageMode::Raw);
        assert_eq!(store.mode, StorageMode::Raw);
        store.set_encoding_mode(StorageMode::Encoded);
        assert_eq!(store.mode, StorageMode::Encoded);
    }
}

// Re-export public types from submodules
pub use self::readonly::{PendingChanges, ReadonlyRefgetStore};
pub use self::core::RefgetStore;
pub use self::alias::{AliasKind, AliasManager};
pub use self::atomic::is_transient_store_file;
pub use self::lock::{LockInfo, LockOptions, StoreLock, force_unlock, lock_status};
pub use self::fhr_metadata::{
    FhrMetadata, FhrAuthor, FhrIdentifier, FhrTaxon, FhrVitalStats,
    // Disk I/O helpers used by persistence and externally
    load_sidecars, write_sidecars, write_sidecar, remove_sidecar, sidecar_path, load_from_json,
};

use serde::{Deserialize, Serialize};
use std::io::{BufReader, Read};
// `BufRead` is only used by the filesystem-only FASTA import path (via `super::*`).
#[cfg(feature = "filesystem")]
#[allow(unused_imports)]
use std::io::BufRead;

pub(crate) use crate::hashkeyable::DigestKey;
use crate::digest::SequenceCollectionMetadata;


// =========================================================================
// Shared constants
// =========================================================================

pub(crate) const DEFAULT_SEQDATA_PATH_TEMPLATE: &str = "sequences/%s2/%s.seq";

// =========================================================================
// Shared types used across multiple submodules
// =========================================================================

/// Paginated result container matching the seqcol spec response format.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PagedResult<T> {
    pub results: Vec<T>,
    pub pagination: Pagination,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pagination {
    pub page: usize,
    pub page_size: usize,
    pub total: usize,
}

/// Enum storing whether sequences will be stored in Raw or Encoded form
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq)]
pub enum StorageMode {
    Raw,
    Encoded,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RetrievedSequence {
    pub sequence: String,
    pub chrom_name: String,
    pub start: u32,
    pub end: u32,
}

/// Options for importing a FASTA file into a RefgetStore.
///
/// ## Single-knob parallelism model
///
/// `jobs` is the number of input FASTA files imported concurrently. `0` = auto
/// (`std::thread::available_parallelism`). `1` = serial. It has no effect with a
/// single input file (each file is processed by the basic read->digest->encode
/// pipeline either way).
///
/// ## Two independent alias indexes
///
/// A store keeps *two* separate alias indexes, and these options feed exactly
/// one each:
///
/// - `namespaces` feeds only the **sequence** alias index. It is a *parsing*
///   directive: it lists which `ns:value` tokens to harvest out of each FASTA
///   header and register per sequence.
/// - `collection_alias` feeds only the **collection** alias index. It is an
///   *assertion* by the caller naming the collection as a whole, because a FASTA
///   says what each sequence is called but never what the assembly is called.
///
/// They are deliberately decoupled: setting `namespaces` never names the
/// collection, and setting `collection_alias` never registers a sequence alias.
#[derive(Clone, Copy)]
pub struct FastaImportOptions<'a> {
    pub(crate) force: bool,
    pub(crate) namespaces: &'a [&'a str],
    /// Number of input FASTA files imported concurrently. Each file gets its own
    /// decoder so gzip decompression is parallelized across files. `0` = auto
    /// (`std::thread::available_parallelism`), `1` = serial. No effect with a
    /// single input file.
    pub(crate) jobs: usize,
    /// Optional `(namespace, alias)` under which the imported collection is
    /// registered in the **collection** alias index (e.g. `("ucsc", "hg38")`).
    ///
    /// Independent of `namespaces`, which only controls per-sequence FASTA
    /// header parsing and never names the collection. Rejected with an error
    /// when more than one FASTA file is imported in a single call, since one
    /// alias cannot name N collections. When unset, nothing is registered and
    /// nothing is derived from the filename.
    pub(crate) collection_alias: Option<(&'a str, &'a str)>,
}

impl<'a> Default for FastaImportOptions<'a> {
    fn default() -> Self {
        Self {
            force: false,
            namespaces: &[],
            jobs: 0,
            collection_alias: None,
        }
    }
}

impl<'a> FastaImportOptions<'a> {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    #[must_use]
    pub fn force(mut self, yes: bool) -> Self {
        self.force = yes;
        self
    }

    #[must_use]
    pub fn namespaces(mut self, ns: &'a [&'a str]) -> Self {
        self.namespaces = ns;
        self
    }

    /// Set the number of input FASTA files imported concurrently.
    /// `0` (the default) means auto via `std::thread::available_parallelism`;
    /// `1` means serial. No effect with a single input file.
    #[must_use]
    pub fn jobs(mut self, n: usize) -> Self {
        self.jobs = n;
        self
    }

    /// Register the imported collection in the **collection** alias index under
    /// `namespace:alias` (e.g. `.collection_alias("ucsc", "hg38")`).
    ///
    /// This is an explicit assertion by the caller — nothing in a FASTA file
    /// names the assembly. It is independent of [`Self::namespaces`], which only
    /// harvests per-sequence aliases from headers.
    ///
    /// Only valid when importing a single FASTA file; importing multiple files
    /// with this set returns an error, because one alias cannot name N
    /// collections. For the multi-file case, import without this option and call
    /// `add_collection_alias` per returned collection (results come back in
    /// input order).
    ///
    /// If the alias already resolves to a *different* collection digest, the
    /// import errors rather than silently remapping it, unless
    /// [`Self::force`] is set. Re-registering the same digest is a no-op.
    #[must_use]
    pub fn collection_alias(mut self, namespace: &'a str, alias: &'a str) -> Self {
        self.collection_alias = Some((namespace, alias));
        self
    }
}

/// Metadata for the entire store.
/// This is used to serialize metadata to `rgstore.json`, which can be loaded by the application.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct StoreMetadata {
    /// Version of the metadata format
    pub(crate) version: u32,
    /// Template for sequence file paths
    pub(crate) seqdata_path_template: String,
    /// Template for collection file paths
    pub(crate) collections_path_template: String,
    /// Path to the sequence metadata index file
    pub(crate) sequence_index: String,
    /// Path to the collection metadata index file (NEW)
    #[serde(default)]
    pub(crate) collection_index: Option<String>,
    /// Storage mode (Raw or Encoded)
    pub(crate) mode: StorageMode,
    /// Creation timestamp
    pub(crate) created_at: String,
    /// Whether ancillary digests are computed and stored
    #[serde(default = "default_true")]
    pub(crate) ancillary_digests: bool,
    /// Whether on-disk attribute index is maintained (Part 2)
    #[serde(default)]
    pub(crate) attribute_index: bool,
    /// Available sequence alias namespaces (for remote discovery)
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub(crate) sequence_alias_namespaces: Vec<String>,
    /// Available collection alias namespaces (for remote discovery)
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub(crate) collection_alias_namespaces: Vec<String>,
    /// Last-modified timestamp (RFC 3339). Updated on every write_index_files().
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) modified: Option<String>,
    /// SHA256 digest of collections.rgci
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) collections_digest: Option<String>,
    /// SHA256 digest of sequences.rgsi
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) sequences_digest: Option<String>,
    /// SHA256 digest of combined alias data
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) aliases_digest: Option<String>,
    /// SHA256 digest of combined FHR sidecar data
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) fhr_digest: Option<String>,
    /// Logical sequence-data size in bytes (sum of length x storage-mode encoding),
    /// computed once at index-write time. Lets a consumer read the store's size from
    /// the manifest without loading the sequence index. Excludes index/sidecar
    /// overhead. `None` for manifests written before this field existed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) logical_sequence_bytes: Option<u64>,
}

pub(crate) fn default_true() -> bool {
    true
}

/// Statistics for a RefgetStore.
///
/// These are a snapshot of the store's CURRENT RAM residency, not a record of
/// what a given import run did. For per-run ingest counts, use [`ImportReport`]
/// returned by `add_sequence_collections_from_fastas`.
#[derive(Debug, Clone)]
pub struct StoreStats {
    /// Total number of sequences (Stub + Full)
    pub n_sequences: usize,
    /// Number of sequences whose bytes are currently held in RAM (`Full`
    /// records). Structurally always 0 for a disk-backed store: importing
    /// downgrades every record to a `Stub` after handing its bytes to the
    /// writer pool (see `add_sequence_record_deferred_write`). Only an
    /// on-demand remote fetch ever makes a record `Full` again.
    pub n_sequences_in_memory: usize,
    /// Total number of collections (Stub + Full)
    pub n_collections: usize,
    /// Number of collections whose sequence list is currently held in RAM
    /// (`Full` records). Collections finalized by this process are `Full`;
    /// collections rehydrated from `collections.rgci` are `Stub`. This is a
    /// residency gauge, NOT a count of collections added by this run — it
    /// resets on process start and also counts collections merely touched by
    /// a read.
    pub n_collections_in_memory: usize,
    /// Storage mode (Raw or Encoded)
    pub storage_mode: String,
    /// Logical size in bytes of all sequence payloads (length x storage-mode
    /// encoding). Excludes index/sidecar/manifest overhead. For the exact full
    /// on-disk footprint use `ReadonlyRefgetStore::actual_disk_usage()`.
    pub logical_sequence_bytes: u64,
}

/// Result of importing one or more FASTA files.
///
/// Unlike [`StoreStats`] (a RAM-residency gauge), these are genuine per-run
/// counters describing what this import actually did.
///
/// `n_sequences_written + n_sequences_deduped` equals the total number of
/// sequence records seen across all *processed* files. Files short-circuited
/// as already-present collections (before their FASTA is ever opened)
/// contribute nothing to either counter.
#[derive(Debug, Clone)]
pub struct ImportReport {
    /// Per-file results in input order: `(metadata, was_new)`.
    pub collections: Vec<(SequenceCollectionMetadata, bool)>,
    /// Sequences whose bytes were dispatched to the writer pool this run
    /// (i.e. genuinely new content). For in-memory stores, sequences newly
    /// inserted into the in-memory map.
    pub n_sequences_written: usize,
    /// Sequences already present (by content digest), so no bytes were written.
    pub n_sequences_deduped: usize,
    /// Number of collections newly added by this run.
    pub n_collections_new: usize,
}

/// Format bytes into human-readable size (KB, MB, GB, etc.)
pub(crate) fn format_bytes(bytes: usize) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    let mut size = bytes as f64;
    let mut unit_idx = 0;

    while size >= 1024.0 && unit_idx < UNITS.len() - 1 {
        size /= 1024.0;
        unit_idx += 1;
    }

    if unit_idx == 0 {
        format!("{} {}", bytes, UNITS[0])
    } else {
        format!("{:.2} {}", size, UNITS[unit_idx])
    }
}

// =========================================================================
// Sidecar sync types
// =========================================================================

/// Conflict resolution strategy for sidecar pull operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SyncStrategy {
    /// Skip fetch if local file exists (default). Local wins.
    KeepOurs,
    /// Always fetch from remote, overwriting local. Remote wins.
    KeepTheirs,
    /// Report what would change without fetching. Returns diff info.
    Notify,
}

/// Result of a sidecar pull operation.
#[derive(Debug, Default)]
pub struct PullResult {
    /// Successfully fetched from remote
    pub pulled: usize,
    /// Skipped (local exists, KeepOurs)
    pub skipped: usize,
    /// Remote 404 (no sidecar exists)
    pub not_found: usize,
    /// For Notify: paths that differ between local and remote
    pub conflicts: Vec<String>,
}

/// Available alias namespaces advertised by a store's manifest.
#[derive(Debug)]
pub struct AvailableAliases<'a> {
    pub sequences: &'a [String],
    pub collections: &'a [String],
}

/// Iterator over BED file regions yielding substrings from a store.
// Constructed only by the filesystem BED-streaming path; fields go unread on
// the WASM (no-filesystem) build.
#[cfg_attr(not(feature = "filesystem"), allow(dead_code))]
pub struct SubstringsFromRegions<'a, K>
where
    K: AsRef<[u8]>,
{
    pub(crate) store: &'a ReadonlyRefgetStore,
    pub(crate) reader: BufReader<Box<dyn Read>>,
    pub(crate) collection_digest: K,
    pub(crate) previous_parsed_chr: String,
    pub(crate) current_seq_digest: String,
    pub(crate) line_num: usize,
}

/// Iterator over in-memory region vectors yielding substrings from a store.
///
/// Unlike `SubstringsFromRegions`, this has no filesystem dependency: it
/// iterates caller-supplied parallel `chroms`/`starts`/`ends` vectors.
// Constructed only by the filesystem-gated `export` module today; fields go
// unread on the WASM (no-filesystem) build.
#[cfg_attr(not(feature = "filesystem"), allow(dead_code))]
pub struct SubstringsFromRegionVectors<'a, K>
where
    K: AsRef<[u8]>,
{
    pub(crate) store: &'a ReadonlyRefgetStore,
    pub(crate) collection_digest: K,
    pub(crate) chroms: Vec<String>,
    pub(crate) starts: Vec<u32>,
    pub(crate) ends: Vec<u32>,
    pub(crate) index: usize,
    pub(crate) previous_parsed_chr: String,
    pub(crate) current_seq_digest: String,
}
