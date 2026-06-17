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

mod readonly;
mod core;
mod alias;
mod fhr_metadata;
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
pub use self::readonly::ReadonlyRefgetStore;
pub use self::core::RefgetStore;
pub use self::alias::{AliasKind, AliasManager};
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
#[derive(Clone, Copy)]
pub struct FastaImportOptions<'a> {
    pub(crate) force: bool,
    pub(crate) namespaces: &'a [&'a str],
    /// Number of input FASTA files imported concurrently. Each file gets its own
    /// decoder so gzip decompression is parallelized across files. `0` = auto
    /// (`std::thread::available_parallelism`), `1` = serial. No effect with a
    /// single input file.
    pub(crate) jobs: usize,
}

impl<'a> Default for FastaImportOptions<'a> {
    fn default() -> Self {
        Self {
            force: false,
            namespaces: &[],
            jobs: 0,
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
}

pub(crate) fn default_true() -> bool {
    true
}

/// Statistics for a RefgetStore
#[derive(Debug, Clone)]
pub struct StoreStats {
    /// Total number of sequences (Stub + Full)
    pub n_sequences: usize,
    /// Number of sequences with data loaded (Full)
    pub n_sequences_loaded: usize,
    /// Total number of collections (Stub + Full)
    pub n_collections: usize,
    /// Number of collections with sequences loaded (Full)
    pub n_collections_loaded: usize,
    /// Storage mode (Raw or Encoded)
    pub storage_mode: String,
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
