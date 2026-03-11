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

mod readonly;
mod core;
mod alias;
mod fhr_metadata;
mod import;
mod persistence;
mod export;

#[cfg(test)]
mod tests;

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
use std::io::{BufRead, BufReader, Read};


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
#[derive(Clone, Copy)]
pub struct FastaImportOptions<'a> {
    pub(crate) force: bool,
    pub(crate) namespaces: &'a [&'a str],
}

impl<'a> Default for FastaImportOptions<'a> {
    fn default() -> Self {
        Self {
            force: false,
            namespaces: &[],
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
