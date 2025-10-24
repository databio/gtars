//! Constants for bbcache configuration and file organization.
//!
//! This module defines environment variable names, file extensions, directory names,
//! and command strings used throughout the bbcache system.

// Environment variable names

/// Environment variable name for setting the cache directory location.
///
/// When set, this overrides the default cache location (`~/.bbcache/`).
///
/// # Example
///
/// ```bash
/// export BBCLIENT_CACHE=/custom/cache/path
/// ```
pub const BBCLIENT_CACHE_ENV: &str = "BBCLIENT_CACHE";

/// Environment variable name for setting the BEDbase API endpoint.
///
/// When set, this overrides the default API endpoint (`https://api.bedbase.org`).
///
/// # Example
///
/// ```bash
/// export BEDBASE_API=https://custom.bedbase.org
/// ```
pub const BEDBASE_API_ENV: &str = "BEDBASE_API";

// Command-line interface command names

/// Main bbcache command name.
pub const BBCACHE_CMD: &str = "bbcache";

/// Subcommand for caching a BED file.
pub const BBCACHE_CACHEBED: &str = "cache-bed";

/// Subcommand for caching a BED set.
pub const BBCACHE_CACHEBEDSET: &str = "cache-bedset";

/// Subcommand for seeking (locating) a cached resource.
pub const BBCACHE_SEEK: &str = "seek";

/// Subcommand for inspecting cached BED files.
pub const BBCACHE_INSPECTBED: &str = "inspect-bedfiles";

/// Subcommand for inspecting cached BED sets.
pub const BBCACHE_INSPECTBEDSET: &str = "inspect-bedsets";

/// Subcommand for removing cached resources.
pub const BBCACHE_REMOVE: &str = "rm";

// Directory structure constants

/// Default subdirectory name for storing individual BED files.
///
/// BED files are stored in `<cache_folder>/bedfiles/`.
pub const DEFAULT_BEDFILE_SUBFOLDER: &str = "bedfiles";

/// Default subdirectory name for storing BED set metadata files.
///
/// BED set files are stored in `<cache_folder>/bedsets/`.
pub const DEFAULT_BEDSET_SUBFOLDER: &str = "bedsets";

// File extension constants

/// Default file extension for cached BED files.
///
/// BED files are stored compressed as gzipped files with the `.bed.gz` extension.
pub const DEFAULT_BEDFILE_EXT: &str = ".bed.gz";

/// Default file extension for BED set metadata files.
///
/// BED sets are stored as plain text files with the `.txt` extension,
/// containing a newline-separated list of BED file identifiers.
pub const DEFAULT_BEDSET_EXT: &str = ".txt";
