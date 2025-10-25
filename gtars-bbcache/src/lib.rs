//! # gtars-bbcache
//!
//! A Rust implementation of bbcache: a caching system for BED files from [BEDbase.org](https://bedbase.org).
//!
//! ## Overview
//!
//! `gtars-bbcache` provides an efficient local caching layer for BED files and BED sets retrieved from
//! the BEDbase API. It handles downloading, storage, and retrieval of genomic region data.
//!
//! ## Features
//!
//! - **Smart Caching**: Automatically downloads and caches BED files from BEDbase API on first access
//! - **BEDset Support**: Manages collections of related BED files as coherent sets
//! - **Local File Import**: Add local BED files to the cache for unified access
//! - **Efficient Storage**: Organizes cached files in a hierarchical directory structure
//! - **SQLite Tracking**: Uses biocache for fast lookups and resource management
//! - **Configurable**: Customize cache location and API endpoints via environment variables
//!
//! ## Quick Start
//!
//! ```rust,no_run
//! use gtars_bbcache::client::BBClient;
//! use std::path::PathBuf;
//!
//! # fn main() -> anyhow::Result<()> {
//! // Create a client with default settings
//! let mut client = BBClient::builder().finish()?;
//!
//! // Load a BED file from BEDbase (downloads and caches if not present)
//! let region_set = client.load_bed("6b2e163a1d4319d99bd465c6c78a9741")?;
//!
//! // Add a local BED file to the cache
//! let bed_id = client.add_local_bed_to_cache(
//!     PathBuf::from("path/to/file.bed.gz"),
//!     None
//! )?;
//!
//! // Check if a file exists in cache
//! let cached_path = client.seek(&bed_id)?;
//! # Ok(())
//! # }
//! ```
//!
//! ## Configuration
//!
//! The cache behavior can be configured through environment variables:
//!
//! - `BBCLIENT_CACHE`: Custom cache directory (default: `~/.bbcache/`)
//! - `BEDBASE_API`: Custom BEDbase API endpoint (default: `https://api.bedbase.org`)
//!
//! Or programmatically via the builder:
//!
//! ```rust,no_run
//! use gtars_bbcache::client::BBClient;
//! use std::path::PathBuf;
//!
//! # fn main() -> anyhow::Result<()> {
//! let client = BBClient::builder()
//!     .with_cache_folder(PathBuf::from("/custom/cache/path"))
//!     .with_bedbase_api("https://api.bedbase.org".to_string())
//!     .finish()?;
//! # Ok(())
//! # }
//! ```
//!
//! ## Cache Structure
//!
//! Cached files are organized hierarchically:
//!
//! ```text
//! <cache_folder>/
//! ├── bedfiles/
//! │   └── <first_char>/
//! │       └── <second_char>/
//! │           └── <identifier>.bed.gz
//! └── bedsets/
//!     └── <first_char>/
//!         └── <second_char>/
//!             └── <identifier>.txt
//! ```

pub mod client;
pub mod consts;
pub mod utils;

#[cfg(test)]
mod tests {
    use super::client::BBClient;
    use rstest::{fixture, rstest};
    use std::fs::read_dir;
    use std::path::PathBuf;

    #[fixture]
    fn path_to_bed_gz_from_bb() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/6b2e163a1d4319d99bd465c6c78a9741.bed.gz")
    }

    #[fixture]
    fn bbid() -> PathBuf {
        "6b2e163a1d4319d99bd465c6c78a9741".into()
    }

    #[fixture]
    fn path_to_bedset() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/bedset")
    }

    #[rstest]
    fn test_bbcache_local(
        path_to_bed_gz_from_bb: PathBuf,
        bbid: PathBuf,
        path_to_bedset: PathBuf,
    ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        fn cleaned_subfolders(subfolder: PathBuf) {
            let subdirs: Vec<_> = read_dir(&subfolder)
                .unwrap_or_else(|e| {
                    panic!("Failed to read directory {}: {}", subfolder.display(), e)
                })
                .filter_map(Result::ok)
                .filter(|entry| entry.path().is_dir())
                .collect();

            // Assert no subdirectories exist
            assert!(
                subdirs.is_empty(),
                "Subfolders found in {}: {:?}",
                subfolder.display(),
                subdirs.iter().map(|e| e.path()).collect::<Vec<_>>()
            );
        }
        let tempdir = tempfile::tempdir()?;
        let cache_folder = PathBuf::from(tempdir.path());

        let mut bbc = BBClient::builder()
            .with_cache_folder(cache_folder.clone())
            .finish()?;

        let bed_id = bbc
            .add_local_bed_to_cache(path_to_bed_gz_from_bb, Some(false))
            .unwrap();
        assert_eq!(&bed_id, &bbid.to_string_lossy());

        let bedset_id = bbc.add_local_folder_as_bedset(path_to_bedset).unwrap();
        assert!(bbc.seek(&bedset_id).is_ok());

        bbc.remove(&bedset_id)
            .expect("Failed to remove bedset file and its bed files");
        let bedset_subfolder = cache_folder.join("bedsets");
        cleaned_subfolders(bedset_subfolder);

        bbc.remove(&bbid.to_string_lossy())
            .expect("Failed to remove cached bed file");
        let bedfile_subfolder = cache_folder.join("bedfiles");
        cleaned_subfolders(bedfile_subfolder);
        Ok(())
    }
}
