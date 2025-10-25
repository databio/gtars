//! Utility functions for bbcache configuration and display.
//!
//! This module provides helper functions for:
//! - Determining default cache locations and API endpoints
//! - Formatting and displaying cached resources in tabular form

use super::consts::{BBCLIENT_CACHE_ENV, BEDBASE_API_ENV};
use biocrs::models::Resource;
use dirs::home_dir;
use std::env;
use std::path::PathBuf;
use tabled::{Table, Tabled};

/// Printable representation of a cached resource for display in tables.
///
/// This struct is used internally by [`print_resources`] to format resource
/// information in a human-readable tabular format.
#[derive(Tabled)]
pub struct ResourcePrint {
    /// Resource identifier (e.g., BED file or BED set ID)
    id: String,
    /// Local filesystem path to the cached resource
    path: String,
}

/// Returns the default cache folder path.
///
/// The cache folder is determined in the following priority order:
/// 1. `BBCLIENT_CACHE` environment variable if set
/// 2. `$HOME/.bbcache/` if home directory is available
/// 3. `/tmp/.bbcache/` as a fallback
///
/// # Returns
///
/// A [`PathBuf`] pointing to the cache folder location.
///
/// # Examples
///
/// ```rust
/// use gtars_bbcache::utils::get_default_cache_folder;
///
/// let cache_path = get_default_cache_folder();
/// println!("Cache will be stored at: {:?}", cache_path);
/// ```
///
/// # Environment Variables
///
/// - `BBCLIENT_CACHE`: Custom cache directory path (highest priority)
/// - `HOME`: User's home directory (used for default location)
pub fn get_default_cache_folder() -> PathBuf {
    if let Ok(val) = env::var(BBCLIENT_CACHE_ENV) {
        PathBuf::from(val)
    } else {
        let home = env::var("HOME")
            .or_else(|_| {
                home_dir()
                    .map(|p| p.to_string_lossy().into_owned())
                    .ok_or(std::env::VarError::NotPresent)
            })
            .unwrap_or_else(|_| "/tmp".to_string());

        let mut path = PathBuf::from(home);
        path.push(".bbcache/");
        path
    }
}

/// Returns the default BEDbase API endpoint URL.
///
/// The API endpoint is determined in the following priority order:
/// 1. `BEDBASE_API` environment variable if set
/// 2. `https://api.bedbase.org` as the default
///
/// # Returns
///
/// A [`String`] containing the BEDbase API endpoint URL.
///
/// # Examples
///
/// ```rust
/// use gtars_bbcache::utils::get_default_bedbase_api;
///
/// let api_url = get_default_bedbase_api();
/// println!("Using BEDbase API at: {}", api_url);
/// ```
///
/// # Environment Variables
///
/// - `BEDBASE_API`: Custom BEDbase API endpoint
pub fn get_default_bedbase_api() -> String {
    env::var(BEDBASE_API_ENV).unwrap_or_else(|_| "https://api.bedbase.org".to_string())
}

/// Prints a list of resources in a formatted table.
///
/// This function takes a vector of [`Resource`] objects and displays them
/// in a tabular format showing identifiers and paths. Useful for displaying
/// the results of [`BBClient::list_beds`] or [`BBClient::list_bedsets`].
///
/// # Arguments
///
/// * `resources` - Vector of resources to display
///
/// # Examples
///
/// ```rust,no_run
/// use gtars_bbcache::client::BBClient;
/// use gtars_bbcache::utils::print_resources;
///
/// # fn main() -> anyhow::Result<()> {
/// let mut client = BBClient::builder().finish()?;
/// let beds = client.list_beds()?;
/// print_resources(beds);
/// # Ok(())
/// # }
/// ```
///
/// # Output Format
///
/// ```text
/// ┌────────────────────────────────┬─────────────────────────┐
/// │ id                             │ path                    │
/// ├────────────────────────────────┼─────────────────────────┤
/// │ 6b2e163a1d4319d99bd465c6c78... │ /path/to/cache/6/b/...  │
/// └────────────────────────────────┴─────────────────────────┘
/// ```
pub fn print_resources(resources: Vec<Resource>) {
    let mut resource_print: Vec<ResourcePrint> = Vec::new();

    for resource in resources {
        resource_print.push(ResourcePrint {
            id: resource.rname,
            path: resource.rpath,
        })
    }

    let table = Table::new(resource_print);

    println!("{}", table);
}
