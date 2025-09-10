use super::consts::{BBCLIENT_CACHE_ENV, BEDBASE_API_ENV};
use biocrs::models::Resource;
use dirs::home_dir;
use std::env;
use std::fs::create_dir_all;
use std::path::PathBuf;
use tabled::{Table, Tabled};

use shellexpand;

#[derive(Tabled)]
pub struct ResourcePrint {
    id: String,
    path: String,
}

/// Get absolute path to the folder and create it if it doesn't exist
/// # Arguments
/// - path: path to the folder
/// - create_folder: create folder if it doesn't exist
///
/// # Returns
/// - absolute path to the folder
pub fn get_abs_path(path: Option<PathBuf>, create_folder: Option<bool>) -> PathBuf {
    let raw_path = path.unwrap_or_else(get_default_cache_folder);

    let raw_str = raw_path.to_string_lossy().into_owned();

    let expanded_str = shellexpand::env(&raw_str)
        .unwrap_or_else(|_| raw_str.clone().into()) // Use clone to satisfy the closure
        .into_owned(); // Result of `shellexpand::env` is a Cow

    let abs_path = PathBuf::from(expanded_str);

    if create_folder.unwrap_or(true) {
        create_dir_all(&abs_path).expect("Failed to create directory");
    }

    abs_path
}

/// Get default cache folder from environment variable, if not available then create it in home folder
///
/// # Returns
/// - path to cache folder
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

/// Get default BEDbase api from environment variable
///
/// # Returns
/// - BEDbase api for url
pub fn get_bedbase_api() -> String {
    env::var(BEDBASE_API_ENV).unwrap_or_else(|_| "https://api.bedbase.org".to_string())
}

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
