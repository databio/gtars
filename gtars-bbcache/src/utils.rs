use super::consts::{BBCLIENT_CACHE_ENV, BEDBASE_API_ENV};
use biocrs::models::Resource;
use dirs::home_dir;
use std::env;
use std::path::PathBuf;
use tabled::{Table, Tabled};

#[derive(Tabled)]
pub struct ResourcePrint {
    id: String,
    path: String,
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
pub fn get_default_bedbase_api() -> String {
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
