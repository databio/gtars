use super::consts::BBCLIENT_CACHE_ENV;
use dirs::home_dir;
use std::env;
use std::fs::create_dir_all;
use std::path::PathBuf;

use shellexpand;

pub fn get_abs_path(path: Option<PathBuf>, create_folder: Option<bool>) -> PathBuf {
    let raw_path = path.unwrap_or_else(get_default_cache_folder);

    // Convert to owned string early to avoid borrowing a temporary
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

pub fn get_default_cache_folder() -> PathBuf {
    if let Ok(val) = env::var(BBCLIENT_CACHE_ENV) {
        PathBuf::from(val)
    } else {
        let home = env::var("HOME")
            .or_else(|_| {
                home_dir()
                    .map(|p| p.to_string_lossy().into_owned())
                    .ok_or_else(|| std::env::VarError::NotPresent)
            })
            .unwrap_or_else(|_| "/tmp".to_string());

        let mut path = PathBuf::from(home);
        path.push(".bbcache/");
        path
    }
}
