use std::fs::create_dir_all;
use std::path::{Path, PathBuf};

use anyhow::Result;
use shellexpand;

use super::consts::DEFAULT_CACHE_FOLDER;

pub fn get_abs_path(path: Option<PathBuf>, create_folder: Option<bool>) -> PathBuf {
    let raw_path = path.unwrap_or_else(|| PathBuf::from(DEFAULT_CACHE_FOLDER));

    let expanded_path = shellexpand::env(&raw_path.to_string_lossy())
        .unwrap_or_else(|_| raw_path.to_string_lossy().into());

    let abs_path = PathBuf::from(expanded_path.into_owned());
    if create_folder.unwrap_or(true) {
        create_dir_all(&abs_path).expect("Failed to create directory");
    }

    abs_path
}


