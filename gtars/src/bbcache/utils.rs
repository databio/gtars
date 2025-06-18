use std::fs::create_dir_all;
use std::path::PathBuf;

use shellexpand;

use super::consts::DEFAULT_CACHE_FOLDER;

pub fn get_abs_path(path: Option<PathBuf>, create_folder: Option<bool>) -> PathBuf {
    let raw_path = path.unwrap_or_else(|| PathBuf::from(DEFAULT_CACHE_FOLDER.clone()));

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

pub fn bb_url_for_regionset(bbid: &str) -> String {
    format!(
        "https://api.bedbase.org/v1/files/files/{}/{}/{}.bed.gz",
        &bbid[0..1],
        &bbid[1..2],
        bbid
    )
}
