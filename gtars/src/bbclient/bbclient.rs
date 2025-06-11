use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use log::info;
use once_cell::sync::Lazy;
use serde::Deserialize;
use std::collections::HashMap;
use std::sync::Mutex;

use crate::utils::{get_abs_path, BedCacheManager};

static MODULE_NAME: &str = "bbclient";

#[derive(Debug, Deserialize)]
struct BedSetEntry {
    id: String,
}

pub struct BBClient {
    pub manager: BedCacheManager,
    pub bedbase_api: String,
}

impl BBClient {
    pub fn new(cache_folder: &str, bedbase_api: &str) -> Self {
        let abs_cache = get_abs_path(cache_folder, true);
        let manager = BedCacheManager::new(&abs_cache);

        BBClient {
            manager,
            bedbase_api: bedbase_api.to_string(),
        }
    }

    pub fn load_bedset(&self, bedset_id: &str) -> Vec<String> {
        let file_path = self._bedset_path(bedset_id);

        if file_path.exists() {
            info!("BED set {} already exists in cache.", bedset_id);
            let content = fs::read_to_string(&file_path).unwrap();
            content.lines().map(|s| s.to_string()).collect()
        } else {
            let extracted_data = self._download_bedset_data(bedset_id);
            let mut file = File::create(&file_path).unwrap();
            for value in &extracted_data {
                writeln!(file, "{}", value).unwrap();
            }
            info!("BED set {} downloaded and cached successfully.", bedset_id);
            extracted_data
        }
    }

    fn _download_bedset_data(&self, bedset_id: &str) -> Vec<String> {
        let url = format!("{}/bedset/{}", self.bedbase_api, bedset_id);
        let response = ureq::get(&url).call().unwrap();
        let json: HashMap<String, Vec<BedSetEntry>> = response.into_json().unwrap();
        json["results"].iter().map(|entry| entry.id.clone()).collect()
    }

    fn _bedset_path(&self, bedset_id: &str) -> PathBuf {
        self._cache_path(bedset_id, "bedsets", ".txt", true)
    }

    fn _cache_path(
        &self,
        identifier: &str,
        subfolder_name: &str,
        file_extension: &str,
        create: bool,
    ) -> PathBuf {
        let filename = format!("{}{}", identifier, file_extension);
        let folder_name = Path::new(&self.manager.cache_folder)
            .join(subfolder_name)
            .join(&identifier[0..1])
            .join(&identifier[1..2]);

        if create {
            self.manager.create_cache_folder(Some(folder_name.to_str().unwrap()));
        }
        folder_name.join(filename)
    }

    pub fn remove_bedfile_from_cache(&self, bedfile_id: &str) {
        let file_path = self._cache_path(bedfile_id, "bedfiles", ".bed.gz", false);
        if file_path.exists() {
            fs::remove_file(&file_path).unwrap();

            let sub_folder_2 = file_path.parent().unwrap();
            let sub_folder_1 = sub_folder_2.parent().unwrap();

            if fs::read_dir(sub_folder_2).unwrap().next().is_none() {
                fs::remove_dir(sub_folder_2).unwrap();
                if fs::read_dir(sub_folder_1).unwrap().next().is_none() {
                    fs::remove_dir(sub_folder_1).unwrap();
                }
            }
        }
    }
}
