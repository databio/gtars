use std::path::PathBuf;

use anyhow::Result;
use reqwest::blocking::Client;

mod consts {
    pub const CACHE_FOLDER: &str = ".bbcache";
    pub const BEDFILE_CACHE_SUBFOLDER: &str = "bedfiles";
    // pub const BEDSET_CACHE_SUBFOLDER: &str = "bedsets";
}

pub struct Bbclient {
    bedbase_api: String,
}

impl Bbclient {
    pub fn new(bedbase_api: &str) -> Self {
        Bbclient {
            bedbase_api: bedbase_api.to_string(),
        }
    }

    pub fn get_bedbase_api(&self) -> &str {
        &self.bedbase_api
    }

    fn get_bed_file_download_link(&self, id: &str) -> Result<String> {
        Ok(format!(
            "{}/v1/objects/bed.{}.bed_file/access/http/bytes",
            self.bedbase_api, id
        ))
    }

    fn get_cache_path_for_bed_file(&self, id: &str) -> Result<PathBuf> {
        let cache_folder = dirs::home_dir().unwrap().join(consts::CACHE_FOLDER);
        let bedfile_cache_folder = cache_folder.join(consts::BEDFILE_CACHE_SUBFOLDER);
        let bedfile_cache_path = bedfile_cache_folder.join(format!("{}.bed", id));

        Ok(bedfile_cache_path)
    }
}

impl Default for Bbclient {
    fn default() -> Self {
        Bbclient {
            bedbase_api: "https://bedbase.org/api".to_string(),
        }
    }
}

pub trait BedbaseFileCache {
    ///
    /// Get the path to a bed file from the cache or download it if it does not exist.
    ///
    /// # Arguments
    /// - `id` - the id of the file to download
    fn get_bed_file(&self, id: &str) -> Result<PathBuf>;
}

impl BedbaseFileCache for Bbclient {
    fn get_bed_file(&self, id: &str) -> Result<PathBuf> {
        let cache_path = self.get_cache_path_for_bed_file(id)?;

        // path exists?
        if cache_path.exists() {
            // read file
            Ok(cache_path)
        } else {
            // download file
            let download_link = self.get_bed_file_download_link(id)?;
            let client = Client::new();
            let response = client.get(download_link).send()?;
            let content = response.bytes()?;

            // write file
            std::fs::create_dir_all(cache_path.parent().unwrap())?;
            std::fs::write(&cache_path, content)?;

            Ok(cache_path)
        }
    }
}
