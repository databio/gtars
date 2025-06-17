use std::fs::{self, File, create_dir_all};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::io::{Error, ErrorKind};
use std::collections::HashMap;
use anyhow::Result; 

use log::info;
use once_cell::sync::Lazy;
use serde::Deserialize;
use reqwest::blocking::get;
use std::sync::Mutex;

use super::utils::get_abs_path;
use super::consts::{
    DEFAULT_CACHE_FOLDER, DEFAULT_BEDBASE_API, BEDBASE_URL_PATTERN,
    DEFAULT_BEDFILE_SUBFOLDER, DEFAULT_BEDFILE_EXT
};
use crate::common::models::region_set::RegionSet;



static MODULE_NAME: &str = "bbcache";


pub struct BBClient {
    pub cache_folder: PathBuf,
    pub bedbase_api: String,
    pub bedfile_cache: HashMap<String, PathBuf>,
}

impl BBClient {
    pub fn new(
        cache_folder: Option<PathBuf>,
        bedbase_api: Option<String>,
    ) -> Result<Self> {
        let cache_folder = get_abs_path(cache_folder, Some(true));
        let bedbase_api = bedbase_api.unwrap_or_else(|| DEFAULT_BEDBASE_API.to_string());
        let bedfile_cache = HashMap::new();

        Ok(BBClient {
            cache_folder,
            bedbase_api,
            bedfile_cache,
        })
    }

    pub fn load_bed(&mut self, bed_id: &str) -> Result<RegionSet> {
        let bedfile_path = self.bedfile_path(bed_id, Some(false));

        if bedfile_path.exists() {
            info!("Loading cached BED file from {:?}", bedfile_path);
            return RegionSet::try_from(bedfile_path.as_path());
        }

        let bed_data = self.download_bed_file_from_bb(bed_id)?;
        fs::write(&bedfile_path, bed_data)?;

        let region_set = RegionSet::try_from(bedfile_path.as_path())?;
        self.bedfile_cache.insert(bed_id.to_string(), bedfile_path.clone());

        Ok(region_set)
    }

    pub fn add_local_bed_to_cache(&mut self, bedfile: PathBuf, force: Option<bool>) -> Result<RegionSet> {
        let region_set = RegionSet::try_from(bedfile.as_path())?;
        self.add_regionset_to_cache(region_set, force)
    }

    pub fn add_regionset_to_cache(
        &mut self,
        regionset: RegionSet,
        force: Option<bool>,
    ) -> Result<RegionSet> {
        let bedfile_id = regionset.identifier();
        let cache_path = self.cache_folder.join(format!("{}{}", bedfile_id, DEFAULT_BEDFILE_EXT));

        if !force.unwrap_or(false) && self.bedfile_cache.contains_key(&bedfile_id) {
            info!("RegionSet {} already cached at {:?}", bedfile_id, cache_path);
            return Ok(regionset);
        }

        let mut file = File::create(&cache_path)?;
        file.write_all(regionset.to_bed_string().as_bytes())?;
        self.bedfile_cache.insert(bedfile_id.clone(), cache_path);

        info!("Cached RegionSet {} at {:?}", bedfile_id, cache_path);
        Ok(regionset)
    }

    fn download_bed_file_from_bb(&self, bedfile: &str) -> Result<Vec<u8>> {
        let bed_url = BEDBASE_URL_PATTERN
            .replace("{bedbase_api}", &self.bedbase_api)
            .replace("{bed_id}", bedfile);

        let response = get(&bed_url)?;
        if !response.status().is_success() {
            anyhow::bail!("Request to {} failed with status {}", bed_url, response.status());
        }

        Ok(response.bytes()?.to_vec())
    }

    fn bedfile_path(&self, bedfile_id: &str, create: Option<bool>) -> PathBuf {
        let subfolder_name = DEFAULT_BEDFILE_SUBFOLDER;
        let file_extension = DEFAULT_BEDFILE_EXT;
        self.cache_path(bedfile_id, subfolder_name, file_extension, create)
    }

    fn cache_path(
        &self,
        identifier: &str,
        subfolder_name: &str,
        file_extension: &str,
        create: Option<bool>,
    ) -> PathBuf {
        let filename = format!("{}{}", identifier, file_extension);
        let folder_path = self
            .cache_folder
            .join(subfolder_name)
            .join(&identifier[0..1])
            .join(&identifier[1..2]);

        if create.unwrap_or(true) {
            self.create_cache_folder(Some(&folder_path));
        }

        folder_path.join(filename)
    }

    pub fn create_cache_folder(&self, subfolder_path: Option<&Path>) {
        let path = match subfolder_path {
            Some(p) => p.to_path_buf(),
            None => self.cache_folder.clone(),
        };

        if !path.exists() {
            create_dir_all(&path).expect("Failed to create cache folder");
        }
    }

    pub fn seek(&self, identifier: &str) -> Result<PathBuf> {
        let file_path = self.bedfile_path(identifier, Some(false));
        if file_path.exists() {
            Ok(file_path)
        } else {
            Err(Error::new(ErrorKind::NotFound, format!("{} does not exist in cache.", identifier)))
        }
    }
}