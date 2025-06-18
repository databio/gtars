use anyhow::Result;
use std::collections::HashMap;
use std::fs::{create_dir_all, read_dir, remove_dir, remove_file, write};
use std::io::{Error, ErrorKind};
use std::path::{Path, PathBuf};

use log::info;
use reqwest::blocking::get;

use super::consts::{
    BEDBASE_URL_PATTERN, DEFAULT_BEDBASE_API, DEFAULT_BEDFILE_EXT, DEFAULT_BEDFILE_SUBFOLDER,
};
use super::utils::{bb_url_for_regionset, get_abs_path};
use crate::common::models::region_set::RegionSet;

static MODULE_NAME: &str = "bbcache";

pub struct BBClient {
    pub cache_folder: PathBuf,
    pub bedbase_api: String,
    pub bedfile_cache: HashMap<String, PathBuf>,
}

impl BBClient {
    pub fn new(cache_folder: Option<PathBuf>, bedbase_api: Option<String>) -> Result<Self> {
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
        let bedfile_path = self.bedfile_path(bed_id, Some(true));

        if bedfile_path.exists() {
            info!("Loading cached BED file from {:?}", bedfile_path);
            return RegionSet::try_from(bedfile_path.as_path());
        }

        let regionset = self.download_bed_file_from_bb(bed_id)?;
        // write(&bedfile_path, bed_data)?;
        regionset.to_bed_gz(bedfile_path.clone().as_path())?;

        // let region_set = RegionSet::try_from(bedfile_path.as_path())?;
        self.bedfile_cache
            .insert(bed_id.to_string(), bedfile_path.clone());

        Ok(regionset)
    }

    pub fn add_local_bed_to_cache(
        &mut self,
        bedfile: PathBuf,
        force: Option<bool>,
    ) -> Result<RegionSet> {
        let regionset = RegionSet::try_from(bedfile.as_path())?;
        self.add_regionset_to_cache(regionset, force)
    }

    pub fn add_regionset_to_cache(
        &mut self,
        regionset: RegionSet,
        force: Option<bool>,
    ) -> Result<RegionSet> {
        let bedfile_id = regionset.identifier();
        let cache_path = self.bedfile_path(&bedfile_id, Some(true));

        if !force.unwrap_or(false) && self.bedfile_cache.contains_key(&bedfile_id) {
            info!("{} already exists in cache", cache_path.clone().display());
            return Ok(regionset);
        }

        regionset.to_bed_gz(cache_path.as_path())?;
        self.bedfile_cache
            .insert(bedfile_id.clone(), cache_path.clone());

        info!(
            "Cached RegionSet {} at {:?}",
            bedfile_id,
            cache_path.display()
        );
        Ok(regionset)
    }

    fn download_bed_file_from_bb(&self, bedfile: &str) -> Result<RegionSet> {
        // let bed_url = BEDBASE_URL_PATTERN
        //     .replace("{bedbase_api}", &self.bedbase_api)
        //     .replace("{bed_id}", bedfile);
        let bed_url = bb_url_for_regionset(bedfile);

        let regionset = RegionSet::try_from(bed_url.clone())
            .expect(&format!("Failed to create RegionSet from URL {}", bed_url));

        Ok(regionset)

        // let response = get(&bed_url)?;
        // if !response.status().is_success() {
        //     anyhow::bail!(
        //         "Request to {} failed with status {}",
        //         bed_url,
        //         response.status()
        //     );
        // }

        // Ok(response.bytes()?.to_vec())
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

    pub fn seek(&self, identifier: &str) -> Result<PathBuf, std::io::Error> {
        self.bedfile_cache.get(identifier).cloned().ok_or_else(|| {
            Error::new(
                ErrorKind::NotFound,
                format!("{} does not exist in cache.", identifier),
            )
        })
    }

    pub fn remove(&mut self, identifier: &str) -> Result<()> {
        let file_path = self.bedfile_path(identifier, Some(false));
        if file_path.exists() {
            // remove file and check if subfolders is cleaned
            let sub_folder_2 = file_path.parent().map(PathBuf::from);
            let sub_folder_1 = sub_folder_2
                .as_ref()
                .and_then(|p| p.parent().map(PathBuf::from));

            remove_file(&file_path)?;

            // Attempt to remove empty subfolders
            if let Some(sub2) = sub_folder_2 {
                if read_dir(&sub2)?.next().is_none() {
                    remove_dir(&sub2)?;
                    if let Some(sub1) = sub_folder_1 {
                        if read_dir(&sub1)?.next().is_none() {
                            remove_dir(&sub1)?;
                        }
                    }
                }
            }

            self.bedfile_cache.remove(identifier);
            info!("{} is removed.", file_path.display());
            Ok(())
        } else {
            Err(Error::new(
                ErrorKind::NotFound,
                format!("{} does not exist in cache.", file_path.display()),
            )
            .into())
        }
    }
}
