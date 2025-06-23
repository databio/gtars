use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::{create_dir_all, read_dir, remove_dir, remove_file};
use std::io::{Error, ErrorKind};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

use log::info;

use super::consts::{DEFAULT_BEDBASE_API, DEFAULT_BEDFILE_EXT, DEFAULT_BEDFILE_SUBFOLDER};
use super::utils::{bb_url_for_regionset, get_abs_path};
use crate::common::models::region_set::RegionSet;

pub struct BBClient {
    pub cache_folder: PathBuf,
    pub bedbase_api: String,
}

impl BBClient {
    pub fn new(cache_folder: Option<PathBuf>, bedbase_api: Option<String>) -> Result<Self> {
        let cache_folder = get_abs_path(cache_folder, Some(true));
        let bedbase_api = bedbase_api.unwrap_or_else(|| DEFAULT_BEDBASE_API.to_string());

        Ok(BBClient {
            cache_folder,
            bedbase_api,
        })
    }

    pub fn load_bed(&mut self, bed_id: &str) -> Result<RegionSet> {
        let bedfile_path = self.bedfile_path(bed_id, Some(false));

        if bedfile_path.exists() {
            info!("Loading cached BED file from {:?}", bedfile_path.display());
            return RegionSet::try_from(bedfile_path.as_path());
        }

        let regionset = self.download_bed_file_from_bb(bed_id)?;
        info!("Downloaded BED file from BEDbase: {}", bed_id);
        // write(&bedfile_path, bed_data)?;
        regionset.to_bed_gz(bedfile_path.clone().as_path())?;

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

        if !force.unwrap_or(false) && cache_path.exists() {
            info!("{} already exists in cache", cache_path.display());
            return Ok(regionset);
        }

        regionset.to_bed_gz(cache_path.as_path())?;
        Ok(regionset)
    }

    fn download_bed_file_from_bb(&self, bedfile: &str) -> Result<RegionSet> {
        let bed_url = bb_url_for_regionset(bedfile);

        // let regionset = RegionSet::try_from(bed_url.clone())
        //     .expect(&format!("Failed to create RegionSet from URL {}", bed_url));

        let regionset = RegionSet::try_from(bed_url.clone())
            .with_context(|| format!("Failed to create RegionSet from URL {}", bed_url))?;

        Ok(regionset)
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

        // Ok((folder_path.join(filename)))
        folder_path.join(filename)
    }

    fn create_cache_folder(&self, subfolder_path: Option<&Path>) {
        let path = match subfolder_path {
            Some(p) => p.to_path_buf(),
            None => self.cache_folder.clone(),
        };

        if !path.exists() {
            create_dir_all(&path).expect("Failed to create cache folder");
        }
    }

    pub fn seek(&self, identifier: &str) -> Result<PathBuf> {
        let cache_path = self.bedfile_path(&identifier, Some(false));
        if cache_path.exists() {
            Ok(cache_path)
        } else {
            Err(anyhow::anyhow!(
                "{} does not exist in cache.",
                cache_path.display()
            ))
        }
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

    pub fn list_beds(&self) -> Result<HashMap<String, PathBuf>> {
        let mut bedfile_map = HashMap::new();

        let bedfile_dir = self.cache_folder.join(DEFAULT_BEDFILE_SUBFOLDER);
        if !bedfile_dir.exists() {
            return Ok(bedfile_map); // return empty map if folder doesn't exist
        }

        for entry in WalkDir::new(&bedfile_dir)
            .into_iter()
            .filter_map(|e| e.ok())
        {
            let path = entry.path();

            if path.is_file() {
                if let Some(_ext) = path.extension().and_then(|s| s.to_str()) {
                    if path
                        .file_name()
                        .and_then(|f| f.to_str())
                        .map(|name| name.ends_with(DEFAULT_BEDFILE_EXT))
                        .unwrap_or(false)
                    {
                        if let Some(file_stem) = path.file_name().and_then(|s| s.to_str()) {
                            let id = file_stem
                                .strip_suffix(DEFAULT_BEDFILE_EXT)
                                .unwrap_or(file_stem)
                                .to_string();
                            bedfile_map.insert(id, path.to_path_buf());
                        }
                    }
                }
            }
        }

        Ok(bedfile_map)
    }
}
