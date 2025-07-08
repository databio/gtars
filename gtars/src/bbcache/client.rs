use anyhow::{anyhow, Context, Ok, Result};
use biocrs::biocache::BioCache;
use biocrs::models::{NewResource, Resource};

use reqwest::blocking::get;
use std::fs::{create_dir_all, read_dir, remove_dir, remove_file, File};
use std::io::{BufRead, BufReader, Error, ErrorKind, Write};
use std::path::{Path, PathBuf};

use super::consts::{
    DEFAULT_BEDFILE_EXT, DEFAULT_BEDFILE_SUBFOLDER, DEFAULT_BEDSET_EXT, DEFAULT_BEDSET_SUBFOLDER,
};
use super::utils::{get_abs_path, get_bedbase_api};
use crate::common::models::bed_set::BedSet;
use crate::common::models::region_set::RegionSet;

pub struct BBClient {
    pub cache_folder: PathBuf,
    pub bedbase_api: String,
    bedfile_cache: BioCache,
    bedset_cache: BioCache,
}

impl BBClient {
    pub fn new(cache_folder: Option<PathBuf>, bedbase_api: Option<String>) -> Result<Self> {
        let cache_folder = get_abs_path(cache_folder, Some(true));
        let bedbase_api = bedbase_api.unwrap_or_else(|| get_bedbase_api());

        let bedfile_subfolder = &cache_folder.join(DEFAULT_BEDFILE_SUBFOLDER);
        create_dir_all(bedfile_subfolder)?;
        let bedfile_cache = BioCache::new(bedfile_subfolder);

        let bedset_subfolder = &cache_folder.join(DEFAULT_BEDSET_SUBFOLDER);
        create_dir_all(bedset_subfolder)?;
        let bedset_cache = BioCache::new(bedset_subfolder);

        Ok(BBClient {
            cache_folder,
            bedbase_api,
            bedfile_cache,
            bedset_cache,
        })
    }

    fn add_resource_to_cache(&mut self, cache_id: &str, cache_path: &str, bedfile: bool) {
        let resource_to_add = NewResource::new(cache_id, cache_path, None, None, None, None);
        if bedfile {
            self.bedfile_cache.add(&resource_to_add);
        } else {
            self.bedset_cache.add(&resource_to_add);
        }
    }

    pub fn load_bed(&mut self, bed_id: &str) -> Result<RegionSet> {
        let bedfile_path = self.bedfile_path(bed_id, Some(false));

        if bedfile_path.exists() {
            println!("Loading cached BED file from {:?}", bedfile_path.display());
            return RegionSet::try_from(bedfile_path);
        }

        let region_set = RegionSet::try_from(bed_id)
            .with_context(|| format!("Failed to create RegionSet from BEDbase id {}", bed_id))?;
        println!("Downloaded BED file from BEDbase: {}", bed_id);

        self.add_resource_to_cache(
            bed_id,
            bedfile_path.to_str().expect("Invalid BED file path"),
            true,
        );
        // self.bedfile_cache.add(bed_resource);
        println!(
            "Downloaded BED file from to path: {}",
            bedfile_path.display()
        );
        region_set.to_bed_gz(bedfile_path)?;
        // println!("Downloaded BED file from to path: {}", bedfile_path.display());
        Ok(region_set)
    }

    pub fn load_bedset(&mut self, bedset_id: &str) -> Result<BedSet> {
        let bedset_path = self.bedset_path(bedset_id, Some(true));

        if bedset_path.exists() {
            println!("Loading cached BED file from {:?}", bedset_path.display());
            return BedSet::try_from(bedset_path);
        }

        let bed_data = self.download_bedset_data(bedset_id).unwrap();
        let mut file = File::create(bedset_path.clone())?;
        let mut region_sets = Vec::new();
        for bbid in bed_data {
            writeln!(file, "{}", bbid)?;
            let rs = self.load_bed(&bbid).unwrap();
            region_sets.push(rs);
        }
        self.add_resource_to_cache(
            bedset_id,
            bedset_path.to_str().expect("Invalid BED set path"),
            false,
        );

        Ok(BedSet::from(region_sets))
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

        let force = force.unwrap_or(false);
        if !force && cache_path.exists() {
            println!("{} already exists in cache", cache_path.display());
            return Ok(regionset);
        }

        regionset.to_bed_gz(cache_path.as_path())?;
        self.add_resource_to_cache(
            &bedfile_id,
            cache_path
                .to_str()
                .expect("cache path cannot be convert to &str"),
            true,
        );

        // println!("{} added to cache", cache_path.display());

        Ok(regionset)
    }

    pub fn add_bedset_to_cache(&mut self, bedset: BedSet) -> Result<String> {
        let bedset_id = bedset.identifier();
        let bedset_path = self.bedset_path(&bedset_id, Some(true));
        if bedset_path.exists() {
            println!("{} already exists in cache", bedset_path.display());
        } else {
            let mut file = File::create(bedset_path.clone())?;
            for rs in bedset.region_sets {
                let bed_id = rs.identifier();
                let _ = self.add_regionset_to_cache(rs, Some(false));
                writeln!(file, "{}", bed_id)?;
            }
        }

        self.add_resource_to_cache(
            &bedset_id,
            bedset_path
                .to_str()
                .expect("cache path cannot be convert to &str"),
            false,
        );

        Ok(bedset_id)
    }

    pub fn add_local_folder_as_bedset(&mut self, folder_path: PathBuf) -> Result<String> {
        let mut region_sets = Vec::new();
        for entry in read_dir(&folder_path).expect("Failed to read directory") {
            let entry = entry.expect("Failed to read directory entry");
            let file_path = entry.path();

            if file_path.is_file() {
                let rs = RegionSet::try_from(file_path).unwrap();
                region_sets.push(rs);
            }
        }
        let bedset = BedSet::from(region_sets);
        Ok(self.add_bedset_to_cache(bedset).unwrap())
    }

    pub fn add_local_file_as_bedset(&mut self, folder_path: PathBuf) -> Result<String> {
        let bedset = BedSet::try_from(folder_path).unwrap();
        Ok(self.add_bedset_to_cache(bedset).unwrap())
    }

    fn download_bedset_data(&self, bedset_id: &str) -> Result<Vec<String>> {
        let bedset_url = format!("{}/v1/bedset/{}/bedfiles", self.bedbase_api, bedset_id);

        let response = get(&bedset_url)?.text()?;

        let json: serde_json::Value = serde_json::from_str(&response)?;

        let results = json["results"]
            .as_array()
            .ok_or_else(|| anyhow!("`results` is not an array"))?;

        let extracted_ids: Vec<String> = results
            .iter()
            .filter_map(|entry| {
                let id_val = entry.get("id");
                println!("[DEBUG] Entry ID field: {:?}", id_val);
                id_val?.as_str().map(|s| s.to_string())
            })
            .collect();

        Ok(extracted_ids)
    }

    fn bedfile_path(&self, bedfile_id: &str, create: Option<bool>) -> PathBuf {
        let subfolder_name = DEFAULT_BEDFILE_SUBFOLDER;
        let file_extension = DEFAULT_BEDFILE_EXT;
        self.cache_path(bedfile_id, subfolder_name, file_extension, create)
    }

    fn bedset_path(&self, bedset_id: &str, create: Option<bool>) -> PathBuf {
        let subfolder_name = DEFAULT_BEDSET_SUBFOLDER;
        let file_extension = DEFAULT_BEDSET_EXT;
        self.cache_path(bedset_id, subfolder_name, file_extension, create)
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
        let file_path = self.bedfile_path(&identifier, Some(false));
        if file_path.exists() {
            Ok(file_path)
        } else {
            let set_path = self.bedset_path(&identifier, Some(false));
            if set_path.exists() {
                Ok(set_path)
            } else {
                Err(anyhow::anyhow!("{} does not exist in cache.", identifier))
            }
        }
    }

    pub fn remove(&mut self, identifier: &str) -> Result<()> {
        let file_path = self.bedfile_path(identifier, Some(false));
        if file_path.exists() {
            // remove file and check if subfolders is cleaned
            let _ = self.local_removal(file_path.clone());
            self.bedfile_cache.remove(identifier);

            println!("{} is removed.", file_path.display());
            Ok(())
        } else {
            let set_path = self.bedset_path(identifier, Some(false));
            if set_path.exists() {
                let bedset_file = File::open(set_path.clone())?;
                let reader = BufReader::new(bedset_file);

                let bed_ids: Vec<String> = reader.lines().collect::<Result<_, _>>()?;

                for bed_id in bed_ids {
                    let _ = self.remove(&bed_id);
                }

                let _ = self.local_removal(set_path.clone());

                self.bedset_cache.remove(identifier);

                println!("{} is removed.", set_path.display());
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

    fn local_removal(&self, file_path: PathBuf) -> Result<()> {
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

        Ok(())
    }

    pub fn list_beds(&mut self) -> Result<Vec<Resource>> {
        let bed_resources = self.bedfile_cache.list_resources(Some(20000 as i64));
        Ok(bed_resources)
        // print_resources(bed_resources);
        // let mut bedfile_map = HashMap::new();

        // let bedfile_dir = self.cache_folder.join(DEFAULT_BEDFILE_SUBFOLDER);
        // if !bedfile_dir.exists() {
        //     return Ok(bedfile_map); // return empty map if folder doesn't exist
        // }

        // for entry in WalkDir::new(&bedfile_dir)
        //     .into_iter()
        //     .filter_map(|e| e.ok())
        // {
        //     let path = entry.path();

        //     if path.is_file() {
        //         if let Some(_ext) = path.extension().and_then(|s| s.to_str()) {
        //             if path
        //                 .file_name()
        //                 .and_then(|f| f.to_str())
        //                 .map(|name| name.ends_with(DEFAULT_BEDFILE_EXT))
        //                 .unwrap_or(false)
        //             {
        //                 if let Some(file_stem) = path.file_name().and_then(|s| s.to_str()) {
        //                     let id = file_stem
        //                         .strip_suffix(DEFAULT_BEDFILE_EXT)
        //                         .unwrap_or(file_stem)
        //                         .to_string();
        //                     bedfile_map.insert(id, path.to_path_buf());
        //                 }
        //             }
        //         }
        //     }
        // }
        // Ok(bedfile_map)
    }

    pub fn list_bedsets(&mut self) -> Result<Vec<Resource>> {
        let bedset_resources = self.bedset_cache.list_resources(Some(20000 as i64));
        Ok(bedset_resources)
    }
}
