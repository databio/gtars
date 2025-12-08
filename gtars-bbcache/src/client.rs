//! BEDbase caching client implementation.
//!
//! This module provides the core [`BBClient`] type and its builder for managing
//! cached BED files and BED sets from the BEDbase API.

use anyhow::{Context, Ok, Result, anyhow};
use biocrs::biocache::BioCache;
use biocrs::models::{NewResource, Resource};

//use reqwest::blocking::get;
use ureq::get;
use std::fs::{File, create_dir_all, read_dir, remove_dir, remove_file};
use std::io::{BufRead, BufReader, Error, ErrorKind, Write};
use std::path::{Path, PathBuf};

use super::consts::{
    DEFAULT_BEDFILE_EXT, DEFAULT_BEDFILE_SUBFOLDER, DEFAULT_BEDSET_EXT, DEFAULT_BEDSET_SUBFOLDER,
};
use super::utils::{get_default_bedbase_api, get_default_cache_folder};
use gtars_core::models::bed_set::BedSet;
use gtars_core::models::region_set::RegionSet;

/// Builder for constructing a [`BBClient`] with custom configuration.
///
/// Use this builder to configure cache location and BEDbase API endpoint
/// before creating a client instance.
///
/// # Examples
///
/// ```rust,no_run
/// use gtars_bbcache::client::BBClient;
/// use std::path::PathBuf;
///
/// # fn main() -> anyhow::Result<()> {
/// let client = BBClient::builder()
///     .with_cache_folder(PathBuf::from("/custom/cache"))
///     .with_bedbase_api("https://api.bedbase.org".to_string())
///     .finish()?;
/// # Ok(())
/// # }
/// ```
#[derive(Default)]
pub struct BBClientBuilder {
    cache_folder: Option<PathBuf>,
    bedbase_api: Option<String>,
}

impl BBClientBuilder {
    /// Creates a new, empty BBClientBuilder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the cache folder for the BBClient.
    pub fn with_cache_folder(mut self, path: PathBuf) -> Self {
        self.cache_folder = Some(path);
        self
    }

    /// Sets the BEDbase API URL for the BBClient.
    pub fn with_bedbase_api(mut self, api: String) -> Self {
        self.bedbase_api = Some(api);
        self
    }

    /// Consumes the builder and creates a BBClient.
    pub fn finish(self) -> Result<BBClient> {
        // handle the cache dir
        let raw_path_to_cache_folder = self.cache_folder.unwrap_or_else(get_default_cache_folder);
        let raw_str_to_cache_folder = raw_path_to_cache_folder.to_string_lossy().into_owned();
        let expanded_str = shellexpand::env(&raw_str_to_cache_folder)
            .unwrap_or_else(|_| raw_str_to_cache_folder.clone().into())
            .into_owned();
        let abs_path_to_cache_folder = PathBuf::from(expanded_str);
        create_dir_all(&abs_path_to_cache_folder)?;

        // handle the bedbase api
        let bedbase_api = self.bedbase_api.unwrap_or_else(get_default_bedbase_api);

        // create sub folders
        let bedfile_subfolder = &abs_path_to_cache_folder.join(DEFAULT_BEDFILE_SUBFOLDER);
        create_dir_all(bedfile_subfolder)?;
        let bedfile_cache = BioCache::new(bedfile_subfolder);

        let bedset_subfolder = &abs_path_to_cache_folder.join(DEFAULT_BEDSET_SUBFOLDER);
        create_dir_all(bedset_subfolder)?;
        let bedset_cache = BioCache::new(bedset_subfolder);

        Ok(BBClient {
            cache_folder: abs_path_to_cache_folder,
            bedbase_api,
            bedfile_cache,
            bedset_cache,
        })
    }
}

/// Client for managing BED file and BED set caching from BEDbase.
///
/// `BBClient` provides a high-level interface for:
/// - Downloading and caching BED files from the BEDbase API
/// - Managing local BED file collections
/// - Organizing BED sets (collections of related BED files)
/// - Querying and removing cached resources
///
/// The client maintains two separate caches:
/// - **bedfile_cache**: Individual BED files (stored as `.bed.gz`)
/// - **bedset_cache**: BED set metadata files (stored as `.txt` with lists of BED IDs)
///
///
/// # Examples
///
/// ```rust,no_run
/// use gtars_bbcache::client::BBClient;
///
/// # fn main() -> anyhow::Result<()> {
/// // Create a client with default settings
/// let mut client = BBClient::builder().finish()?;
///
/// // Download and cache a BED file from BEDbase
/// let region_set = client.load_bed("6b2e163a1d4319d99bd465c6c78a9741")?;
///
/// // Check if file exists in cache
/// let path = client.seek("6b2e163a1d4319d99bd465c6c78a9741")?;
/// println!("Cached at: {:?}", path);
///
/// // List all cached BED files
/// let beds = client.list_beds()?;
/// println!("Found {} cached BED files", beds.len());
///
/// // Remove from cache
/// client.remove("6b2e163a1d4319d99bd465c6c78a9741")?;
/// # Ok(())
/// # }
/// ```
pub struct BBClient {
    /// Path to the root cache directory
    pub cache_folder: PathBuf,
    /// BEDbase API endpoint URL
    pub bedbase_api: String,
    /// Internal cache manager for BED files
    bedfile_cache: BioCache,
    /// Internal cache manager for BED sets
    bedset_cache: BioCache,
}

impl BBClient {
    /// Creates a new builder for constructing a [`BBClient`].
    ///
    /// The builder pattern allows you to configure the cache folder and API endpoint
    /// before creating the client instance.
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use gtars_bbcache::client::BBClient;
    /// use std::path::PathBuf;
    ///
    /// # fn main() -> anyhow::Result<()> {
    /// let client = BBClient::builder()
    ///     .with_cache_folder(PathBuf::from("/tmp/bbcache"))
    ///     .finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn builder() -> BBClientBuilder {
        BBClientBuilder::default()
    }

    /// Add id and path of cached bed file or bed set into file cacher (SQLite).
    /// # Arguments
    /// - cache_id: id of bed file or bed set
    /// - cache_path: path to the cached bed file or bed set
    /// - bedfile: if the the id and file is a bed file
    ///
    fn add_resource_to_cache(&mut self, cache_id: &str, cache_path: &str, bedfile: bool) {
        let resource_to_add = NewResource::new(cache_id, cache_path).set_fpath(cache_path);
        if bedfile {
            self.bedfile_cache.add(&resource_to_add);
        } else {
            self.bedset_cache.add(&resource_to_add);
        }
    }

    /// Loads a BED file from cache, or downloads and caches it if it doesn't exist
    /// # Arguments
    /// - bed_id: unique identifier of a BED file
    ///
    /// # Returns
    /// - the RegionSet object of the loaded bed file
    pub fn load_bed(&mut self, bed_id: &str) -> Result<RegionSet> {
        let bedfile_path = self.bedfile_path(bed_id, Some(false));

        if bedfile_path.exists() {
            println!("Loading cached BED file from {:?}", bedfile_path.display());
            return RegionSet::try_from(bedfile_path);
        }

        let region_set = RegionSet::try_from(bed_id)
            .with_context(|| format!("Failed to create RegionSet from BEDbase id {}", bed_id))?;

        self.add_resource_to_cache(
            bed_id,
            bedfile_path.to_str().expect("Invalid BED file path"),
            true,
        );

        region_set.to_bed_gz(bedfile_path.clone())?;
        println!(
            "Downloaded BED file {} from BEDbase to path: {}",
            bed_id,
            bedfile_path.display()
        );
        Ok(region_set)
    }

    /// Load a BEDset from cache, or download and add it to the cache with its BED files
    /// # Arguments
    /// - bedset_id: unique identifier of a BED set
    ///
    /// # Returns
    /// - the BedSet object of the loaded bed set
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

    ///  Add a BED file to the cache
    /// # Arguments
    /// - bedfile: a path or url to the BED file
    /// - force: whether to overwrite the existing file in cache
    ///
    /// # Returns
    /// - the RegionSet identifier
    pub fn add_local_bed_to_cache(
        &mut self,
        bedfile: PathBuf,
        force: Option<bool>,
    ) -> Result<String> {
        let regionset = RegionSet::try_from(bedfile.as_path())?;
        self.add_regionset_to_cache(regionset, force)
    }

    ///  Add a RegionSet object to the cache
    /// # Arguments
    /// - regionset:  a RegionSet object
    /// - force: whether to overwrite the existing file in cache
    ///
    /// # Returns
    /// - the RegionSet identifier
    pub fn add_regionset_to_cache(
        &mut self,
        regionset: RegionSet,
        force: Option<bool>,
    ) -> Result<String> {
        let bedfile_id = regionset.identifier();
        let cache_path = self.bedfile_path(&bedfile_id, Some(true));

        let force = force.unwrap_or(false);
        if !force && cache_path.exists() {
            println!("{} already exists in cache", cache_path.display());
            return Ok(bedfile_id);
        }

        regionset.to_bed_gz(cache_path.as_path())?;
        self.add_resource_to_cache(
            &bedfile_id,
            cache_path
                .to_str()
                .expect("cache path cannot be convert to &str"),
            true,
        );
        println!("BED file cached to {}", cache_path.display());

        Ok(bedfile_id)
    }

    ///  Add a BED set to the cache
    /// # Arguments
    /// - bedset: the BED set to be added, a BedSet class
    ///
    /// # Returns
    /// - the identifier if the BedSet object
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
        println!("BED set cached to {}", bedset_path.display());

        Ok(bedset_id)
    }

    ///  Add a folder of bed files to the cache as a bedset
    /// # Arguments
    /// - folder_path: path to the folder where bed files are stored
    ///
    /// # Returns
    /// - the identifier if the BedSet object
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

    ///  Add a local file that contains bed file paths as a bed set
    /// # Arguments
    /// - file_path: path to the file of bedset info
    ///
    /// # Returns
    /// - the identifier if the BedSet object
    pub fn add_local_file_as_bedset(&mut self, file_path: PathBuf) -> Result<String> {
        let bedset = BedSet::try_from(file_path).unwrap();
        Ok(self.add_bedset_to_cache(bedset).unwrap())
    }

    ///  Download BED set from BEDbase API and return the list of identifiers of BED files in the set
    /// # Arguments
    /// - bedset_id: unique identifier of a BED set
    ///
    /// # Returns
    /// - the list of identifiers of BED files in the set
    fn download_bedset_data(&self, bedset_id: &str) -> Result<Vec<String>> {
        let bedset_url = format!("{}/v1/bedset/{}/bedfiles", self.bedbase_api, bedset_id);

        let response = get(&bedset_url)
            .call()
            .map_err(|e| anyhow!("Failed to GET {}: {}", bedset_url, e))?
            .into_string()
            .map_err(|e| anyhow!("Failed to read response body for {}: {}", bedset_url, e))?;

        let json: serde_json::Value = serde_json::from_str(&response)?;

        let results = json["results"]
            .as_array()
            .ok_or_else(|| anyhow!("`results` is not an array"))?;

        let extracted_ids: Vec<String> = results
            .iter()
            .filter_map(|entry| {
                let id_val = entry.get("id");
                id_val?.as_str().map(|s| s.to_string())
            })
            .collect();

        Ok(extracted_ids)
    }

    ///  Get the path of a BED file's .bed.gz file with given identifier
    /// # Arguments
    /// - bedfile_id: the identifier of BED file
    /// - create: whether the cache path needs creating
    ///
    /// # Returns
    /// - the path to the .bed.gz file
    fn bedfile_path(&self, bedfile_id: &str, create: Option<bool>) -> PathBuf {
        let subfolder_name = DEFAULT_BEDFILE_SUBFOLDER;
        let file_extension = DEFAULT_BEDFILE_EXT;
        self.cache_path(bedfile_id, subfolder_name, file_extension, create)
    }

    ///  Get the path of a BED set's .txt file with given identifier
    /// # Arguments
    /// - bedset_id: the identifier of BED set
    /// - create: whether the cache path needs creating
    ///
    /// # Returns
    /// - the path to the .txt file
    fn bedset_path(&self, bedset_id: &str, create: Option<bool>) -> PathBuf {
        let subfolder_name = DEFAULT_BEDSET_SUBFOLDER;
        let file_extension = DEFAULT_BEDSET_EXT;
        self.cache_path(bedset_id, subfolder_name, file_extension, create)
    }

    ///  Get the path of a file in cache folder
    /// # Arguments
    /// - identifier: the identifier of BED set or BED file
    /// - subfolder_name: "bedsets" or "bedfiles"
    /// - file_extension: ".txt" or ".bed.gz"
    /// - create: whether the cache path needs creating
    ///
    /// # Returns
    /// - the path to the file
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

    ///  Create cache folder if it doesn't exist
    /// # Arguments
    /// - subfolder_path: path to the subfolder
    fn create_cache_folder(&self, subfolder_path: Option<&Path>) {
        let path = match subfolder_path {
            Some(p) => p.to_path_buf(),
            None => self.cache_folder.clone(),
        };

        if !path.exists() {
            create_dir_all(&path).expect("Failed to create cache folder");
        }
    }

    /// Get local path to BED file or BED set with specific identifier
    /// # Arguments
    /// - identifier: the unique identifier
    ///
    /// # Returns
    /// - the local path of the file
    pub fn seek(&self, identifier: &str) -> Result<PathBuf> {
        let file_path = self.bedfile_path(identifier, Some(false));
        if file_path.exists() {
            Ok(file_path)
        } else {
            let set_path = self.bedset_path(identifier, Some(false));
            if set_path.exists() {
                Ok(set_path)
            } else {
                Err(anyhow::anyhow!("{} does not exist in cache.", identifier))
            }
        }
    }

    /// Remove a BED file or BED set from the cache folder as well as biocfilcache (SQLite)
    /// # Arguments
    /// - identifier: the identifier of BED file / BED set to be removed
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

    /// Remove a BED file or BED set from the cache folder, and make sure the removal won't cause empty subfolders
    /// # Arguments
    /// - identifier: the identifier of BED file / BED set to be removed
    fn local_removal(&self, file_path: PathBuf) -> Result<()> {
        let sub_folder_2 = file_path.parent().map(PathBuf::from);
        let sub_folder_1 = sub_folder_2
            .as_ref()
            .and_then(|p| p.parent().map(PathBuf::from));

        remove_file(&file_path)?;

        // Attempt to remove empty subfolders
        if let Some(sub2) = sub_folder_2
            && read_dir(&sub2)?.next().is_none()
        {
            remove_dir(&sub2)?;
            if let Some(sub1) = sub_folder_1
                && read_dir(&sub1)?.next().is_none()
            {
                remove_dir(&sub1)?;
            }
        }

        Ok(())
    }

    /// List identifiers and paths of all BED files in cache
    /// # Returns
    /// - the list of resource stored in biocfilecache (identifiers & paths)
    pub fn list_beds(&mut self) -> Result<Vec<Resource>> {
        let bed_resources = self.bedfile_cache.list_resources(Some(20_000));
        Ok(bed_resources)
    }

    /// List identifiers and paths of all BED sets in cache
    /// # Returns
    /// - the list of resource stored in biocfilecache (identifiers & paths)
    pub fn list_bedsets(&mut self) -> Result<Vec<Resource>> {
        let bedset_resources = self.bedset_cache.list_resources(Some(20_000));
        Ok(bedset_resources)
    }
}
