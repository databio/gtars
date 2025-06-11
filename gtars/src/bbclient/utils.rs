use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use once_cell::sync::Lazy;
use std::sync::Mutex;

static DEFAULT_CACHE_FOLDER: &str = "./cache";
static CACHE_MUTEX: Lazy<Mutex<()>> = Lazy::new(|| Mutex::new(()));

pub struct BedCacheManager {
    pub cache_folder: String,
}

impl BedCacheManager {
    pub fn new(cache_folder: &str) -> Self {
        let manager = BedCacheManager {
            cache_folder: cache_folder.to_string(),
        };
        manager.create_cache_folder(None);
        manager
    }

    pub fn create_cache_folder(&self, subfolder_path: Option<&str>) {
        let path = match subfolder_path {
            Some(p) => PathBuf::from(p),
            None => PathBuf::from(&self.cache_folder),
        };
        if !path.exists() {
            if let Err(e) = fs::create_dir_all(&path) {
                eprintln!("Failed to create cache folder: {}", e);
            }
        }
    }

    pub fn process_local_bed_data(file_path: &str) -> io::Result<Vec<u8>> {
        let mut file = File::open(file_path)?;
        let mut content = Vec::new();
        file.read_to_end(&mut content)?;
        Ok(Self::decompress_bed_bytes(content))
    }

    pub fn decompress_bed_bytes(content: Vec<u8>) -> Vec<u8> {
        if content.starts_with(&[0x1F, 0x8B]) {
            let mut d = GzDecoder::new(&content[..]);
            let mut out = Vec::new();
            if let Err(e) = d.read_to_end(&mut out) {
                eprintln!("Failed to decompress gzipped content: {}", e);
                return Vec::new();
            }
            out
        } else {
            content
        }
    }
}

pub fn get_abs_path(path: &str, create_folder: bool) -> String {
    let expanded_path = shellexpand::env(path).unwrap_or_else(|_| path.into());
    let abs_path = Path::new(&expanded_path).canonicalize().unwrap_or_else(|_| PathBuf::from(&expanded_path));

    if create_folder {
        let _lock = CACHE_MUTEX.lock().unwrap();
        if let Err(e) = fs::create_dir_all(&abs_path) {
            eprintln!("Failed to create folder: {}", e);
        }
    }

    abs_path.to_string_lossy().to_string()
}
