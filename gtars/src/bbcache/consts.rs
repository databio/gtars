use once_cell::sync::Lazy;
use std::env;
use std::path::PathBuf;
use dirs::home_dir;

/// Equivalent of Python's BBCLIENT_CACHE_ENV constant
pub const BBCLIENT_CACHE_ENV: &str = "BBCLIENT_CACHE";

pub const BEDBASE_URL_PATTERN: &str = "{bedbase_api}/v1/objects/bed.{bed_id}.bed_file/access/http/bytes";

pub const DEFAULT_BEDFILE_SUBFOLDER: &str = "bedfiles";

pub const DEFAULT_BEDFILE_EXT: &str = ".bed.gz";

/// Equivalent of Python's DEFAULT_CACHE_FOLDER, but lazy-evaluated once at runtime
pub static DEFAULT_CACHE_FOLDER: Lazy<PathBuf> = Lazy::new(|| {
    match env::var(BBCLIENT_CACHE_ENV) {
        Ok(val) => PathBuf::from(val),
        Err(_) => {
            let home = env::var("HOME")
                .or_else(|_| {
                    home_dir()
                        .map(|p| p.to_string_lossy().into_owned())
                        .ok_or_else(|| std::env::VarError::NotPresent)
                })
                .unwrap_or_else(|_| "/tmp".to_string());

            let mut path = PathBuf::from(home);
            path.push(".bbcache/");
            path
        }
    }
});


pub static DEFAULT_BEDBASE_API: Lazy<&'static str> = Lazy::new(|| {
    Box::leak(Box::new(
        env::var("BEDBASE_API").unwrap_or_else(|_| "https://api.bedbase.org".to_string())
    ))
});
