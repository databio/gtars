use once_cell::sync::Lazy;
use std::env;

pub const BBCLIENT_CACHE_ENV: &str = "BBCLIENT_CACHE";
pub const BEDBASE_API_ENV: &str = "BEDBASE_API";

pub const BBCACHE_CMD: &str = "bbcache";
pub const BBCACHE_CACHEBED: &str = "cache-bed";
pub const BBCACHE_SEEK: &str = "seek";
pub const BBCACHE_INSPECTBED: &str = "inspect-bedfiles";
pub const BBCACHE_REMOVE: &str = "rm";

pub const DEFAULT_BEDFILE_SUBFOLDER: &str = "bedfiles";
pub const DEFAULT_BEDSET_SUBFOLDER: &str = "bedsets";

pub const DEFAULT_BEDFILE_EXT: &str = ".bed.gz";
pub const DEFAULT_BEDSET_EXT: &str = ".txt";
