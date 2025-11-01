use clap::{Arg, Command};

pub const BBCACHE_CMD: &str = "bbcache";
pub const BBCACHE_CACHEBED: &str = "cache-bed";
pub const BBCACHE_CACHEBEDSET: &str = "cache-bedset";
pub const BBCACHE_SEEK: &str = "seek";
pub const BBCACHE_INSPECTBED: &str = "inspect-bedfiles";
pub const BBCACHE_INSPECTBEDSET: &str = "inspect-bedsets";
pub const BBCACHE_REMOVE: &str = "rm";

pub fn create_bbcache_cli() -> Command {
    Command::new(BBCACHE_CMD)
        .author("ZH")
        .about("Downloads, processes, and caches BED files from the BEDbase API")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new(BBCACHE_CACHEBED)
                .about("Cache a BED file from local file or BEDbase")
                .arg(
                    Arg::new("identifier")
                        .long("identifier")
                        .short('i')
                        .required(true)
                        .help("BED file identifier, url, or file path"),
                )
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .short('f')
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new(BBCACHE_CACHEBEDSET)
                .about("Cache a BED set from local file or BEDbase")
                .arg(
                    Arg::new("identifier")
                        .long("identifier")
                        .short('i')
                        .required(true)
                        .help("BED set identifier, url, or file path"),
                )
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .short('f')
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new(BBCACHE_SEEK)
                .about("Seek the BED file path by giving identifier")
                .arg(
                    Arg::new("identifier")
                        .long("identifier")
                        .short('i')
                        .required(true)
                        .help("BED file identifier, url, or file path"),
                )
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .short('f')
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new(BBCACHE_INSPECTBED)
                .about("Inspect the contents of bedfile cache folder")
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .short('f')
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new(BBCACHE_INSPECTBEDSET)
                .about("Inspect the contents of bedsets cache folder")
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .short('f')
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new(BBCACHE_REMOVE)
                .about("Remove the BED file or BED set from cache with given identifier")
                .arg(
                    Arg::new("identifier")
                        .long("identifier")
                        .short('i')
                        .required(true)
                        .help("BED file identifier, url, or file path"),
                )
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .short('f')
                        .help("Cache folder path"),
                ),
        )
}
