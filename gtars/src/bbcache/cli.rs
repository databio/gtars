use clap::{Arg, Command};

use crate::bbcache::consts::{BBCACHE_CMD, DEFAULT_CACHE_FOLDER};

pub fn create_bbcache_cli() -> Command {
    Command::new(BBCACHE_CMD)
        .author("DRC")
        .about("downloads, processes, and caches BED files from the BEDbase API")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("cache-bed")
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
                        .default_value(DEFAULT_CACHE_FOLDER.to_str().expect("Invalid UTF-8 path"))
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("seek")
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
                        .default_value(DEFAULT_CACHE_FOLDER.to_str().expect("Invalid UTF-8 path"))
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("inspect-bedfiles")
                .about("Inspect the contents of bedfile cache folder")
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .short('f')
                        .default_value(DEFAULT_CACHE_FOLDER.to_str().expect("Invalid UTF-8 path"))
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("rm")
                .about("Remove the BED file from cache with given identifier")
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
                        .default_value(DEFAULT_CACHE_FOLDER.to_str().expect("Invalid UTF-8 path"))
                        .help("Cache folder path"),
                ),
        )
}
