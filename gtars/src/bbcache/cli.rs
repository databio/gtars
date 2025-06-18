use clap::{Arg, ArgAction, Command};

use crate::bbclient::consts::{BBCLIENT_CMD, DEFAULT_CACHE_FOLDER};

static DEFAULT_CACHE_FOLDER: Lazy<String> = Lazy::new(|| {
    env::var(BBCLIENT_CACHE_ENV).unwrap_or_else(|_| {
        let home = env::var("HOME").unwrap_or_else(|_| String::from("/tmp"));
        format!("{}/.bbcache", home)
    })
});

pub fn create_bbclient_cli() -> Command {
    Command::new(BBCLIENT_CMD)
        .author("DRC")
        .about("downloads, processes, and caches BED files and BED sets from the BEDbase API")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("cache-bed")
                .about("Cache a BED file from local file or BEDbase")
                .arg(
                    Arg::new("identifier")
                        .required(true)
                        .num_args(1)
                        .help("BED file identifier, url, or file path"),
                )
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .default_value(DEFAULT_CACHE_FOLDER)
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("cache-tokens")
                .about("Cache tokens from local file or BEDbase")
                .arg(
                    Arg::new("bed-id")
                        .long("bed-id")
                        .required(true)
                        .num_args(1)
                        .help("Token file identifier, url, or file path"),
                )
                .arg(
                    Arg::new("universe-id")
                        .long("universe-id")
                        .required(true)
                        .num_args(1)
                        .help("Unique identifier for the universe of the token file"),
                )
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .default_value(DEFAULT_CACHE_FOLDER)
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("seek")
                .about("Seek the BED / BEDset path by giving identifier")
                .arg(
                    Arg::new("identifier")
                        .required(true)
                        .num_args(1)
                        .help("BED file identifier, url, or file path"),
                )
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .default_value(DEFAULT_CACHE_FOLDER)
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("inspect-bedfiles")
                .about("Inspect the contents of bedfile cache folder")
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .default_value(DEFAULT_CACHE_FOLDER)
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("inspect-bedsets")
                .about("Inspect the contents of bedset cache folder")
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .default_value(DEFAULT_CACHE_FOLDER)
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("rm")
                .about("Remove the BED/BEDset from cache with given identifier")
                .arg(
                    Arg::new("identifier")
                        .required(true)
                        .num_args(1)
                        .help("BED file identifier, url, or file path"),
                )
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .default_value(DEFAULT_CACHE_FOLDER)
                        .help("Cache folder path"),
                ),
        )
}
