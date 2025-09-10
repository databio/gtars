pub use bbcache::consts::BBCACHE_CMD;
use clap::{Arg, Command};

pub fn create_bbcache_cli() -> Command {
    Command::new(BBCACHE_CMD)
        .author("ZH")
        .about("Downloads, processes, and caches BED files from the BEDbase API")
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
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("cache-bedset")
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
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("inspect-bedsets")
                .about("Inspect the contents of bedsets cache folder")
                .arg(
                    Arg::new("cache-folder")
                        .long("cache-folder")
                        .short('f')
                        .help("Cache folder path"),
                ),
        )
        .subcommand(
            Command::new("rm")
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
