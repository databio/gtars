use biocrs::common::print_resources;
use clap::ArgMatches;
use client::BBClient;
use std::path::PathBuf;

pub mod cli;
pub mod client;
pub mod consts;
pub mod utils;

use crate::bbcache::utils::get_default_cache_folder;

/// Excute the input commands from CLI
/// # Arguments
/// - subcmd: the subcommand under bbcache
/// - matches: matched items from CLAP args
pub fn run_bbcache(subcmd: &str, matches: &ArgMatches) {
    let cache_folder = matches
        .get_one::<String>("cache-folder")
        .map(PathBuf::from)
        .unwrap_or_else(get_default_cache_folder);
    let mut bbc = BBClient::new(Some(cache_folder), None).expect("Failed to create BBClient");

    match subcmd {
        consts::BBCACHE_INSPECTBED => {
            let bed_resources = bbc.list_beds().unwrap();
            let n = bed_resources.len();
            print_resources(bed_resources);
            println!("Number of BED files: {}", n);
        }
        consts::BBCACHE_INSPECTBEDSET => {
            let bedset_resources = bbc.list_bedsets().unwrap();
            let n = bedset_resources.len();
            print_resources(bedset_resources);
            println!("Number of BED sets: {}", n);
        }
        consts::BBCACHE_SEEK => {
            let bed_id = matches
                .get_one::<String>("identifier")
                .expect("BED file identifier is required");

            let path = bbc.seek(bed_id).expect("Failed to seek BED file in cache");
            println!("{}", path.display());
        }
        consts::BBCACHE_CACHEBED => {
            let bed_id = matches
                .get_one::<String>("identifier")
                .expect("BED file identifier is required");

            let path = PathBuf::from(bed_id);

            if path.is_dir() {
                println!(
                    "Detected '{}' as a directory. Please use `cache-bedset` command to cache all files as a BedSet",
                    path.display()
                );
            } else if path.is_file() {
                println!(
                    "Detected '{}' as a local file. Adding to cache...",
                    path.display()
                );
                bbc.add_local_bed_to_cache(path, None)
                    .expect("Failed to add local BED file to cache");
            } else {
                println!(
                    "'{}' not found locally. Attempting to load from BEDbase...",
                    path.display()
                );
                bbc.load_bed(bed_id)
                    .expect("Failed to load BED file from BEDbase");
            }
        }
        consts::BBCACHE_CACHEBEDSET => {
            let bedset_id = matches
                .get_one::<String>("identifier")
                .expect("BED file identifier is required");

            let path = PathBuf::from(bedset_id);

            if path.is_dir() {
                bbc.add_local_folder_as_bedset(path)
                    .expect("Failed to cache BED set from folder");
            } else if path.is_file() {
                println!(
                    "Detected '{}' as a local file. Adding to cache...",
                    path.display()
                );
                bbc.add_local_file_as_bedset(path)
                    .expect("Failed to add local file as a BED set to cache");
            } else {
                println!(
                    "'{}' not found locally. Attempting to load from BEDbase...",
                    path.display()
                );
                bbc.load_bedset(bedset_id)
                    .expect("Failed to load BED file from BEDbase");
            }
        }
        consts::BBCACHE_REMOVE => {
            let bed_id = matches
                .get_one::<String>("identifier")
                .expect("BED file identifier is required");

            bbc.remove(bed_id)
                .expect("Failed to remove BED file from cache");
        }
        _ => unreachable!("Unknown BBCACHE subcommand: {subcmd}"),
    }
}
