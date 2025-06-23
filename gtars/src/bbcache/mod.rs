use clap::ArgMatches;
use client::BBClient;
use log::info;
use std::path::{Path, PathBuf};
use tabled::settings::Style;
use tabled::Table;

pub mod cli;
pub mod client;
pub mod consts;
pub mod utils;

pub fn run_bbcache(matches: &ArgMatches) {
    let cache_folder_str = matches
        .get_one::<String>("cache-folder")
        .expect("cache folder path is required");

    let cache_folder = PathBuf::from(cache_folder_str);

    let mut bbc = BBClient::new(Some(cache_folder), None).expect("Failed to create BBClient");

    match matches.subcommand() {
        Some((consts::BBCACHE_INSPECTBED, matches)) => {
            let rows: Vec<(String, String)> = match bbc.list_beds() {
                Ok(map) => map
                    .iter()
                    .filter_map(|(id, path)| path.to_str().map(|p| (id.clone(), p.to_string())))
                    .collect(),
                Err(e) => {
                    eprintln!("Error reading cache: {e}");
                    return;
                }
            };

            if rows.is_empty() {
                println!("No cached bed files found.");
            } else {
                let styled_table = {
                    let mut table = Table::new(&rows);
                    table.with(Style::rounded());
                    table
                };

                println!("{styled_table}");
            }
        }
        Some((consts::BBCACHE_SEEK, matches)) => {
            let bed_id = matches
                .get_one::<String>("identifier")
                .expect("BED file identifier, url, or file path is required");
            info!(
                "{}",
                bbc.seek(bed_id.as_str())
                    .expect("Failed to seek BED file in cache")
                    .display()
            );
        }
        Some((consts::BBCACHE_CACHEBED, matches)) => {
            //cache bed file
            let bed_id = matches
                .get_one::<String>("identifier")
                .expect("BED file identifier, url, or file path is required");

            if Path::new(bed_id).exists() {
                // If the identifier is a local file path
                let bedfile_path = PathBuf::from(bed_id);
                bbc.add_local_bed_to_cache(bedfile_path, None)
                    .expect("Failed to add local BED file to cache");
            } else {
                // If the identifier is a BEDbase ID or URL
                let _ = bbc
                    .load_bed(bed_id.as_str())
                    .expect("Failed to load BED file from BEDbase");
            }
        }
        Some((consts::BBCACHE_REMOVE, matches)) => {
            let bed_id = matches
                .get_one::<String>("identifier")
                .expect("BED file identifier, url, or file path is required");

            bbc.remove(bed_id.as_str())
                .expect("Failed to remove BED file from cache");
        }
        _ => unreachable!("BBCACHE Subcommand not found"),
    };
}
