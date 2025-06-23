use clap::ArgMatches;
use client::BBClient;
use std::path::{Path, PathBuf};
use tabled::settings::Style;
use tabled::{Table, Tabled};

pub mod cli;
pub mod client;
pub mod consts;
pub mod utils;

#[derive(Tabled)]
struct BedEntry {
    #[tabled(rename = "ID")]
    id: String,
    #[tabled(rename = "Path")]
    path: String,
}

pub fn run_bbcache(subcmd: &str, matches: &ArgMatches) {
    let cache_folder_str = matches
        .get_one::<String>("cache-folder")
        .expect("cache folder path is required");

    let cache_folder = PathBuf::from(cache_folder_str);
    let mut bbc = BBClient::new(Some(cache_folder), None).expect("Failed to create BBClient");

    match subcmd {
        consts::BBCACHE_INSPECTBED => {
            let rows: Vec<BedEntry> = match bbc.list_beds() {
                Ok(map) => map
                    .iter()
                    .filter_map(|(id, path)| {
                        path.to_str().map(|p| BedEntry {
                            id: id.clone(),
                            path: p.to_string(),
                        })
                    })
                    .collect(),
                Err(e) => {
                    eprintln!("Error reading cache: {e}");
                    return;
                }
            };

            if rows.is_empty() {
                println!("Number of bed files: 0");
            } else {
                let styled_table = {
                    let mut table = Table::new(&rows);
                    table.with(Style::rounded());
                    table
                };

                println!("{styled_table}");
                println!("Number of bed files: {}", rows.len())
            }
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

            if Path::new(bed_id).exists() {
                let bedfile_path = PathBuf::from(bed_id);
                bbc.add_local_bed_to_cache(bedfile_path, None)
                    .expect("Failed to add local BED file to cache");
            } else {
                bbc.load_bed(bed_id)
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
