mod igd;
mod overlaprs;
mod uniwig;
mod bbcache;

use anyhow::Result;
use clap::Command;

pub mod consts {
    pub const VERSION: &str = env!("CARGO_PKG_VERSION");
    pub const PKG_NAME: &str = "gtars";
    pub const BIN_NAME: &str = "gtars";
}

fn build_parser() -> Command {
    Command::new(consts::BIN_NAME)
        .bin_name(consts::BIN_NAME)
        .version(consts::VERSION)
        .author("Databio")
        .about("Performance critical tools for working with genomic interval data with an emphasis on preprocessing for machine learning pipelines.")
        .subcommand_required(true)
        .subcommand(uniwig::cli::create_uniwig_cli())
        .subcommand(overlaprs::cli::create_overlap_cli())
        .subcommand(bbcache::cli::create_bbcache_cli())
}

fn main() -> Result<()> {
    let app = build_parser();
    let matches = app.get_matches();

    match matches.subcommand() {
        //
        // UNIWIG
        //
        Some((uniwig::cli::UNIWIG_CMD, matches)) => {
            uniwig::handlers::run_uniwig(matches);
        }

        //
        // OVERLAPPES
        //
        Some((overlaprs::cli::OVERLAP_CMD, matches)) => {
            let _ = overlaprs::handlers::run_overlaprs(matches);
        }

        //
        // BBCACHE
        //
        Some((bbcache::cli::BBCACHE_CMD, matches)) => {
            bbcache::handlers::run_bbcache(matches);
        }

        //
        // IGD
        //
        Some((igd::cli::IGD_CMD, matches)) => match matches.subcommand() {
            Some((igd::cli::IGD_CREATE, matches)) => {
                igd::handlers::igd_get_create_matches(matches);
            }
            Some((igd::cli::IGD_SEARCH, matches)) => {
                igd::handlers::igd_get_search_matches(matches);
            }
            _ => unreachable!("IGD Subcommand not found"),
        },

        _ => unreachable!("Subcommand not found"),
    };

    Ok(())
}
