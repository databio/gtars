use anyhow::Result;
use clap::Command;
use std::env;
use once_cell::sync::Lazy;

// go through the library crate to get the interfaces
use gtars::fragsplit;
use gtars::igd;
use gtars::scoring;
use gtars::uniwig;
use gtars::bbclient;

pub mod consts {
    pub const VERSION: &str = env!("CARGO_PKG_VERSION");
    pub const PKG_NAME: &str = env!("CARGO_PKG_NAME");
    pub const BIN_NAME: &str = env!("CARGO_PKG_NAME");
    pub const UNIWIG_CMD: &str = "uniwig";
    pub const BBCLIENT_CMD: &str = "bbclient";

    pub static DEFAULT_CACHE_FOLDER: Lazy<&'static str> = Lazy::new(|| {
        Box::leak(Box::new(env::var("BBCLIENT_CACHE")
            .unwrap_or_else(|_| format!("{}/.bbcache", env::var("HOME").unwrap_or_else(|_| String::from("/tmp"))))))
    });
}

fn build_parser() -> Command {
    Command::new(consts::BIN_NAME)
        .bin_name(consts::BIN_NAME)
        .version(consts::VERSION)
        .author("Databio")
        .about("Performance critical tools for working with genomic interval data with an emphasis on preprocessing for machine learning pipelines.")
        .subcommand_required(true)
        .subcommand(fragsplit::cli::make_fragsplit_cli())
        .subcommand(uniwig::cli::create_uniwig_cli())
        .subcommand(bbclient::cli::create_bbclient_cli())
        .subcommand(igd::cli::create_igd_cli())
        .subcommand(scoring::cli::make_fscoring_cli())
}

fn main() -> Result<()> {
    let app = build_parser();
    let matches = app.get_matches();

    match matches.subcommand() {
        Some((fragsplit::consts::FRAGSPLIT_CMD, matches)) => {
            fragsplit::cli::handlers::split_fragment_files(matches)?;
        }
        Some((uniwig::consts::UNIWIG_CMD, matches)) => {
            uniwig::run_uniwig(matches);
        }
        Some((bbclient::consts::BBCLIENT_CMD, matches)) => {
            uniwig::run_bbclient(matches);
        }
        Some((igd::consts::IGD_CMD, matches)) => match matches.subcommand() {
            Some((igd::consts::IGD_CREATE, matches)) => {
                igd::create::igd_get_create_matches(matches);
            }
            Some((igd::consts::IGD_SEARCH, matches)) => {
                igd::search::igd_get_search_matches(matches);
            }
            _ => unreachable!("IGD Subcommand not found"),
        },
        Some((scoring::consts::FSCORING_CMD, matches)) => {
            scoring::cli::handlers::region_fragment_scoring(matches)?;
        }
        _ => unreachable!("Subcommand not found"),
    };

    Ok(())
}
