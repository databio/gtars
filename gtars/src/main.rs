use anyhow::Result;
use clap::Command;
// go through the library crate to get the interfaces
use gtars::overlap;
// use gtars::bbcache;
use gtars::fragsplit;
use gtars::igd;
use gtars::scoring;
use gtars::uniwig;
use gtars::tokenizers;

pub mod consts {
    pub const VERSION: &str = env!("CARGO_PKG_VERSION");
    pub const PKG_NAME: &str = env!("CARGO_PKG_NAME");
    pub const BIN_NAME: &str = env!("CARGO_PKG_NAME");
}

fn build_parser() -> Command {
    Command::new(consts::BIN_NAME)
        .bin_name(consts::BIN_NAME)
        .version(consts::VERSION)
        .author("Databio")
        .about("Performance critical tools for working with genomic interval data with an emphasis on preprocessing for machine learning pipelines.")
        .subcommand_required(true)
        .subcommand(overlap::cli::create_overlap_cli())
        .subcommand(fragsplit::cli::make_fragsplit_cli())
        .subcommand(uniwig::cli::create_uniwig_cli())
        .subcommand(igd::cli::create_igd_cli())
        .subcommand(scoring::cli::make_fscoring_cli())
        // .subcommand(bbcache::cli::create_bbcache_cli())
        .subcommand(tokenizers::cli::create_tokenizer_cli())
}

fn main() -> Result<()> {
    let app = build_parser();
    let matches = app.get_matches();

    match matches.subcommand() {
        //
        // OVERLAP
        //
        Some((overlap::consts::OVERLAP_CMD, matches)) => {
            overlap::cli::handlers::overlap_query_with_universe(matches)?;
        }

        //
        // FRAGMENT SPLITTING UTIL
        //
        Some((fragsplit::consts::FRAGSPLIT_CMD, matches)) => {
            fragsplit::cli::handlers::split_fragment_files(matches)?;
        }

        //
        // UNIWIG
        //
        Some((uniwig::consts::UNIWIG_CMD, matches)) => {
            uniwig::run_uniwig(matches);
        }

        //
        // BBCACHE
        //
        // Some((bbcache::consts::BBCACHE_CMD, sub_m)) => {
        //     if let Some((subcmd, sub_matches)) = sub_m.subcommand() {
        //         bbcache::run_bbcache(subcmd, sub_matches);
        //     }
        // }

        //
        // IGD
        //
        Some((igd::consts::IGD_CMD, matches)) => match matches.subcommand() {
            Some((igd::consts::IGD_CREATE, matches)) => {
                igd::create::igd_get_create_matches(matches);
            }
            Some((igd::consts::IGD_SEARCH, matches)) => {
                igd::search::igd_get_search_matches(matches);
            }
            _ => unreachable!("IGD Subcommand not found"),
        },

        //
        // REGION SCORING MATRIX
        //
        Some((scoring::consts::FSCORING_CMD, matches)) => {
            scoring::cli::handlers::region_fragment_scoring(matches)?;
        }

        //
        // INTERVAL TOKENIZERS
        //
        Some((tokenizers::consts::TOKENIZERS_CMD, matches)) => {
            tokenizers::cli::handlers::tokenize_query_into_universe(matches)?;
        }

        _ => unreachable!("Subcommand not found"),
    };

    Ok(())
}
