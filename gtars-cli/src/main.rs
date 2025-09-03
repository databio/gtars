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
            let _ = overlaprs::handlers::overlap_query_with_universe(matches);
        }

        _ => unreachable!("Subcommand not found"),
    };

    Ok(())
}
