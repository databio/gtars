mod uniwig;

use anyhow::Result;
use clap::Command;

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
        .subcommand(uniwig::cli::create_uniwig_cli())
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

        _ => unreachable!("Subcommand not found"),
    };

    Ok(())
}
