use anyhow::Result;
use clap::Command;

// go through the library crate to get the interfaces
use gtars::igd;
use gtars::tokenizers;
use gtars::uniwig;

pub mod consts {
    pub const VERSION: &str = env!("CARGO_PKG_VERSION");
    pub const PKG_NAME: &str = env!("CARGO_PKG_NAME");
    pub const BIN_NAME: &str = env!("CARGO_PKG_NAME");
    pub const UNIWIG_CMD: &str = "uniwig";

}

fn build_parser() -> Command {
    Command::new(consts::BIN_NAME)
        .bin_name(consts::BIN_NAME)
        .version(consts::VERSION)
        .author("Databio")
        .about("Performance critical tools for working with genomic interval data with an emphasis on preprocessing for machine learning pipelines.")
        .subcommand_required(true)
        .subcommand(tokenizers::cli::make_tokenization_cli())
        .subcommand(uniwig::cli::create_uniwig_cli())
        .subcommand(igd::cli::create_igd_cli())
}

fn main() -> Result<()> {
    let app = build_parser();
    let matches = app.get_matches();

    match matches.subcommand() {
        Some((tokenizers::consts::TOKENIZE_CMD, matches)) => {
            tokenizers::cli::handlers::tokenize_bed_file(matches)?;
        }
        Some((uniwig::consts::UNIWIG_CMD, matches)) => {

            uniwig::run_uniwig(matches);
        }
        Some((igd::consts::IGD_CMD, matches)) => {

            match  matches.subcommand() {
                Some((igd::consts::IGD_CREATE, matches)) =>{

                    igd::create::igd_get_create_matches(matches);
                }
                Some((igd::consts::IGD_SEARCH, matches)) =>{

                    igd::search::search_igd(matches);
                }
                _ => unreachable!("IGD Subcommand not found"),
            }
        }

        _ => unreachable!("Subcommand not found"),
    };

    Ok(())
}
