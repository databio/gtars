use anyhow::Result;
use clap::Command;

// go through the library crate to get the interfaces
use genimtools::tokenizers;
use genimtools::tools;
use genimtools::vocab;
// use genimtools::uniwig;

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
        .subcommand(vocab::cli::make_prune_cli())
        .subcommand(tokenizers::cli::make_tokenization_cli())
        .subcommand(tools::cli::make_tools_cli())
}

fn main() -> Result<()> {
    let app = build_parser();
    let matches = app.get_matches();

    match matches.subcommand() {
        Some((vocab::consts::PRUNE_CMD, matches)) => {
            vocab::cli::handlers::prune_universe(matches);
        }
        Some((tokenizers::consts::TOKENIZE_CMD, matches)) => {
            tokenizers::cli::handlers::tokenize_bed_file(matches)?;
        }
        Some((tools::consts::TOOLS_CMD, matches)) => {
            let _ = tools::cli::handlers::tools_handler(matches);
        }

        _ => unreachable!("Subcommand not found"),
    };

    Ok(())
}
