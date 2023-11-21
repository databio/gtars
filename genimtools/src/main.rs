use clap::Command;

// go through the library crate to get the interfaces
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
        .subcommand(
            vocab::cli::make_prune_cli()
        )
    // .subcommand(
    //     uniwig::cli::interfaces::make_uniwig_cli()
    // )
}

fn main() {
    let app = build_parser();
    let matches = app.get_matches();

    match matches.subcommand() {
        Some((vocab::consts::PRUNE_CMD, matches)) => {
            vocab::cli::handlers::prune_universe(matches);
        }
        Some((consts::UNIWIG_CMD, _matches)) => {
            println!("uniwig");
        }

        _ => unreachable!("Subcommand not found"),
    };
}
