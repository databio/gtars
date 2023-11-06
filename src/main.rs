use clap::Command;

pub mod cli;

pub mod consts {
    pub const VERSION: &str = env!("CARGO_PKG_VERSION");
    pub const BIN_NAME: &str = "genim";
    pub const PRUNE_CMD: &str = "prune";
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
            cli::interfaces::make_prune_cli()
        )
        .subcommand(
            cli::interfaces::make_uniwig_cli()
        )
}

fn main() {
    let app = build_parser();
    let matches = app.get_matches();

    match matches.subcommand() {
        Some((consts::PRUNE_CMD, matches)) => {
            cli::functions::prune_universe(matches);
        }
        Some((consts::UNIWIG_CMD, matches)) => {
            cli::functions::uniwig(matches);
        }

        _ => unreachable!("Subcommand not found"),
    };
}
