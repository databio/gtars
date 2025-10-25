use clap::{Arg, Command, arg};

pub use gtars_scoring::consts::*;

pub fn create_scoring_cli() -> Command {
    Command::new(FSCORING_CMD)
        .author("Nathan LeRoy")
        .about("Create a scoring matrix for a set of fragment files over a consensus peak set.")
        .arg(Arg::new("fragments"))
        .arg(Arg::new("consensus"))
        .arg(arg!(--mode <mode>))
        .arg(arg!(--output <output>))
        .arg(
            arg!(--barcode)
                .help("Enable barcode-based scoring for single-cell data (outputs sparse Matrix Market format)")
                .action(clap::ArgAction::SetTrue)
        )
}
