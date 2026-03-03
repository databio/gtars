use clap::{Arg, Command};

pub const PREP_CMD: &str = "prep";

pub fn create_prep_cli() -> Command {
    Command::new(PREP_CMD)
        .about("Pre-serialize GTF gene models or signal matrices to bincode for fast loading.")
        .arg(
            Arg::new("gtf")
                .long("gtf")
                .required(false)
                .help("Path to GTF/GTF.gz gene model to serialize"),
        )
        .arg(
            Arg::new("signal-matrix")
                .long("signal-matrix")
                .required(false)
                .help("Path to signal matrix TSV/TSV.gz to serialize"),
        )
        .arg(
            Arg::new("output")
                .long("output")
                .short('o')
                .required(false)
                .help("Output path (default: input path with .bin extension, stripping .gz first)"),
        )
}
