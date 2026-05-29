use clap::{Arg, ArgAction, Command};

pub const REFGET_CMD: &str = "refget";
pub const REFGET_BUILD: &str = "build";

pub fn create_refget_cli() -> Command {
    Command::new(REFGET_CMD)
        .about("Build and manage GA4GH refget sequence stores.")
        .subcommand_required(true)
        .subcommand(
            Command::new(REFGET_BUILD)
                .about("Build a RefgetStore on disk from one or more FASTA files.")
                .arg(
                    Arg::new("fasta")
                        .required(true)
                        .num_args(1..)
                        .help("Path(s) to FASTA file(s) (.fa or .fa.gz) to import"),
                )
                .arg(
                    Arg::new("output")
                        .long("output")
                        .short('o')
                        .required(true)
                        .help("Output directory for the RefgetStore"),
                )
                .arg(
                    Arg::new("threads")
                        .long("threads")
                        .short('t')
                        .value_parser(clap::value_parser!(usize))
                        .default_value("0")
                        .help("Total worker threads for the digest+encode stage (0 = auto, 1 = serial)"),
                )
                .arg(
                    Arg::new("file-jobs")
                        .long("file-jobs")
                        .short('j')
                        .value_parser(clap::value_parser!(usize))
                        .default_value("0")
                        .help("Input FASTA files decoded concurrently, parallelizing gzip decode across files (0 = auto)"),
                )
                .arg(
                    Arg::new("raw")
                        .long("raw")
                        .action(ArgAction::SetTrue)
                        .help("Use Raw storage mode instead of the default Encoded (2-bit) mode"),
                )
                .arg(
                    Arg::new("force")
                        .long("force")
                        .action(ArgAction::SetTrue)
                        .help("Overwrite existing collections/sequences in the store"),
                ),
        )
}
