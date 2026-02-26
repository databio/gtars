use clap::{Arg, Command};

pub const CONSENSUS_CMD: &str = "consensus";

pub fn create_consensus_cli() -> Command {
    Command::new(CONSENSUS_CMD)
        .about("Compute consensus regions across multiple BED files. Outputs BED4 (chr, start, end, count).")
        .arg(
            Arg::new("beds")
                .long("beds")
                .required(true)
                .num_args(2..)
                .help("Two or more input BED files"),
        )
        .arg(
            Arg::new("min-count")
                .long("min-count")
                .required(false)
                .default_value("1")
                .help("Minimum overlap count to include a region in output"),
        )
        .arg(
            Arg::new("output")
                .long("output")
                .required(false)
                .help("Output BED file (default: stdout)"),
        )
}
