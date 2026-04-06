use clap::{Arg, Command, arg};

pub const OVERLAP_CMD: &str = "overlaprs";

pub fn create_overlap_cli() -> Command {
    Command::new(OVERLAP_CMD)
        .author("NJL")
        .about("Tokenize a BED file against a universe of regions (overlap-based encoding).")
        .arg_required_else_help(true)
        .arg(
            Arg::new("query")
                .short('q')
                .long("query")
                .required(true)
                .help("Path to the BED file to tokenize"),
        )
        .arg(
            Arg::new("universe")
                .short('u')
                .long("universe")
                .required(true)
                .help("Path to the universe BED file to tokenize against"),
        )
        .arg(arg!(-e --backend <backend> "Which backend to use (ailist or bits)"))
        .arg(arg!(--streaming "Use streaming mode for very large universes (lower memory usage)"))
}
