use clap::{Arg, Command, arg};

pub const RANGES_CMD: &str = "ranges";

pub fn create_ranges_cli() -> Command {
    Command::new(RANGES_CMD)
        .about("Interval set algebra operations on BED files.")
        .subcommand_required(true)
        .subcommand(
            Command::new("reduce")
                .about("Merge overlapping and adjacent intervals.")
                .arg(arg!(--input <BED> "Input BED file").required(true))
                .arg(arg!(--output <OUTPUT> "Output BED file (default: stdout)").required(false)),
        )
        .subcommand(
            Command::new("trim")
                .about("Trim regions to chromosome boundaries.")
                .arg(arg!(--input <BED> "Input BED file").required(true))
                .arg(
                    Arg::new("chrom-sizes")
                        .long("chrom-sizes")
                        .required(true)
                        .help("Path to chrom.sizes file"),
                )
                .arg(arg!(--output <OUTPUT> "Output BED file (default: stdout)").required(false)),
        )
        .subcommand(
            Command::new("promoters")
                .about("Generate promoter regions from region starts.")
                .arg(arg!(--input <BED> "Input BED file").required(true))
                .arg(
                    arg!(--upstream <UPSTREAM> "Bases upstream of start")
                        .required(false)
                        .default_value("2000"),
                )
                .arg(
                    arg!(--downstream <DOWNSTREAM> "Bases downstream of start")
                        .required(false)
                        .default_value("200"),
                )
                .arg(arg!(--output <OUTPUT> "Output BED file (default: stdout)").required(false)),
        )
        .subcommand(
            Command::new("setdiff")
                .about("Subtract regions in B from regions in A.")
                .arg(arg!(-a <BED_A> "Input BED file A").required(true))
                .arg(arg!(-b <BED_B> "Input BED file B to subtract").required(true))
                .arg(arg!(--output <OUTPUT> "Output BED file (default: stdout)").required(false)),
        )
        .subcommand(
            Command::new("pintersect")
                .about("Pairwise intersection by index position.")
                .arg(arg!(-a <BED_A> "Input BED file A").required(true))
                .arg(arg!(-b <BED_B> "Input BED file B").required(true))
                .arg(arg!(--output <OUTPUT> "Output BED file (default: stdout)").required(false)),
        )
        .subcommand(
            Command::new("concat")
                .about("Concatenate two region sets without merging.")
                .arg(arg!(-a <BED_A> "Input BED file A").required(true))
                .arg(arg!(-b <BED_B> "Input BED file B").required(true))
                .arg(arg!(--output <OUTPUT> "Output BED file (default: stdout)").required(false)),
        )
        .subcommand(
            Command::new("union")
                .about("Merge two region sets into a minimal non-overlapping result.")
                .arg(arg!(-a <BED_A> "Input BED file A").required(true))
                .arg(arg!(-b <BED_B> "Input BED file B").required(true))
                .arg(arg!(--output <OUTPUT> "Output BED file (default: stdout)").required(false)),
        )
        .subcommand(
            Command::new("jaccard")
                .about("Compute nucleotide-level Jaccard similarity between two BED files.")
                .arg(arg!(-a <BED_A> "Input BED file A").required(true))
                .arg(arg!(-b <BED_B> "Input BED file B").required(true)),
        )
}
