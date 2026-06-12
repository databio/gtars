use clap::{Arg, Command};

pub const PREP_CMD: &str = "prep";

pub fn create_prep_cli() -> Command {
    Command::new(PREP_CMD)
        .about("Pre-serialize GTF gene models, signal matrices, or FASTA files to binary for fast loading.")
        .arg(
            Arg::new("gtf")
                .long("gtf")
                .required(false)
                .help("Path to GTF/GTF.gz gene model to serialize to GDA binary"),
        )
        .arg(
            Arg::new("signal-matrix")
                .long("signal-matrix")
                .required(false)
                .help("Path to signal matrix TSV/TSV.gz to serialize"),
        )
        .arg(
            Arg::new("fasta")
                .long("fasta")
                .required(false)
                .help("Path to FASTA file to convert to .fab binary (zero-copy mmap format)"),
        )
        .arg(
            Arg::new("output")
                .long("output")
                .short('o')
                .required(false)
                .help("Output path (default: input path with .bin/.fab extension)"),
        )
}
