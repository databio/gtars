use clap::{Arg, ArgAction, Command};

pub const REFGET_CMD: &str = "refget";
pub const REFGET_BUILD: &str = "build";
pub const REFGET_EXPORT: &str = "export";

pub fn create_refget_cli() -> Command {
    Command::new(REFGET_CMD)
        .about("Build and manage GA4GH refget sequence stores.")
        .subcommand_required(true)
        .subcommand(
            Command::new(REFGET_BUILD)
                .about("Build a RefgetStore on disk from one or more FASTA files.")
                .arg(
                    Arg::new("fasta")
                        .required(false)
                        .num_args(1..)
                        .help("Path(s), glob pattern(s) (e.g. 'fasta/*.fa.gz'), or directory of FASTA files to import"),
                )
                .arg(
                    Arg::new("file_list")
                        .long("file-list")
                        .short('f')
                        .help("Path to a file listing FASTA paths, one per line (blank lines and #-comments ignored). May contain globs/directories."),
                )
                .arg(
                    Arg::new("output")
                        .long("output")
                        .short('o')
                        .required(true)
                        .help("Output directory for the RefgetStore"),
                )
                .arg(
                    Arg::new("jobs")
                        .long("jobs")
                        .short('j')
                        .value_parser(clap::value_parser!(usize))
                        .default_value("0")
                        .help("Number of FASTA files imported concurrently (0 = auto, 1 = serial). Does not affect input ordering or store output."),
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
                )
                .arg(
                    Arg::new("collection_alias")
                        .long("collection-alias")
                        .value_name("NAMESPACE:ALIAS")
                        .help("Register the imported collection under this collection alias, e.g. 'ucsc:hg38'. This names the collection as a whole and is distinct from per-sequence aliases parsed from FASTA headers. Only valid with a single input FASTA."),
                ),
        )
        .subcommand(
            Command::new(REFGET_EXPORT)
                .about("Export sequences from a RefgetStore to a FASTA file.")
                .arg(
                    Arg::new("store")
                        .long("store")
                        .short('s')
                        .required(true)
                        .help("Path to an existing RefgetStore directory"),
                )
                .arg(
                    Arg::new("output")
                        .long("output")
                        .short('o')
                        .required(true)
                        .help("Output FASTA path ('.gz' extension triggers gzip compression)"),
                )
                .arg(
                    Arg::new("collection")
                        .long("collection")
                        .short('c')
                        .help("Collection digest to export. If omitted and the store holds exactly one collection, that collection is used; otherwise the available digests are listed."),
                )
                .arg(
                    Arg::new("names")
                        .long("names")
                        .num_args(1..)
                        .help("Restrict export to these sequence names (default: all sequences in the collection)"),
                )
                .arg(
                    Arg::new("line_width")
                        .long("line-width")
                        .short('w')
                        .value_parser(clap::value_parser!(usize))
                        .default_value("80")
                        .help("Sequence line width. 0 = unwrapped, one sequence per line (what GGCAT/SSHash expect)"),
                ),
        )
}
