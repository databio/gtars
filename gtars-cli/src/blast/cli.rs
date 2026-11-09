use clap::Subcommand;

pub const BLAST_CMD: &str = "blast";

/// `gtars blast <verb>`: the same parsers as the `rgmakeblastdb` and `rgblast` binaries.
pub fn create_blast_cli() -> clap::Command {
    gtars_rgblast::cli::BlastCommand::augment_subcommands(
        clap::Command::new(BLAST_CMD)
            .about("BLAST over a RefgetStore: build databases, search, and report where hits come from")
            .long_about(
                "BLAST over a RefgetStore: build databases, search, and report where hits come from.\n\
                 The verbs and flags are those of the rgmakeblastdb (build, verify, view, nstats)\n\
                 and rgblast (search, report, origin, reps) binaries.",
            )
            .subcommand_required(true)
            .arg_required_else_help(true),
    )
}
