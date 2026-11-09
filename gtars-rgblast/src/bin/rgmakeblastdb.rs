//! `rgmakeblastdb`: BLAST nucleotide databases straight from a RefgetStore.
fn main() {
    std::process::exit(gtars_rgblast::cli::main_with_verbs(
        "rgmakeblastdb",
        "BLAST nucleotide databases built straight from a RefgetStore",
        &["build", "verify", "view", "nstats"],
    ));
}
