//! `rgblast`: search a RefgetStore with BLAST, and report where each query comes from.
fn main() {
    std::process::exit(gtars_rgblast::cli::main_with_verbs(
        "rgblast",
        "Search a RefgetStore with BLAST and report where each query comes from",
        &["search", "report", "origin", "reps"],
    ));
}
