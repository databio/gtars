//! `rgblastn`: nucleotide BLAST (NCBI BLASTN 2.17.0+ algorithms) that searches a RefgetStore
//! directly, with no database files.
fn main() {
    std::process::exit(gtars_rgblast::cli::rgblastn_main());
}
