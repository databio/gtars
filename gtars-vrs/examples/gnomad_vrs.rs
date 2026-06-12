//! Compute VRS Allele IDs from a gnomAD VCF using a real reference genome.
//!
//! Usage:
//!   cargo run --release --example gnomad_vrs -- <fasta_path> <vcf_path>

use std::time::Instant;

use gtars_refget::store::{FastaImportOptions, RefgetStore};
use gtars_vrs::vcf::compute_vrs_ids_streaming;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <fasta_path> <vcf_path>", args[0]);
        std::process::exit(1);
    }
    let fasta_path = &args[1];
    let vcf_path = &args[2];

    // Step 1: Load reference genome into RefgetStore
    eprintln!("Loading FASTA: {}", fasta_path);
    let t0 = Instant::now();
    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(fasta_path, FastaImportOptions::new())
        .expect("Failed to load FASTA");
    let load_time = t0.elapsed();
    eprintln!("Loaded in {:.1}s", load_time.as_secs_f64());

    // Print collection info
    let paged = store.list_collections(0, 10, &[]).expect("Failed to list collections");
    let collections = &paged.results;
    assert!(!collections.is_empty(), "No collections found in store");
    let collection_digest = &collections[0].digest;
    eprintln!(
        "Collection: {} ({} sequences)",
        collection_digest, collections[0].n_sequences
    );

    // Step 2: Stream VRS IDs from VCF directly to stdout
    eprintln!("\nProcessing VCF: {}", vcf_path);
    println!("chrom\tpos\tref\talt\tvrs_id");
    let t1 = Instant::now();
    let count = compute_vrs_ids_streaming(&mut store, collection_digest, vcf_path, |r| {
        println!(
            "{}\t{}\t{}\t{}\t{}",
            r.chrom, r.pos, r.ref_allele, r.alt_allele, r.vrs_id
        );
    })
    .expect("Failed to compute VRS IDs");
    let vrs_time = t1.elapsed();

    eprintln!(
        "\nComputed {} VRS IDs in {:.3}s ({:.0} variants/sec)",
        count,
        vrs_time.as_secs_f64(),
        count as f64 / vrs_time.as_secs_f64()
    );
    eprintln!("Total time: {:.1}s", t0.elapsed().as_secs_f64());
}
