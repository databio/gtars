//! Full HGVS-to-VRS pipeline benchmark (parsing + normalization + digest)

use std::io::Write;
use std::time::Instant;

use gtars_refget::store::{FastaImportOptions, RefgetStore};
use gtars_vrs::hgvs::bridge::hgvs_str_to_vrs_id;
use gtars_vrs::NoTranscriptProvider;

fn main() {
    // Create temp FASTA with a synthetic sequence
    let tmpdir = tempfile::tempdir().expect("create tmpdir");
    let fasta_path = tmpdir.path().join("bench.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).expect("create fasta");
        writeln!(f, ">NC_000001.11").unwrap();
        let seq = "ACGTACGTACGT".repeat(1000); // 12kb
        writeln!(f, "{}", seq).unwrap();
    }

    // Create RefgetStore from FASTA
    let mut store = RefgetStore::in_memory();
    store.disable_encoding();
    store.set_quiet(true);
    store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
        .expect("import fasta");

    let collection_digest = store
        .iter_collections()
        .next()
        .map(|c| c.metadata.digest.clone())
        .expect("collection");

    // Create test variants (genomic only - no transcript lookup needed)
    let test_cases: Vec<String> = (1..=1000)
        .map(|i| format!("NC_000001.11:g.{}A>T", i * 10))
        .collect();

    let provider = NoTranscriptProvider;

    println!("Test cases: {}", test_cases.len());

    // Warmup
    for case in &test_cases {
        let _ = hgvs_str_to_vrs_id(case, &provider, &mut store, &collection_digest, &mut Vec::new());
    }

    // Benchmark
    let iterations = 1000;
    let total = test_cases.len() * iterations;

    let start = Instant::now();
    for _ in 0..iterations {
        for case in &test_cases {
            let _ = hgvs_str_to_vrs_id(case, &provider, &mut store, &collection_digest, &mut Vec::new());
        }
    }
    let elapsed = start.elapsed();

    let secs = elapsed.as_secs_f64();
    let throughput = total as f64 / secs;

    println!("Processed {} variants in {:.3}s", total, secs);
    println!("Full VRS throughput: {:.0} variants/sec", throughput);
}
