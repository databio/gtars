//! Quick HGVS parsing benchmark
//!
//! Usage (run from the gtars-vrs/ package directory):
//!   cargo bench --bench hgvs_bench

use std::fs;
use std::time::Instant;

use gtars_vrs::hgvs::parse;

fn main() {
    let corpus_path = "tests/fixtures/external_hgvs_corpus.json";
    let contents = fs::read_to_string(corpus_path).expect("read corpus");
    let entries: Vec<serde_json::Value> = serde_json::from_str(&contents).expect("parse json");
    
    let inputs: Vec<&str> = entries
        .iter()
        .map(|e| e["input"].as_str().unwrap())
        .collect();
    
    println!("Corpus size: {} variants", inputs.len());
    
    // Warmup
    for input in &inputs {
        let _ = parse(input);
    }
    
    // Benchmark: run many iterations
    let iterations = 10_000;
    let total_variants = inputs.len() * iterations;
    
    let start = Instant::now();
    for _ in 0..iterations {
        for input in &inputs {
            let _ = parse(input);
        }
    }
    let elapsed = start.elapsed();
    
    let secs = elapsed.as_secs_f64();
    let throughput = total_variants as f64 / secs;
    
    println!("Parsed {} variants in {:.3}s", total_variants, secs);
    println!("Throughput: {:.0} variants/sec", throughput);
}
