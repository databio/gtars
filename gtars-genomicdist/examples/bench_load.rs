use std::env;
use std::io::Read;
use std::time::Instant;

use gtars_genomicdist::asset::GenomicDistAnnotation;
use gtars_genomicdist::partitions::GeneModel;

fn bench_json_gz(path: &str, iterations: u32) {
    // Warm cache
    let _ = std::fs::read(path).unwrap();

    let start = Instant::now();
    for _ in 0..iterations {
        let file = std::fs::File::open(path).unwrap();
        let mut decoder = flate2::read::GzDecoder::new(file);
        let mut json_str = String::new();
        decoder.read_to_string(&mut json_str).unwrap();
        let _parsed: serde_json::Value = serde_json::from_str(&json_str).unwrap();
    }
    let elapsed = start.elapsed();
    let per_iter = elapsed / iterations;
    println!("json.gz decompress+parse: {iterations} iterations in {elapsed:?}");
    println!("  per iteration: {per_iter:?}");
}

fn main() {
    let bin_path = env::args().nth(1).expect("Usage: bench_load <path.bin> [path.json.gz] [path.gtf.gz]");
    let json_gz_path = env::args().nth(2);
    let gtf_path = env::args().nth(3);

    // Warm filesystem cache
    let _ = std::fs::read(&bin_path).unwrap();

    let iterations = 100;
    let start = Instant::now();
    for _ in 0..iterations {
        let _asset = GenomicDistAnnotation::load_bin(&bin_path).unwrap();
    }
    let elapsed = start.elapsed();
    let per_iter = elapsed / iterations;
    println!("load_bin (GDA): {iterations} iterations in {elapsed:?}");
    println!("  per iteration: {per_iter:?}");

    if let Some(json_gz) = json_gz_path {
        println!();
        bench_json_gz(&json_gz, iterations);
    }

    if let Some(gtf) = gtf_path {
        let iters = 5;
        println!();
        let start = Instant::now();
        for _ in 0..iters {
            let _model = GeneModel::from_gtf(&gtf, true, true).unwrap();
        }
        let elapsed = start.elapsed();
        let per_iter = elapsed / iters;
        println!("from_gtf: {iters} iterations in {elapsed:?}");
        println!("  per iteration: {per_iter:?}");
    }
}
