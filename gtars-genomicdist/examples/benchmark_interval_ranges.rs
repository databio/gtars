//! Benchmark for interval range operations.
//!
//! Usage:
//!   cargo run -p gtars-genomicdist --example benchmark_interval_ranges -- --scale 1000
//!   cargo run -p gtars-genomicdist --example benchmark_interval_ranges -- --scale 1000 --output results.csv
//!   cargo run -p gtars-genomicdist --example benchmark_interval_ranges -- --input path.bed
//!   cargo run -p gtars-genomicdist --release --example benchmark_interval_ranges -- --scale 1000000

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::Write;
use std::time::Instant;

use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::IntervalRanges;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

const CHROMS: &[&str] = &[
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22",
];

const CHROM_SIZES: &[(& str, u32)] = &[
    ("chr1", 248956422),
    ("chr2", 242193529),
    ("chr3", 198295559),
    ("chr4", 190214555),
    ("chr5", 181538259),
    ("chr6", 170805979),
    ("chr7", 159345973),
    ("chr8", 145138636),
    ("chr9", 138394717),
    ("chr10", 133797422),
    ("chr11", 135086622),
    ("chr12", 133275309),
    ("chr13", 114364328),
    ("chr14", 107043718),
    ("chr15", 101991189),
    ("chr16", 90338345),
    ("chr17", 83257441),
    ("chr18", 80373285),
    ("chr19", 58617616),
    ("chr20", 64444167),
    ("chr21", 46709983),
    ("chr22", 50818468),
];

fn generate_synthetic_regionset(n: usize, seed: u64) -> RegionSet {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut regions = Vec::with_capacity(n);

    for _ in 0..n {
        let chr_idx = rng.gen_range(0..CHROMS.len());
        let chr = CHROMS[chr_idx];
        let chrom_size = CHROM_SIZES[chr_idx].1;

        // Lognormal-ish length: exp(uniform(4.6, 9.2)) ~= 100 to 10000
        let log_len: f64 = rng.gen_range(4.6..9.2);
        let length = (log_len.exp() as u32).min(10000).max(100);

        let max_start = chrom_size.saturating_sub(length);
        let start = rng.gen_range(0..=max_start);
        let end = start + length;

        regions.push(Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        });
    }

    let mut rs = RegionSet::from(regions);
    rs.sort();
    rs
}

fn chrom_sizes_map() -> HashMap<String, u32> {
    CHROM_SIZES
        .iter()
        .map(|(k, v)| (k.to_string(), *v))
        .collect()
}

struct TimingResult {
    function: String,
    scale: String,
    iterations: u32,
    min_us: f64,
    median_us: f64,
    mean_us: f64,
    max_us: f64,
}

fn benchmark_fn<F: Fn()>(f: F, iterations: u32) -> Vec<f64> {
    let mut times = Vec::with_capacity(iterations as usize);
    for _ in 0..iterations {
        let start = Instant::now();
        f();
        let elapsed = start.elapsed();
        times.push(elapsed.as_micros() as f64);
    }
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    times
}

fn stats_from_times(times: &[f64]) -> (f64, f64, f64, f64) {
    let min = times[0];
    let max = times[times.len() - 1];
    let mean = times.iter().sum::<f64>() / times.len() as f64;
    let median = if times.len() % 2 == 0 {
        (times[times.len() / 2 - 1] + times[times.len() / 2]) / 2.0
    } else {
        times[times.len() / 2]
    };
    (min, median, mean, max)
}

fn run_benchmarks(rs_a: &RegionSet, rs_b: &RegionSet, scale_label: &str, iterations: u32) -> Vec<TimingResult> {
    let chrom_sizes = chrom_sizes_map();
    let mut results = Vec::new();

    // trim
    let times = benchmark_fn(|| { rs_a.trim(&chrom_sizes); }, iterations);
    let (min, median, mean, max) = stats_from_times(&times);
    results.push(TimingResult {
        function: "trim".to_string(),
        scale: scale_label.to_string(),
        iterations,
        min_us: min, median_us: median, mean_us: mean, max_us: max,
    });

    // promoters
    let times = benchmark_fn(|| { rs_a.promoters(2000, 200); }, iterations);
    let (min, median, mean, max) = stats_from_times(&times);
    results.push(TimingResult {
        function: "promoters".to_string(),
        scale: scale_label.to_string(),
        iterations,
        min_us: min, median_us: median, mean_us: mean, max_us: max,
    });

    // reduce
    let times = benchmark_fn(|| { rs_a.reduce(); }, iterations);
    let (min, median, mean, max) = stats_from_times(&times);
    results.push(TimingResult {
        function: "reduce".to_string(),
        scale: scale_label.to_string(),
        iterations,
        min_us: min, median_us: median, mean_us: mean, max_us: max,
    });

    // setdiff
    let times = benchmark_fn(|| { rs_a.setdiff(rs_b); }, iterations);
    let (min, median, mean, max) = stats_from_times(&times);
    results.push(TimingResult {
        function: "setdiff".to_string(),
        scale: scale_label.to_string(),
        iterations,
        min_us: min, median_us: median, mean_us: mean, max_us: max,
    });

    // pintersect
    let times = benchmark_fn(|| { rs_a.pintersect(rs_b); }, iterations);
    let (min, median, mean, max) = stats_from_times(&times);
    results.push(TimingResult {
        function: "pintersect".to_string(),
        scale: scale_label.to_string(),
        iterations,
        min_us: min, median_us: median, mean_us: mean, max_us: max,
    });

    results
}

fn write_output_beds(rs_a: &RegionSet, rs_b: &RegionSet, scale_label: &str, output_dir: &str) {
    let chrom_sizes = chrom_sizes_map();
    std::fs::create_dir_all(output_dir).unwrap();

    // Write input files for cross-language validation
    rs_a.to_bed(format!("{}/input_a_{}.bed", output_dir, scale_label)).unwrap();
    rs_b.to_bed(format!("{}/input_b_{}.bed", output_dir, scale_label)).unwrap();

    rs_a.trim(&chrom_sizes)
        .to_bed(format!("{}/rust_trim_{}.bed", output_dir, scale_label))
        .unwrap();
    rs_a.promoters(2000, 200)
        .to_bed(format!("{}/rust_promoters_{}.bed", output_dir, scale_label))
        .unwrap();
    rs_a.reduce()
        .to_bed(format!("{}/rust_reduce_{}.bed", output_dir, scale_label))
        .unwrap();
    rs_a.setdiff(rs_b)
        .to_bed(format!("{}/rust_setdiff_{}.bed", output_dir, scale_label))
        .unwrap();
    rs_a.pintersect(rs_b)
        .to_bed(format!("{}/rust_pintersect_{}.bed", output_dir, scale_label))
        .unwrap();
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut scale: Option<usize> = None;
    let mut input_path: Option<String> = None;
    let mut output_csv: Option<String> = None;
    let mut output_beds_dir: Option<String> = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--scale" => {
                i += 1;
                scale = Some(args[i].parse().expect("Invalid scale value"));
            }
            "--input" => {
                i += 1;
                input_path = Some(args[i].clone());
            }
            "--output" => {
                i += 1;
                output_csv = Some(args[i].clone());
            }
            "--output-beds" => {
                i += 1;
                output_beds_dir = Some(args[i].clone());
            }
            _ => {
                eprintln!("Unknown argument: {}", args[i]);
                std::process::exit(1);
            }
        }
        i += 1;
    }

    let mut all_results: Vec<TimingResult> = Vec::new();

    if let Some(n) = scale {
        let iterations = if n >= 100_000 { 10 } else { 100 };
        let label = format!("{}K", n / 1000);
        eprintln!("Generating synthetic data: {} regions (seed=42,43)...", n);
        let rs_a = generate_synthetic_regionset(n, 42);
        let rs_b = generate_synthetic_regionset(n, 43);
        eprintln!("Running benchmarks ({} iterations)...", iterations);
        let results = run_benchmarks(&rs_a, &rs_b, &label, iterations);
        for r in &results {
            eprintln!(
                "  {:<12} min={:>10.0}us  median={:>10.0}us  mean={:>10.0}us  max={:>10.0}us",
                r.function, r.min_us, r.median_us, r.mean_us, r.max_us
            );
        }
        if let Some(ref dir) = output_beds_dir {
            eprintln!("Writing output BED files to {}...", dir);
            write_output_beds(&rs_a, &rs_b, &label, dir);
        }
        all_results.extend(results);
    } else if let Some(ref path) = input_path {
        eprintln!("Loading BED file: {}...", path);
        let rs_a = RegionSet::try_from(path.as_str()).expect("Failed to load BED file");
        let n = rs_a.regions.len();
        // Generate a second set of same size for pairwise ops
        let rs_b = generate_synthetic_regionset(n, 43);
        let iterations = if n >= 100_000 { 10 } else { 100 };
        let label = format!("real_{}K", n / 1000);
        eprintln!("Running benchmarks on {} regions ({} iterations)...", n, iterations);
        let results = run_benchmarks(&rs_a, &rs_b, &label, iterations);
        for r in &results {
            eprintln!(
                "  {:<12} min={:>10.0}us  median={:>10.0}us  mean={:>10.0}us  max={:>10.0}us",
                r.function, r.min_us, r.median_us, r.mean_us, r.max_us
            );
        }
        if let Some(ref dir) = output_beds_dir {
            eprintln!("Writing output BED files to {}...", dir);
            write_output_beds(&rs_a, &rs_b, &label, dir);
        }
        all_results.extend(results);
    } else {
        // Default: run all 4 scales
        for &n in &[1_000, 10_000, 100_000, 1_000_000] {
            let iterations = if n >= 100_000 { 10 } else { 100 };
            let label = format!("{}K", n / 1000);
            eprintln!("=== Scale: {} ({} regions) ===", label, n);
            eprintln!("Generating synthetic data (seed=42,43)...");
            let rs_a = generate_synthetic_regionset(n, 42);
            let rs_b = generate_synthetic_regionset(n, 43);
            eprintln!("Running benchmarks ({} iterations)...", iterations);
            let results = run_benchmarks(&rs_a, &rs_b, &label, iterations);
            for r in &results {
                eprintln!(
                    "  {:<12} min={:>10.0}us  median={:>10.0}us  mean={:>10.0}us  max={:>10.0}us",
                    r.function, r.min_us, r.median_us, r.mean_us, r.max_us
                );
            }
            if let Some(ref dir) = output_beds_dir {
                eprintln!("Writing output BED files to {}...", dir);
                write_output_beds(&rs_a, &rs_b, &label, dir);
            }
            all_results.extend(results);
        }
    }

    // Write CSV output
    if let Some(ref csv_path) = output_csv {
        let mut f = File::create(csv_path).expect("Failed to create output CSV");
        writeln!(f, "language,function,scale,iterations,min_us,median_us,mean_us,max_us").unwrap();
        for r in &all_results {
            writeln!(
                f,
                "rust,{},{},{},{:.0},{:.0},{:.0},{:.0}",
                r.function, r.scale, r.iterations, r.min_us, r.median_us, r.mean_us, r.max_us
            )
            .unwrap();
        }
        eprintln!("Results written to {}", csv_path);
    } else {
        // Print CSV to stdout
        println!("language,function,scale,iterations,min_us,median_us,mean_us,max_us");
        for r in &all_results {
            println!(
                "rust,{},{},{},{:.0},{:.0},{:.0},{:.0}",
                r.function, r.scale, r.iterations, r.min_us, r.median_us, r.mean_us, r.max_us
            );
        }
    }
}
