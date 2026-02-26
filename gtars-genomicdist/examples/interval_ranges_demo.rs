//! Demo: load a BED file and run the 5 IntervalRanges operations.
//!
//! Usage:
//!   cargo run -p gtars-genomicdist --example interval_ranges_demo -- <bed_file> [bed_file_b]
//!
//! With one file:  runs trim, promoters, reduce (prints first 20 lines each)
//! With two files: also runs setdiff(a, b) and pintersect(a, b)

use std::collections::HashMap;
use std::env;

use gtars_core::models::RegionSet;
use gtars_genomicdist::interval_ranges::IntervalRanges;

fn print_regionset(label: &str, rs: &RegionSet, max: usize) {
    println!("\n=== {} ({} regions) ===", label, rs.regions.len());
    for r in rs.regions.iter().take(max) {
        println!("{}\t{}\t{}", r.chr, r.start, r.end);
    }
    if rs.regions.len() > max {
        println!("... ({} more)", rs.regions.len() - max);
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: interval_ranges_demo <bed_file> [bed_file_b]");
        std::process::exit(1);
    }

    let a = RegionSet::try_from(args[1].as_str()).expect("Failed to load BED file A");
    println!("Loaded A: {} regions from {}", a.regions.len(), &args[1]);

    // -- trim (use a simple chrom_sizes: hg38 main chroms) --
    let chrom_sizes: HashMap<String, u32> = [
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
        ("chrX", 156040895),
        ("chrY", 57227415),
        ("chrM", 16569),
    ]
    .iter()
    .map(|(k, v)| (k.to_string(), *v))
    .collect();

    print_regionset("trim(A)", &a.trim(&chrom_sizes), 20);
    print_regionset("promoters(A, 2000, 200)", &a.promoters(2000, 200), 20);
    print_regionset("reduce(A)", &a.reduce(), 20);

    // -- two-file operations --
    if args.len() >= 3 {
        let b = RegionSet::try_from(args[2].as_str()).expect("Failed to load BED file B");
        println!("\nLoaded B: {} regions from {}", b.regions.len(), &args[2]);

        print_regionset("setdiff(A, B)", &a.setdiff(&b), 20);
        print_regionset("pintersect(A, B)", &a.pintersect(&b), 20);
    } else {
        println!("\n(Pass a second BED file to also run setdiff and pintersect)");
    }
}
