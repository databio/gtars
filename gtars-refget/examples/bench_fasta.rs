//! Benchmark: load a FASTA file into a RefgetStore and report peak RSS.
//! Usage: cargo run --release --example bench_fasta -- <fasta_path>

use std::time::Instant;

fn peak_rss_mb() -> f64 {
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    for line in status.lines() {
        if line.starts_with("VmHWM:") {
            let kb: f64 = line
                .split_whitespace()
                .nth(1)
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0);
            return kb / 1024.0;
        }
    }
    0.0
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: bench_fasta <fasta_path>");
        std::process::exit(1);
    }
    let fasta_path = &args[1];

    println!("Loading: {}", fasta_path);
    let start = Instant::now();

    let tmp = tempfile::tempdir()?;
    let mut store = gtars_refget::store::RefgetStore::on_disk(tmp.path())?;
    store.set_quiet(true);

    let opts = gtars_refget::store::FastaImportOptions::default();
    let (metadata, _) = store.add_sequence_collection_from_fasta(fasta_path, opts)?;

    let elapsed = start.elapsed();
    let peak = peak_rss_mb();

    println!("Digest:   {}", metadata.digest);
    println!("Seqs:     {}", metadata.n_sequences);
    println!("Time:     {:.1}s", elapsed.as_secs_f64());
    println!("Peak RSS: {:.0} MB", peak);

    Ok(())
}
