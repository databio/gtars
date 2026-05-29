use std::time::Instant;

use anyhow::Result;
use clap::ArgMatches;

use gtars_refget::store::{FastaImportOptions, RefgetStore, StorageMode};

pub fn run_refget(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some((super::cli::REFGET_BUILD, sub)) => run_build(sub),
        _ => unreachable!("refget subcommand not found"),
    }
}

fn run_build(matches: &ArgMatches) -> Result<()> {
    let fastas: Vec<&String> = matches
        .get_many::<String>("fasta")
        .expect("fasta is required")
        .collect();
    let output = matches
        .get_one::<String>("output")
        .expect("output is required");
    let threads = *matches.get_one::<usize>("threads").unwrap_or(&0);
    let raw = matches.get_flag("raw");
    let force = matches.get_flag("force");

    let mut store = RefgetStore::on_disk(output)
        .map_err(|e| anyhow::anyhow!("Failed to create store at {}: {}", output, e))?;
    if raw {
        store.set_encoding_mode(StorageMode::Raw);
    } else {
        store.set_encoding_mode(StorageMode::Encoded);
    }

    let mode = if raw { "Raw" } else { "Encoded" };
    eprintln!(
        "Building RefgetStore at {} (mode={}, threads={})",
        output,
        mode,
        if threads == 0 {
            "auto".to_string()
        } else {
            threads.to_string()
        }
    );

    let mut total_bases: u64 = 0;
    let mut total_seqs: usize = 0;
    let start = Instant::now();

    for fa in &fastas {
        let opts = FastaImportOptions::new().force(force).threads(threads);
        let (metadata, was_new) = store
            .add_sequence_collection_from_fasta(fa.as_str(), opts)
            .map_err(|e| anyhow::anyhow!("Failed to import {}: {}", fa, e))?;
        total_seqs += metadata.n_sequences;
        eprintln!(
            "  {} {}: {} ({} sequences)",
            if was_new { "added" } else { "skipped" },
            fa,
            metadata.digest,
            metadata.n_sequences
        );
    }

    store
        .write()
        .map_err(|e| anyhow::anyhow!("Failed to write store: {}", e))?;

    let elapsed = start.elapsed().as_secs_f64();

    // Best-effort base count from the loaded sequence index (for throughput).
    for meta in store.list_sequences() {
        total_bases += meta.length as u64;
    }

    let mbps = if elapsed > 0.0 {
        (total_bases as f64 / 1_000_000.0) / elapsed
    } else {
        0.0
    };
    eprintln!(
        "Done: {} sequences, {} bases in {:.3}s ({:.1} Mbase/s, threads={})",
        total_seqs,
        total_bases,
        elapsed,
        mbps,
        if threads == 0 {
            "auto".to_string()
        } else {
            threads.to_string()
        }
    );

    Ok(())
}
