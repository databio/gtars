use std::time::Instant;

use anyhow::Result;
use clap::ArgMatches;

use gtars_refget::store::{FastaImportOptions, RefgetStore, StorageMode};
use gtars_refget::{expand_fasta_inputs, FastaInputs};

pub fn run_refget(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some((super::cli::REFGET_BUILD, sub)) => run_build(sub),
        _ => unreachable!("refget subcommand not found"),
    }
}

fn run_build(matches: &ArgMatches) -> Result<()> {
    let paths: Vec<std::path::PathBuf> = matches
        .get_many::<String>("fasta")
        .into_iter()
        .flatten()
        .map(std::path::PathBuf::from)
        .collect();
    let file_list = matches
        .get_one::<String>("file_list")
        .map(std::path::PathBuf::from);
    let inputs = FastaInputs { paths, file_list };
    let fastas = expand_fasta_inputs(&inputs)
        .map_err(|e| anyhow::anyhow!("Failed to expand FASTA inputs: {}", e))?;
    let output = matches
        .get_one::<String>("output")
        .expect("output is required");
    let jobs = *matches.get_one::<usize>("jobs").unwrap_or(&0);
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
    let fmt_auto = |n: usize| if n == 0 { "auto".to_string() } else { n.to_string() };
    eprintln!(
        "Building RefgetStore at {} (mode={}, jobs={})",
        output,
        mode,
        fmt_auto(jobs),
    );

    let mut total_bases: u64 = 0;
    let mut total_seqs: usize = 0;
    let start = Instant::now();

    let opts = FastaImportOptions::new()
        .force(force)
        .jobs(jobs);
    let results = store
        .add_sequence_collections_from_fastas(&fastas, opts)
        .map_err(|e| anyhow::anyhow!("Failed to import FASTA files: {}", e))?;

    for (fa, (metadata, was_new)) in fastas.iter().zip(results.iter()) {
        total_seqs += metadata.n_sequences;
        eprintln!(
            "  {} {}: {} ({} sequences)",
            if *was_new { "added" } else { "skipped" },
            fa.display(),
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
        "Done: {} sequences, {} bases in {:.3}s ({:.1} Mbase/s, jobs={})",
        total_seqs,
        total_bases,
        elapsed,
        mbps,
        fmt_auto(jobs),
    );

    Ok(())
}
