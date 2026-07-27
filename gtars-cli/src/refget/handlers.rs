use std::time::Instant;

use anyhow::Result;
use clap::ArgMatches;

use gtars_refget::store::{FastaImportOptions, RefgetStore, StorageMode};
use gtars_refget::{expand_fasta_inputs, FastaInputs};

pub fn run_refget(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some((super::cli::REFGET_BUILD, sub)) => run_build(sub),
        Some((super::cli::REFGET_EXPORT, sub)) => run_export(sub),
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

    // Parse --collection-alias NAMESPACE:ALIAS. Owned up front so the borrows
    // in FastaImportOptions outlive the builder chain.
    let collection_alias: Option<(String, String)> = matches
        .get_one::<String>("collection_alias")
        .map(|raw| {
            let (ns, alias) = raw.split_once(':').ok_or_else(|| {
                anyhow::anyhow!(
                    "--collection-alias expects NAMESPACE:ALIAS (e.g. 'ucsc:hg38'), got '{}'",
                    raw
                )
            })?;
            if ns.is_empty() || alias.is_empty() {
                return Err(anyhow::anyhow!(
                    "--collection-alias expects a non-empty namespace and alias in \
                     NAMESPACE:ALIAS (e.g. 'ucsc:hg38'), got '{}'",
                    raw
                ));
            }
            Ok((ns.to_string(), alias.to_string()))
        })
        .transpose()?;

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

    let mut opts = FastaImportOptions::new()
        .force(force)
        .jobs(jobs);
    if let Some((ns, alias)) = collection_alias.as_ref() {
        opts = opts.collection_alias(ns, alias);
    }
    let report = store
        .add_sequence_collections_from_fastas(&fastas, opts)
        .map_err(|e| anyhow::anyhow!("Failed to import FASTA files: {}", e))?;

    for (fa, (metadata, was_new)) in fastas.iter().zip(report.collections.iter()) {
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
    // Per-run ingest counters: what THIS run actually added, as opposed to the
    // store-wide residency numbers reported by `store stats`.
    eprintln!(
        "Ingested this run: {} collection(s) new, {} sequence(s) written, {} sequence(s) deduped",
        report.n_collections_new, report.n_sequences_written, report.n_sequences_deduped,
    );

    Ok(())
}

fn run_export(matches: &ArgMatches) -> Result<()> {
    let store_path = matches
        .get_one::<String>("store")
        .expect("store is required");
    let output = matches
        .get_one::<String>("output")
        .expect("output is required");
    let requested_collection = matches.get_one::<String>("collection");
    let names: Option<Vec<&str>> = matches
        .get_many::<String>("names")
        .map(|vals| vals.map(|s| s.as_str()).collect());
    let line_width = *matches.get_one::<usize>("line_width").unwrap_or(&80);

    // Open the store and load collection metadata (stub records + name_lookup).
    // Sequence BYTES are loaded further down, after the digest is validated and
    // only for what this export actually needs.
    let mut store = RefgetStore::open_local(store_path)
        .map_err(|e| anyhow::anyhow!("Failed to open store at {}: {}", store_path, e))?;
    store
        .load_all_collections()
        .map_err(|e| anyhow::anyhow!("Failed to load collections: {}", e))?;

    // Resolve the collection digest.
    let collections = store
        .list_collections(0, usize::MAX, &[])
        .map_err(|e| anyhow::anyhow!("Failed to list collections: {}", e))?;
    let digest = match requested_collection {
        Some(c) => {
            if !collections.results.iter().any(|m| &m.digest == c) {
                let available: Vec<String> =
                    collections.results.iter().map(|m| m.digest.clone()).collect();
                return Err(anyhow::anyhow!(
                    "Collection '{}' not found in store. Available: {}",
                    c,
                    available.join(", ")
                ));
            }
            c.clone()
        }
        None => match collections.results.len() {
            0 => return Err(anyhow::anyhow!("Store contains no collections to export")),
            1 => collections.results[0].digest.clone(),
            _ => {
                let available: Vec<String> =
                    collections.results.iter().map(|m| m.digest.clone()).collect();
                return Err(anyhow::anyhow!(
                    "Store contains multiple collections; specify one with --collection. Available: {}",
                    available.join(", ")
                ));
            }
        },
    };

    // Load ONLY the sequence bytes this export needs. `load_all_sequences()`
    // would pull EVERY sequence in the store into RAM, not just this
    // collection's -- fatal on a large store (the vgp store holds ~384k
    // sequences / hundreds of GB) and wasteful even when it fits. With
    // `--names`, narrow further to just the requested sequences.
    let collection = store
        .get_collection(&digest)
        .map_err(|e| anyhow::anyhow!("Failed to load collection {}: {}", digest, e))?;
    let wanted: Option<std::collections::HashSet<&str>> =
        names.as_ref().map(|v| v.iter().copied().collect());
    for record in &collection.sequences {
        let meta = record.metadata();
        if wanted
            .as_ref()
            .is_some_and(|wanted| !wanted.contains(meta.name.as_str()))
        {
            continue;
        }
        store.load_sequence(&meta.sha512t24u).map_err(|e| {
            anyhow::anyhow!("Failed to load sequence '{}': {}", meta.name, e)
        })?;
    }
    let store = store.into_readonly();

    let n_names = names.as_ref().map(|v| v.len());
    store
        .export_fasta(&digest, output, names, Some(line_width))
        .map_err(|e| anyhow::anyhow!("Failed to export FASTA: {}", e))?;

    let wrap_desc = if line_width == 0 {
        "unwrapped (one sequence per line)".to_string()
    } else {
        format!("wrapped at {} bases/line", line_width)
    };
    let seq_desc = match n_names {
        Some(n) => format!("{} named sequence(s)", n),
        None => "all sequences".to_string(),
    };
    eprintln!(
        "Exported collection {} ({}) to {} [{}]",
        digest, seq_desc, output, wrap_desc
    );

    Ok(())
}
