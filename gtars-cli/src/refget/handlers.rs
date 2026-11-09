use std::time::Instant;

use anyhow::Result;
use clap::ArgMatches;

use gtars_refget::store::{
    ExportOptions, ExportSelection, FastaImportOptions, FetchConfig, LabelTemplate, LockOptions,
    Molecule, RefgetStore, StorageMode, VerifyOptions, VerifyReport, default_cache_root,
    force_unlock, lock_status,
};
use gtars_refget::{FastaInputs, expand_fasta_inputs};

pub fn run_refget(matches: &ArgMatches) -> Result<()> {
    match matches.subcommand() {
        Some((super::cli::REFGET_BUILD, sub)) => run_build(sub),
        Some((super::cli::REFGET_EXPORT, sub)) => run_export(sub),
        Some((super::cli::REFGET_LOCK_STATUS, sub)) => run_lock_status(sub),
        Some((super::cli::REFGET_VERIFY, sub)) => run_verify(sub),
        _ => unreachable!("refget subcommand not found"),
    }
}

fn run_lock_status(matches: &ArgMatches) -> Result<()> {
    let store = matches
        .get_one::<String>("store")
        .expect("store is required");
    let path = std::path::Path::new(store);

    match lock_status(path)? {
        Some(info) => {
            println!("RefgetStore write lock HELD at {}", store);
            println!("  holder:        {}", info.describe());
            println!("  last heartbeat: {}", info.heartbeat_at);
        }
        None => println!("RefgetStore write lock is free at {}", store),
    }

    if matches.get_flag("force_unlock") {
        if force_unlock(path)? {
            println!("Lock forcibly cleared.");
        } else {
            println!("Nothing to clear.");
        }
    }

    Ok(())
}

/// Build [`LockOptions`] from the shared `--lock-timeout` / `--force-unlock`
/// flags, falling back to the environment-driven defaults.
fn lock_options_from(matches: &ArgMatches) -> LockOptions {
    let mut opts = LockOptions::default();
    if let Some(secs) = matches.get_one::<u64>("lock_timeout") {
        opts = opts.timeout_secs(*secs);
    }
    opts.force(matches.get_flag("force_unlock"))
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
    store.set_lock_options(lock_options_from(matches));
    store.set_force_alias(matches.get_flag("force_alias"));
    if raw {
        store.set_encoding_mode(StorageMode::Raw)?;
    } else {
        store.set_encoding_mode(StorageMode::Encoded)?;
    }

    let mode = if raw { "Raw" } else { "Encoded" };
    let fmt_auto = |n: usize| {
        if n == 0 {
            "auto".to_string()
        } else {
            n.to_string()
        }
    };
    eprintln!(
        "Building RefgetStore at {} (mode={}, jobs={})",
        output,
        mode,
        fmt_auto(jobs),
    );

    let mut total_bases: u64 = 0;
    let mut total_seqs: usize = 0;
    let start = Instant::now();

    let mut opts = FastaImportOptions::new().force(force).jobs(jobs);
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
    for meta in store.list_sequences()? {
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

/// Resolve `-c` to a collection digest. Accepts a bare digest or
/// `NAMESPACE:ALIAS` (same syntax as `build --collection-alias`).
fn resolve_collection_arg(store: &RefgetStore, raw: &str) -> Result<String> {
    if let Some((ns, alias)) = raw.split_once(':') {
        if ns.is_empty() || alias.is_empty() {
            return Err(anyhow::anyhow!(
                "--collection expects a digest or NAMESPACE:ALIAS (e.g. 'ucsc:hg38'), got '{}'",
                raw
            ));
        }
        return store
            .get_collection_metadata_by_alias(ns, alias)
            .map(|m| m.digest.clone())
            .ok_or_else(|| {
                let known = known_collection_aliases(store);
                anyhow::anyhow!(
                    "Collection alias '{}' not found in store. Known aliases: {}",
                    raw,
                    if known.is_empty() {
                        "(none)".to_string()
                    } else {
                        known.join(", ")
                    }
                )
            });
    }
    Ok(raw.to_string())
}

/// All registered collection aliases as `NAMESPACE:ALIAS` strings, for error messages.
fn known_collection_aliases(store: &RefgetStore) -> Vec<String> {
    store
        .list_collection_alias_namespaces()
        .into_iter()
        .flat_map(|ns| {
            let aliases = store.list_collection_aliases(&ns).unwrap_or_default();
            aliases
                .into_iter()
                .map(move |alias| format!("{}:{}", ns, alias))
        })
        .collect()
}

/// Format a collection digest for an "available collections" listing, appending
/// its aliases (if any) so the alias path is discoverable from the error.
fn describe_collection(store: &RefgetStore, digest: &str) -> String {
    let aliases = store.get_aliases_for_collection(digest);
    if aliases.is_empty() {
        digest.to_string()
    } else {
        let alias_list: Vec<String> = aliases
            .iter()
            .map(|(ns, alias)| format!("{}:{}", ns, alias))
            .collect();
        format!("{} ({})", digest, alias_list.join(", "))
    }
}

fn run_export(matches: &ArgMatches) -> Result<()> {
    let output = matches
        .get_one::<String>("output")
        .map(String::as_str)
        .unwrap_or("-");
    let requested: Vec<&String> = matches
        .get_many::<String>("collection")
        .map(|v| v.collect())
        .unwrap_or_default();
    let names: Option<Vec<String>> = matches
        .get_many::<String>("names")
        .map(|vals| vals.cloned().collect());
    let label = LabelTemplate::parse(matches.get_one::<String>("label").expect("has default"))?;
    let molecule: Option<Molecule> = matches
        .get_one::<String>("molecule")
        .map(|m| m.parse())
        .transpose()?;
    let line_width = *matches.get_one::<usize>("line_width").unwrap_or(&80);

    let remote = matches.get_one::<String>("remote");
    let cache_dir = matches.get_one::<String>("cache_dir");
    let resume = matches.get_flag("resume");
    let mut store = match remote {
        Some(url) => {
            let fetch = FetchConfig {
                jobs: *matches.get_one::<usize>("fetch_jobs").expect("has default"),
                ..FetchConfig::default()
            };
            let root = cache_dir
                .map(std::path::PathBuf::from)
                .unwrap_or_else(default_cache_root);
            RefgetStore::open_remote(&root, url, &fetch)
                .map_err(|e| anyhow::anyhow!("Failed to open remote store {}: {:#}", url, e))?
        }
        None => {
            let store_path = matches
                .get_one::<String>("store")
                .expect("store or remote is required");
            RefgetStore::open_local(store_path)
                .map_err(|e| anyhow::anyhow!("Failed to open store at {}: {}", store_path, e))?
        }
    };

    let selection = if matches.get_flag("all") {
        ExportSelection::All
    } else {
        let collections = store
            .list_collections(0, usize::MAX, &[])
            .map_err(|e| anyhow::anyhow!("Failed to list collections: {}", e))?;
        let available = || -> String {
            collections
                .results
                .iter()
                .map(|m| describe_collection(&store, &m.digest))
                .collect::<Vec<_>>()
                .join(", ")
        };
        let digests: Vec<String> = if requested.is_empty() {
            match collections.results.len() {
                0 => return Err(anyhow::anyhow!("Store contains no collections to export")),
                1 => vec![collections.results[0].digest.clone()],
                _ => {
                    return Err(anyhow::anyhow!(
                        "Store contains multiple collections; specify one with --collection, or use --all. Available: {}",
                        available()
                    ));
                }
            }
        } else {
            let mut out = Vec::new();
            for c in requested {
                let resolved = resolve_collection_arg(&store, c)?;
                if !collections.results.iter().any(|m| m.digest == resolved) {
                    return Err(anyhow::anyhow!(
                        "Collection '{}' not found in store. Available: {}",
                        c,
                        available()
                    ));
                }
                out.push(resolved);
            }
            out
        };
        ExportSelection::Collections { digests, names }
    };

    let opts = ExportOptions {
        label,
        molecule,
        line_width: Some(line_width),
        write_cache: resume || cache_dir.is_some(),
    };
    let summary = match store.export_fasta(&selection, output, &opts) {
        Ok(s) => s,
        Err(e) if is_broken_pipe(&e) => return Ok(()),
        Err(e) => return Err(anyhow::anyhow!("Failed to export FASTA: {}", e)),
    };

    let mut notes = Vec::new();
    if summary.n_fallbacks > 0 {
        notes.push(format!(
            "{} label(s) fell back to a digest",
            summary.n_fallbacks
        ));
    }
    if summary.n_deduplicated > 0 {
        notes.push(format!(
            "{} repeated record(s) written once",
            summary.n_deduplicated
        ));
    }
    let dest = if output == "-" { "stdout" } else { output };
    eprintln!(
        "Exported {} record(s), {} bases, to {}{}",
        summary.n_records,
        summary.n_bases,
        dest,
        if notes.is_empty() {
            String::new()
        } else {
            format!(" ({})", notes.join("; "))
        }
    );

    Ok(())
}

/// Did the export fail because the reader closed the pipe (as `| head` does)?
fn is_broken_pipe(e: &anyhow::Error) -> bool {
    e.chain().any(|c| {
        c.downcast_ref::<std::io::Error>()
            .is_some_and(|io| io.kind() == std::io::ErrorKind::BrokenPipe)
    })
}

/// Strip an optional `SQ.` prefix, matching the convention used elsewhere for
/// user-supplied sequence digests.
fn strip_sq_prefix(digest: &str) -> &str {
    digest.strip_prefix("SQ.").unwrap_or(digest)
}

fn write_verify_tsv<W: std::io::Write>(mut out: W, report: &VerifyReport) -> Result<()> {
    writeln!(
        out,
        "digest\tname\talphabet\tlength\tkind\tcomputed_sha512t24u\tstored_md5\tcomputed_md5\terror"
    )?;
    for f in &report.failures {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            f.digest,
            f.name,
            f.alphabet,
            f.length,
            f.kind.as_str(),
            f.computed_sha512t24u.as_deref().unwrap_or(""),
            f.stored_md5,
            f.computed_md5.as_deref().unwrap_or(""),
            f.error.as_deref().unwrap_or(""),
        )?;
    }
    Ok(())
}

fn run_verify(matches: &ArgMatches) -> Result<()> {
    let store_path = matches
        .get_one::<String>("store")
        .expect("store is required");
    let requested_collection = matches.get_one::<String>("collection");
    let digests: Option<Vec<String>> = matches
        .get_many::<String>("digest")
        .map(|vals| vals.map(|s| strip_sq_prefix(s).to_string()).collect());
    let jobs = *matches.get_one::<usize>("jobs").unwrap_or(&0);
    let check_md5 = !matches.get_flag("no_md5");
    let as_json = matches.get_flag("json");
    let output = matches.get_one::<String>("output");

    let mut store = RefgetStore::open_local(store_path)
        .map_err(|e| anyhow::anyhow!("Failed to open store at {}: {}", store_path, e))?;

    let opts = VerifyOptions { jobs, check_md5 };

    eprintln!(
        "Verifying store at {} (jobs={}, check_md5={})",
        store_path,
        if jobs == 0 {
            "auto".to_string()
        } else {
            jobs.to_string()
        },
        check_md5,
    );

    let report = if let Some(c) = requested_collection {
        let resolved = resolve_collection_arg(&store, c)?;
        store
            .verify_collection(&resolved, &opts)
            .map_err(|e| anyhow::anyhow!("Failed to verify collection {}: {}", resolved, e))?
    } else {
        store
            .verify_sequences(digests.as_deref(), &opts)
            .map_err(|e| anyhow::anyhow!("Failed to verify store: {}", e))?
    };

    // Write the report before returning any error, so a failing exit status
    // still leaves the failure listing on disk / on stdout for scripts to read.
    let rendered = if as_json {
        serde_json::to_string_pretty(&report)
            .map_err(|e| anyhow::anyhow!("Failed to serialize verify report as JSON: {}", e))?
    } else {
        let mut buf = Vec::new();
        write_verify_tsv(&mut buf, &report)?;
        String::from_utf8(buf).expect("TSV output is always valid UTF-8")
    };

    match output {
        Some(path) => {
            std::fs::write(path, &rendered)
                .map_err(|e| anyhow::anyhow!("Failed to write report to {}: {}", path, e))?;
        }
        None => {
            print!("{}", rendered);
            if !rendered.ends_with('\n') {
                println!();
            }
        }
    }

    eprintln!(
        "Verified {} sequences: OK {}, FAILED {}",
        report.n_checked, report.n_ok, report.n_failed
    );
    for (alphabet, counts) in &report.by_alphabet {
        eprintln!(
            "  {}: checked {}, failed {}",
            alphabet, counts.checked, counts.failed
        );
    }

    if report.n_failed > 0 {
        eprintln!(
            "Failing sequences cannot be repaired in place; re-import the affected \
             collection(s) from the source FASTA (gtars refget build --force ...)."
        );
        return Err(anyhow::anyhow!(
            "{} of {} sequences failed verification",
            report.n_failed,
            report.n_checked
        ));
    }

    Ok(())
}
