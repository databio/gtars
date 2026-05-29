//! FASTA import pipeline for RefgetStore.
//!
//! Contains the multithreaded FASTA import pipeline and cached metadata fast path.
//!
//! ## Two-axis parallelism model
//!
//! ### Per-sequence parallelism (within one file)
//!
//! The CPU-bound per-sequence work (SHA-512 + MD5 + alphabet guess + 2-bit
//! encode) is embarrassingly parallel across sequences and touches no shared
//! state. We fan that work out to N scoped worker threads. A reader thread
//! streams FASTA records and dispatches them round-robin to per-worker channels,
//! tagging each item with a monotonic `input_index`. The collector reorders
//! worker output by `input_index` so records are observed in exact FASTA order;
//! this is what makes `name_lookup` (and therefore the collection digest)
//! deterministic and identical to a single-threaded build. The pattern mirrors
//! `gtars-vrs/src/vcf.rs::compute_vrs_ids_parallel_encoded`.
//!
//! ### File-level parallelism (across files)
//!
//! A single reader thread doing gzip DECOMPRESSION is the bottleneck on
//! compressed input (the per-sequence workers starve). To parallelize
//! decompression we run K file pipelines concurrently, each with its OWN
//! decoder (so K gzip streams decompress in parallel). The build half (decode +
//! parse + digest + encode -> in-memory `BuiltCollection`) touches NO shared
//! store state and runs concurrently across files; the insert half (`&mut self`
//! registration of collections/sequences/aliases/name_lookup) runs on a single
//! owner in FIXED input-file order, so the resulting store is byte-identical to
//! a serial build regardless of which builder finishes first.
//!
//! `threads` is the TOTAL encode-worker budget; `file_jobs` is the number of
//! files in flight. They are reconciled by `resolve_threads` /
//! `resolve_file_jobs` so `file_jobs * per_file_workers <= threads`.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{anyhow, Result};

use crate::collection::SequenceCollectionExt;
use crate::digest::{
    SequenceCollection, SequenceCollectionMetadata,
    SequenceEncoder, SequenceMetadata, SequenceRecord,
};
use crate::hashkeyable::HashKeyable;

/// Threshold above which the reader processes sequences inline (to bound peak
/// memory: the few huge contigs are hashed+encoded on the reader instead of
/// fanning N copies into the worker pool).
const LARGE_SEQ_THRESHOLD: usize = 500 * 1024 * 1024; // 500 MB

/// Resolve a requested thread count, mapping `0` to the machine parallelism.
fn resolve_threads(requested: usize) -> usize {
    if requested == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        requested
    }
    .max(1)
}

/// Resolve how many files to decode/build concurrently (`file_jobs`).
///
/// `requested == 0` means auto. The auto heuristic biases toward more files in
/// flight (more concurrent gzip decoders, since decode is the bottleneck on
/// compressed input) while giving each file at least `PER_FILE_MIN_WORKERS`
/// encode workers, and never exceeding the number of input files or the total
/// thread budget. With a single file we always return 1 (all workers go to that
/// file, preserving the prior single-file behavior).
fn resolve_file_jobs(requested: usize, num_files: usize, n_workers: usize) -> usize {
    /// Each file gets at least this many encode workers under the auto heuristic.
    const PER_FILE_MIN_WORKERS: usize = 2;

    let num_files = num_files.max(1);
    if num_files == 1 {
        return 1;
    }
    let resolved = if requested == 0 {
        // ceil(n_workers / PER_FILE_MIN_WORKERS), capped to num_files.
        let by_budget = n_workers.div_ceil(PER_FILE_MIN_WORKERS);
        by_budget.clamp(1, num_files)
    } else {
        requested.min(num_files)
    };
    resolved.max(1)
}

/// Derive per-file encode workers from the total budget and files-in-flight,
/// so `file_jobs * per_file <= n_workers` (each file gets at least 1).
fn per_file_workers(n_workers: usize, file_jobs: usize) -> usize {
    (n_workers / file_jobs.max(1)).max(1)
}

/// A fully-built sequence collection produced by the build half (no `&mut self`).
///
/// Holds everything needed to insert the collection into the store: the encoded
/// sequences in FASTA order, the (namespace, alias, sha512) triples, and the
/// computed collection metadata. The insert half consumes this on a single
/// owner thread in fixed file order.
struct BuiltCollection {
    /// Collection metadata (digest + ancillary digests already computed).
    metadata: SequenceCollectionMetadata,
    /// Sequence stub records in FASTA order (used to register the collection).
    stub_records: Vec<SequenceRecord>,
    /// Fully-encoded sequences in FASTA order, ready to write.
    sequences: Vec<ReadyOrdered>,
    /// (namespace, alias_value, sha512t24u) alias triples in FASTA order.
    aliases: Vec<(String, String, String)>,
    /// Source FASTA path (for metadata / RGSI cache).
    source_path: PathBuf,
    /// RGSI cache path to (re)write, if caching is enabled and this was built
    /// from raw FASTA (not from an existing cache).
    rgsi_cache_path: Option<PathBuf>,
    /// Number of sequences (for reporting).
    seq_count: usize,
}

/// A fully-built sequence ready for the single inserter, tagged with its input
/// (FASTA) index so the collector can restore order.
struct ReadyOrdered {
    input_index: u64,
    metadata: SequenceMetadata,
    sequence_data: Vec<u8>,
    aliases: Vec<(String, String)>,
}

/// Run the reader -> N workers -> in-order collector pipeline.
///
/// - `reader` runs on its own scoped thread. It is handed N round-robin work
///   senders (one per worker) and must tag each emitted item with a monotonic
///   `input_index` starting at 0 with no gaps. It returns the number of items
///   dispatched (unused by callers but useful for assertions).
/// - `work_fn` is the per-item CPU work, run on N scoped worker threads. It maps
///   one input item `T` (already carrying its `input_index`) to a `ReadyOrdered`.
/// - `on_ready` runs on the **main** thread, receiving `ReadyOrdered` items in
///   exact `input_index` order. This is where `&mut self` insertion happens.
///
/// `T` must be `Send` and carry its own `input_index` (the reader assigns it).
fn run_encode_pipeline<T, R, W, C>(
    n_workers: usize,
    reader: R,
    work_fn: W,
    mut on_ready: C,
) -> Result<()>
where
    T: Send,
    R: FnOnce(&[crossbeam_channel::Sender<T>]) -> Result<u64> + Send,
    W: Fn(T) -> Result<ReadyOrdered> + Send + Sync,
    C: FnMut(ReadyOrdered) -> Result<()>,
{
    use crossbeam_channel::{bounded, unbounded};

    let n_workers = n_workers.max(1);
    let work_fn = &work_fn;

    // Output channel: workers push ReadyOrdered, the collector (main thread)
    // drains it. Unbounded so workers never block on the collector; in-flight
    // raw bytes are bounded by the small per-worker input channels instead.
    let (out_tx, out_rx) = unbounded::<Result<ReadyOrdered>>();

    std::thread::scope(|scope| -> Result<()> {
        // One bounded work channel per worker. bounded(2) keeps the reader from
        // pulling the whole FASTA into RAM (a few contigs in flight per worker).
        let mut work_txs: Vec<crossbeam_channel::Sender<T>> = Vec::with_capacity(n_workers);
        let mut worker_handles = Vec::with_capacity(n_workers);

        for _ in 0..n_workers {
            let (work_tx, work_rx) = bounded::<T>(2);
            work_txs.push(work_tx);
            let out_tx = out_tx.clone();
            let handle = scope.spawn(move || {
                for item in work_rx.iter() {
                    let result = work_fn(item);
                    let is_err = result.is_err();
                    if out_tx.send(result).is_err() {
                        break;
                    }
                    if is_err {
                        break;
                    }
                }
            });
            worker_handles.push(handle);
        }
        drop(out_tx); // workers hold the only remaining output senders

        // Reader thread: dispatches work round-robin and assigns input indices.
        let reader_handle = scope.spawn(move || -> Result<u64> {
            let n = reader(&work_txs)?;
            // Dropping all work senders signals workers to finish.
            drop(work_txs);
            Ok(n)
        });

        // Collector (main thread): reorder by input_index and emit in order.
        let mut pending: HashMap<u64, ReadyOrdered> = HashMap::new();
        let mut next_emit: u64 = 0;
        let mut first_err: Option<anyhow::Error> = None;

        for msg in out_rx.iter() {
            match msg {
                Ok(ready) => {
                    pending.insert(ready.input_index, ready);
                    while let Some(item) = pending.remove(&next_emit) {
                        if let Err(e) = on_ready(item) {
                            first_err = Some(e);
                            break;
                        }
                        next_emit += 1;
                    }
                    if first_err.is_some() {
                        break;
                    }
                }
                Err(e) => {
                    first_err = Some(e);
                    break;
                }
            }
        }

        // Drain any remaining in-order items (none if a worker errored).
        if first_err.is_none() {
            while let Some(item) = pending.remove(&next_emit) {
                if let Err(e) = on_ready(item) {
                    first_err = Some(e);
                    break;
                }
                next_emit += 1;
            }
        }

        // Join reader and workers, surfacing panics and the reader error.
        let reader_result = reader_handle
            .join()
            .map_err(|_| anyhow!("FASTA reader thread panicked"))?;
        for h in worker_handles {
            h.join()
                .map_err(|_| anyhow!("FASTA worker thread panicked"))?;
        }

        if let Some(e) = first_err {
            return Err(e);
        }
        reader_result?;
        Ok(())
    })
}

// ============================================================================
// Build half (no &mut self): decode + parse + digest + encode -> BuiltCollection
//
// These free functions run concurrently across files. They touch NO shared
// store state. The returned `BuiltCollection` is consumed by the insert half
// on a single owner thread in fixed file order.
// ============================================================================

/// Configuration needed by the build half, extracted from the store so the
/// build can run without borrowing `&self` across threads.
#[derive(Clone, Copy)]
struct BuildConfig<'a> {
    mode: StorageMode,
    quiet: bool,
    ancillary_digests: bool,
    /// Whether the store is disk-backed (controls RGSI cache read/write).
    use_cache: bool,
    namespaces: &'a [&'a str],
    n_workers: usize,
}

/// Build a `BuiltCollection` from a FASTA file (raw path: digest + encode).
///
/// Checks for an RGSI metadata cache first and dispatches to the cached fast
/// path if present and valid; otherwise runs the full digest+encode pipeline.
/// Does NOT touch `&mut self`.
fn build_collection_from_fasta(
    file_path: &Path,
    cfg: BuildConfig<'_>,
) -> Result<BuiltCollection> {
    use crate::utils::PathExtension;

    let rgsi_path = file_path.replace_exts_with("rgsi");
    let have_rgsi = cfg.use_cache && rgsi_path.exists();

    if have_rgsi {
        if let Ok(built) = build_collection_from_cached_metadata(file_path, &rgsi_path, cfg) {
            return Ok(built);
        }
        // Cache was stale/empty; fall through to the full pipeline.
    }

    build_collection_full(file_path, cfg)
}

/// Full build half: digest + encode every sequence and compute collection
/// metadata. Returns a `BuiltCollection` ready for the insert half.
fn build_collection_full(file_path: &Path, cfg: BuildConfig<'_>) -> Result<BuiltCollection> {
    use crate::utils::PathExtension;

    // Reader item type: a raw record needing work, or a large record already
    // processed inline by the reader. Both carry their input_index.
    struct RawItem {
        input_index: u64,
        name: String,
        description: Option<String>,
        raw_header: String,
        raw_bytes: Vec<u8>,
    }
    enum WorkUnit {
        NeedsWork(RawItem),
        AlreadyDone(ReadyOrdered),
    }

    let file_path_buf = file_path.to_path_buf();
    let namespaces: Vec<String> = cfg.namespaces.iter().map(|s| s.to_string()).collect();
    let mode = cfg.mode;
    let quiet = cfg.quiet;

    let mut sequences: Vec<ReadyOrdered> = Vec::new();

    // --- Reader closure (runs on a scoped thread) ---
    let reader = {
        let namespaces = &namespaces;
        move |work_txs: &[crossbeam_channel::Sender<WorkUnit>]| -> Result<u64> {
            use md5::Md5;
            use sha2::{Digest, Sha512};

            let mut fasta_reader = crate::fasta::FastaReader::from_path(&file_path_buf)?;
            let n = work_txs.len();
            let mut input_index: u64 = 0;
            let mut next_worker = 0usize;

            while let Some(record) = fasta_reader.next_record()? {
                let unit = if record.raw_bytes.len() > LARGE_SEQ_THRESHOLD {
                    if !quiet {
                        println!(
                            "  Large sequence '{}' ({} MB) -- processing inline to reduce memory",
                            record.name,
                            record.raw_bytes.len() / (1024 * 1024),
                        );
                    }
                    let crate::fasta::FastaRecord { name, description, raw_header, raw_bytes } = record;

                    let mut sha512_hasher = Sha512::new();
                    sha512_hasher.update(&raw_bytes);
                    let sha512 = base64_url::encode(&sha512_hasher.finalize()[0..24]);

                    let mut md5_hasher = Md5::new();
                    md5_hasher.update(&raw_bytes);
                    let md5 = format!("{:x}", md5_hasher.finalize());

                    let mut guesser = crate::digest::AlphabetGuesser::new();
                    guesser.update(&raw_bytes);
                    let alphabet = guesser.guess();

                    let length = raw_bytes.len();
                    let aliases = if !namespaces.is_empty() {
                        let ns_refs: Vec<&str> = namespaces.iter().map(|s| s.as_str()).collect();
                        crate::digest::fasta::extract_aliases_from_header(&raw_header, &ns_refs)
                    } else {
                        vec![]
                    };
                    let sequence_data = match mode {
                        StorageMode::Encoded => {
                            let mut encoder = SequenceEncoder::new(alphabet, length);
                            encoder.update(&raw_bytes);
                            drop(raw_bytes);
                            encoder.finalize()
                        }
                        StorageMode::Raw => raw_bytes,
                    };

                    let metadata = SequenceMetadata {
                        name,
                        description,
                        length,
                        sha512t24u: sha512,
                        md5,
                        alphabet,
                        fai: None,
                    };

                    WorkUnit::AlreadyDone(ReadyOrdered {
                        input_index,
                        metadata,
                        sequence_data,
                        aliases,
                    })
                } else {
                    WorkUnit::NeedsWork(RawItem {
                        input_index,
                        name: record.name,
                        description: record.description,
                        raw_header: record.raw_header,
                        raw_bytes: record.raw_bytes,
                    })
                };

                if work_txs[next_worker].send(unit).is_err() {
                    break; // worker hung up
                }
                input_index += 1;
                next_worker = (next_worker + 1) % n;
            }
            Ok(input_index)
        }
    };

    // --- Worker closure (runs on N scoped threads) ---
    let work_fn = {
        let namespaces = &namespaces;
        move |unit: WorkUnit| -> Result<ReadyOrdered> {
            use md5::Md5;
            use sha2::{Digest, Sha512};

            match unit {
                WorkUnit::AlreadyDone(ready) => Ok(ready),
                WorkUnit::NeedsWork(item) => {
                    let mut sha512_hasher = Sha512::new();
                    sha512_hasher.update(&item.raw_bytes);
                    let sha512 = base64_url::encode(&sha512_hasher.finalize()[0..24]);

                    let mut md5_hasher = Md5::new();
                    md5_hasher.update(&item.raw_bytes);
                    let md5 = format!("{:x}", md5_hasher.finalize());

                    let mut guesser = crate::digest::AlphabetGuesser::new();
                    guesser.update(&item.raw_bytes);
                    let alphabet = guesser.guess();

                    let length = item.raw_bytes.len();
                    let aliases = if !namespaces.is_empty() {
                        let ns_refs: Vec<&str> = namespaces.iter().map(|s| s.as_str()).collect();
                        crate::digest::fasta::extract_aliases_from_header(&item.raw_header, &ns_refs)
                    } else {
                        vec![]
                    };

                    let sequence_data = match mode {
                        StorageMode::Encoded => {
                            let mut encoder = SequenceEncoder::new(alphabet, length);
                            encoder.update(&item.raw_bytes);
                            drop(item.raw_bytes);
                            encoder.finalize()
                        }
                        StorageMode::Raw => item.raw_bytes,
                    };

                    let metadata = SequenceMetadata {
                        name: item.name,
                        description: item.description,
                        length,
                        sha512t24u: sha512,
                        md5,
                        alphabet,
                        fai: None,
                    };

                    Ok(ReadyOrdered {
                        input_index: item.input_index,
                        metadata,
                        sequence_data,
                        aliases,
                    })
                }
            }
        }
    };

    // Collector: gather encoded sequences in FASTA order (no &mut self).
    run_encode_pipeline(
        cfg.n_workers,
        reader,
        work_fn,
        |ready: ReadyOrdered| -> Result<()> {
            sequences.push(ready);
            Ok(())
        },
    )?;

    // Compute collection metadata from the sequence stubs in FASTA order.
    let stub_records: Vec<SequenceRecord> = sequences
        .iter()
        .map(|s| SequenceRecord::Stub(s.metadata.clone()))
        .collect();

    let mut seqcol_metadata = SequenceCollectionMetadata::from_sequences(
        &stub_records,
        Some(file_path.to_path_buf()),
    );
    if cfg.ancillary_digests {
        seqcol_metadata.compute_ancillary_digests(&stub_records);
    }

    // Flatten alias triples in FASTA order.
    let mut aliases: Vec<(String, String, String)> = Vec::new();
    for ready in &sequences {
        for (ns, alias_value) in &ready.aliases {
            aliases.push((ns.clone(), alias_value.clone(), ready.metadata.sha512t24u.clone()));
        }
    }

    let seq_count = sequences.len();
    let rgsi_cache_path = if cfg.use_cache {
        Some(file_path.replace_exts_with("rgsi"))
    } else {
        None
    };

    Ok(BuiltCollection {
        metadata: seqcol_metadata,
        stub_records,
        sequences,
        aliases,
        source_path: file_path.to_path_buf(),
        rgsi_cache_path,
        seq_count,
    })
}

/// Build half for the RGSI cached-metadata fast path: skip digesting, only
/// decompress + encode. Does NOT touch `&mut self`.
fn build_collection_from_cached_metadata(
    file_path: &Path,
    rgsi_path: &Path,
    cfg: BuildConfig<'_>,
) -> Result<BuiltCollection> {
    let mut seqcol = crate::collection::read_rgsi_file(rgsi_path)?;

    if seqcol.sequences.is_empty() {
        let _ = std::fs::remove_file(rgsi_path);
        return Err(anyhow!("Empty RGSI cache"));
    }

    if cfg.ancillary_digests {
        seqcol.metadata.compute_ancillary_digests(&seqcol.sequences);
    }

    let seqmeta_hashmap: HashMap<String, SequenceMetadata> = seqcol
        .sequences
        .iter()
        .map(|r| {
            let meta = r.metadata().clone();
            (meta.name.clone(), meta)
        })
        .collect();

    let mode = cfg.mode;
    let namespaces: Vec<String> = cfg.namespaces.iter().map(|s| s.to_string()).collect();
    let file_path_buf = file_path.to_path_buf();

    // Reader item: a raw record plus its already-known cached metadata.
    struct CachedRaw {
        input_index: u64,
        metadata: SequenceMetadata,
        raw_bytes: Vec<u8>,
        aliases: Vec<(String, String)>,
    }

    let mut sequences: Vec<ReadyOrdered> = Vec::new();

    let reader = {
        let namespaces = &namespaces;
        let seqmeta_hashmap = &seqmeta_hashmap;
        move |work_txs: &[crossbeam_channel::Sender<CachedRaw>]| -> Result<u64> {
            let mut fasta_reader = crate::fasta::FastaReader::from_path(&file_path_buf)?;
            let n = work_txs.len();
            let mut input_index: u64 = 0;
            let mut next_worker = 0usize;

            while let Some(record) = fasta_reader.next_record()? {
                let (name, _) = crate::fasta::parse_fasta_header(&record.raw_header);

                let aliases = if !namespaces.is_empty() {
                    let ns_refs: Vec<&str> = namespaces.iter().map(|s| s.as_str()).collect();
                    crate::digest::fasta::extract_aliases_from_header(&record.raw_header, &ns_refs)
                } else {
                    vec![]
                };

                let dr = seqmeta_hashmap
                    .get(&name)
                    .ok_or_else(|| anyhow!("Sequence '{}' not found in cached metadata", name))?
                    .clone();

                let unit = CachedRaw {
                    input_index,
                    metadata: dr,
                    raw_bytes: record.raw_bytes,
                    aliases,
                };
                if work_txs[next_worker].send(unit).is_err() {
                    break;
                }
                input_index += 1;
                next_worker = (next_worker + 1) % n;
            }
            Ok(input_index)
        }
    };

    let work_fn = move |unit: CachedRaw| -> Result<ReadyOrdered> {
        let sequence_data = match mode {
            StorageMode::Encoded => {
                let mut encoder =
                    SequenceEncoder::new(unit.metadata.alphabet, unit.metadata.length);
                encoder.update(&unit.raw_bytes);
                drop(unit.raw_bytes);
                encoder.finalize()
            }
            StorageMode::Raw => unit.raw_bytes,
        };
        Ok(ReadyOrdered {
            input_index: unit.input_index,
            metadata: unit.metadata,
            sequence_data,
            aliases: unit.aliases,
        })
    };

    run_encode_pipeline(
        cfg.n_workers,
        reader,
        work_fn,
        |ready: ReadyOrdered| -> Result<()> {
            sequences.push(ready);
            Ok(())
        },
    )?;

    let stub_records = seqcol.sequences.clone();
    let metadata = seqcol.metadata.clone();

    let mut aliases: Vec<(String, String, String)> = Vec::new();
    for ready in &sequences {
        for (ns, alias_value) in &ready.aliases {
            aliases.push((ns.clone(), alias_value.clone(), ready.metadata.sha512t24u.clone()));
        }
    }

    let seq_count = sequences.len();

    Ok(BuiltCollection {
        metadata,
        stub_records,
        sequences,
        aliases,
        source_path: file_path.to_path_buf(),
        // Cache already exists and is valid; no need to rewrite it.
        rgsi_cache_path: None,
        seq_count,
    })
}

// ============================================================================
// ReadonlyRefgetStore import methods
// ============================================================================

impl ReadonlyRefgetStore {
    /// Import multiple FASTA files into the store with file-level parallelism.
    ///
    /// Up to `file_jobs` files are decoded/built concurrently (each with its own
    /// decoder, so gzip decompression is parallelized across files). The insert
    /// half (`&mut self`) is applied on this single owner thread in FIXED input
    /// file order, guaranteeing a byte-identical store to a serial build.
    ///
    /// Returns per-file `(collection_metadata, was_new)` results in input order.
    pub fn add_sequence_collections_from_fastas(
        &mut self,
        files: &[PathBuf],
        opts: FastaImportOptions<'_>,
    ) -> Result<Vec<(SequenceCollectionMetadata, bool)>> {
        if files.is_empty() {
            return Ok(Vec::new());
        }

        let n_workers = resolve_threads(opts.threads);
        let file_jobs = resolve_file_jobs(opts.file_jobs, files.len(), n_workers);
        let workers_per_file = per_file_workers(n_workers, file_jobs);

        let cfg = BuildConfig {
            mode: self.mode,
            quiet: self.quiet,
            ancillary_digests: self.ancillary_digests,
            use_cache: self.local_path.is_some(),
            namespaces: opts.namespaces,
            n_workers: workers_per_file,
        };

        let overall_start = Instant::now();

        // --- Build half: run up to `file_jobs` builders concurrently. ---
        // Results are streamed back tagged with the file index; the owner
        // collects them and inserts in fixed file index order.
        let mut built: Vec<Option<Result<BuiltCollection>>> =
            (0..files.len()).map(|_| None).collect();

        if file_jobs <= 1 {
            // Serial build (preserves single-file behavior; also used when the
            // budget only allows one file at a time).
            for (idx, file) in files.iter().enumerate() {
                if !self.quiet {
                    println!("Processing {}...", file.display());
                }
                built[idx] = Some(build_collection_from_fasta(file, cfg));
            }
        } else {
            use crossbeam_channel::{bounded, unbounded};
            // Bounded job queue feeds a fixed pool of `file_jobs` builder
            // threads; results stream back unbounded (each result is one
            // BuiltCollection, bounded in count by the number of files).
            let (job_tx, job_rx) = bounded::<(usize, PathBuf)>(file_jobs);
            let (res_tx, res_rx) = unbounded::<(usize, Result<BuiltCollection>)>();

            std::thread::scope(|scope| {
                for _ in 0..file_jobs {
                    let job_rx = job_rx.clone();
                    let res_tx = res_tx.clone();
                    scope.spawn(move || {
                        for (idx, path) in job_rx.iter() {
                            let result = build_collection_from_fasta(&path, cfg);
                            if res_tx.send((idx, result)).is_err() {
                                break;
                            }
                        }
                    });
                }
                drop(res_tx);

                // Feed jobs from this thread, draining results as we go so the
                // result channel never grows beyond the in-flight builders.
                for (idx, file) in files.iter().enumerate() {
                    if !self.quiet {
                        println!("Processing {}...", file.display());
                    }
                    job_tx.send((idx, file.clone())).expect("builder pool dropped");
                }
                drop(job_tx);

                for (idx, result) in res_rx.iter() {
                    built[idx] = Some(result);
                }
            });
        }

        // --- Insert half: single owner, FIXED file order. ---
        let mut results: Vec<(SequenceCollectionMetadata, bool)> = Vec::with_capacity(files.len());
        for slot in built.into_iter() {
            let built = slot.expect("every file index must be built")?;
            let (metadata, was_new) = self.insert_built_collection(built, opts.force)?;
            results.push((metadata, was_new));
        }

        if !self.quiet {
            println!(
                "Imported {} file(s) in {:.1}s (threads={}, file_jobs={}, workers/file={})",
                files.len(),
                overall_start.elapsed().as_secs_f64(),
                n_workers,
                file_jobs,
                workers_per_file,
            );
        }

        Ok(results)
    }

    /// Import a single FASTA file. Thin wrapper over the multi-file path.
    pub fn add_sequence_collection_from_fasta<P: AsRef<Path>>(
        &mut self,
        file_path: P,
        opts: FastaImportOptions<'_>,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
        let files = [file_path.as_ref().to_path_buf()];
        let mut results = self.add_sequence_collections_from_fastas(&files, opts)?;
        Ok(results.pop().expect("one file yields one result"))
    }

    /// Insert half: register a fully-built collection on `&mut self` in fixed
    /// order. Honors `force` and the already-exists skip. Dedup-by-digest is
    /// handled by `add_sequence_record` (keyed on content digest).
    fn insert_built_collection(
        &mut self,
        built: BuiltCollection,
        force: bool,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
        let BuiltCollection {
            metadata,
            stub_records,
            sequences,
            aliases,
            source_path,
            rgsi_cache_path,
            seq_count,
        } = built;

        let coll_key = metadata.digest.to_key();
        let coll_digest_display = metadata.digest.clone();

        if !force && self.collections.contains_key(&coll_key) {
            if !self.quiet {
                println!("Skipped {} (already exists)", coll_digest_display);
            }
            return Ok((metadata, false));
        }

        let insert_start = Instant::now();

        // Write the RGSI cache for next time (only for raw-FASTA builds).
        if let Some(rgsi_path) = &rgsi_cache_path {
            let seqcol_for_cache = SequenceCollection {
                metadata: metadata.clone(),
                sequences: stub_records.clone(),
            };
            let _ = seqcol_for_cache.write_collection_rgsi(rgsi_path);
        }

        // Register the collection.
        let seqcol = SequenceCollection {
            metadata: metadata.clone(),
            sequences: stub_records,
        };
        self.add_sequence_collection_internal(seqcol, force)?;

        // Insert encoded sequences in FASTA order and populate name_lookup.
        // `add_sequence_collection_internal` above already inserted stub records
        // (no sequence data) for these digests, so we MUST force here to write
        // the actual sequence bytes through to disk. This matches the prior
        // behavior (both the old raw and cached paths wrote sequences with
        // force=true). Dedup-by-digest is preserved: a shared sequence is keyed
        // on its content digest, so a later collection writes byte-identical
        // bytes to the same path.
        for ready in sequences {
            let ReadyOrdered { metadata: seq_meta, sequence_data, .. } = ready;
            self.add_sequence(
                SequenceRecord::Full {
                    metadata: seq_meta,
                    sequence: sequence_data,
                },
                coll_key,
                true,
            )?;
        }

        // Register aliases in FASTA order.
        for (ns, alias_value, sha512t24u) in &aliases {
            self.add_sequence_alias(ns, alias_value, sha512t24u)?;
        }

        if !self.quiet {
            println!(
                "Added {} ({} seqs) from {} in {:.1}s",
                coll_digest_display,
                seq_count,
                source_path.display(),
                insert_start.elapsed().as_secs_f64(),
            );
        }

        Ok((metadata, true))
    }
}
