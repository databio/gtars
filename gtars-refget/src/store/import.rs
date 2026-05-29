//! FASTA import pipeline for RefgetStore.
//!
//! Contains the multithreaded FASTA import pipeline and cached metadata fast path.
//!
//! ## Parallelism model (Option A from `parallel_store_encoding.md`)
//!
//! The CPU-bound per-sequence work (SHA-512 + MD5 + alphabet guess + 2-bit
//! encode) is embarrassingly parallel across sequences and touches no shared
//! state. We fan that work out to N scoped worker threads while keeping a single
//! inserter on the main thread (the only `&mut self` holder) calling the
//! existing `add_sequence_record` / `add_sequence`. A reader thread streams
//! FASTA records and dispatches them round-robin to per-worker channels, tagging
//! each item with a monotonic `input_index`. The collector reorders worker
//! output by `input_index` so the main thread observes records in exact FASTA
//! order; this is what makes `name_lookup` (and therefore the collection digest)
//! deterministic and identical to a single-threaded build. The pattern mirrors
//! `gtars-vrs/src/vcf.rs::compute_vrs_ids_parallel_encoded`.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::collections::HashMap;
use std::path::Path;
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
// ReadonlyRefgetStore import methods
// ============================================================================

impl ReadonlyRefgetStore {
    /// Import a FASTA file into the store using a multithreaded pipeline.
    ///
    /// After the pipeline finishes, computes collection metadata, registers the
    /// collection, and inserts all sequences.
    pub fn add_sequence_collection_from_fasta<P: AsRef<Path>>(
        &mut self,
        file_path: P,
        opts: FastaImportOptions<'_>,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
        if !self.quiet {
            println!("Processing {}...", file_path.as_ref().display());
        }

        let pipeline_start = Instant::now();

        // Check for RGSI cache
        let use_cache = self.local_path.is_some();
        use crate::utils::PathExtension;
        let rgsi_path = file_path.as_ref().replace_exts_with("rgsi");
        let have_rgsi = use_cache && rgsi_path.exists();

        if have_rgsi {
            match self.add_fasta_with_cached_metadata(&file_path, &rgsi_path, opts, pipeline_start) {
                Ok(result) => return Ok(result),
                Err(_) => {
                    // Cache was stale/empty, fall through to pipeline
                }
            }
        }

        let n_workers = resolve_threads(opts.threads);

        // --- Reader item type: a raw record needing work, or a large record
        // already processed inline by the reader. Both carry their input_index.
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

        let file_path_buf = file_path.as_ref().to_path_buf();
        let namespaces: Vec<String> = opts.namespaces.iter().map(|s| s.to_string()).collect();

        let mode = self.mode;
        let quiet = self.quiet;

        let mut sequence_metadata: Vec<SequenceMetadata> = Vec::new();
        let mut all_aliases: Vec<(String, String, String)> = Vec::new();

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

        // --- In-order inserter (main thread): single &mut self writer ---
        run_encode_pipeline(
            n_workers,
            reader,
            work_fn,
            |ready: ReadyOrdered| -> Result<()> {
                for (ns, alias_value) in &ready.aliases {
                    all_aliases.push((
                        ns.clone(),
                        alias_value.clone(),
                        ready.metadata.sha512t24u.clone(),
                    ));
                }
                self.add_sequence_record(
                    SequenceRecord::Full {
                        metadata: ready.metadata.clone(),
                        sequence: ready.sequence_data,
                    },
                    true,
                )?;
                sequence_metadata.push(ready.metadata);
                Ok(())
            },
        )?;

        let pipeline_elapsed = pipeline_start.elapsed();

        // --- Post-pipeline: compute collection metadata, register ---
        let metadata_start = std::time::Instant::now();
        let seq_count = sequence_metadata.len();

        let stub_records: Vec<SequenceRecord> = sequence_metadata
            .iter()
            .map(|meta| SequenceRecord::Stub(meta.clone()))
            .collect();

        let mut seqcol_metadata = SequenceCollectionMetadata::from_sequences(
            &stub_records,
            Some(file_path.as_ref().to_path_buf()),
        );
        if self.ancillary_digests {
            seqcol_metadata.compute_ancillary_digests(&stub_records);
        }

        let coll_key = seqcol_metadata.digest.to_key();
        let coll_digest_display = seqcol_metadata.digest.clone();
        let metadata = seqcol_metadata.clone();
        let metadata_elapsed = metadata_start.elapsed();

        if !opts.force && self.collections.contains_key(&coll_key) {
            if !self.quiet {
                println!("Skipped {} (already exists)", coll_digest_display);
            }
            return Ok((metadata, false));
        }

        let index_start = std::time::Instant::now();

        // Write RGSI cache for next time
        if use_cache {
            let seqcol_for_cache = SequenceCollection {
                metadata: seqcol_metadata.clone(),
                sequences: stub_records.clone(),
            };
            let _ = seqcol_for_cache.write_collection_rgsi(&rgsi_path);
        }

        // Register the collection
        let seqcol = SequenceCollection {
            metadata: seqcol_metadata,
            sequences: stub_records,
        };
        self.add_sequence_collection_internal(seqcol, opts.force)?;

        // Register aliases
        for (ns, alias_value, sha512t24u) in &all_aliases {
            self.add_sequence_alias(ns, alias_value, sha512t24u)?;
        }

        // Populate name_lookup for already-written sequences (FASTA order).
        for meta in &sequence_metadata {
            self.name_lookup
                .entry(coll_key)
                .or_default()
                .insert(meta.name.clone(), meta.sha512t24u.to_key());
        }
        let index_elapsed = index_start.elapsed();

        if !self.quiet {
            let total = pipeline_elapsed + metadata_elapsed + index_elapsed;
            println!(
                "Added {} {} seqs in {:.1}s ({:.1} proc | {:.1} meta | {:.1} index)",
                coll_digest_display,
                seq_count,
                total.as_secs_f64(),
                pipeline_elapsed.as_secs_f64(),
                metadata_elapsed.as_secs_f64(),
                index_elapsed.as_secs_f64(),
            );
        }

        Ok((metadata, true))
    }

    /// Fast path when RGSI cache exists: skip digesting, only decompress for
    /// encoding. Encoding fans out to the same worker pool; the main thread
    /// inserts in FASTA order so `name_lookup` matches a serial build.
    fn add_fasta_with_cached_metadata<P: AsRef<Path>>(
        &mut self,
        file_path: P,
        rgsi_path: &Path,
        opts: FastaImportOptions<'_>,
        start_time: Instant,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
        let mut seqcol = crate::collection::read_rgsi_file(rgsi_path)?;

        if seqcol.sequences.is_empty() {
            let _ = std::fs::remove_file(rgsi_path);
            return Err(anyhow!("Empty RGSI cache"));
        }
        let coll_key = seqcol.metadata.digest.to_key();
        let coll_digest_display = seqcol.metadata.digest.clone();

        if !opts.force && self.collections.contains_key(&coll_key) {
            if !self.quiet {
                println!("Skipped {} (already exists)", coll_digest_display);
            }
            return Ok((seqcol.metadata.clone(), false));
        }

        if self.ancillary_digests {
            seqcol.metadata.compute_ancillary_digests(&seqcol.sequences);
        }
        let metadata = seqcol.metadata.clone();

        let seqmeta_hashmap: HashMap<String, SequenceMetadata> = seqcol
            .sequences
            .iter()
            .map(|r| {
                let meta = r.metadata().clone();
                (meta.name.clone(), meta)
            })
            .collect();

        self.add_sequence_collection_internal(seqcol, opts.force)?;

        let n_workers = resolve_threads(opts.threads);
        let mode = self.mode;
        let namespaces: Vec<String> = opts.namespaces.iter().map(|s| s.to_string()).collect();
        let file_path_buf = file_path.as_ref().to_path_buf();

        // Reader item: a raw record plus its already-known cached metadata.
        struct CachedRaw {
            input_index: u64,
            metadata: SequenceMetadata,
            raw_bytes: Vec<u8>,
            aliases: Vec<(String, String)>,
        }

        let mut seq_count: usize = 0;

        // --- Reader closure ---
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

        // --- Worker closure: encode only (metadata already known) ---
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

        // --- In-order inserter (main thread). Uses add_sequence so each call
        // updates name_lookup in FASTA order. ---
        run_encode_pipeline(
            n_workers,
            reader,
            work_fn,
            |ready: ReadyOrdered| -> Result<()> {
                for (ns, alias_value) in &ready.aliases {
                    self.add_sequence_alias(ns, alias_value, &ready.metadata.sha512t24u)?;
                }
                self.add_sequence(
                    SequenceRecord::Full {
                        metadata: ready.metadata,
                        sequence: ready.sequence_data,
                    },
                    coll_key,
                    true,
                )?;
                seq_count += 1;
                Ok(())
            },
        )?;

        let elapsed = start_time.elapsed();
        if !self.quiet {
            println!(
                "Added {} ({} seqs) in {:.1}s [cached metadata]",
                coll_digest_display, seq_count, elapsed.as_secs_f64()
            );
        }

        Ok((metadata, true))
    }
}
