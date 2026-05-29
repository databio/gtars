//! FASTA import pipeline for RefgetStore.
//!
//! Contains the FASTA import pipeline and cached metadata fast path.
//!
//! ## Single-knob parallelism model
//!
//! There is ONE concurrency knob: `jobs` = the number of input FASTA files
//! imported concurrently. With multiple files, up to `jobs` files are
//! decoded/built concurrently, each with its OWN decoder (so K gzip streams
//! decompress in parallel). With a single file `jobs` has no effect.
//!
//! ### Per-file pipeline (basic chain)
//!
//! Each file is processed by a simple linear pipeline of ~3 threads with
//! `bounded(1)` channels providing back-pressure (a couple of contigs in flight):
//! a decompress/read thread, a digest thread, and an encode step on the build
//! thread. Records flow through the chain in FASTA order, so output is naturally
//! in FASTA order without any reordering machinery. This makes `name_lookup`
//! (and therefore the collection digest) deterministic and identical to a
//! single-threaded build.
//!
//! ### File-level parallelism (across files)
//!
//! The build half (decode + parse + digest + encode -> in-memory
//! `BuiltCollection`) touches NO shared store state and runs concurrently across
//! files; the insert half (`&mut self` registration of
//! collections/sequences/aliases/name_lookup) runs on a single owner in FIXED
//! input-file order, so the resulting store is byte-identical to a serial build
//! regardless of which builder finishes first.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::thread::available_parallelism;
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
/// being passed down the pipeline holding extra copies in flight).
const LARGE_SEQ_THRESHOLD: usize = 500 * 1024 * 1024; // 500 MB

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
    sequences: Vec<ReadySequence>,
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

/// A fully-built sequence ready for the single inserter, in FASTA order.
struct ReadySequence {
    metadata: SequenceMetadata,
    sequence_data: Vec<u8>,
    aliases: Vec<(String, String)>,
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
///
/// Uses the basic 3-thread chain: a reader/decompress thread, a digest thread,
/// and the encode step on this (build) thread. `bounded(1)` channels keep at
/// most a couple of contigs in flight, capping RAM. Records flow in FASTA order.
fn build_collection_full(file_path: &Path, cfg: BuildConfig<'_>) -> Result<BuiltCollection> {
    use crate::utils::PathExtension;
    use crossbeam_channel::bounded;
    use md5::Md5;
    use sha2::{Digest, Sha512};

    // --- Channel message types ---
    struct DecompressedSequence {
        name: String,
        description: Option<String>,
        raw_header: String,
        raw_bytes: Vec<u8>,
    }
    struct DigestedSequence {
        metadata: SequenceMetadata,
        raw_bytes: Vec<u8>,
        aliases: Vec<(String, String)>,
    }
    enum ToDigest {
        NeedsWork(DecompressedSequence),
        AlreadyDone(ReadySequence),
    }
    enum ToEncode {
        NeedsWork(DigestedSequence),
        AlreadyDone(ReadySequence),
    }

    let (decompress_tx, decompress_rx) = bounded::<ToDigest>(1);
    let (digest_tx, digest_rx) = bounded::<ToEncode>(1);

    let file_path_buf = file_path.to_path_buf();
    let namespaces: Vec<String> = cfg.namespaces.iter().map(|s| s.to_string()).collect();
    let ns_for_digest = namespaces.clone();

    let quiet = cfg.quiet;
    let mode = cfg.mode;

    // --- Thread 1: read FASTA (decompress) ---
    let decompress_handle = std::thread::spawn(move || -> Result<()> {
        let mut fasta_reader = crate::fasta::FastaReader::from_path(&file_path_buf)?;

        while let Some(record) = fasta_reader.next_record()? {
            if record.raw_bytes.len() > LARGE_SEQ_THRESHOLD {
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

                decompress_tx
                    .send(ToDigest::AlreadyDone(ReadySequence {
                        metadata,
                        sequence_data,
                        aliases,
                    }))
                    .map_err(|_| anyhow!("Digest thread stopped receiving"))?;
            } else {
                decompress_tx
                    .send(ToDigest::NeedsWork(DecompressedSequence {
                        name: record.name,
                        description: record.description,
                        raw_header: record.raw_header,
                        raw_bytes: record.raw_bytes,
                    }))
                    .map_err(|_| anyhow!("Digest thread stopped receiving"))?;
            }
        }
        Ok(())
    });

    // --- Thread 2: digest ---
    let digest_handle = std::thread::spawn(move || -> Result<()> {
        let ns_refs: Vec<&str> = ns_for_digest.iter().map(|s| s.as_str()).collect();

        for msg in decompress_rx {
            match msg {
                ToDigest::NeedsWork(seq) => {
                    let mut sha512_hasher = Sha512::new();
                    sha512_hasher.update(&seq.raw_bytes);
                    let sha512 = base64_url::encode(&sha512_hasher.finalize()[0..24]);

                    let mut md5_hasher = Md5::new();
                    md5_hasher.update(&seq.raw_bytes);
                    let md5 = format!("{:x}", md5_hasher.finalize());

                    let mut guesser = crate::digest::AlphabetGuesser::new();
                    guesser.update(&seq.raw_bytes);
                    let alphabet = guesser.guess();

                    let metadata = SequenceMetadata {
                        name: seq.name,
                        description: seq.description,
                        length: seq.raw_bytes.len(),
                        sha512t24u: sha512,
                        md5,
                        alphabet,
                        fai: None,
                    };

                    let aliases = if !ns_refs.is_empty() {
                        crate::digest::fasta::extract_aliases_from_header(&seq.raw_header, &ns_refs)
                    } else {
                        vec![]
                    };

                    digest_tx
                        .send(ToEncode::NeedsWork(DigestedSequence {
                            metadata,
                            raw_bytes: seq.raw_bytes,
                            aliases,
                        }))
                        .map_err(|_| anyhow!("Encode thread stopped receiving"))?;
                }
                ToDigest::AlreadyDone(ready) => {
                    digest_tx
                        .send(ToEncode::AlreadyDone(ready))
                        .map_err(|_| anyhow!("Encode thread stopped receiving"))?;
                }
            }
        }
        drop(digest_tx);
        Ok(())
    });

    // --- Thread 3 (this thread): encode + collect in FASTA order ---
    let mut sequences: Vec<ReadySequence> = Vec::new();

    for msg in digest_rx {
        match msg {
            ToEncode::NeedsWork(digested) => {
                let sequence_data = match mode {
                    StorageMode::Encoded => {
                        let mut encoder =
                            SequenceEncoder::new(digested.metadata.alphabet, digested.metadata.length);
                        encoder.update(&digested.raw_bytes);
                        drop(digested.raw_bytes);
                        encoder.finalize()
                    }
                    StorageMode::Raw => digested.raw_bytes,
                };
                sequences.push(ReadySequence {
                    metadata: digested.metadata,
                    sequence_data,
                    aliases: digested.aliases,
                });
            }
            ToEncode::AlreadyDone(ready) => {
                sequences.push(ready);
            }
        }
    }

    // Join threads and propagate errors.
    decompress_handle
        .join()
        .map_err(|e| anyhow!("Decompress thread panicked: {:?}", e))??;
    digest_handle
        .join()
        .map_err(|e| anyhow!("Digest thread panicked: {:?}", e))??;

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
///
/// Uses a basic 2-thread split: a reader/decompress thread feeds raw records
/// (with their already-known cached metadata) to this thread, which encodes them
/// in FASTA order.
fn build_collection_from_cached_metadata(
    file_path: &Path,
    rgsi_path: &Path,
    cfg: BuildConfig<'_>,
) -> Result<BuiltCollection> {
    use crossbeam_channel::bounded;

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
        metadata: SequenceMetadata,
        raw_bytes: Vec<u8>,
        aliases: Vec<(String, String)>,
    }

    let (read_tx, read_rx) = bounded::<CachedRaw>(1);

    // --- Thread 1: read FASTA + look up cached metadata ---
    let read_handle = std::thread::spawn(move || -> Result<()> {
        let ns_refs: Vec<&str> = namespaces.iter().map(|s| s.as_str()).collect();
        let mut fasta_reader = crate::fasta::FastaReader::from_path(&file_path_buf)?;

        while let Some(record) = fasta_reader.next_record()? {
            let (name, _) = crate::fasta::parse_fasta_header(&record.raw_header);

            let aliases = if !ns_refs.is_empty() {
                crate::digest::fasta::extract_aliases_from_header(&record.raw_header, &ns_refs)
            } else {
                vec![]
            };

            let dr = seqmeta_hashmap
                .get(&name)
                .ok_or_else(|| anyhow!("Sequence '{}' not found in cached metadata", name))?
                .clone();

            read_tx
                .send(CachedRaw {
                    metadata: dr,
                    raw_bytes: record.raw_bytes,
                    aliases,
                })
                .map_err(|_| anyhow!("Encode thread stopped receiving"))?;
        }
        Ok(())
    });

    // --- Thread 2 (this thread): encode + collect in FASTA order ---
    let mut sequences: Vec<ReadySequence> = Vec::new();
    for unit in read_rx {
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
        sequences.push(ReadySequence {
            metadata: unit.metadata,
            sequence_data,
            aliases: unit.aliases,
        });
    }

    read_handle
        .join()
        .map_err(|e| anyhow!("Decompress thread panicked: {:?}", e))??;

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
    /// Up to `jobs` files are decoded/built concurrently (each with its own
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

        // Resolve concurrency: 0 = auto (available_parallelism). Has no effect
        // with a single file (clamped to 1).
        let jobs = match opts.jobs {
            0 => available_parallelism().map(|n| n.get()).unwrap_or(1),
            n => n,
        };
        let jobs = jobs.min(files.len()).max(1);

        let cfg = BuildConfig {
            mode: self.mode,
            quiet: self.quiet,
            ancillary_digests: self.ancillary_digests,
            use_cache: self.local_path.is_some(),
            namespaces: opts.namespaces,
        };

        let overall_start = Instant::now();

        // --- Build half: run up to `jobs` builders concurrently. ---
        // Results are streamed back tagged with the file index; the owner
        // collects them and inserts in fixed file index order.
        let mut built: Vec<Option<Result<BuiltCollection>>> =
            (0..files.len()).map(|_| None).collect();

        if jobs <= 1 {
            // Serial build (one file at a time, basic pipeline per file).
            for (idx, file) in files.iter().enumerate() {
                if !self.quiet {
                    println!("Processing {}...", file.display());
                }
                built[idx] = Some(build_collection_from_fasta(file, cfg));
            }
        } else {
            use crossbeam_channel::{bounded, unbounded};
            // Bounded job queue feeds a fixed pool of `jobs` builder threads;
            // results stream back unbounded (each result is one BuiltCollection,
            // bounded in count by the number of files).
            let (job_tx, job_rx) = bounded::<(usize, PathBuf)>(jobs);
            let (res_tx, res_rx) = unbounded::<(usize, Result<BuiltCollection>)>();

            std::thread::scope(|scope| {
                for _ in 0..jobs {
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
                "Imported {} file(s) in {:.1}s (jobs={})",
                files.len(),
                overall_start.elapsed().as_secs_f64(),
                jobs,
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
            let ReadySequence { metadata: seq_meta, sequence_data, .. } = ready;
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
