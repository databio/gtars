//! FASTA import pipeline for RefgetStore.
//!
//! Contains the streaming FASTA import pipeline and cached metadata fast path.
//!
//! ## Single-knob parallelism model
//!
//! There is ONE concurrency knob: `jobs` = the number of input FASTA files
//! imported concurrently. With multiple files, up to `jobs` files are
//! decoded/built concurrently, each with its OWN decoder (so K gzip streams
//! decompress in parallel). With a single file `jobs` has no effect.
//!
//! ## Streaming, bounded-memory design
//!
//! Each builder emits sequences ONE AT A TIME over a bounded cross-file channel
//! rather than accumulating a whole collection in RAM. The single inserter
//! (owner, `&mut self`) writes each `.seq` to disk immediately as it arrives and
//! retains only lightweight per-sequence METADATA (name + digest), so peak
//! memory is bounded by `jobs x channel_depth x avg_encoded_seq_size` plus the
//! on-disk index metadata, INDEPENDENT of how many sequences a collection holds.
//!
//! ### Per-file pipeline (basic chain)
//!
//! Each file is processed by a simple linear pipeline of ~3 threads with
//! `bounded(1)` channels providing back-pressure (a couple of contigs in flight):
//! a decompress/read thread, a digest thread, and an encode step on the build
//! thread. Records flow through the chain in FASTA order; the build thread
//! forwards each finished sequence onto the cross-file output channel instead of
//! pushing into a `Vec`. This keeps `name_lookup` (and therefore the collection
//! digest) deterministic and identical to a single-threaded build.
//!
//! ### File-level parallelism (across files)
//!
//! The build half (decode + parse + digest + encode) touches NO shared store
//! state and runs concurrently across files; the insert half (`&mut self`
//! registration of sequences/collections/aliases/name_lookup) runs on a single
//! owner. Because the on-disk indexes are written once at the end, sorted by
//! content/collection digest, the resulting store is byte-identical to a serial
//! build regardless of which builder finishes first or how `Seq`/`End` messages
//! from different files interleave.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::thread::available_parallelism;
use std::time::Instant;

use anyhow::{anyhow, Result};
use crossbeam_channel::{bounded, Sender};
use indexmap::IndexMap;

use crate::collection::SequenceCollectionExt;
use crate::digest::{
    SequenceCollection, SequenceCollectionMetadata,
    SequenceCollectionRecord, SequenceEncoder, SequenceMetadata, SequenceRecord,
};
use crate::hashkeyable::{DigestKey, HashKeyable};

/// Threshold above which the reader processes sequences inline (to bound peak
/// memory: the few huge contigs are hashed+encoded on the reader instead of
/// being passed down the pipeline holding extra copies in flight).
const LARGE_SEQ_THRESHOLD: usize = 500 * 1024 * 1024; // 500 MB

/// Per-builder output channel depth (multiplied by `jobs` for the shared
/// cross-file channel). Small constant: a couple sequences in flight per builder
/// is enough to keep builders and the inserter busy without unbounded buffering.
const CHANNEL_DEPTH: usize = 4;

/// A fully-built sequence ready for the single inserter, in FASTA order.
///
/// Carries the encoded bytes ONLY from emission until the inserter writes them
/// to disk; the inserter then drops the bytes and keeps only metadata.
struct ReadySequence {
    metadata: SequenceMetadata,
    sequence_data: Vec<u8>,
    aliases: Vec<(String, String)>,
}

/// One message in the streaming import protocol, consumed by the single inserter.
///
/// A builder emits `Begin`, then one `Seq` per record in FASTA order, then `End`.
/// Messages from different files may interleave on the shared channel; the
/// `file_idx` tag disambiguates which in-flight collection each `Seq` belongs to.
enum InsertMsg {
    /// A builder started processing this file.
    Begin { file_idx: usize },
    /// The collection was found already present in the store (cached `.rgsi`
    /// short-circuit). No `Begin`/`Seq`/`End` is emitted for this file: the FASTA
    /// is never opened or decoded. The inserter records `was_new = false`.
    Skip {
        file_idx: usize,
        metadata: SequenceCollectionMetadata,
    },
    /// One fully digested+encoded sequence, in FASTA order within its file.
    Seq {
        file_idx: usize,
        ready: ReadySequence,
    },
    /// The collection finished: metadata computed from the retained per-sequence
    /// stubs (NOT from bytes), plus the RGSI cache path (raw builds only).
    End {
        file_idx: usize,
        metadata: SequenceCollectionMetadata,
        stub_records: Vec<SequenceRecord>,
        rgsi_cache_path: Option<PathBuf>,
        source_path: PathBuf,
        seq_count: usize,
    },
}

// ============================================================================
// Build half (no &mut self): decode + parse + digest + encode -> stream of msgs
//
// These free functions run concurrently across files. They touch NO shared
// store state; they emit `InsertMsg`s onto a bounded channel. The inserter on a
// single owner thread persists them.
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
    /// Re-process collections even if already present (disables the cached
    /// `.rgsi` pre-decode short-circuit).
    force: bool,
    /// Read-only snapshot of collection digest keys already present in the store,
    /// taken once before builders spawn. Lets the build half (which has no
    /// `&self`) short-circuit an already-present collection BEFORE opening the
    /// FASTA, avoiding a wasteful decode+encode on resume.
    present_collections: &'a HashSet<DigestKey>,
}

/// Build half for a single FASTA file: digest + encode, streaming each sequence
/// onto `out` as `InsertMsg::Seq`, framed by `Begin`/`End`.
///
/// Checks for an RGSI metadata cache first and dispatches to the cached fast
/// path if present and valid; otherwise runs the full digest+encode pipeline.
/// Does NOT touch `&mut self`.
fn build_collection_from_fasta(
    file_idx: usize,
    file_path: &Path,
    cfg: BuildConfig<'_>,
    out: &Sender<InsertMsg>,
) -> Result<()> {
    use crate::utils::PathExtension;

    let rgsi_path = file_path.replace_exts_with("rgsi");
    let have_rgsi = cfg.use_cache && rgsi_path.exists();

    if have_rgsi {
        match build_collection_from_cached_metadata(file_idx, file_path, &rgsi_path, cfg, out) {
            Ok(()) => return Ok(()),
            Err(BuildError::Channel(e)) => return Err(BuildError::Channel(e).into()),
            // Cache was stale/empty; fall through to the full pipeline. Note: a
            // partial cached run cannot have emitted `Seq`s (the cache validity
            // is checked before any `Begin`/`Seq` is sent), so re-running the
            // full pipeline is safe.
            Err(BuildError::Other(_)) => {}
        }
    }

    build_collection_full(file_idx, file_path, cfg, out)
}

/// Errors from the build half. `Channel` means the inserter hung up (fatal);
/// `Other` means a recoverable build error (e.g. stale cache) the caller may
/// choose to handle by falling through to another path.
enum BuildError {
    Channel(anyhow::Error),
    Other(anyhow::Error),
}

impl From<BuildError> for anyhow::Error {
    fn from(e: BuildError) -> Self {
        match e {
            BuildError::Channel(e) | BuildError::Other(e) => e,
        }
    }
}

/// Full build half: digest + encode every sequence and stream it to the inserter.
///
/// Uses the basic 3-thread chain: a reader/decompress thread, a digest thread,
/// and the encode step on this (build) thread. `bounded(1)` channels keep at
/// most a couple of contigs in flight, capping per-builder RAM. Records flow in
/// FASTA order, and each finished sequence is forwarded onto `out` (the shared
/// cross-file channel) instead of being collected into a Vec. The builder
/// retains ONLY the lightweight stub metadata to compute the collection digest
/// at `End` -- never the encoded bytes.
fn build_collection_full(
    file_idx: usize,
    file_path: &Path,
    cfg: BuildConfig<'_>,
    out: &Sender<InsertMsg>,
) -> Result<()> {
    use crate::utils::PathExtension;
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

    // --- Thread 3 (this thread): encode + stream in FASTA order ---
    // Emit Begin, then a Seq per record. Retain ONLY stub metadata + alias
    // triples (no bytes) to compute the collection digest at End.
    out.send(InsertMsg::Begin { file_idx })
        .map_err(|_| anyhow!("Inserter stopped receiving"))?;

    let mut stub_records: Vec<SequenceRecord> = Vec::new();

    for msg in digest_rx {
        let ready = match msg {
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
                ReadySequence {
                    metadata: digested.metadata,
                    sequence_data,
                    aliases: digested.aliases,
                }
            }
            ToEncode::AlreadyDone(ready) => ready,
        };

        stub_records.push(SequenceRecord::Stub(ready.metadata.clone()));
        out.send(InsertMsg::Seq { file_idx, ready })
            .map_err(|_| anyhow!("Inserter stopped receiving"))?;
    }

    // Join threads and propagate errors.
    decompress_handle
        .join()
        .map_err(|e| anyhow!("Decompress thread panicked: {:?}", e))??;
    digest_handle
        .join()
        .map_err(|e| anyhow!("Digest thread panicked: {:?}", e))??;

    // Compute collection metadata from the retained sequence stubs (FASTA order).
    let mut seqcol_metadata = SequenceCollectionMetadata::from_sequences(
        &stub_records,
        Some(file_path.to_path_buf()),
    );
    if cfg.ancillary_digests {
        seqcol_metadata.compute_ancillary_digests(&stub_records);
    }

    let seq_count = stub_records.len();
    let rgsi_cache_path = if cfg.use_cache {
        Some(file_path.replace_exts_with("rgsi"))
    } else {
        None
    };

    out.send(InsertMsg::End {
        file_idx,
        metadata: seqcol_metadata,
        stub_records,
        rgsi_cache_path,
        source_path: file_path.to_path_buf(),
        seq_count,
    })
    .map_err(|_| anyhow!("Inserter stopped receiving"))?;

    Ok(())
}

/// Build half for the RGSI cached-metadata fast path: skip digesting, only
/// decompress + encode, streaming each sequence to the inserter.
///
/// Uses a basic 2-thread split: a reader/decompress thread feeds raw records
/// (with their already-known cached metadata) to this thread, which encodes them
/// in FASTA order and forwards them. Returns `BuildError::Other` if the cache is
/// invalid (caller falls through to the full pipeline) -- this is checked BEFORE
/// any `Begin`/`Seq` is emitted.
fn build_collection_from_cached_metadata(
    file_idx: usize,
    file_path: &Path,
    rgsi_path: &Path,
    cfg: BuildConfig<'_>,
    out: &Sender<InsertMsg>,
) -> Result<(), BuildError> {
    let mut seqcol = crate::collection::read_rgsi_file(rgsi_path)
        .map_err(BuildError::Other)?;

    if seqcol.sequences.is_empty() {
        let _ = std::fs::remove_file(rgsi_path);
        return Err(BuildError::Other(anyhow!("Empty RGSI cache")));
    }

    // Pre-decode short-circuit: if this collection is already present in the
    // store (and we're not forcing), skip it WITHOUT opening/decoding the FASTA.
    // The `.rgsi` sidecar gave us the collection digest for free; re-decoding +
    // re-encoding every record only to discover the duplicate at `End` wastes
    // enormous CPU on resume. Emit a `Skip` so the inserter records
    // `was_new = false`, then return before spawning the reader thread.
    if !cfg.force && cfg.present_collections.contains(&seqcol.metadata.digest.to_key()) {
        if !cfg.quiet {
            println!(
                "Skipped {} (already exists) [cached metadata]",
                seqcol.metadata.digest
            );
        }
        out.send(InsertMsg::Skip {
            file_idx,
            metadata: seqcol.metadata.clone(),
        })
        .map_err(|_| BuildError::Channel(anyhow!("Inserter stopped receiving")))?;
        return Ok(());
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

    // The cache is valid; from here on we emit Begin/Seq/End. A channel error is
    // fatal (inserter hung up), not a recoverable cache miss.
    out.send(InsertMsg::Begin { file_idx })
        .map_err(|_| BuildError::Channel(anyhow!("Inserter stopped receiving")))?;

    // --- Thread 2 (this thread): encode + stream in FASTA order ---
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
        let ready = ReadySequence {
            metadata: unit.metadata,
            sequence_data,
            aliases: unit.aliases,
        };
        out.send(InsertMsg::Seq { file_idx, ready })
            .map_err(|_| BuildError::Channel(anyhow!("Inserter stopped receiving")))?;
    }

    read_handle
        .join()
        .map_err(|e| BuildError::Other(anyhow!("Decompress thread panicked: {:?}", e)))?
        .map_err(BuildError::Other)?;

    let stub_records = seqcol.sequences.clone();
    let metadata = seqcol.metadata.clone();
    let seq_count = stub_records.len();

    out.send(InsertMsg::End {
        file_idx,
        metadata,
        stub_records,
        // Cache already exists and is valid; no need to rewrite it.
        rgsi_cache_path: None,
        source_path: file_path.to_path_buf(),
        seq_count,
    })
    .map_err(|_| BuildError::Channel(anyhow!("Inserter stopped receiving")))?;

    Ok(())
}

// ============================================================================
// Per-in-flight collection scratch held by the inserter
// ============================================================================

/// Lightweight per-in-flight-collection state held by the single inserter
/// between `Begin` and `End`. Metadata-sized only -- never holds encoded bytes
/// (those are written to disk as each `Seq` arrives, then dropped). At most
/// `jobs` of these are live at once.
struct InFlight {
    /// Ordered (name -> sha512 digest key) in FASTA order; becomes name_lookup.
    name_to_digest: IndexMap<String, DigestKey>,
    /// (namespace, alias_value, sha512t24u) alias triples in FASTA order.
    aliases: Vec<(String, String, String)>,
}

impl InFlight {
    fn new() -> Self {
        InFlight {
            name_to_digest: IndexMap::new(),
            aliases: Vec::new(),
        }
    }
}

// ============================================================================
// ReadonlyRefgetStore import methods
// ============================================================================

impl ReadonlyRefgetStore {
    /// Import multiple FASTA files into the store with file-level parallelism.
    ///
    /// Up to `jobs` files are decoded/built concurrently (each with its own
    /// decoder, so gzip decompression is parallelized across files). Each builder
    /// STREAMS its sequences one at a time over a bounded cross-file channel; the
    /// single inserter (`&mut self`) on this owner thread persists each `.seq`
    /// immediately and retains only metadata, so peak memory is bounded by
    /// `jobs x channel_depth x avg_encoded_seq_size`, independent of collection
    /// size. The on-disk indexes are written once at the very end (sorted by
    /// content/collection digest), guaranteeing a byte-identical store to a serial
    /// build regardless of build/arrival order.
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

        // Snapshot the digest keys of collections already present in the store,
        // taken ONCE here (before any builder spawns). The build half has no
        // `&self`, so this read-only set is how it learns which collections can
        // be skipped before opening their FASTAs. New collections added during
        // this run are NOT in the snapshot, so they are never wrongly skipped;
        // duplicates within the same batch are still deduped by the inserter.
        let present_collections: HashSet<DigestKey> = self.collections.keys().copied().collect();

        let cfg = BuildConfig {
            mode: self.mode,
            quiet: self.quiet,
            ancillary_digests: self.ancillary_digests,
            use_cache: self.local_path.is_some(),
            namespaces: opts.namespaces,
            force: opts.force,
            present_collections: &present_collections,
        };

        let overall_start = Instant::now();

        // Per-file results, filled in by the inserter as each `End` arrives, then
        // returned in input order.
        let mut results: Vec<Option<(SequenceCollectionMetadata, bool)>> =
            (0..files.len()).map(|_| None).collect();

        // Bounded cross-file output channel: builders block on `send` when it is
        // full, throttling fast builders to the inserter's pace (back-pressure)
        // and capping in-flight encoded bytes.
        let (out_tx, out_rx) = bounded::<InsertMsg>(jobs * CHANNEL_DEPTH);

        // Bounded job queue feeds a fixed pool of `jobs` builder threads.
        let (job_tx, job_rx) = bounded::<(usize, PathBuf)>(jobs);

        // The first builder error encountered (build threads stash it here).
        let build_err: std::sync::Mutex<Option<anyhow::Error>> = std::sync::Mutex::new(None);
        let build_err_ref = &build_err;

        std::thread::scope(|scope| -> Result<()> {
            // Spawn the builder pool.
            for _ in 0..jobs {
                let job_rx = job_rx.clone();
                let out_tx = out_tx.clone();
                scope.spawn(move || {
                    for (idx, path) in job_rx.iter() {
                        if let Err(e) = build_collection_from_fasta(idx, &path, cfg, &out_tx) {
                            let mut slot = build_err_ref.lock().unwrap();
                            if slot.is_none() {
                                *slot = Some(e);
                            }
                            // Stop pulling more jobs; let the inserter drain and
                            // the scope tear down.
                            break;
                        }
                    }
                });
            }
            // Drop our extra senders so the inserter's loop terminates once all
            // builders finish.
            drop(out_tx);
            drop(job_rx);

            // Feeder thread: push all jobs, then drop the sender. Runs in the
            // scope so it can proceed while the inserter (this thread) drains the
            // output channel -- otherwise a bounded job queue + bounded output
            // channel could deadlock.
            let feeder = {
                let job_tx = job_tx.clone();
                let files = files.to_vec();
                let quiet = self.quiet;
                scope.spawn(move || {
                    for (idx, file) in files.into_iter().enumerate() {
                        if !quiet {
                            println!("Processing {}...", file.display());
                        }
                        if job_tx.send((idx, file)).is_err() {
                            break;
                        }
                    }
                })
            };
            drop(job_tx);

            // --- Single inserter: drain the streaming channel. ---
            // Per-in-flight-collection scratch keyed by file_idx (metadata-sized,
            // at most `jobs` live).
            let mut in_flight: HashMap<usize, InFlight> = HashMap::new();

            // For each unique sequence digest, the input file_idx whose name
            // currently "owns" the global `sequence_store` metadata entry. The
            // global RGSI row for a shared digest carries one (arbitrary but
            // DETERMINISTIC) name; we make the highest-file-index name win,
            // matching the serial build's fixed-file-order, last-writer-wins
            // behavior so parallel == serial byte-identical. (Per-collection
            // names are always correct via the per-collection name_lookup.)
            let mut seq_name_owner: HashMap<DigestKey, usize> = HashMap::new();

            for msg in out_rx.iter() {
                match msg {
                    InsertMsg::Begin { file_idx } => {
                        in_flight.insert(file_idx, InFlight::new());
                    }
                    InsertMsg::Skip { file_idx, metadata } => {
                        // Already-present collection skipped before any decode.
                        // No scratch was created (no Begin), nothing to persist.
                        results[file_idx] = Some((metadata, false));
                    }
                    InsertMsg::Seq { file_idx, ready } => {
                        let ReadySequence { metadata, sequence_data, aliases } = ready;
                        let seq_key = metadata.sha512t24u.to_key();

                        // Record name -> digest in FASTA order for this file's
                        // name_lookup, and stash alias triples.
                        let scratch = in_flight
                            .get_mut(&file_idx)
                            .expect("Seq before Begin");
                        scratch
                            .name_to_digest
                            .insert(metadata.name.clone(), seq_key);
                        for (ns, alias_value) in &aliases {
                            scratch.aliases.push((
                                ns.clone(),
                                alias_value.clone(),
                                metadata.sha512t24u.clone(),
                            ));
                        }

                        // Persist the bytes immediately (write-through + dedup),
                        // then drop them. Dedup is keyed on the content digest.
                        // The first occurrence of a digest always writes the
                        // `.seq` and its metadata. For a sequence shared across
                        // collections we write again ONLY when this file's index
                        // is higher than the current name owner, so the
                        // highest-file-index name wins deterministically
                        // (matching the serial fixed-file-order build). The bytes
                        // are identical (same digest), so the re-write is a no-op
                        // for `.seq`; only the stored name metadata changes. For
                        // disk-backed stores the record is downgraded to a Stub
                        // (no bytes retained); for in-memory stores the Full
                        // record is retained (memory is inherently O(total bytes)
                        // there).
                        let force = match seq_name_owner.get(&seq_key) {
                            None => {
                                seq_name_owner.insert(seq_key, file_idx);
                                false // first occurrence: insert unconditionally
                            }
                            Some(&owner) if file_idx > owner => {
                                seq_name_owner.insert(seq_key, file_idx);
                                true // higher file index: overwrite stored name
                            }
                            Some(_) => {
                                // Lower-or-equal file index already owns the name;
                                // skip (still deduped, no re-write).
                                continue;
                            }
                        };
                        self.add_sequence_record(
                            SequenceRecord::Full {
                                metadata,
                                sequence: sequence_data,
                            },
                            force,
                        )?;
                    }
                    InsertMsg::End {
                        file_idx,
                        metadata,
                        stub_records,
                        rgsi_cache_path,
                        source_path,
                        seq_count,
                    } => {
                        let scratch = in_flight
                            .remove(&file_idx)
                            .expect("End before Begin");
                        let (meta, was_new) = self.finalize_collection(
                            metadata,
                            stub_records,
                            scratch,
                            rgsi_cache_path,
                            &source_path,
                            seq_count,
                            opts.force,
                        )?;
                        results[file_idx] = Some((meta, was_new));
                    }
                }
            }

            feeder.join().map_err(|e| anyhow!("Feeder thread panicked: {:?}", e))?;
            Ok(())
        })?;

        // Propagate the first builder error, if any.
        if let Some(e) = build_err.into_inner().unwrap() {
            return Err(e);
        }

        // Finalize ALL indexes ONCE, after every collection is persisted. This
        // is what makes parallel == serial byte-identical: sequences.rgsi sorts
        // by sha512t24u and collections.rgci sorts by collection digest, so the
        // arrival/build order is irrelevant.
        if self.persist_to_disk && self.local_path.is_some() {
            self.write_index_files()?;
        }

        if !self.quiet {
            println!(
                "Imported {} file(s) in {:.1}s (jobs={})",
                files.len(),
                overall_start.elapsed().as_secs_f64(),
                jobs,
            );
        }

        let results: Vec<(SequenceCollectionMetadata, bool)> = results
            .into_iter()
            .map(|slot| slot.expect("every file index must be finalized"))
            .collect();

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

    /// Finalize one collection at `End`: register the collection record, install
    /// the buffered name_lookup (FASTA order) and aliases. Sequence bytes were
    /// already written to disk as each `Seq` arrived. Does NOT write the global
    /// index files (the import owner does that ONCE at the end).
    #[allow(clippy::too_many_arguments)]
    fn finalize_collection(
        &mut self,
        metadata: SequenceCollectionMetadata,
        stub_records: Vec<SequenceRecord>,
        scratch: InFlight,
        rgsi_cache_path: Option<PathBuf>,
        source_path: &Path,
        seq_count: usize,
        force: bool,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
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

        // Register the collection record (stub sequences only).
        let record = SequenceCollectionRecord::Full {
            metadata: metadata.clone(),
            sequences: stub_records,
        };
        if self.persist_to_disk && self.local_path.is_some() {
            self.write_collection_to_disk_single(&record)?;
        }
        self.collections.insert(coll_key, record);

        // Install the buffered name_lookup in FASTA order.
        self.name_lookup.insert(coll_key, scratch.name_to_digest);

        // Register aliases in FASTA order.
        for (ns, alias_value, sha512t24u) in &scratch.aliases {
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
