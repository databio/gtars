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
use std::sync::Mutex;
use std::thread::available_parallelism;
use std::time::Instant;

use anyhow::{anyhow, Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
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

/// Bounded depth of the inserter -> writer-pool channel, per writer thread. The
/// channel capacity is `K * WRITER_CHANNEL_DEPTH` so a few `.seq` writes can be
/// queued ahead of the writers without unbounded buffering. Keeping this a small
/// constant preserves the bounded-memory property (peak RAM stays O(K x channel
/// depth x avg encoded seq size), independent of collection size).
const WRITER_CHANNEL_DEPTH: usize = 8;

/// A unit of work for the writer pool. The inserter has already committed the
/// dedup decision + in-memory metadata (and, for packs, the deterministic
/// `(ordinal, offset)` assignment); only the byte write is offloaded here.
///
/// `Standalone` writes a per-digest `.seq` file (race-free by content address).
/// `PackChunk` does a positional write into a disjoint range of a temp pack
/// file — also race-free, since the inserter chose non-overlapping ranges.
enum WriteJob {
    Standalone { full_path: PathBuf, bytes: Vec<u8> },
    PackChunk { ordinal: usize, offset: u64, bytes: Vec<u8> },
}

/// A pending sidecar row recorded by the inserter for one new packed sequence,
/// carrying the pack ORDINAL (mapped to the final content-hash name at seal).
struct PendingPackRow {
    key: DigestKey,
    ordinal: usize,
    offset: u64,
    len: u64,
}

/// The currently-open pack being filled by the inserter.
struct OpenPack {
    ordinal: usize,
    len: u64,
    hasher: sha2::Sha256,
}

/// A sealed pack: its final 32-hex content-hash name and byte length.
#[derive(Clone)]
struct SealedPack {
    ordinal: usize,
    #[allow(dead_code)]
    len: u64,
    hash32: String,
}

/// Deterministic pack assignment, run ONLY on the single inserter thread.
///
/// `assign` returns `(ordinal, offset)` for a sequence's bytes in arrival order
/// — the same order the existing machinery already makes deterministic — so
/// parallel and serial imports produce identical packs (names + content),
/// identical offsets, and an identical sidecar. Writers only ever do positional
/// writes into ranges this chooses; writer completion order affects no byte.
struct PackAssigner {
    cap: u64,
    next_ordinal: usize,
    open: Option<OpenPack>,
    sealed: Vec<SealedPack>,
}

impl PackAssigner {
    fn new(cap: u64) -> Self {
        PackAssigner {
            cap: cap.max(1),
            next_ordinal: 0,
            open: None,
            sealed: Vec::new(),
        }
    }

    /// Seal-before-place: if the open pack is non-empty and cannot fit `bytes`,
    /// seal it first. A sequence never straddles two packs; one larger than the
    /// cap lands alone at offset 0 of its own pack.
    fn assign(&mut self, bytes: &[u8]) -> (usize, u64) {
        use sha2::Digest;
        let len = bytes.len() as u64;
        if matches!(&self.open, Some(o) if o.len > 0 && o.len + len > self.cap) {
            self.seal_open();
        }
        if self.open.is_none() {
            self.open = Some(OpenPack {
                ordinal: self.next_ordinal,
                len: 0,
                hasher: sha2::Sha256::new(),
            });
            self.next_ordinal += 1;
        }
        let open = self.open.as_mut().unwrap();
        let offset = open.len;
        // Assignment order == byte order, so the content hash costs no extra I/O.
        open.hasher.update(bytes);
        open.len += len;
        (open.ordinal, offset)
    }

    fn seal_open(&mut self) {
        use sha2::Digest;
        if let Some(open) = self.open.take() {
            let digest = open.hasher.finalize();
            // 32-hex-char prefix of the sha256 (16 bytes).
            let hash32 = digest[..16].iter().map(|b| format!("{:02x}", b)).collect();
            self.sealed.push(SealedPack {
                ordinal: open.ordinal,
                len: open.len,
                hash32,
            });
        }
    }

    fn finish(&mut self) -> Vec<SealedPack> {
        self.seal_open();
        std::mem::take(&mut self.sealed)
    }
}

/// Shared temp-pack file handles for the writer pool, keyed by pack ordinal.
/// Positional writes into disjoint ranges of the same file are race-free, so a
/// single shared handle per ordinal is enough.
struct PackWriters {
    packs_dir: PathBuf,
    pid: u32,
    files: Mutex<HashMap<usize, std::sync::Arc<std::fs::File>>>,
}

impl PackWriters {
    fn new(packs_dir: PathBuf) -> Self {
        PackWriters {
            packs_dir,
            pid: std::process::id(),
            files: Mutex::new(HashMap::new()),
        }
    }

    fn tmp_path(&self, ordinal: usize) -> PathBuf {
        self.packs_dir
            .join(format!(".rgstore.tmp.pack.{}.{}", self.pid, ordinal))
    }

    fn get_or_create(&self, ordinal: usize) -> Result<std::sync::Arc<std::fs::File>> {
        {
            let g = self.files.lock().unwrap();
            if let Some(f) = g.get(&ordinal) {
                return Ok(std::sync::Arc::clone(f));
            }
        }
        std::fs::create_dir_all(&self.packs_dir)?;
        let path = self.tmp_path(ordinal);
        let file = std::sync::Arc::new(
            std::fs::OpenOptions::new()
                .create(true)
                .truncate(false)
                .read(true)
                .write(true)
                .open(&path)?,
        );
        let mut g = self.files.lock().unwrap();
        if let Some(existing) = g.get(&ordinal) {
            return Ok(std::sync::Arc::clone(existing));
        }
        g.insert(ordinal, std::sync::Arc::clone(&file));
        Ok(file)
    }

    /// Positional write of `bytes` at `offset` into the temp pack for `ordinal`.
    fn write_at(&self, ordinal: usize, offset: u64, bytes: &[u8]) -> Result<()> {
        let file = self.get_or_create(ordinal)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::FileExt;
            file.write_all_at(bytes, offset)?;
        }
        #[cfg(not(unix))]
        {
            use std::io::{Seek, SeekFrom, Write};
            // No pwrite in stable std off-unix: clone for an independent cursor.
            let mut handle = file.try_clone()?;
            handle.seek(SeekFrom::Start(offset))?;
            handle.write_all(bytes)?;
        }
        Ok(())
    }
}

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

    let mut seqmeta_hashmap: HashMap<String, SequenceMetadata> =
        HashMap::with_capacity(seqcol.sequences.len());
    for r in &seqcol.sequences {
        let meta = r.metadata().clone();
        if seqmeta_hashmap.insert(meta.name.clone(), meta).is_some() {
            return Err(BuildError::Other(anyhow!(
                "RGSI cache has duplicate sequence names; rebuilding from FASTA"
            )));
        }
    }

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

/// Count FASTA records (`>` header lines) in a possibly-gzipped file, without
/// digesting or encoding. Used by the standalone pre-flight count gate.
fn count_fasta_records(path: &Path) -> Result<usize> {
    use std::io::BufRead;
    let mut reader = gtars_core::utils::get_dynamic_reader(path)
        .with_context(|| format!("opening {} for a pre-flight count", path.display()))?;
    let mut n = 0usize;
    let mut buf = Vec::new();
    loop {
        buf.clear();
        let read = reader.read_until(b'\n', &mut buf)?;
        if read == 0 {
            break;
        }
        if buf.first() == Some(&b'>') {
            n += 1;
        }
    }
    Ok(n)
}

/// Count data rows (non-comment, non-blank) in an `.rgsi` sidecar — the cached
/// per-file sequence count, reused for free by the pre-flight count gate.
fn count_rgsi_rows(path: &Path) -> Result<usize> {
    use std::io::BufRead;
    let file = std::fs::File::open(path)?;
    let mut n = 0usize;
    for line in std::io::BufReader::new(file).lines() {
        let line = line?;
        if !line.starts_with('#') && !line.trim().is_empty() {
            n += 1;
        }
    }
    Ok(n)
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
    /// Returns an [`ImportReport`]: per-file `(collection_metadata, was_new)`
    /// results in input order, plus genuine per-run ingest counters
    /// (`n_sequences_written`, `n_sequences_deduped`, `n_collections_new`).
    /// Note these are per-RUN counts, unlike [`StoreStats`], which is a
    /// RAM-residency snapshot.
    ///
    /// ## Error handling / non-transactional semantics
    ///
    /// This function is NOT fully transactional. On a build or writer error, any
    /// sequences already emitted (as `InsertMsg::Seq`) before the failure may have
    /// been written to disk (as `.seq` files) and inserted into the in-memory
    /// `sequence_store`, but the corresponding collection will NOT be in
    /// `self.collections` (it never received an `End`). These are orphan sequences:
    /// referenced by no surviving collection.
    ///
    /// On error this function calls [`Self::remove_orphan_seq_files`] which
    /// performs a best-effort cleanup: it removes any sequence from `sequence_store`
    /// (and its `.seq` file) that is not referenced by any surviving collection.
    /// The in-memory store is left consistent (no dangling references), but if
    /// you require a fully clean store on error, discard the store and rebuild.
    pub fn add_sequence_collections_from_fastas(
        &mut self,
        files: &[PathBuf],
        opts: FastaImportOptions<'_>,
    ) -> Result<ImportReport> {
        if files.is_empty() {
            return Ok(ImportReport {
                collections: Vec::new(),
                n_sequences_written: 0,
                n_sequences_deduped: 0,
                n_collections_new: 0,
            });
        }

        // One collection alias cannot name N collections. Reject BEFORE any
        // thread spawns or any file is opened, so a rejected call leaves the
        // store completely untouched. Applying the alias to every collection
        // would let an arbitrary one win (finalize order across builder threads
        // is non-deterministic), which is a worse version of the silent
        // misnaming this option exists to eliminate.
        if opts.collection_alias.is_some() && files.len() > 1 {
            return Err(anyhow!(
                "collection_alias names a single collection but {} FASTA files were given; \
                 import them without collection_alias and call add_collection_alias() per \
                 returned collection metadata (results are in input order)",
                files.len(),
            ));
        }

        // Resolve concurrency: 0 = auto (available_parallelism). Has no effect
        // with a single file (clamped to 1).
        let jobs = match opts.jobs {
            0 => available_parallelism().map(|n| n.get()).unwrap_or(1),
            n => n,
        };
        let jobs = jobs.min(files.len()).max(1);

        // Pre-flight count gate (STANDALONE disk imports only). Packed stores
        // bypass it; `force` overrides it. Checked BEFORE any bytes are written,
        // via a cheap counting pass (no digesting, no encoding).
        if self.persist_to_disk
            && self.layout == Layout::Standalone
            && !opts.force
        {
            use crate::utils::PathExtension;
            let existing = self.sequence_store.len();
            let limit = opts.standalone_file_limit;
            let mut counted = 0usize;
            for file in files {
                let rgsi = file.replace_exts_with("rgsi");
                let n = if self.local_path.is_some() && rgsi.exists() {
                    count_rgsi_rows(&rgsi).unwrap_or(0)
                } else {
                    count_fasta_records(file)?
                };
                counted += n;
                if existing + counted > limit {
                    return Err(anyhow!(
                        "importing {} sequences would put the store at {} standalone .seq \
                         files, over the limit of {}. Re-create the store with --packed \
                         (RefgetStore::on_disk_packed) for a packed layout, or pass force \
                         to override.",
                        counted,
                        existing + counted,
                        limit
                    ));
                }
            }
        }

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

        // --- Writer pool (disk-backed stores only) ---------------------------
        // The single inserter does the FAST in-memory work (dedup decision,
        // name_lookup/alias buffering, collection finalization) and DISPATCHES
        // the actual `.seq` byte write for each NEW sequence to a bounded pool
        // of `jobs` writer threads. Per-digest `.seq` files are independent, so
        // these writes are race-free. A BOUNDED channel provides back-pressure
        // so memory stays O(jobs x channel depth x avg seq size). In-memory mode
        // never dispatches (the writer pool is effectively unused). The shard-dir
        // cache lets writers skip a `create_dir_all` syscall per (tiny) sequence.
        let writer_disk_backed = self.persist_to_disk && self.local_path.is_some();
        let packed = self.layout == Layout::Packed;
        let (write_tx, write_rx) = bounded::<WriteJob>(jobs * WRITER_CHANNEL_DEPTH);
        let created_shards: Mutex<HashSet<PathBuf>> = Mutex::new(HashSet::new());
        let created_shards_ref = &created_shards;
        let write_err: Mutex<Option<anyhow::Error>> = Mutex::new(None);
        let write_err_ref = &write_err;

        // Packed layout: shared temp-pack writer handles, plus the inserter-side
        // deterministic assigner and the sidecar rows it records.
        let pack_writers: Option<PackWriters> = if packed && writer_disk_backed {
            Some(PackWriters::new(
                self.local_path.as_ref().unwrap().join("packs"),
            ))
        } else {
            None
        };
        let pack_writers_ref = pack_writers.as_ref();
        let mut assigner = PackAssigner::new(self.pack_cap_bytes);
        let mut pending_pack_rows: Vec<PendingPackRow> = Vec::new();
        // Filled after the writer barrier (below), read after the scope closes.
        let mut sealed_packs: Vec<SealedPack> = Vec::new();

        // --- Per-run ingest counters ----------------------------------------
        // Plain locals, not atomics: the inserter loop runs on THIS thread
        // (inside the scope closure), so these are ordinary register/stack
        // increments with no synchronization in the per-sequence hot path.
        let mut n_seqs_written: usize = 0;
        let mut n_seqs_deduped: usize = 0;

        let scope_result = std::thread::scope(|scope| -> Result<()> {
            // Spawn the writer pool. Each writer drains `WriteJob`s and writes
            // the `.seq` bytes; the first error is stashed and stops the writer.
            let mut writer_handles = Vec::with_capacity(jobs);
            if writer_disk_backed {
                for _ in 0..jobs {
                    let write_rx: Receiver<WriteJob> = write_rx.clone();
                    let handle = scope.spawn(move || {
                        for job in write_rx.iter() {
                            let result = match job {
                                WriteJob::Standalone { full_path, bytes } => {
                                    ReadonlyRefgetStore::write_seq_bytes_to_full_path(
                                        &full_path,
                                        &bytes,
                                        created_shards_ref,
                                    )
                                }
                                WriteJob::PackChunk {
                                    ordinal,
                                    offset,
                                    bytes,
                                } => match pack_writers_ref {
                                    Some(pw) => pw.write_at(ordinal, offset, &bytes),
                                    None => Err(anyhow!(
                                        "internal: PackChunk dispatched without a pack writer"
                                    )),
                                },
                            };
                            if let Err(e) = result {
                                let mut slot = write_err_ref.lock().unwrap();
                                if slot.is_none() {
                                    *slot = Some(e);
                                }
                                break;
                            }
                        }
                    });
                    writer_handles.push(handle);
                }
            }
            // The inserter keeps its own clone of `write_tx`; drop this extra one
            // so the writers' channel closes once the inserter drops its clone.
            drop(write_rx);

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

            // Move ownership of out_rx into the inserter block so we can
            // explicitly drop it after the loop (before joining the feeder).
            let out_rx = out_rx;
            for msg in out_rx.iter() {
                match msg {
                    InsertMsg::Begin { file_idx } => {
                        in_flight.insert(file_idx, InFlight::new());
                    }
                    InsertMsg::Skip { file_idx, metadata } => {
                        // Already-present collection skipped before any decode.
                        // No scratch was created (no Begin), nothing to persist.
                        // This path never reaches finalize_collection, so the
                        // requested collection alias must be registered HERE --
                        // otherwise re-importing an already-imported FASTA would
                        // silently drop the name the caller asked for.
                        if let Err(e) = self.register_import_collection_alias(
                            opts.collection_alias,
                            &metadata.digest,
                            opts.force,
                        ) {
                            let mut slot = build_err_ref.lock().unwrap();
                            if slot.is_none() { *slot = Some(e); }
                            break;
                        }
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

                        // Commit the dedup decision + in-memory metadata on this
                        // (single) inserter thread, then DISPATCH the `.seq` byte
                        // write to the bounded writer pool. Dedup is keyed on the
                        // content digest. The first occurrence of a digest always
                        // dispatches a write. For a sequence shared across
                        // collections we re-do the in-memory metadata ONLY when
                        // this file's index is higher than the current name owner,
                        // so the highest-file-index name wins deterministically
                        // (matching the serial fixed-file-order build). The bytes
                        // are identical (same digest) on such an overwrite, so we
                        // do NOT re-dispatch the disk write -- the on-disk `.seq`
                        // is already correct (byte-identical). For disk-backed
                        // stores the in-memory record is a Stub (no bytes
                        // retained); for in-memory stores the Full record is
                        // retained and no write is dispatched.
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
                                n_seqs_deduped += 1;
                                continue;
                            }
                        };
                        // A forced overwrite only changes in-memory name metadata
                        // (identical bytes); skip dispatching a redundant write.
                        let already_written = force;
                        // Per-run counter classification costs NOTHING on the
                        // disk-backed (throughput-critical) path: it is derived
                        // from `force` plus whether a write was dispatched, both
                        // of which we already have. Only an in-memory store is
                        // ambiguous -- there a `None` return means either a dedup
                        // hit or a fresh insert -- so it alone pays for a lookup,
                        // and the `!writer_disk_backed` short-circuit keeps that
                        // lookup out of the disk path entirely.
                        let mem_is_new =
                            !writer_disk_backed && !self.sequence_store.contains_key(&seq_key);
                        match self.add_sequence_record_deferred_write(
                            SequenceRecord::Full {
                                metadata,
                                sequence: std::sync::Arc::new(sequence_data),
                            },
                            force,
                        )? {
                            Some((full_path, bytes)) => {
                                if !already_written {
                                    n_seqs_written += 1;
                                    let job = if packed {
                                        // Deterministic pack assignment on THIS
                                        // (single) inserter thread, in arrival
                                        // order, then a positional-write job.
                                        let len = bytes.len() as u64;
                                        let (ordinal, offset) = assigner.assign(&bytes);
                                        pending_pack_rows.push(PendingPackRow {
                                            key: seq_key,
                                            ordinal,
                                            offset,
                                            len,
                                        });
                                        WriteJob::PackChunk { ordinal, offset, bytes }
                                    } else {
                                        WriteJob::Standalone { full_path, bytes }
                                    };
                                    // BOUNDED send: blocks (back-pressure) when the
                                    // writer pool is saturated, capping in-flight RAM.
                                    if write_tx.send(job).is_err() {
                                        // Writers all hung up (an earlier write error);
                                        // stop feeding and let the error propagate.
                                        break;
                                    }
                                } else {
                                    // `force` name-overwrite: identical bytes are
                                    // already on disk, so nothing was written.
                                    n_seqs_deduped += 1;
                                }
                            }
                            // Disk-backed: `None` is always a dedup hit. In-memory:
                            // the pre-check above decided it.
                            None => {
                                if !writer_disk_backed && mem_is_new {
                                    n_seqs_written += 1;
                                } else {
                                    n_seqs_deduped += 1;
                                }
                            }
                        }
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
                        match self.finalize_collection(
                            metadata,
                            stub_records,
                            scratch,
                            rgsi_cache_path,
                            &source_path,
                            seq_count,
                            opts.force,
                            opts.collection_alias,
                        ) {
                            Ok((meta, was_new)) => results[file_idx] = Some((meta, was_new)),
                            Err(e) => {
                                let mut slot = build_err_ref.lock().unwrap();
                                if slot.is_none() { *slot = Some(e); }
                                break;
                            }
                        }
                    }
                }
            }
            // Disconnect the cross-file channel BEFORE joining: if we exited the
            // loop early (writer hang-up `break`, or finalize error `break`),
            // builders may be blocked in `out_tx.send()`. Dropping the receiver
            // makes those sends return Err so builders break and the scope can join.
            drop(out_rx);

            feeder.join().map_err(|e| anyhow!("Feeder thread panicked: {:?}", e))?;

            // END-OF-RUN BARRIER: close the writer channel and JOIN every writer
            // so all dispatched `.seq` files are on disk BEFORE the global index
            // files are written. Writer-completion order does not affect any
            // persisted bytes (per-digest files are independent) nor index
            // ordering (indexes are sorted and written after this barrier), so
            // determinism is preserved.
            drop(write_tx);
            for handle in writer_handles {
                handle
                    .join()
                    .map_err(|e| anyhow!("Writer thread panicked: {:?}", e))?;
            }

            // SEAL the packs (packed layout): every dispatched positional write
            // has landed. Seal each open/filled pack to its final content-hash
            // name. Writer-completion order affected no byte, so this is
            // deterministic.
            if packed {
                sealed_packs = assigner.finish();
                if let Some(pw) = pack_writers_ref {
                    for sp in &sealed_packs {
                        let tmp = pw.tmp_path(sp.ordinal);
                        let final_path = pw.packs_dir.join(format!("{}.pack", sp.hash32));
                        if final_path.exists() {
                            // Content is identical by construction (content-hash
                            // name); drop the redundant temp.
                            let _ = std::fs::remove_file(&tmp);
                        } else {
                            if let Some(f) = pw.files.lock().unwrap().get(&sp.ordinal) {
                                let _ = f.sync_all();
                            }
                            std::fs::rename(&tmp, &final_path).with_context(|| {
                                format!("sealing pack {}", final_path.display())
                            })?;
                        }
                    }
                    if let Ok(dir) = std::fs::File::open(&pw.packs_dir) {
                        let _ = dir.sync_all();
                    }
                }
            }
            Ok(())
        });

        // Helper: clean up orphan sequences (and temp packs) on any failure path.
        let cleanup_on_err = |store: &mut ReadonlyRefgetStore| {
            store.remove_orphan_seq_files();
            if packed {
                store.remove_orphan_temp_packs();
            }
        };

        // Surface scope errors (e.g. feeder/writer panics), then builder/writer errors.
        if let Err(e) = scope_result {
            cleanup_on_err(self);
            return Err(e);
        }

        // Propagate the first builder error, if any.
        if let Some(e) = build_err.into_inner().unwrap() {
            cleanup_on_err(self);
            return Err(e);
        }

        // Propagate the first writer error, if any (every `.seq` must be on disk
        // before the index files reference it).
        if let Some(e) = write_err.into_inner().unwrap() {
            cleanup_on_err(self);
            return Err(e);
        }

        // Finalize ALL indexes ONCE, after every collection is persisted -- and
        // therefore take the store write lock ONCE, for the length of a delta
        // commit rather than the length of the import.
        //
        // WITHIN ONE PROCESS this is byte-identical to a serial import:
        // sequences.rgsi sorts by sha512t24u and collections.rgci sorts by
        // collection digest, so arrival/build order is irrelevant.
        //
        // ACROSS PROCESSES that guarantee does not hold, and the reason is worth
        // knowing. The commit starts from whatever is on disk and does not
        // overwrite an already-published row for a sequence digest, so the `name`
        // column is first-committer-wins: if another process publishes the same
        // sequence as `1` before we publish it as `chr1`, `1` is what stays in
        // the index. Every content-derived column (length, alphabet, md5, and all
        // the collection digests) is unaffected -- only `name`/`description`,
        // which are advisory in this file. The authoritative per-collection names
        // live in each collections/<digest>.rgsi.
        // Packed layout: convert the inserter's recorded pack rows into sidecar
        // rows (ordinal -> final content-hash name) and stage them so the commit
        // publishes `sequences.pack.idx` alongside the other indexes.
        if packed && !pending_pack_rows.is_empty() {
            let mut ord_to_name: HashMap<usize, std::sync::Arc<str>> = HashMap::new();
            for sp in &sealed_packs {
                ord_to_name.insert(sp.ordinal, std::sync::Arc::from(sp.hash32.as_str()));
            }
            let mut pending = self.pending.lock().unwrap();
            for row in &pending_pack_rows {
                let name = ord_to_name
                    .get(&row.ordinal)
                    .ok_or_else(|| anyhow!("internal: pack ordinal {} not sealed", row.ordinal))?;
                pending.pack_rows.insert(
                    row.key,
                    super::packidx::PackRow {
                        pack: std::sync::Arc::clone(name),
                        offset: row.offset,
                        len: row.len,
                    },
                );
            }
        }

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

        let collections: Vec<(SequenceCollectionMetadata, bool)> = results
            .into_iter()
            .map(|slot| {
                slot.ok_or_else(|| anyhow!("internal error: a file index was not finalized"))
            })
            .collect::<Result<_>>()?;

        let n_collections_new = collections.iter().filter(|(_, was_new)| *was_new).count();

        Ok(ImportReport {
            collections,
            n_sequences_written: n_seqs_written,
            n_sequences_deduped: n_seqs_deduped,
            n_collections_new,
        })
    }

    /// Import a single FASTA file. Thin wrapper over the multi-file path.
    pub fn add_sequence_collection_from_fasta<P: AsRef<Path>>(
        &mut self,
        file_path: P,
        opts: FastaImportOptions<'_>,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
        let files = [file_path.as_ref().to_path_buf()];
        let mut report = self.add_sequence_collections_from_fastas(&files, opts)?;
        report
            .collections
            .pop()
            .ok_or_else(|| anyhow!("internal error: importing one file yielded no result"))
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
        collection_alias: Option<(&str, &str)>,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
        let coll_key = metadata.digest.to_key();
        let coll_digest_display = metadata.digest.clone();

        // Validate the requested collection alias BEFORE touching any state.
        // A conflicting alias is a hard error, and it must be raised while the
        // store is still untouched: everything below this point (the collection
        // record, its on-disk `.rgsi`, name_lookup, sequence aliases) is
        // committed immediately, but the top-level index is only written once
        // at the very end of the import. Erroring after those mutations would
        // report a failed import while leaving the collection in the store --
        // and, on disk, an orphaned `collections/<digest>.rgsi` that no index
        // ever references. The actual alias registration still happens below,
        // once the collection is registered, and is published with it.
        self.check_import_collection_alias(collection_alias, &metadata.digest, force)?;

        if !force && self.collections.contains_key(&coll_key) {
            // Register the name even on the in-run duplicate path: the
            // collection IS in the store and the caller DID ask for it to be
            // named. Skipping here would mean a re-import silently drops the
            // alias. The helper is idempotent, so a second registration of the
            // same digest is a no-op.
            self.register_import_collection_alias(collection_alias, &metadata.digest, force)?;
            if !self.quiet {
                println!("Skipped {} (already exists)", coll_digest_display);
            }
            return Ok((metadata, false));
        }

        let insert_start = Instant::now();

        // Write the RGSI cache for next time (only for raw-FASTA builds).
        // Empty collections are never cached (a delete at the read site is
        // defensive only — empty caches are never written here).
        if let Some(rgsi_path) = &rgsi_cache_path {
            if !stub_records.is_empty() {
                let seqcol_for_cache = SequenceCollection {
                    metadata: metadata.clone(),
                    sequences: stub_records.clone(),
                };
                let _ = seqcol_for_cache.write_collection_rgsi(rgsi_path);
            }
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
        self.record(|p| p.add_collection(coll_key));

        // Install the buffered name_lookup in FASTA order.
        self.name_lookup.insert(coll_key, scratch.name_to_digest);

        // Register aliases in FASTA order. Pending only: they are published
        // with the indexes in the final commit, not one TSV rewrite per header.
        for (ns, alias_value, sha512t24u) in &scratch.aliases {
            self.add_sequence_alias_pending(ns, alias_value, sha512t24u);
        }

        // Register the collection alias (if requested). Conflicts were already
        // rejected up front, while nothing had been mutated; this records the
        // alias as pending so it lands on disk in the same commit as the
        // collection index that makes its target resolvable.
        self.register_import_collection_alias(collection_alias, &metadata.digest, force)?;

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
