//! VCF parsing and VRS ID computation.
//!
//! Reads a VCF file (plain text or gzipped/bgzf), normalizes each variant,
//! and computes GA4GH VRS Allele identifiers using a fast zero-allocation
//! digest path (no serde_json per variant).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::sync::mpsc;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore};

use crate::digest::DigestWriter;
use crate::normalize::{RefView, normalize_ref};
// WASM-safe VCF primitives live in `vcf_core` so they compile for the wasm
// build too; this filesystem-only module reuses them.
// Re-export the public VCF primitives at `crate::vcf::` so existing callers
// (examples, downstream crates) keep their import paths after the wasm carve-out.
pub use crate::vcf_core::{VrsResult, is_real_alt, parse_vcf_record};
use crate::vcf_core::{compute_vrs_ids_streaming_readonly_from_reader, read_vcf_line, ref_view_for};

/// Detect a BGZF stream by inspecting the first gzip member's header for the
/// BGZF `BC` extra subfield. Reads 18 bytes and rewinds to the start. If the
/// read is short or the bytes don't match, we assume plain gzip.
fn is_bgzf(file: &mut File) -> Result<bool> {
    let mut buf = [0u8; 18];
    let n = file.read(&mut buf)?;
    file.seek(SeekFrom::Start(0))?;
    if n < 18 {
        return Ok(false);
    }
    // gzip magic
    if buf[0] != 0x1f || buf[1] != 0x8b {
        return Ok(false);
    }
    // FEXTRA flag + BC subfield id
    let fextra = (buf[3] & 0x04) != 0;
    Ok(fextra && buf[12] == b'B' && buf[13] == b'C')
}

/// Open a VCF file, auto-detecting gzip/bgzf compression.
fn open_vcf(path: &str) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).context(format!("Failed to open VCF: {}", path))?;
    let capacity = 256 * 1024; // 256KB buffer for large VCF files
    if path.ends_with(".gz") || path.ends_with(".bgz") {
        Ok(Box::new(BufReader::with_capacity(
            capacity,
            MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(capacity, file)))
    }
}

/// Build the name→digest mapping for a collection.
pub(crate) fn build_name_to_digest(
    store: &mut RefgetStore,
    collection_digest: &str,
) -> Result<HashMap<String, String>> {
    let collection = store
        .get_collection(collection_digest)
        .context("Failed to get sequence collection")?;

    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    for seq_record in &collection.sequences {
        let meta = seq_record.metadata();
        name_to_digest.insert(meta.name.clone(), meta.sha512t24u.clone());
    }
    Ok(name_to_digest)
}

/// Ensure a sequence's bytes are resident (Full) without decoding. For on-disk
/// stores this loads the encoded `.seq`; for in-memory stores it is a no-op.
pub(crate) fn ensure_resident(store: &mut RefgetStore, raw_digest: &str) -> Result<()> {
    let resident = store
        .get_sequence(raw_digest)
        .map(|r| r.sequence().is_some())
        .unwrap_or(false);
    if !resident {
        store.load_sequence(raw_digest)?;
    }
    Ok(())
}

// ── Streaming APIs (primary) ────────────────────────────────────────────

/// Stream VRS results via callback. Lazily decodes chromosomes, clears cache when done.
/// Returns the number of results processed.
pub fn compute_vrs_ids_streaming(
    store: &mut RefgetStore,
    collection_digest: &str,
    vcf_path: &str,
    mut on_result: impl FnMut(VrsResult),
) -> Result<usize> {
    let name_to_digest = build_name_to_digest(store, collection_digest)?;
    let mut chrom_accessions: HashMap<String, String> = HashMap::new();
    let mut digest_writer = DigestWriter::new();
    let mut line_buf = String::new();
    let mut count = 0;

    let mut reader = open_vcf(vcf_path)?;

    while read_vcf_line(&mut *reader, &mut line_buf)? {
        let Some(rec) = parse_vcf_record(&line_buf) else {
            continue;
        };
        let chrom = rec.chrom;
        let pos = rec.pos;
        let ref_allele = rec.ref_allele;

        if !chrom_accessions.contains_key(chrom) {
            if let Some(digest) = name_to_digest.get(chrom) {
                // Ensure the encoded bytes are resident (no decode); in-memory is a no-op.
                ensure_resident(store, digest.as_str())?;
                chrom_accessions.insert(chrom.to_string(), format!("SQ.{}", digest));
            }
        }

        let seq_accession = match chrom_accessions.get(chrom) {
            Some(acc) => acc.clone(),
            None => continue,
        };

        let raw_digest = &name_to_digest[chrom];
        let view = ref_view_for(store, raw_digest.as_str())
            .context(format!("Chromosome {} not available in store", chrom))?;

        for alt in rec.real_alts() {
            let norm = normalize_ref(&view, pos, ref_allele.as_bytes(), alt.as_bytes())
                .context(format!("Failed to normalize variant at {}:{}", chrom, pos + 1))?;
            let norm_seq = std::str::from_utf8(&norm.allele)
                .context(format!("Normalized allele is not valid UTF-8 at {}:{}", chrom, pos + 1))?;

            let vrs_id = digest_writer.allele_identifier_literal(
                &seq_accession,
                norm.start,
                norm.end,
                norm_seq,
            );

            on_result(VrsResult {
                chrom: chrom.to_string(),
                pos,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt.to_string(),
                vrs_id,
            });
            count += 1;
        }
    }

    Ok(count)
}

/// Stream VRS results via callback using a pre-loaded read-only store, reading
/// the VCF from a filesystem path (plain/gzip/bgzf). All referenced sequences
/// must already be decodable from the store. Returns the number of results.
///
/// This is a thin filesystem wrapper around the reader-generic core
/// [`compute_vrs_ids_streaming_readonly_from_reader`]; it opens the path and
/// delegates the per-record loop. The wasm build calls the core directly with a
/// `std::io::Cursor` over the dropped VCF bytes.
pub fn compute_vrs_ids_streaming_readonly(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    on_result: impl FnMut(VrsResult),
) -> Result<usize> {
    let reader = open_vcf(vcf_path)?;
    compute_vrs_ids_streaming_readonly_from_reader(store, name_to_digest, reader, on_result)
}

// ── Parallel, encoded-on-the-fly API ────────────────────────────────────

/// A unit of work handed to a worker: one VCF record's parsed fields plus the
/// index of the sequence (digest) it maps to. `alts` is split out by the reader
/// so workers do no string scanning beyond UTF-8 of the normalized allele.
struct WorkItem {
    /// Monotonic index of this record among accepted records (for reordering).
    record_index: u64,
    /// Index into the shared `seqs` table identifying the reference sequence.
    seq_index: usize,
    /// 0-based interbase start.
    pos: u64,
    /// Original chromosome name (carried through to the result).
    chrom: String,
    /// Reference allele bytes.
    ref_allele: String,
    /// One or more alternate alleles already filtered of symbolic/`*`/`.` entries.
    alts: Vec<String>,
}

/// Per-record output: the results for one VCF record, tagged with its index so
/// the collector can emit records in input order regardless of worker timing.
struct OrderedOutput {
    record_index: u64,
    results: Vec<VrsResult>,
}

/// Immutable per-sequence entry shared (by `&`) with every worker. Pairs the
/// `SQ.<digest>` accession with the canonical [`RefView`] built once up front via
/// [`crate::normalize::ref_view_for`]. The `RefView` borrows the readonly store's
/// bytes and decodes on the fly when the store is encoded; nothing is decoded up
/// front. `RefView<'a>: Sync`, so `&[SeqEntry]` is freely shareable across the
/// scoped worker threads.
struct SeqEntry<'a> {
    /// `SQ.<digest>` accession used in the VRS identifier.
    accession: String,
    /// Canonical reference view (decoded or encoded), built once via the shared builder.
    view: RefView<'a>,
}

/// Compute VRS Allele identifiers for every variant in a VCF using an in-process
/// pool of `threads` workers, reading the reference directly from the immutable
/// 2-bit-encoded readonly store (decode-on-the-fly, no decoded cache, no mmap).
///
/// One reader thread parses the VCF (plain or gzip/bgzf via [`open_vcf`]) and
/// dispatches per-record [`WorkItem`]s round-robin to the workers. Each worker
/// owns its own [`DigestWriter`] (which is not `Sync`) and borrows the shared
/// [`RefView`] for its record's sequence (built once up front). Mode-awareness
/// (decoded vs 2-bit-encoded) derives from the shared length-based builder
/// [`crate::normalize::ref_view_for`], so a Raw-mode store works here too (the
/// dedicated equivalence test is owned by a separate test-driven plan). Results
/// are tagged with the record's input index and reordered by a collector so the
/// callback observes records in exact input order. Within a record, alts are
/// emitted in VCF order.
///
/// The caller must make all referenced sequences resident first (e.g. via
/// `RefgetStore::load_sequence` for each needed digest, then `into_readonly`).
/// Workers only read; the `RefView` borrows are kept valid by scoped threads.
///
/// Output is identical to [`compute_vrs_ids_streaming_readonly`] on well-formed
/// input; on malformed input both paths return Err and the set of results
/// delivered to the callback before the Err is unspecified.
/// Returns the number of results emitted.
pub fn compute_vrs_ids_parallel_encoded(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    threads: usize,
    mut on_result_ordered: impl FnMut(VrsResult),
) -> Result<usize> {
    let n_workers = threads.max(1);

    // Build the shared per-sequence entries once via the canonical builder,
    // borrowing the store's bytes. The builder picks decoded vs 2-bit-encoded by
    // length, so the view is mode-correct without the worker knowing the storage
    // mode. We index records by sequence position to avoid per-record hashing in
    // the workers and to give each entry a stable index.
    let mut seq_index_of_digest: HashMap<String, usize> = HashMap::new();
    let mut seqs: Vec<SeqEntry> = Vec::new();
    // Map chrom name -> seq_index, lazily filled as the reader sees new chroms.
    let mut chrom_to_seq_index: HashMap<String, usize> = HashMap::new();
    for (name, digest) in name_to_digest.iter() {
        let idx = if let Some(&i) = seq_index_of_digest.get(digest) {
            i
        } else {
            let view = ref_view_for(store, digest.as_str()).with_context(|| {
                format!("sequence {digest} not resident/encoded in readonly store")
            })?;
            let rec = store
                .get_sequence(digest.as_str())
                .with_context(|| format!("get_sequence {digest}"))?;
            let accession = format!("SQ.{}", rec.metadata().sha512t24u);
            let i = seqs.len();
            seqs.push(SeqEntry { accession, view });
            seq_index_of_digest.insert(digest.clone(), i);
            i
        };
        chrom_to_seq_index.insert(name.clone(), idx);
    }

    let seqs = &seqs; // share by reference with the scoped workers

    let mut count = 0usize;

    // Shared first-error slot: must outlive the scope so its borrow is valid
    // for the full lifetime of the scope (scoped threads borrow it).
    let worker_err: std::sync::Mutex<Option<anyhow::Error>> = std::sync::Mutex::new(None);
    let worker_err_ref = &worker_err;

    std::thread::scope(|scope| -> Result<()> {

        // Work channels: one per worker (round-robin keeps per-worker order, and
        // the collector restores global order via record_index).
        let mut work_txs: Vec<mpsc::Sender<WorkItem>> = Vec::with_capacity(n_workers);
        let (out_tx, out_rx) = mpsc::channel::<OrderedOutput>();

        let mut worker_handles = Vec::with_capacity(n_workers);
        for _ in 0..n_workers {
            let (work_tx, work_rx) = mpsc::channel::<WorkItem>();
            work_txs.push(work_tx);
            let out_tx = out_tx.clone();
            let handle = scope.spawn(move || {
                let mut writer = DigestWriter::new();
                for item in work_rx.iter() {
                    let entry = &seqs[item.seq_index];
                    let mut results = Vec::with_capacity(item.alts.len());
                    for alt in &item.alts {
                        let norm = match normalize_ref(
                            &entry.view,
                            item.pos,
                            item.ref_allele.as_bytes(),
                            alt.as_bytes(),
                        ) {
                            Ok(n) => n,
                            Err(e) => {
                                let mut s = worker_err_ref.lock().unwrap();
                                if s.is_none() {
                                    *s = Some(anyhow::anyhow!(
                                        "Failed to normalize variant at {}:{}: {e}",
                                        item.chrom,
                                        item.pos + 1
                                    ));
                                }
                                continue;
                            }
                        };
                        let norm_seq = match std::str::from_utf8(&norm.allele) {
                            Ok(s) => s,
                            Err(e) => {
                                let mut s = worker_err_ref.lock().unwrap();
                                if s.is_none() {
                                    *s = Some(anyhow::anyhow!(
                                        "Normalized allele not valid UTF-8 at {}:{}: {e}",
                                        item.chrom,
                                        item.pos + 1
                                    ));
                                }
                                continue;
                            }
                        };
                        let vrs_id = writer.allele_identifier_literal(
                            &entry.accession,
                            norm.start,
                            norm.end,
                            norm_seq,
                        );
                        results.push(VrsResult {
                            chrom: item.chrom.clone(),
                            pos: item.pos,
                            ref_allele: item.ref_allele.clone(),
                            alt_allele: alt.clone(),
                            vrs_id,
                        });
                    }
                    // Send even empty result sets so the collector can advance the
                    // expected index without stalling.
                    if out_tx.send(OrderedOutput {
                        record_index: item.record_index,
                        results,
                    }).is_err()
                    {
                        break;
                    }
                }
            });
            worker_handles.push(handle);
        }
        drop(out_tx); // workers hold the only remaining senders

        // Reader thread: parse the VCF and dispatch work round-robin.
        let reader_handle = scope.spawn(move || -> Result<u64> {
            let mut reader = open_vcf(vcf_path)?;
            let mut line_buf = String::new();
            let mut record_index: u64 = 0;
            let mut next_worker = 0usize;

            while read_vcf_line(&mut *reader, &mut line_buf)? {
                let Some(rec) = parse_vcf_record(&line_buf) else {
                    continue;
                };
                let seq_index = match chrom_to_seq_index.get(rec.chrom) {
                    Some(&i) => i,
                    None => continue, // unknown chrom: skip (matches serial)
                };
                let alts: Vec<String> = rec.real_alts().map(|s| s.to_string()).collect();
                if alts.is_empty() {
                    continue;
                }

                let item = WorkItem {
                    record_index,
                    seq_index,
                    pos: rec.pos,
                    chrom: rec.chrom.to_string(),
                    ref_allele: rec.ref_allele.to_string(),
                    alts,
                };
                record_index += 1;
                // If a worker has hung up, stop early.
                if work_txs[next_worker].send(item).is_err() {
                    break;
                }
                next_worker = (next_worker + 1) % n_workers;
            }
            // Dropping all work senders signals workers to finish.
            drop(work_txs);
            Ok(record_index)
        });

        // Collector: reorder OrderedOutput by record_index and emit in order.
        // Buffer out-of-order records in a HashMap keyed by index.
        let mut pending: HashMap<u64, Vec<VrsResult>> = HashMap::new();
        let mut next_emit: u64 = 0;
        for out in out_rx.iter() {
            pending.insert(out.record_index, out.results);
            while let Some(results) = pending.remove(&next_emit) {
                for r in results {
                    on_result_ordered(r);
                    count += 1;
                }
                next_emit += 1;
            }
        }
        // Drain any stragglers (should be none if all workers completed).
        while let Some(results) = pending.remove(&next_emit) {
            for r in results {
                on_result_ordered(r);
                count += 1;
            }
            next_emit += 1;
        }

        // Surface reader errors (e.g. bad POS) and worker panics.
        reader_handle
            .join()
            .map_err(|_| anyhow::anyhow!("VCF reader thread panicked"))??;
        for h in worker_handles {
            h.join()
                .map_err(|_| anyhow::anyhow!("VRS worker thread panicked"))?;
        }
        // Surface the first worker error (if any), matching the serial path's
        // error behaviour on malformed input.
        if let Some(e) = worker_err_ref.lock().unwrap().take() {
            return Err(e);
        }
        Ok(())
    })?;

    Ok(count)
}

// ── BGZF-block-parallel, encoded-on-the-fly API ─────────────────────────
//
// The single-reader path above (`compute_vrs_ids_parallel_encoded`) is bound
// by one thread that decompresses + parses the whole (B)GZF stream and feeds
// workers; it does not scale past the reader. This path goes one level deeper:
// the reader thread does *only* raw block I/O (no decompression). It reads raw
// compressed BGZF *blocks* (~64 KB each), tags each with a monotonic batch_id,
// and hands them to workers. Each worker decompresses its own block, finds line
// boundaries by byte search, and parses + normalizes + digests the complete
// lines inside the block, returning per-block results plus the two "fragment"
// byte slices at each boundary (bytes before the first `\n` and after the last
// `\n`). The collector sorts blocks by id, stitches each block's tail with the
// next block's head to form the cross-boundary lines, and emits final VRS
// results in VCF order. All decompression + parsing happens on workers, so the
// expensive work is fully parallel.
//
// Crucially, the workers compute VRS via the canonical ENCODED `RefView` path
// (decode-on-the-fly against the 2-bit store) -- NOT an eager full-collection
// decode. The store stays resident as encoded bytes, so peak RSS is in the tens
// of MB while throughput scales with worker count.

/// Bounded-channel depth between the reader, workers, and collector. Caps the
/// number of in-flight raw blocks / per-block results so a slow stage applies
/// backpressure instead of letting the queue grow without bound.
const PARALLEL_CHANNEL_CAPACITY: usize = 64;

/// Per-block processing output from a worker.
struct BgzfBlockResult {
    batch_id: usize,
    /// Bytes before the first `\n` in the decompressed block. Only meaningful
    /// when `had_newline == true` -- in that case these bytes belong to the line
    /// that started in a previous block and finishes here. Empty otherwise.
    head_fragment: Vec<u8>,
    /// Bytes after the last `\n` in the decompressed block if `had_newline ==
    /// true`. If `had_newline == false` the entire decompressed block is stored
    /// here; the collector accumulates it into `prev_tail` verbatim.
    tail_fragment: Vec<u8>,
    /// VRS results (or errors) for every complete line fully inside this block,
    /// in line order. Empty when `had_newline == false`.
    results: Vec<Result<VrsResult>>,
    /// Whether this block's decompressed payload contained at least one `\n`.
    /// `false` means the block is entirely interior to a single VCF line that
    /// spans multiple BGZF blocks (e.g. a very long INFO field in gnomAD).
    had_newline: bool,
}

/// Read the next raw BGZF block from `reader`. Returns `Ok(None)` at EOF,
/// `Ok(Some(bytes))` with the full block (header + deflate + trailer) --
/// possibly the zero-length EOF block at the end of a BGZF stream. No
/// decompression happens here; this is raw I/O only.
fn read_raw_bgzf_block<R: std::io::Read>(reader: &mut R) -> Result<Option<Vec<u8>>> {
    // BGZF header is 18 bytes: 12 bytes of fixed gzip/FEXTRA, then a 6-byte
    // "extra subfield" where bytes 12-13 are the `BC` subfield id and bytes
    // 16-17 are BSIZE (u16 little-endian). Total block size is `BSIZE + 1`.
    let mut header = [0u8; 18];
    match reader.read_exact(&mut header) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e.into()),
    }
    // Some BGZF files (notably gnomAD VCFs) append non-BGZF trailing bytes after
    // the stream's EOF marker. Match flate2::MultiGzDecoder's behavior (used by
    // open_vcf) and treat non-BGZF content as end-of-stream.
    if header[0] != 0x1f
        || header[1] != 0x8b
        || header[2] != 0x08
        || (header[3] & 0x04) == 0
        || header[12] != b'B'
        || header[13] != b'C'
    {
        return Ok(None);
    }
    let bsize = u16::from_le_bytes([header[16], header[17]]) as usize;
    let total = bsize + 1;
    // A valid BGZF block has at minimum: 18-byte header + 8-byte gzip trailer
    // (CRC32 + ISIZE) = 26 bytes. Anything shorter is malformed.
    if total < 26 {
        return Err(anyhow::anyhow!("Invalid BGZF block size: {}", bsize));
    }
    let mut block = Vec::with_capacity(total);
    block.extend_from_slice(&header);
    block.resize(total, 0);
    reader.read_exact(&mut block[18..])?;
    Ok(Some(block))
}

/// Decompress one complete BGZF block into its payload bytes. Each BGZF block
/// is a complete gzip member, so `MultiGzDecoder` handles it directly.
fn decompress_bgzf_block(block: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(64 * 1024);
    MultiGzDecoder::new(block).read_to_end(&mut out)?;
    Ok(out)
}

/// Parse one complete VCF data line (bytes, no trailing newline) and, on
/// success, append one `VrsResult` per non-symbolic ALT allele to `out`.
/// Errors are pushed as `Err` entries to preserve ordering.
///
/// Computes VRS via the canonical encoded `RefView` path: the chromosome is
/// looked up in `chrom_to_seq_index` to find its shared [`SeqEntry`], and each
/// allele is normalized against that entry's [`RefView`] (decode-on-the-fly).
fn process_vcf_line_bytes(
    line: &[u8],
    seqs: &[SeqEntry<'_>],
    chrom_to_seq_index: &HashMap<String, usize>,
    digest_writer: &mut DigestWriter,
    out: &mut Vec<Result<VrsResult>>,
) {
    // Strip trailing \r (Windows-style line endings inside a block).
    let line = if line.last() == Some(&b'\r') {
        &line[..line.len() - 1]
    } else {
        line
    };
    if line.is_empty() || line[0] == b'#' {
        return;
    }
    let mut fields = line.splitn(6, |&b| b == b'\t');
    let chrom_b = match fields.next() { Some(f) => f, None => return };
    let pos_b = match fields.next() { Some(f) => f, None => return };
    let _id_b = match fields.next() { Some(f) => f, None => return };
    let ref_b = match fields.next() { Some(f) => f, None => return };
    let alt_b = match fields.next() { Some(f) => f, None => return };

    let chrom = match std::str::from_utf8(chrom_b) {
        Ok(s) => s,
        Err(_) => return,
    };
    let pos: u64 = match std::str::from_utf8(pos_b)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
    {
        Some(p) => p.saturating_sub(1),
        None => return, // matches serial parse_vcf_record: unparseable POS -> skip
    };
    let ref_str = match std::str::from_utf8(ref_b) {
        Ok(s) => s,
        Err(_) => return,
    };
    let alt_field = match std::str::from_utf8(alt_b) {
        Ok(s) => s,
        Err(_) => return,
    };

    // Unknown chromosome: skip (matches serial).
    let &seq_index = match chrom_to_seq_index.get(chrom) {
        Some(i) => i,
        None => return,
    };
    let entry = &seqs[seq_index];

    for alt in alt_field.split(',') {
        if !is_real_alt(alt) {
            continue;
        }
        let norm = match normalize_ref(&entry.view, pos, ref_str.as_bytes(), alt.as_bytes()) {
            Ok(n) => n,
            Err(e) => {
                out.push(Err(anyhow::anyhow!(
                    "Failed to normalize variant at {}:{}: {}",
                    chrom,
                    pos + 1,
                    e
                )));
                continue;
            }
        };
        let norm_seq = match std::str::from_utf8(&norm.allele) {
            Ok(s) => s,
            Err(e) => {
                out.push(Err(anyhow::anyhow!(
                    "Normalized allele is not valid UTF-8 at {}:{}: {}",
                    chrom,
                    pos + 1,
                    e
                )));
                continue;
            }
        };
        let vrs_id =
            digest_writer.allele_identifier_literal(&entry.accession, norm.start, norm.end, norm_seq);
        out.push(Ok(VrsResult {
            chrom: chrom.to_string(),
            pos,
            ref_allele: ref_str.to_string(),
            alt_allele: alt.to_string(),
            vrs_id,
        }));
    }
}

/// Decompress one block, split into head_fragment / complete_lines /
/// tail_fragment, and process the complete lines via the encoded `RefView`
/// path. All CPU-heavy work (inflate + normalize + digest) happens here, on a
/// worker thread.
fn process_bgzf_block(
    batch_id: usize,
    raw_block: &[u8],
    seqs: &[SeqEntry<'_>],
    chrom_to_seq_index: &HashMap<String, usize>,
    digest_writer: &mut DigestWriter,
) -> BgzfBlockResult {
    let decompressed = match decompress_bgzf_block(raw_block) {
        Ok(d) => d,
        Err(e) => {
            return BgzfBlockResult {
                batch_id,
                head_fragment: Vec::new(),
                tail_fragment: Vec::new(),
                results: vec![Err(e)],
                had_newline: false,
            };
        }
    };
    if decompressed.is_empty() {
        // BGZF EOF block or empty block -- contributes nothing.
        return BgzfBlockResult {
            batch_id,
            head_fragment: Vec::new(),
            tail_fragment: Vec::new(),
            results: Vec::new(),
            had_newline: false,
        };
    }

    let first_nl = decompressed.iter().position(|&b| b == b'\n');
    let last_nl = decompressed.iter().rposition(|&b| b == b'\n');

    // If the block has no newline at all, it sits entirely inside a single VCF
    // line that spans multiple BGZF blocks (common with very long INFO fields).
    // Emit nothing; stash the whole decompressed payload in `tail_fragment` so
    // the collector accumulates it into `prev_tail` without treating it as a
    // complete line.
    let (first, last) = match (first_nl, last_nl) {
        (Some(f), Some(l)) => (f, l),
        _ => {
            return BgzfBlockResult {
                batch_id,
                head_fragment: Vec::new(),
                tail_fragment: decompressed,
                results: Vec::new(),
                had_newline: false,
            };
        }
    };

    // Block has >=1 newline: split head / middle / tail. The middle slice is
    // empty when `first == last` (single newline); `split(b'\n')` on an empty
    // slice yields one empty chunk which is filtered below.
    let head_fragment = decompressed[..first].to_vec();
    let tail_fragment = decompressed[last + 1..].to_vec();
    let middle_slice: &[u8] = &decompressed[first + 1..last + 1];

    let mut results: Vec<Result<VrsResult>> = Vec::new();
    for line in middle_slice.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        process_vcf_line_bytes(line, seqs, chrom_to_seq_index, digest_writer, &mut results);
    }

    BgzfBlockResult {
        batch_id,
        head_fragment,
        tail_fragment,
        results,
        had_newline: true,
    }
}

/// Build the shared per-sequence table (accession + canonical [`RefView`]) and
/// the chrom-name -> table-index map, once, borrowing the readonly store's
/// bytes. The [`RefView`] decodes on the fly when the store is encoded; nothing
/// is decoded up front. Shared by `&` with every worker (RefView is Sync).
fn build_seq_table<'a>(
    store: &'a ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
) -> Result<(Vec<SeqEntry<'a>>, HashMap<String, usize>)> {
    let mut seq_index_of_digest: HashMap<String, usize> = HashMap::new();
    let mut seqs: Vec<SeqEntry<'a>> = Vec::new();
    let mut chrom_to_seq_index: HashMap<String, usize> = HashMap::new();
    for (name, digest) in name_to_digest.iter() {
        let idx = if let Some(&i) = seq_index_of_digest.get(digest) {
            i
        } else {
            let view = ref_view_for(store, digest.as_str()).with_context(|| {
                format!("sequence {digest} not resident/encoded in readonly store")
            })?;
            let rec = store
                .get_sequence(digest.as_str())
                .with_context(|| format!("get_sequence {digest}"))?;
            let accession = format!("SQ.{}", rec.metadata().sha512t24u);
            let i = seqs.len();
            seqs.push(SeqEntry { accession, view });
            seq_index_of_digest.insert(digest.clone(), i);
            i
        };
        chrom_to_seq_index.insert(name.clone(), idx);
    }
    Ok((seqs, chrom_to_seq_index))
}

/// BGZF-block-parallel VRS computation with encoded-on-the-fly workers.
///
/// Requires BGZF input (gnomAD-style `.bgz`). The reader thread reads raw BGZF
/// blocks (no decompression); `num_workers` worker threads each decompress +
/// parse + normalize + digest their own blocks against the shared encoded
/// [`RefView`]s; a streaming collector reorders blocks by id, stitches
/// cross-block boundary lines, and invokes `on_result` for each `VrsResult` in
/// exact VCF order. Only one block's tail fragment plus the current reorder gap
/// sit in memory, so peak RSS stays in the tens of MB.
///
/// Output is byte-identical to [`compute_vrs_ids_streaming_readonly`] on
/// well-formed input. On malformed input (e.g. a REF extending past the
/// sequence end) it returns the first `Err`, matching the serial path.
///
/// The caller must make all referenced sequences resident (encoded) first, then
/// `into_readonly()`. Returns the number of results emitted.
pub fn compute_vrs_ids_parallel_bgzf_encoded_with_sink<F: FnMut(VrsResult)>(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    num_workers: usize,
    mut on_result: F,
) -> Result<usize> {
    let mut file = File::open(vcf_path).context(format!("Failed to open VCF: {}", vcf_path))?;
    if !is_bgzf(&mut file)? {
        return Err(anyhow::anyhow!(
            "compute_vrs_ids_parallel_bgzf_encoded_with_sink requires BGZF input; \
             use compute_vrs_ids_parallel_encoded for plain gzip / uncompressed"
        ));
    }

    let num_workers = num_workers.max(1);

    // Build the shared encoded-view table once; workers borrow it by `&`.
    let (seqs, chrom_to_seq_index) = build_seq_table(store, name_to_digest)?;
    let seqs = &seqs;
    let chrom_to_seq_index = &chrom_to_seq_index;

    std::thread::scope(|s| -> Result<usize> {
        let (block_tx, block_rx) =
            crossbeam_channel::bounded::<(usize, Vec<u8>)>(PARALLEL_CHANNEL_CAPACITY);
        let (result_tx, result_rx) =
            crossbeam_channel::bounded::<BgzfBlockResult>(PARALLEL_CHANNEL_CAPACITY);

        // Reader thread: raw block I/O only. No decompression, no parsing.
        let reader_handle = s.spawn(move || -> Result<()> {
            let mut reader = BufReader::with_capacity(1 << 20, file);
            let mut batch_id: usize = 0;
            loop {
                match read_raw_bgzf_block(&mut reader)? {
                    Some(block) => {
                        if block_tx.send((batch_id, block)).is_err() {
                            return Ok(());
                        }
                        batch_id += 1;
                    }
                    None => break,
                }
            }
            Ok(())
        });

        // Workers decompress + parse + normalize + digest their blocks.
        let mut worker_handles = Vec::with_capacity(num_workers);
        for _ in 0..num_workers {
            let block_rx = block_rx.clone();
            let result_tx = result_tx.clone();
            let handle = s.spawn(move || {
                let mut digest_writer = DigestWriter::new();
                while let Ok((batch_id, raw_block)) = block_rx.recv() {
                    let br = process_bgzf_block(
                        batch_id,
                        &raw_block,
                        seqs,
                        chrom_to_seq_index,
                        &mut digest_writer,
                    );
                    if result_tx.send(br).is_err() {
                        break;
                    }
                }
            });
            worker_handles.push(handle);
        }

        drop(block_rx);
        drop(result_tx);

        // Streaming collector. Blocks arrive in arbitrary order; buffer
        // out-of-order ones in a HashMap by batch_id and process the contiguous
        // run starting at `next_batch_id`. Each processed block's results are
        // emitted through the sink and dropped, so only one block's
        // tail_fragment and the current reorder gap sit in memory.
        let mut first_err: Option<anyhow::Error> = None;
        let mut count: usize = 0;
        let mut prev_tail: Vec<u8> = Vec::new();
        let mut stitch_writer = DigestWriter::new();
        let mut pending: HashMap<usize, BgzfBlockResult> = HashMap::new();
        let mut next_batch_id: usize = 0;

        #[allow(clippy::too_many_arguments)]
        fn process_block<F: FnMut(VrsResult)>(
            block: BgzfBlockResult,
            next_batch_id: usize,
            prev_tail: &mut Vec<u8>,
            stitch_writer: &mut DigestWriter,
            seqs: &[SeqEntry<'_>],
            chrom_to_seq_index: &HashMap<String, usize>,
            on_result: &mut F,
            count: &mut usize,
            first_err: &mut Option<anyhow::Error>,
        ) {
            if !block.had_newline {
                // Block sits entirely inside a single VCF line that spans
                // multiple BGZF blocks. Accumulate bytes into prev_tail; the
                // line completes when a later block produces its first `\n`.
                for r in block.results {
                    if let Err(e) = r {
                        if first_err.is_none() {
                            *first_err = Some(e);
                        }
                    }
                }
                prev_tail.extend_from_slice(&block.tail_fragment);
                return;
            }

            // Block contains >=1 newline. Stitch prev.tail + this.head into the
            // line that straddles the boundary and process it. For
            // next_batch_id == 0 with empty prev_tail the "stitched" bytes are
            // typically the header line.
            let mut stitched = std::mem::take(prev_tail);
            stitched.extend_from_slice(&block.head_fragment);
            if !stitched.is_empty() || next_batch_id > 0 {
                let mut tmp: Vec<Result<VrsResult>> = Vec::new();
                process_vcf_line_bytes(
                    &stitched,
                    seqs,
                    chrom_to_seq_index,
                    stitch_writer,
                    &mut tmp,
                );
                for r in tmp {
                    match r {
                        Ok(v) => {
                            on_result(v);
                            *count += 1;
                        }
                        Err(e) => {
                            if first_err.is_none() {
                                *first_err = Some(e);
                            }
                        }
                    }
                }
            }
            for r in block.results {
                match r {
                    Ok(v) => {
                        on_result(v);
                        *count += 1;
                    }
                    Err(e) => {
                        if first_err.is_none() {
                            *first_err = Some(e);
                        }
                    }
                }
            }
            *prev_tail = block.tail_fragment;
        }

        while let Ok(br) = result_rx.recv() {
            pending.insert(br.batch_id, br);
            while let Some(block) = pending.remove(&next_batch_id) {
                process_block(
                    block,
                    next_batch_id,
                    &mut prev_tail,
                    &mut stitch_writer,
                    seqs,
                    chrom_to_seq_index,
                    &mut on_result,
                    &mut count,
                    &mut first_err,
                );
                next_batch_id += 1;
            }
        }
        // Drain anything left (only if the channel closed with gaps, which
        // shouldn't happen in practice).
        let mut remaining_ids: Vec<usize> = pending.keys().copied().collect();
        remaining_ids.sort_unstable();
        for id in remaining_ids {
            if let Some(block) = pending.remove(&id) {
                process_block(
                    block,
                    id,
                    &mut prev_tail,
                    &mut stitch_writer,
                    seqs,
                    chrom_to_seq_index,
                    &mut on_result,
                    &mut count,
                    &mut first_err,
                );
            }
        }

        // Final unterminated line (if any).
        if !prev_tail.is_empty() {
            let mut tmp: Vec<Result<VrsResult>> = Vec::new();
            process_vcf_line_bytes(
                &prev_tail,
                seqs,
                chrom_to_seq_index,
                &mut stitch_writer,
                &mut tmp,
            );
            for r in tmp {
                match r {
                    Ok(v) => {
                        on_result(v);
                        count += 1;
                    }
                    Err(e) => {
                        if first_err.is_none() {
                            first_err = Some(e);
                        }
                    }
                }
            }
        }

        // Surface reader / worker panics.
        match reader_handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => {
                if first_err.is_none() {
                    first_err = Some(e);
                }
            }
            Err(_) => {
                if first_err.is_none() {
                    first_err = Some(anyhow::anyhow!("BGZF reader thread panicked"));
                }
            }
        }
        for h in worker_handles {
            if h.join().is_err() && first_err.is_none() {
                first_err = Some(anyhow::anyhow!("VRS worker thread panicked"));
            }
        }
        if let Some(e) = first_err {
            return Err(e);
        }
        Ok(count)
    })
}

// ── Vec-collecting APIs (convenience wrappers) ──────────────────────────

/// Convenience wrapper: collects all VRS results into a Vec.
pub fn compute_vrs_ids_from_vcf(
    store: &mut RefgetStore,
    collection_digest: &str,
    vcf_path: &str,
) -> Result<Vec<VrsResult>> {
    let mut results = Vec::new();
    compute_vrs_ids_streaming(store, collection_digest, vcf_path, |r| results.push(r))?;
    Ok(results)
}

/// Collect all VRS results from a pre-loaded read-only store.
pub fn compute_vrs_ids_from_vcf_readonly(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
) -> Result<Vec<VrsResult>> {
    let mut results = Vec::new();
    compute_vrs_ids_streaming_readonly(store, name_to_digest, vcf_path, |r| results.push(r))?;
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_line_is_none() {
        assert!(parse_vcf_record("##fileformat=VCFv4.2").is_none());
        assert!(parse_vcf_record("#CHROM\tPOS\tID\tREF\tALT").is_none());
    }

    #[test]
    fn blank_line_is_none() {
        assert!(parse_vcf_record("").is_none());
        assert!(parse_vcf_record("\n").is_none());
        assert!(parse_vcf_record("\r\n").is_none());
    }

    #[test]
    fn too_few_fields_is_none() {
        assert!(parse_vcf_record("chr1\t100\t.\tA").is_none());
        assert!(parse_vcf_record("chr1\t100").is_none());
    }

    #[test]
    fn non_numeric_pos_is_none() {
        assert!(parse_vcf_record("chr1\tNOPE\t.\tA\tG").is_none());
    }

    #[test]
    fn normal_line_parses() {
        let rec = parse_vcf_record("chr1\t100\trs1\tA\tG\t.\tPASS\tINFO").unwrap();
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.pos, 99); // 1-based 100 -> 0-based interbase 99
        assert_eq!(rec.ref_allele, "A");
        let alts: Vec<&str> = rec.real_alts().collect();
        assert_eq!(alts, vec!["G"]);
    }

    #[test]
    fn pos_one_saturates() {
        let rec = parse_vcf_record("chr1\t1\t.\tA\tG").unwrap();
        assert_eq!(rec.pos, 0);
    }

    #[test]
    fn multi_alt_filters_symbolic_and_missing() {
        let rec = parse_vcf_record("chr1\t100\t.\tA\tG,<DEL>,*,.,,T").unwrap();
        let alts: Vec<&str> = rec.real_alts().collect();
        assert_eq!(alts, vec!["G", "T"]);
    }

    #[test]
    fn is_real_alt_truth_table() {
        assert!(is_real_alt("A"));
        assert!(!is_real_alt("<DEL>"));
        assert!(!is_real_alt("*"));
        assert!(!is_real_alt("."));
        assert!(!is_real_alt(""));
    }
}
