//! VCF parsing and VRS ID computation.
//!
//! Reads a VCF file (plain text or gzipped/bgzf), normalizes each variant,
//! and computes GA4GH VRS Allele identifiers using a fast zero-allocation
//! digest path (no serde_json per variant).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore};

use arrayvec::ArrayString;
use compact_str::CompactString;

use crate::digest::DigestWriter;
use crate::normalize::normalize;

/// Result of computing a VRS identifier for a single VCF variant.
///
/// Fields use short-string-optimized types so emitted results do not
/// heap-allocate. On the gnomAD chr22 workload (~11.6M emitted results)
/// this saves ~46M small `String` allocations relative to using `String`
/// for every field.
///
/// - `chrom`: `CompactString` — SSO inline up to 24 bytes. Chromosome
///   names (`chr1`..`chrY`, 4-5 bytes) always inline.
/// - `ref_allele` / `alt_allele`: `CompactString` — SSO inline up to 24
///   bytes. ~99% SNVs (1 byte), remainder short indels. Inlines in all
///   realistic cases; allocations only for pathological inputs.
/// - `vrs_id`: `ArrayString<48>` — fixed-length 41-byte `ga4gh:VA.<32-byte
///   base64url digest>`. Stack-allocated; no heap touch.
#[derive(Debug, Clone)]
pub struct VrsResult {
    pub chrom: CompactString,
    pub pos: u64,
    pub ref_allele: CompactString,
    pub alt_allele: CompactString,
    pub vrs_id: ArrayString<48>,
}

/// Peek the first 18 bytes to determine if a gzip-framed file is actually BGZF.
/// BGZF files carry a `BC` (0x42 0x43) extra subfield at offset 12–13 of the
/// first gzip header. If either the read fails or the bytes don't match, we
/// assume plain gzip. The stream is rewound to the start before returning.
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

/// Open a VCF file, auto-detecting plain/gzip/BGZF compression.
///
/// All gzip-framed input (plain `.gz` and BGZF `.bgz`) goes through
/// `flate2::MultiGzDecoder`, which handles single-member gzip, multi-member
/// gzip, and concatenated BGZF streams uniformly. We deliberately avoid
/// `noodles_bgzf` here: its 0.46 `Reader` / `MultithreadedReader` both run
/// every frame through a strict header validator, which rejects gnomAD-style
/// `.bgz` files that have non-BGZF bytes appended after the BGZF EOF block
/// ("invalid BGZF header"). `MultiGzDecoder` stops at the first non-gzip
/// byte and surfaces that as `ErrorKind::InvalidInput`, which
/// `read_vcf_line` below already treats as clean EOF. Block-parallel
/// decompression is still available via `compute_vrs_ids_parallel_bgzf_with_sink`,
/// which does its own raw block I/O and does not route through `open_vcf`.
/// Uncompressed input is passed through a `BufReader` directly.
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
fn build_name_to_digest(
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

/// Read the next non-empty line from a VCF reader. Returns false at EOF.
fn read_vcf_line(reader: &mut dyn BufRead, buf: &mut String) -> Result<bool> {
    buf.clear();
    match reader.read_line(buf) {
        Ok(0) => Ok(false),
        Ok(_) => Ok(true),
        Err(e) => {
            // BGZF files end with an empty gzip block that MultiGzDecoder
            // may interpret as an invalid header. Treat as EOF.
            if e.kind() == std::io::ErrorKind::InvalidInput {
                Ok(false)
            } else {
                Err(e.into())
            }
        }
    }
}

// ── Streaming APIs (primary) ────────────────────────────────────────────

/// Stream VRS results via callback. Lazily decodes chromosomes as it
/// encounters them. Decoded sequence data is left in the store so the
/// caller can run additional queries or decide when to reclaim memory
/// via `store.clear_decoded_cache()` (which downgrades cached records
/// to `Stub`, so only call it when you are done with the store or have
/// a disk-backed source to reload from).
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
        let line = line_buf.trim_end_matches('\n').trim_end_matches('\r');
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.splitn(10, '\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let chrom = fields[0];
        let pos: u64 = fields[1]
            .parse::<u64>()
            .context("Invalid POS field")?
            .saturating_sub(1);
        let ref_allele = fields[3];
        let alt_field = fields[4];

        if !chrom_accessions.contains_key(chrom) {
            if let Some(digest) = name_to_digest.get(chrom) {
                store.ensure_decoded(digest.as_str())?;
                chrom_accessions.insert(chrom.to_string(), format!("SQ.{}", digest));
            }
        }

        let seq_accession = match chrom_accessions.get(chrom) {
            Some(acc) => acc,
            None => continue,
        };

        let raw_digest = &name_to_digest[chrom];
        let sequence = store
            .sequence_bytes(raw_digest.as_str())
            .context(format!("Chromosome {} not decoded in store", chrom))?;

        for alt in alt_field.split(',') {
            if alt.starts_with('<') || alt == "*" || alt == "." {
                continue;
            }

            let norm = normalize(sequence, pos, ref_allele.as_bytes(), alt.as_bytes())
                .context(format!("Failed to normalize variant at {}:{}", chrom, pos + 1))?;
            let norm_seq = std::str::from_utf8(&norm.allele)
                .context(format!("Normalized allele is not valid UTF-8 at {}:{}", chrom, pos + 1))?;

            let vrs_id =
                digest_writer.allele_identifier_literal(seq_accession, norm.start, norm.end, norm_seq);

            on_result(VrsResult {
                chrom: CompactString::from(chrom),
                pos,
                ref_allele: CompactString::from(ref_allele),
                alt_allele: CompactString::from(alt),
                vrs_id,
            });
            count += 1;
        }
    }

    // Intentionally do not call store.clear_decoded_cache() here:
    // since commit 95a7ca5, clear_decoded_cache downgrades decoded
    // records to Stub, which permanently destroys data for in-memory
    // stores (no disk backing to reload from). Leaving the cache
    // populated also keeps subsequent calls on the same store cheap.
    // Callers that need to reclaim memory can call it explicitly.
    Ok(count)
}

/// Stream VRS results via callback using a pre-loaded read-only store.
/// All referenced sequences must already be decoded.
/// Returns the number of results processed.
pub fn compute_vrs_ids_streaming_readonly(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    mut on_result: impl FnMut(VrsResult),
) -> Result<usize> {
    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(name, digest)| (name.clone(), format!("SQ.{}", digest)))
        .collect();

    let mut digest_writer = DigestWriter::new();
    let mut line_buf = String::new();
    let mut count = 0;

    let mut reader = open_vcf(vcf_path)?;

    while read_vcf_line(&mut *reader, &mut line_buf)? {
        let line = line_buf.trim_end_matches('\n').trim_end_matches('\r');
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.splitn(10, '\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let chrom = fields[0];
        let pos: u64 = fields[1]
            .parse::<u64>()
            .context("Invalid POS field")?
            .saturating_sub(1);
        let ref_allele = fields[3];
        let alt_field = fields[4];

        let seq_accession = match chrom_accessions.get(chrom) {
            Some(acc) => acc,
            None => continue,
        };

        let raw_digest = match name_to_digest.get(chrom) {
            Some(d) => d,
            None => continue,
        };
        let sequence = store
            .sequence_bytes(raw_digest.as_str())
            .context(format!("Chromosome {} not decoded in store", chrom))?;

        for alt in alt_field.split(',') {
            if alt.starts_with('<') || alt == "*" || alt == "." {
                continue;
            }

            let norm = normalize(sequence, pos, ref_allele.as_bytes(), alt.as_bytes())
                .context(format!("Failed to normalize variant at {}:{}", chrom, pos + 1))?;
            let norm_seq = std::str::from_utf8(&norm.allele)
                .context(format!("Normalized allele is not valid UTF-8 at {}:{}", chrom, pos + 1))?;

            let vrs_id =
                digest_writer.allele_identifier_literal(seq_accession, norm.start, norm.end, norm_seq);

            on_result(VrsResult {
                chrom: CompactString::from(chrom),
                pos,
                ref_allele: CompactString::from(ref_allele),
                alt_allele: CompactString::from(alt),
                vrs_id,
            });
            count += 1;
        }
    }

    Ok(count)
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

// ── Parallel API shared constants ───────────────────────────────────────

const PARALLEL_CHANNEL_CAPACITY: usize = 64;
const PARALLEL_BATCH_SIZE: usize = 1024;

/// Build name_to_digest from a ReadonlyRefgetStore and a collection digest.
///
/// The store must already be preloaded (via `load_all_collections` /
/// `load_all_sequences` on the mutable RefgetStore before `into_readonly`).
pub fn build_name_to_digest_readonly(
    store: &ReadonlyRefgetStore,
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

// ── Parallel API (workers parse) ────────────────────────────────────────
//
// Reader thread batches *raw* VCF lines; the worker pool does all parsing
// + normalization + digest. The reader's per-line work collapses to:
// read_line into a String, skip headers, push into a batch Vec. That's ~1
// allocation per line, and every hot-path operation (field split,
// parse::<u64>, HashMap lookup, multi-alt expansion, normalize, digest)
// runs across the worker pool. See compute_vrs_ids_parallel_with_sink.

/// A batch of raw VCF data lines, tagged with a monotonic batch id so the
/// collector can restore VCF order across workers.
struct LineBatch {
    batch_id: usize,
    lines: Vec<String>,
}

/// A batch of per-line VRS results (or errors). Preserves input line order
/// within the batch; batches themselves are reordered by `batch_id`.
struct LineResultBatch {
    batch_id: usize,
    results: Vec<Result<VrsResult>>,
}

/// Streaming core for the parallel path. Invokes `on_result` for
/// each `VrsResult` in VCF order (via a per-batch reorder buffer), holding
/// at most one batch worth of results at a time — not the full result set.
///
/// Preconditions: every chromosome referenced by the VCF must already be
/// decoded in `store`.
pub fn compute_vrs_ids_parallel_with_sink<F: FnMut(VrsResult)>(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    num_workers: usize,
    mut on_result: F,
) -> Result<usize> {
    let num_workers = num_workers.max(1);

    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(name, digest)| (name.clone(), format!("SQ.{}", digest)))
        .collect();

    std::thread::scope(|s| -> Result<usize> {
        let (line_tx, line_rx) =
            crossbeam_channel::bounded::<LineBatch>(PARALLEL_CHANNEL_CAPACITY);
        let (result_tx, result_rx) =
            crossbeam_channel::bounded::<LineResultBatch>(PARALLEL_CHANNEL_CAPACITY);

        // Reader thread: read lines, skip headers, batch, send.
        let reader_handle = s.spawn(move || -> Result<()> {
            let mut reader = open_vcf(vcf_path)?;
            let mut batch: Vec<String> = Vec::with_capacity(PARALLEL_BATCH_SIZE);
            let mut batch_id: usize = 0;

            let mut scratch = String::new();
            loop {
                scratch.clear();
                // Use read_vcf_line (not raw read_line) so a BGZF empty
                // trailer block's ErrorKind::InvalidInput is treated as
                // EOF rather than a hard failure, matching the streaming
                // paths.
                if !read_vcf_line(&mut *reader, &mut scratch)? {
                    break;
                }
                let trimmed = scratch.trim_end_matches('\n').trim_end_matches('\r');
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    continue;
                }
                // `scratch` now owns this line; move it into the batch
                // and allocate a fresh String for the next iteration.
                let line = std::mem::take(&mut scratch);
                batch.push(line);
                if batch.len() >= PARALLEL_BATCH_SIZE {
                    let full = std::mem::replace(
                        &mut batch,
                        Vec::with_capacity(PARALLEL_BATCH_SIZE),
                    );
                    if line_tx.send(LineBatch { batch_id, lines: full }).is_err() {
                        return Ok(());
                    }
                    batch_id += 1;
                }
            }
            if !batch.is_empty() {
                let _ = line_tx.send(LineBatch { batch_id, lines: batch });
            }
            Ok(())
        });

        // Worker threads: parse + normalize + digest.
        let mut worker_handles = Vec::with_capacity(num_workers);
        for _ in 0..num_workers {
            let line_rx = line_rx.clone();
            let result_tx = result_tx.clone();
            let store_ref: &ReadonlyRefgetStore = store;
            let name_to_digest_ref = name_to_digest;
            let chrom_accessions_ref = &chrom_accessions;
            let handle = s.spawn(move || {
                let mut digest_writer = DigestWriter::new();
                while let Ok(line_batch) = line_rx.recv() {
                    let mut out: Vec<Result<VrsResult>> =
                        Vec::with_capacity(line_batch.lines.len());
                    for raw in &line_batch.lines {
                        let line = raw.trim_end_matches('\n').trim_end_matches('\r');
                        let fields: Vec<&str> = line.splitn(10, '\t').collect();
                        if fields.len() < 5 {
                            continue;
                        }

                        let chrom = fields[0];
                        let pos: u64 = match fields[1].parse::<u64>() {
                            Ok(p) => p.saturating_sub(1),
                            Err(_) => {
                                out.push(Err(anyhow::anyhow!(
                                    "Invalid POS field at {}:{}",
                                    chrom,
                                    fields[1]
                                )));
                                continue;
                            }
                        };
                        let ref_allele = fields[3];
                        let alt_field = fields[4];

                        let seq_accession = match chrom_accessions_ref.get(chrom) {
                            Some(acc) => acc.as_str(),
                            None => continue,
                        };
                        let raw_digest = match name_to_digest_ref.get(chrom) {
                            Some(d) => d,
                            None => continue,
                        };
                        let sequence = match store_ref.sequence_bytes(raw_digest.as_str()) {
                            Some(b) => b,
                            None => {
                                out.push(Err(anyhow::anyhow!(
                                    "Chromosome {} not decoded in store",
                                    chrom
                                )));
                                continue;
                            }
                        };

                        for alt in alt_field.split(',') {
                            if alt.starts_with('<') || alt == "*" || alt == "." {
                                continue;
                            }
                            let norm = match normalize(
                                sequence,
                                pos,
                                ref_allele.as_bytes(),
                                alt.as_bytes(),
                            ) {
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
                            let vrs_id = digest_writer.allele_identifier_literal(
                                seq_accession,
                                norm.start,
                                norm.end,
                                norm_seq,
                            );
                            out.push(Ok(VrsResult {
                                chrom: CompactString::from(chrom),
                                pos,
                                ref_allele: CompactString::from(ref_allele),
                                alt_allele: CompactString::from(alt),
                                vrs_id,
                            }));
                        }
                    }
                    if result_tx
                        .send(LineResultBatch {
                            batch_id: line_batch.batch_id,
                            results: out,
                        })
                        .is_err()
                    {
                        break;
                    }
                }
            });
            worker_handles.push(handle);
        }

        drop(line_rx);
        drop(result_tx);

        // Streaming reorder: as batches arrive from workers (in arbitrary
        // order), buffer out-of-order ones in a HashMap keyed by batch_id.
        // Emit the contiguous run starting at `next_batch_id` through the
        // sink, dropping each batch after emission. Memory usage is
        // bounded by the reorder gap (≤ channel capacity × batch size),
        // not total variants.
        let mut first_err: Option<anyhow::Error> = None;
        let mut pending: HashMap<usize, Vec<Result<VrsResult>>> = HashMap::new();
        let mut next_batch_id: usize = 0;
        let mut count: usize = 0;

        let emit = |results: Vec<Result<VrsResult>>,
                    count: &mut usize,
                    first_err: &mut Option<anyhow::Error>,
                    on_result: &mut F| {
            for r in results {
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
        };

        while let Ok(rb) = result_rx.recv() {
            pending.insert(rb.batch_id, rb.results);
            while let Some(results) = pending.remove(&next_batch_id) {
                emit(results, &mut count, &mut first_err, &mut on_result);
                next_batch_id += 1;
            }
        }
        // Drain any remaining (shouldn't happen if all batches were
        // received and emitted in order, but handle gracefully).
        let mut remaining_ids: Vec<usize> = pending.keys().copied().collect();
        remaining_ids.sort_unstable();
        for id in remaining_ids {
            if let Some(results) = pending.remove(&id) {
                emit(results, &mut count, &mut first_err, &mut on_result);
            }
        }

        // Surface reader errors / panics.
        match reader_handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => {
                if first_err.is_none() {
                    first_err = Some(e);
                }
            }
            Err(_) => {
                if first_err.is_none() {
                    first_err = Some(anyhow::anyhow!("VCF reader thread panicked"));
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

// ── BGZF-block parallel API ─────────────────────────────────────────────
//
// The parallel variant above still routes every VCF byte through one
// reader thread: it decompresses (via MultiGzDecoder) and line-batches
// sequentially before dispatching. That reader is the ceiling (~3.7× on
// 14 cores).
//
// This variant goes one level deeper: the reader thread never decompresses
// anything. It reads raw compressed BGZF *blocks* (~64 KB each), tagged
// with a monotonic batch_id, and hands them to workers. Each worker
// decompresses its own block, finds line boundaries via byte search,
// parses the complete lines inside the block, and returns per-block
// results plus the two "fragment" byte slices at each boundary (bytes
// before the first `\n` and after the last `\n`). The collector sorts
// blocks by id, stitches each block's tail with the next block's head to
// form the cross-boundary lines, and emits final VRS results in VCF
// order. All decompression and parsing happens on workers.

/// Per-block processing output from a worker.
struct BgzfBlockResult {
    batch_id: usize,
    /// Bytes before the first `\n` in the decompressed block. Only
    /// meaningful when `had_newline == true` — in that case, these bytes
    /// belong to the line that started in a previous block and finishes
    /// here. Empty when the block has no newlines.
    head_fragment: Vec<u8>,
    /// Bytes after the last `\n` in the decompressed block if
    /// `had_newline == true`. If `had_newline == false` the entire
    /// decompressed block is stored here; the collector accumulates it
    /// into `prev_tail` verbatim without emitting.
    tail_fragment: Vec<u8>,
    /// VRS results (or errors) for every complete line fully inside this
    /// block, in line order. Empty when `had_newline == false`.
    results: Vec<Result<VrsResult>>,
    /// Whether this block's decompressed payload contained at least one
    /// `\n`. `false` means the block is entirely interior to a single
    /// VCF line that spans multiple BGZF blocks (e.g. a very long INFO
    /// field in gnomAD).
    had_newline: bool,
}

/// Read the next raw BGZF block from `reader`. Returns `Ok(None)` at EOF,
/// `Ok(Some(bytes))` with the full block (header + deflate + trailer) —
/// possibly the zero-length EOF block at the end of a BGZF stream.
fn read_raw_bgzf_block<R: std::io::Read>(reader: &mut R) -> Result<Option<Vec<u8>>> {
    // BGZF header is 18 bytes: 12 bytes of fixed gzip/FEXTRA, then a
    // 6-byte "extra subfield" where bytes 12-13 are the `BC` subfield id
    // and bytes 16-17 are BSIZE (u16 little-endian). Total block size is
    // `BSIZE + 1` bytes.
    let mut header = [0u8; 18];
    match reader.read_exact(&mut header) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e.into()),
    }
    // Some BGZF files (notably gnomAD VCFs) append non-BGZF trailing bytes
    // after the stream's EOF marker. Match flate2::MultiGzDecoder's behavior
    // (used by open_vcf) and treat non-BGZF content as end-of-stream.
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
    // A valid BGZF block has at minimum: 18-byte header + 8-byte gzip
    // trailer (CRC32 + ISIZE) = 26 bytes. Anything shorter is malformed.
    if total < 26 {
        return Err(anyhow::anyhow!("Invalid BGZF block size: {}", bsize));
    }
    let mut block = Vec::with_capacity(total);
    block.extend_from_slice(&header);
    block.resize(total, 0);
    reader.read_exact(&mut block[18..])?;
    Ok(Some(block))
}

/// Decompress one complete BGZF block into its payload bytes. Uses
/// `flate2::read::MultiGzDecoder` because each BGZF block is a complete
/// gzip member.
fn decompress_bgzf_block(block: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(64 * 1024);
    MultiGzDecoder::new(block).read_to_end(&mut out)?;
    Ok(out)
}

/// Parse one complete VCF data line (bytes, no trailing newline) and, on
/// success, append one `VrsResult` per non-symbolic ALT allele to `out`.
/// Errors are pushed as `Err` entries to preserve ordering.
fn process_vcf_line_bytes(
    line: &[u8],
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    chrom_accessions: &HashMap<String, String>,
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
    let mut fields = line.splitn(10, |&b| b == b'\t');
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
        None => {
            out.push(Err(anyhow::anyhow!("Invalid POS at chrom {}", chrom)));
            return;
        }
    };
    let ref_str = match std::str::from_utf8(ref_b) {
        Ok(s) => s,
        Err(_) => return,
    };
    let alt_field = match std::str::from_utf8(alt_b) {
        Ok(s) => s,
        Err(_) => return,
    };

    let seq_accession = match chrom_accessions.get(chrom) {
        Some(a) => a.as_str(),
        None => return,
    };
    let raw_digest = match name_to_digest.get(chrom) {
        Some(d) => d,
        None => return,
    };
    let sequence = match store.sequence_bytes(raw_digest.as_str()) {
        Some(b) => b,
        None => {
            out.push(Err(anyhow::anyhow!(
                "Chromosome {} not decoded in store",
                chrom
            )));
            return;
        }
    };

    for alt in alt_field.split(',') {
        if alt.starts_with('<') || alt == "*" || alt == "." {
            continue;
        }
        let norm = match normalize(sequence, pos, ref_str.as_bytes(), alt.as_bytes()) {
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
        let vrs_id = digest_writer.allele_identifier_literal(
            seq_accession,
            norm.start,
            norm.end,
            norm_seq,
        );
        out.push(Ok(VrsResult {
            chrom: CompactString::from(chrom),
            pos,
            ref_allele: CompactString::from(ref_str),
            alt_allele: CompactString::from(alt),
            vrs_id,
        }));
    }
}

/// Decompress one block, split into head_fragment / complete_lines /
/// tail_fragment, and process the complete lines.
fn process_bgzf_block(
    batch_id: usize,
    raw_block: &[u8],
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    chrom_accessions: &HashMap<String, String>,
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
        // BGZF EOF block or empty block — contributes nothing.
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

    // If the block has no newline at all, it sits entirely inside a
    // single VCF line that spans multiple BGZF blocks (common with very
    // long INFO fields). Emit nothing; stash the whole decompressed
    // payload in `tail_fragment` so the collector accumulates it into
    // `prev_tail` without treating it as a complete line.
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

    // Block has ≥1 newline: split head / middle / tail. The middle slice
    // is empty when `first == last` (single newline); `split(b'\n')` on
    // an empty slice yields one empty chunk which is filtered below.
    let head_fragment = decompressed[..first].to_vec();
    let tail_fragment = decompressed[last + 1..].to_vec();
    let middle_slice: &[u8] = &decompressed[first + 1..last + 1];

    let mut results: Vec<Result<VrsResult>> = Vec::new();
    for line in middle_slice.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        process_vcf_line_bytes(
            line,
            store,
            name_to_digest,
            chrom_accessions,
            digest_writer,
            &mut results,
        );
    }

    BgzfBlockResult {
        batch_id,
        head_fragment,
        tail_fragment,
        results,
        had_newline: true,
    }
}

/// Streaming core for the BGZF-block parallel path. Invokes `on_result`
/// for each `VrsResult` in VCF order via a per-block reorder buffer. The
/// collector still needs to stitch cross-block boundary lines, so it
/// holds one block's tail_fragment at a time — constant extra memory.
pub fn compute_vrs_ids_parallel_bgzf_with_sink<F: FnMut(VrsResult)>(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    num_workers: usize,
    mut on_result: F,
) -> Result<usize> {
    let mut file = File::open(vcf_path).context(format!("Failed to open VCF: {}", vcf_path))?;
    if !is_bgzf(&mut file)? {
        return Err(anyhow::anyhow!(
            "compute_vrs_ids_parallel_bgzf_with_sink requires BGZF input; use compute_vrs_ids_parallel_with_sink for plain gzip"
        ));
    }

    let num_workers = num_workers.max(1);
    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(name, digest)| (name.clone(), format!("SQ.{}", digest)))
        .collect();

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
            let store_ref: &ReadonlyRefgetStore = store;
            let name_to_digest_ref = name_to_digest;
            let chrom_accessions_ref = &chrom_accessions;
            let handle = s.spawn(move || {
                let mut digest_writer = DigestWriter::new();
                while let Ok((batch_id, raw_block)) = block_rx.recv() {
                    let br = process_bgzf_block(
                        batch_id,
                        &raw_block,
                        store_ref,
                        name_to_digest_ref,
                        chrom_accessions_ref,
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

        // Streaming collector. Blocks arrive from workers in arbitrary
        // order; buffer out-of-order ones in a HashMap by batch_id and
        // process the contiguous run starting at `next_batch_id`. Each
        // processed block's results are emitted through the sink and
        // dropped, so only one block's tail_fragment and the current
        // reorder gap sit in memory.
        let mut first_err: Option<anyhow::Error> = None;
        let mut count: usize = 0;
        let mut prev_tail: Vec<u8> = Vec::new();
        let mut stitch_writer = DigestWriter::new();
        let mut pending: HashMap<usize, BgzfBlockResult> = HashMap::new();
        let mut next_batch_id: usize = 0;

        fn process_block<F: FnMut(VrsResult)>(
            block: BgzfBlockResult,
            next_batch_id: usize,
            prev_tail: &mut Vec<u8>,
            stitch_writer: &mut DigestWriter,
            store: &ReadonlyRefgetStore,
            name_to_digest: &HashMap<String, String>,
            chrom_accessions: &HashMap<String, String>,
            on_result: &mut F,
            count: &mut usize,
            first_err: &mut Option<anyhow::Error>,
        ) {
            if !block.had_newline {
                // Block sits entirely inside a single VCF line that
                // spans multiple BGZF blocks. Accumulate bytes into
                // prev_tail; the line will complete when a later block
                // produces its first `\n`. Surface any worker errors.
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

            // Block contains ≥1 newline. Stitch prev.tail + this.head
            // into the line that straddles the boundary and process it.
            // For next_batch_id == 0 with empty prev_tail the "stitched"
            // bytes are typically the header line.
            let mut stitched = std::mem::take(prev_tail);
            stitched.extend_from_slice(&block.head_fragment);
            if !stitched.is_empty() || next_batch_id > 0 {
                let mut tmp: Vec<Result<VrsResult>> = Vec::new();
                process_vcf_line_bytes(
                    &stitched,
                    store,
                    name_to_digest,
                    chrom_accessions,
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
                    store,
                    name_to_digest,
                    &chrom_accessions,
                    &mut on_result,
                    &mut count,
                    &mut first_err,
                );
                next_batch_id += 1;
            }
        }
        // Drain anything left (only if channel closed with gaps, which
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
                    store,
                    name_to_digest,
                    &chrom_accessions,
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
                store,
                name_to_digest,
                &chrom_accessions,
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

// ── Streaming TSV output entrypoints ────────────────────────────────────
//
// These pair the parallel compute paths with a direct-to-disk TSV writer
// so the peak result set never materializes as a `Vec<VrsResult>`. For a
// run with 11.6M chr22 variants this drops ~2–3 GB of peak RSS.

/// Shared TSV-writer helper. Opens the output file, writes the header row,
/// returns a `BufWriter` the caller can push lines into.
fn open_tsv_writer(output_path: &str) -> Result<std::io::BufWriter<std::fs::File>> {
    use std::io::Write;
    let file = std::fs::File::create(output_path)
        .context(format!("Failed to create TSV output {}", output_path))?;
    let mut w = std::io::BufWriter::with_capacity(1 << 20, file);
    writeln!(w, "chrom\tpos\tref\talt\tvrs_id")?;
    Ok(w)
}

/// Format a single `VrsResult` as a TSV line (chrom\tpos\tref\talt\tvrs_id\n)
/// appended to `buf`. Uses `itoa` for the integer `pos` and byte-level
/// `extend_from_slice` instead of the `writeln!` formatter machinery.
/// Called from worker threads in the `_to_tsv` fast paths so the
/// formatting cost (~20 s on 11.6M chr22 variants) is parallelized
/// instead of starving workers on single-threaded collector format.
#[inline]
fn append_tsv_line(buf: &mut Vec<u8>, r: &VrsResult) {
    let mut itoa_buf = itoa::Buffer::new();
    buf.extend_from_slice(r.chrom.as_bytes());
    buf.push(b'\t');
    buf.extend_from_slice(itoa_buf.format(r.pos).as_bytes());
    buf.push(b'\t');
    buf.extend_from_slice(r.ref_allele.as_bytes());
    buf.push(b'\t');
    buf.extend_from_slice(r.alt_allele.as_bytes());
    buf.push(b'\t');
    buf.extend_from_slice(r.vrs_id.as_bytes());
    buf.push(b'\n');
}

/// Per-batch TSV output from a parallel-TSV worker. `bytes` is a
/// ready-to-write TSV chunk (lines already terminated with `\n`).
/// Errors are collected separately so the collector can surface the
/// first one without needing to touch `bytes`.
struct LineResultTsvBatch {
    batch_id: usize,
    bytes: Vec<u8>,
    errors: Vec<anyhow::Error>,
}

/// Per-block TSV output from a BGZF-TSV worker. Mirrors `BgzfBlockResult`
/// but with pre-formatted TSV bytes instead of `Vec<Result<VrsResult>>`.
/// The collector still needs the structured `boundary_results` field to
/// handle the stitched cross-block line (worker can't format that line
/// itself — it spans two blocks).
struct BgzfBlockTsvResult {
    batch_id: usize,
    head_fragment: Vec<u8>,
    tail_fragment: Vec<u8>,
    bytes: Vec<u8>,
    errors: Vec<anyhow::Error>,
    had_newline: bool,
}

/// Compute VRS IDs via the BGZF-block parallel path and stream them
/// directly to a TSV (chrom, pos, ref, alt, vrs_id). Returns the number
/// of variants written. Requires BGZF input.
///
/// Unlike the `_with_sink` variant, formatting of TSV bytes happens on
/// worker threads (in parallel), so the collector only reassembles the
/// pre-formatted bytes in VCF order and writes them. This is critical
/// for throughput on large VCFs (gnomAD chr22: ~500k v/s vs. ~180k v/s
/// for the serial-format sink path).
pub fn compute_vrs_ids_parallel_bgzf_to_tsv(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    output_path: &str,
    num_workers: usize,
) -> Result<usize> {
    use std::io::Write;
    let mut file = File::open(vcf_path).context(format!("Failed to open VCF: {}", vcf_path))?;
    if !is_bgzf(&mut file)? {
        return Err(anyhow::anyhow!(
            "compute_vrs_ids_parallel_bgzf_to_tsv requires BGZF input; use compute_vrs_ids_parallel_to_tsv for plain gzip"
        ));
    }

    let num_workers = num_workers.max(1);
    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(name, digest)| (name.clone(), format!("SQ.{}", digest)))
        .collect();

    let mut w = open_tsv_writer(output_path)?;

    let count = std::thread::scope(|s| -> Result<usize> {
        let (block_tx, block_rx) =
            crossbeam_channel::bounded::<(usize, Vec<u8>)>(PARALLEL_CHANNEL_CAPACITY);
        let (result_tx, result_rx) =
            crossbeam_channel::bounded::<BgzfBlockTsvResult>(PARALLEL_CHANNEL_CAPACITY);

        // Reader thread: raw block I/O only.
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

        // Workers: decompress + parse + normalize + digest + format-to-TSV.
        let mut worker_handles = Vec::with_capacity(num_workers);
        for _ in 0..num_workers {
            let block_rx = block_rx.clone();
            let result_tx = result_tx.clone();
            let store_ref: &ReadonlyRefgetStore = store;
            let name_to_digest_ref = name_to_digest;
            let chrom_accessions_ref = &chrom_accessions;
            let handle = s.spawn(move || {
                let mut digest_writer = DigestWriter::new();
                let mut scratch_results: Vec<Result<VrsResult>> = Vec::new();
                while let Ok((batch_id, raw_block)) = block_rx.recv() {
                    let br = process_bgzf_block(
                        batch_id,
                        &raw_block,
                        store_ref,
                        name_to_digest_ref,
                        chrom_accessions_ref,
                        &mut digest_writer,
                    );
                    // Reformat the block's structured results into TSV
                    // bytes here (on the worker) so the collector never
                    // has to do per-result formatting serially.
                    let mut bytes: Vec<u8> = Vec::with_capacity(br.results.len() * 80);
                    let mut errors: Vec<anyhow::Error> = Vec::new();
                    scratch_results.clear();
                    scratch_results.extend(br.results);
                    for r in scratch_results.drain(..) {
                        match r {
                            Ok(v) => append_tsv_line(&mut bytes, &v),
                            Err(e) => errors.push(e),
                        }
                    }
                    let tsv_br = BgzfBlockTsvResult {
                        batch_id: br.batch_id,
                        head_fragment: br.head_fragment,
                        tail_fragment: br.tail_fragment,
                        bytes,
                        errors,
                        had_newline: br.had_newline,
                    };
                    if result_tx.send(tsv_br).is_err() {
                        break;
                    }
                }
            });
            worker_handles.push(handle);
        }

        drop(block_rx);
        drop(result_tx);

        // Collector: reassemble blocks in order, stitch boundary lines,
        // write pre-formatted bytes to disk.
        let mut first_err: Option<anyhow::Error> = None;
        let mut count: usize = 0;
        let mut prev_tail: Vec<u8> = Vec::new();
        let mut stitch_writer = DigestWriter::new();
        let mut pending: HashMap<usize, BgzfBlockTsvResult> = HashMap::new();
        let mut next_batch_id: usize = 0;
        let mut io_err: Option<std::io::Error> = None;

        let mut process_tsv_block = |block: BgzfBlockTsvResult,
                                     next_batch_id: usize,
                                     prev_tail: &mut Vec<u8>,
                                     stitch_writer: &mut DigestWriter,
                                     count: &mut usize,
                                     first_err: &mut Option<anyhow::Error>,
                                     io_err: &mut Option<std::io::Error>|
         -> () {
            for e in block.errors {
                if first_err.is_none() {
                    *first_err = Some(e);
                }
            }
            if !block.had_newline {
                // Interior of a long VCF line: just accumulate bytes.
                prev_tail.extend_from_slice(&block.tail_fragment);
                return;
            }
            // Stitch prev_tail + head_fragment, process as one line.
            let mut stitched = std::mem::take(prev_tail);
            stitched.extend_from_slice(&block.head_fragment);
            if !stitched.is_empty() || next_batch_id > 0 {
                let mut tmp: Vec<Result<VrsResult>> = Vec::new();
                process_vcf_line_bytes(
                    &stitched,
                    store,
                    name_to_digest,
                    &chrom_accessions,
                    stitch_writer,
                    &mut tmp,
                );
                // Format the stitched line's results into a local buffer
                // and write it, then write the worker's pre-formatted
                // block bytes. Keeping order is critical: the stitched
                // line is the *first* line of this block, so it writes
                // before `block.bytes`.
                let mut stitched_bytes: Vec<u8> = Vec::new();
                for r in tmp {
                    match r {
                        Ok(v) => {
                            append_tsv_line(&mut stitched_bytes, &v);
                            *count += 1;
                        }
                        Err(e) => {
                            if first_err.is_none() {
                                *first_err = Some(e);
                            }
                        }
                    }
                }
                if io_err.is_none() {
                    if let Err(e) = w.write_all(&stitched_bytes) {
                        *io_err = Some(e);
                    }
                }
            }
            if io_err.is_none() {
                if let Err(e) = w.write_all(&block.bytes) {
                    *io_err = Some(e);
                }
            }
            // bytes.len() / avg_line_len isn't exact — count worker
            // results accurately by counting `\n` in the bytes buffer.
            *count += memchr_count_newlines(&block.bytes);
            *prev_tail = block.tail_fragment;
        };

        while let Ok(br) = result_rx.recv() {
            pending.insert(br.batch_id, br);
            while let Some(block) = pending.remove(&next_batch_id) {
                process_tsv_block(
                    block,
                    next_batch_id,
                    &mut prev_tail,
                    &mut stitch_writer,
                    &mut count,
                    &mut first_err,
                    &mut io_err,
                );
                next_batch_id += 1;
            }
        }
        // Drain anything left (shouldn't happen in practice).
        let mut remaining_ids: Vec<usize> = pending.keys().copied().collect();
        remaining_ids.sort_unstable();
        for id in remaining_ids {
            if let Some(block) = pending.remove(&id) {
                process_tsv_block(
                    block,
                    id,
                    &mut prev_tail,
                    &mut stitch_writer,
                    &mut count,
                    &mut first_err,
                    &mut io_err,
                );
            }
        }

        // Final unterminated line.
        if !prev_tail.is_empty() {
            let mut tmp: Vec<Result<VrsResult>> = Vec::new();
            process_vcf_line_bytes(
                &prev_tail,
                store,
                name_to_digest,
                &chrom_accessions,
                &mut stitch_writer,
                &mut tmp,
            );
            let mut tail_bytes: Vec<u8> = Vec::new();
            for r in tmp {
                match r {
                    Ok(v) => {
                        append_tsv_line(&mut tail_bytes, &v);
                        count += 1;
                    }
                    Err(e) => {
                        if first_err.is_none() {
                            first_err = Some(e);
                        }
                    }
                }
            }
            if io_err.is_none() {
                if let Err(e) = w.write_all(&tail_bytes) {
                    io_err = Some(e);
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

        if let Some(e) = io_err {
            return Err(e.into());
        }
        if let Some(e) = first_err {
            return Err(e);
        }
        Ok(count)
    })?;

    w.flush()?;
    Ok(count)
}

/// Count of `\n` bytes in `buf`. Used to tally TSV lines produced by a
/// worker when the collector only sees pre-formatted bytes.
#[inline]
fn memchr_count_newlines(buf: &[u8]) -> usize {
    buf.iter().filter(|&&b| b == b'\n').count()
}

/// Compute VRS IDs via the parallel path and stream them
/// directly to a TSV (chrom, pos, ref, alt, vrs_id). Returns the number
/// of variants written.
///
/// Like `compute_vrs_ids_parallel_bgzf_to_tsv`, TSV formatting happens
/// on worker threads (parallel) so the collector just reassembles bytes
/// in VCF order.
pub fn compute_vrs_ids_parallel_to_tsv(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    output_path: &str,
    num_workers: usize,
) -> Result<usize> {
    use std::io::Write;
    let num_workers = num_workers.max(1);
    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(name, digest)| (name.clone(), format!("SQ.{}", digest)))
        .collect();

    let mut w = open_tsv_writer(output_path)?;

    let count = std::thread::scope(|s| -> Result<usize> {
        let (line_tx, line_rx) =
            crossbeam_channel::bounded::<LineBatch>(PARALLEL_CHANNEL_CAPACITY);
        let (result_tx, result_rx) =
            crossbeam_channel::bounded::<LineResultTsvBatch>(PARALLEL_CHANNEL_CAPACITY);

        let reader_handle = s.spawn(move || -> Result<()> {
            let mut reader = open_vcf(vcf_path)?;
            let mut batch: Vec<String> = Vec::with_capacity(PARALLEL_BATCH_SIZE);
            let mut batch_id: usize = 0;
            let mut scratch = String::new();
            loop {
                scratch.clear();
                if !read_vcf_line(&mut *reader, &mut scratch)? {
                    break;
                }
                let trimmed = scratch.trim_end_matches('\n').trim_end_matches('\r');
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    continue;
                }
                let line = std::mem::take(&mut scratch);
                batch.push(line);
                if batch.len() >= PARALLEL_BATCH_SIZE {
                    let full = std::mem::replace(
                        &mut batch,
                        Vec::with_capacity(PARALLEL_BATCH_SIZE),
                    );
                    if line_tx.send(LineBatch { batch_id, lines: full }).is_err() {
                        return Ok(());
                    }
                    batch_id += 1;
                }
            }
            if !batch.is_empty() {
                let _ = line_tx.send(LineBatch { batch_id, lines: batch });
            }
            Ok(())
        });

        let mut worker_handles = Vec::with_capacity(num_workers);
        for _ in 0..num_workers {
            let line_rx = line_rx.clone();
            let result_tx = result_tx.clone();
            let store_ref: &ReadonlyRefgetStore = store;
            let name_to_digest_ref = name_to_digest;
            let chrom_accessions_ref = &chrom_accessions;
            let handle = s.spawn(move || {
                let mut digest_writer = DigestWriter::new();
                while let Ok(line_batch) = line_rx.recv() {
                    let mut results: Vec<Result<VrsResult>> =
                        Vec::with_capacity(line_batch.lines.len());
                    for raw in &line_batch.lines {
                        let line = raw.trim_end_matches('\n').trim_end_matches('\r');
                        process_vcf_line_bytes(
                            line.as_bytes(),
                            store_ref,
                            name_to_digest_ref,
                            chrom_accessions_ref,
                            &mut digest_writer,
                            &mut results,
                        );
                    }
                    // Parallel-format into bytes on this worker.
                    let mut bytes: Vec<u8> = Vec::with_capacity(results.len() * 80);
                    let mut errors: Vec<anyhow::Error> = Vec::new();
                    for r in results {
                        match r {
                            Ok(v) => append_tsv_line(&mut bytes, &v),
                            Err(e) => errors.push(e),
                        }
                    }
                    if result_tx
                        .send(LineResultTsvBatch {
                            batch_id: line_batch.batch_id,
                            bytes,
                            errors,
                        })
                        .is_err()
                    {
                        break;
                    }
                }
            });
            worker_handles.push(handle);
        }

        drop(line_rx);
        drop(result_tx);

        let mut first_err: Option<anyhow::Error> = None;
        let mut io_err: Option<std::io::Error> = None;
        let mut pending: HashMap<usize, LineResultTsvBatch> = HashMap::new();
        let mut next_batch_id: usize = 0;
        let mut count: usize = 0;

        let mut emit = |batch: LineResultTsvBatch,
                        count: &mut usize,
                        first_err: &mut Option<anyhow::Error>,
                        io_err: &mut Option<std::io::Error>| {
            for e in batch.errors {
                if first_err.is_none() {
                    *first_err = Some(e);
                }
            }
            *count += memchr_count_newlines(&batch.bytes);
            if io_err.is_none() {
                if let Err(e) = w.write_all(&batch.bytes) {
                    *io_err = Some(e);
                }
            }
        };

        while let Ok(rb) = result_rx.recv() {
            pending.insert(rb.batch_id, rb);
            while let Some(batch) = pending.remove(&next_batch_id) {
                emit(batch, &mut count, &mut first_err, &mut io_err);
                next_batch_id += 1;
            }
        }
        let mut remaining_ids: Vec<usize> = pending.keys().copied().collect();
        remaining_ids.sort_unstable();
        for id in remaining_ids {
            if let Some(batch) = pending.remove(&id) {
                emit(batch, &mut count, &mut first_err, &mut io_err);
            }
        }

        match reader_handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => {
                if first_err.is_none() {
                    first_err = Some(e);
                }
            }
            Err(_) => {
                if first_err.is_none() {
                    first_err = Some(anyhow::anyhow!("Reader thread panicked"));
                }
            }
        }
        for h in worker_handles {
            if h.join().is_err() && first_err.is_none() {
                first_err = Some(anyhow::anyhow!("VRS worker thread panicked"));
            }
        }

        if let Some(e) = io_err {
            return Err(e.into());
        }
        if let Some(e) = first_err {
            return Err(e);
        }
        Ok(count)
    })?;

    w.flush()?;
    Ok(count)
}

// ── Lazy per-chromosome decode helper ───────────────────────────────────

/// Pre-scan a VCF and return the unique chromosome names observed in
/// column 1 of data lines. Used by the Python bindings to decode only the
/// chromosomes the VCF actually references instead of the whole
/// collection (for a single-chromosome VCF this saves ~9 GB of RAM on
/// GRCh38). Handles plain/gzip/BGZF input via `open_vcf`.
pub fn unique_chroms_from_vcf(vcf_path: &str) -> Result<Vec<String>> {
    use std::collections::HashSet;
    let mut reader = open_vcf(vcf_path)?;
    let mut line_buf = String::new();
    let mut seen: HashSet<String> = HashSet::new();
    while read_vcf_line(&mut *reader, &mut line_buf)? {
        let line = line_buf.trim_end_matches('\n').trim_end_matches('\r');
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // First field up to TAB is the chromosome.
        let chrom = match line.split_once('\t') {
            Some((c, _)) => c,
            None => continue,
        };
        if !seen.contains(chrom) {
            seen.insert(chrom.to_string());
        }
    }
    Ok(seen.into_iter().collect())
}

/// Decode only the collection sequences referenced by the chromosome
/// column of a VCF. For a single-chromosome VCF against GRCh38 this
/// avoids ~3-9 GB of RSS that eager `load_all_sequences` +
/// full-collection `ensure_decoded` would use.
///
/// Three paths:
///
/// 1. **Fast path 1 (all decoded):** if every sequence in the collection
///    is already decoded, return `Ok(())` without scanning the VCF.
/// 2. **Fast path 2 (all loaded, not decoded):** if every sequence is
///    loaded but not yet decoded (e.g. the caller just called
///    `load_all_sequences`), decode them all without scanning the VCF.
///    Scanning a BGZF-compressed VCF to extract the chromosome column
///    can take 30-60 s, so this fast path is substantially cheaper than
///    the slow path when it applies.
/// 3. **Slow path:** scan the VCF for unique chromosome names via
///    [`unique_chroms_from_vcf`], then `load_sequence` + `ensure_decoded`
///    for each referenced chrom. Chroms not in the collection are
///    silently skipped (matches the compute path behavior — see
///    `build_name_to_digest_readonly` and
///    `compute_vrs_ids_parallel_bgzf_with_sink`).
pub fn decode_vcf_chroms(
    store: &mut RefgetStore,
    collection_digest: &str,
    vcf_path: &str,
) -> Result<()> {
    // Use the mutable wrapper's lazy-loading get_collection so the
    // collection data is auto-loaded if it was only a Stub. Build
    // name→digest from its sequences.
    let collection = store
        .get_collection(collection_digest)
        .context("Failed to get sequence collection")?;
    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    let mut collection_digests: Vec<String> = Vec::with_capacity(collection.sequences.len());
    for seq_record in &collection.sequences {
        let meta = seq_record.metadata();
        name_to_digest.insert(meta.name.clone(), meta.sha512t24u.clone());
        collection_digests.push(meta.sha512t24u.clone());
    }

    // Fast path 1: if every sequence in the collection is already
    // decoded, skip both the VCF scan and the decode. Cheap HashSet
    // lookup per sequence via `is_sequence_decoded`.
    let all_decoded = !collection_digests.is_empty()
        && collection_digests
            .iter()
            .all(|d| store.is_sequence_decoded(d.as_str()));
    if all_decoded {
        return Ok(());
    }

    // Fast path 2: if every sequence is loaded (Full, just not
    // decoded), skip the VCF scan and decode them all in place. The
    // `load_all_sequences()` preload puts the store in this state,
    // and decoding everything is cheap (~1 s for GRCh38) compared
    // to scanning a BGZF-compressed VCF (~30-50 s). The memory cost
    // of decoding every chromosome instead of only the VCF-referenced
    // ones is only ~3 GB of extra RAM for GRCh38 — a deliberate
    // trade the caller opted into by calling `load_all_sequences()`.
    let all_loaded = !collection_digests.is_empty()
        && collection_digests
            .iter()
            .all(|d| store.is_sequence_loaded(d.as_str()));
    if all_loaded {
        for digest in &collection_digests {
            store
                .ensure_decoded(digest.as_str())
                .with_context(|| format!("Failed to decode sequence {}", digest))?;
        }
        return Ok(());
    }

    // Slow path: scan the VCF for unique chromosome names and
    // load+decode only those sequences.
    let chroms =
        unique_chroms_from_vcf(vcf_path).context("Failed to scan VCF for chromosomes")?;

    for chrom in &chroms {
        if let Some(digest) = name_to_digest.get(chrom) {
            let digest = digest.clone();
            // Lazy load the sequence bytes (Stub → Full) from disk
            // first if not already loaded. Then decode in place.
            store
                .load_sequence(digest.as_str())
                .with_context(|| format!("Failed to load sequence {}", digest))?;
            store
                .ensure_decoded(digest.as_str())
                .with_context(|| format!("Failed to decode sequence {}", digest))?;
        }
    }
    Ok(())
}
