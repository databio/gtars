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

use crate::digest::DigestWriter;
use crate::normalize::normalize;

/// Result of computing a VRS identifier for a single VCF variant.
#[derive(Debug, Clone)]
pub struct VrsResult {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub vrs_id: String,
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
/// decompression is still available via `compute_vrs_ids_parallel_bgzf`,
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
                chrom: chrom.to_string(),
                pos,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt.to_string(),
                vrs_id,
            });
            count += 1;
        }
    }

    store.clear_decoded_cache();
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

// ── Parallel API (reader thread + worker pool) ──────────────────────────

/// A parsed VCF variant ready for normalization and digest computation.
/// Sent from the reader thread to worker threads via channel.
struct ParsedVariant {
    /// Monotonic 0-based index assigned by the reader (one per emitted ALT).
    /// Used to restore VCF order in the final output.
    index: usize,
    chrom: String,
    pos: u64,
    ref_allele: String,
    alt_allele: String,
    /// Pre-resolved "SQ.xxx" accession for this chromosome.
    seq_accession: String,
    /// SHA512t24u digest key for sequence_bytes() lookup.
    raw_digest: String,
}

/// A computed VRS result tagged with its original input index for reordering.
struct IndexedVrsResult {
    index: usize,
    result: VrsResult,
}

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

/// Compute VRS IDs for all variants in a VCF using a reader thread plus a
/// pool of worker threads, all operating on a preloaded `ReadonlyRefgetStore`.
///
/// The reader thread parses the VCF sequentially (decompression and line
/// parsing are inherently serial) and forwards `ParsedVariant`s to workers via
/// a bounded crossbeam channel. Each worker holds its own `DigestWriter` and
/// computes VRS IDs independently. Results carry input indices so the final
/// `Vec<VrsResult>` is reordered to match VCF line order for determinism.
///
/// Preconditions: every chromosome referenced by the VCF must already be
/// decoded in `store`. Callers typically preload by calling
/// `load_all_sequences` on the underlying mutable store before converting
/// it with `into_readonly`.
pub fn compute_vrs_ids_parallel(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    num_workers: usize,
) -> Result<Vec<VrsResult>> {
    let num_workers = num_workers.max(1);

    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(name, digest)| (name.clone(), format!("SQ.{}", digest)))
        .collect();

    std::thread::scope(|s| -> Result<Vec<VrsResult>> {
        // Batching: channels carry Vec<ParsedVariant> / Vec<Result<IndexedVrsResult>>
        // to amortize channel overhead. With per-variant messaging the channel
        // send/recv cost dominates the actual digest work (~5 us/variant),
        // killing parallel scaling. Batches of 1024 drop channel overhead to
        // <0.1% of per-variant work.
        let (variant_tx, variant_rx) =
            crossbeam_channel::bounded::<Vec<ParsedVariant>>(PARALLEL_CHANNEL_CAPACITY);
        let (result_tx, result_rx) = crossbeam_channel::bounded::<
            Vec<Result<IndexedVrsResult>>,
        >(PARALLEL_CHANNEL_CAPACITY);

        // Reader thread: sequential VCF parse, batch ParsedVariants.
        let reader_handle = s.spawn({
            let name_to_digest = name_to_digest;
            let chrom_accessions = &chrom_accessions;
            move || -> Result<()> {
                let mut reader = open_vcf(vcf_path)?;
                let mut line_buf = String::new();
                let mut index: usize = 0;
                let mut batch: Vec<ParsedVariant> = Vec::with_capacity(PARALLEL_BATCH_SIZE);

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

                    for alt in alt_field.split(',') {
                        if alt.starts_with('<') || alt == "*" || alt == "." {
                            continue;
                        }
                        batch.push(ParsedVariant {
                            index,
                            chrom: chrom.to_string(),
                            pos,
                            ref_allele: ref_allele.to_string(),
                            alt_allele: alt.to_string(),
                            seq_accession: seq_accession.clone(),
                            raw_digest: raw_digest.clone(),
                        });
                        index += 1;

                        if batch.len() >= PARALLEL_BATCH_SIZE {
                            let full = std::mem::replace(
                                &mut batch,
                                Vec::with_capacity(PARALLEL_BATCH_SIZE),
                            );
                            if variant_tx.send(full).is_err() {
                                return Ok(());
                            }
                        }
                    }
                }

                if !batch.is_empty() && variant_tx.send(batch).is_err() {
                    return Ok(());
                }
                Ok(())
            }
        });

        // Worker threads: process batches, return a batch of results.
        let mut worker_handles = Vec::with_capacity(num_workers);
        for _ in 0..num_workers {
            let variant_rx = variant_rx.clone();
            let result_tx = result_tx.clone();
            let store_ref: &ReadonlyRefgetStore = store;
            let handle = s.spawn(move || {
                let mut digest_writer = DigestWriter::new();
                while let Ok(pv_batch) = variant_rx.recv() {
                    let mut out_batch: Vec<Result<IndexedVrsResult>> =
                        Vec::with_capacity(pv_batch.len());
                    for pv in pv_batch {
                        let r: Result<IndexedVrsResult> = (|| {
                            let sequence =
                                store_ref.sequence_bytes(pv.raw_digest.as_str()).context(
                                    format!("Chromosome {} not decoded in store", pv.chrom),
                                )?;
                            let norm = normalize(
                                sequence,
                                pv.pos,
                                pv.ref_allele.as_bytes(),
                                pv.alt_allele.as_bytes(),
                            )
                            .context(format!(
                                "Failed to normalize variant at {}:{}",
                                pv.chrom,
                                pv.pos + 1
                            ))?;
                            let norm_seq = std::str::from_utf8(&norm.allele).context(format!(
                                "Normalized allele is not valid UTF-8 at {}:{}",
                                pv.chrom,
                                pv.pos + 1
                            ))?;
                            let vrs_id = digest_writer.allele_identifier_literal(
                                &pv.seq_accession,
                                norm.start,
                                norm.end,
                                norm_seq,
                            );
                            Ok(IndexedVrsResult {
                                index: pv.index,
                                result: VrsResult {
                                    chrom: pv.chrom,
                                    pos: pv.pos,
                                    ref_allele: pv.ref_allele,
                                    alt_allele: pv.alt_allele,
                                    vrs_id,
                                },
                            })
                        })();
                        out_batch.push(r);
                    }
                    if result_tx.send(out_batch).is_err() {
                        break;
                    }
                }
            });
            worker_handles.push(handle);
        }

        // Drop our local handles so the channels close once reader and all
        // workers finish.
        drop(variant_rx);
        drop(result_tx);

        // Collect result batches as they arrive.
        let mut indexed: Vec<IndexedVrsResult> = Vec::new();
        let mut first_err: Option<anyhow::Error> = None;
        while let Ok(out_batch) = result_rx.recv() {
            for r in out_batch {
                match r {
                    Ok(ir) => indexed.push(ir),
                    Err(e) => {
                        if first_err.is_none() {
                            first_err = Some(e);
                        }
                    }
                }
            }
        }

        // Surface reader errors.
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
        // Workers have no Result return; join to surface panics.
        for h in worker_handles {
            if h.join().is_err() && first_err.is_none() {
                first_err = Some(anyhow::anyhow!("VRS worker thread panicked"));
            }
        }

        if let Some(e) = first_err {
            return Err(e);
        }

        indexed.sort_by_key(|i| i.index);
        Ok(indexed.into_iter().map(|i| i.result).collect())
    })
}

/// Convenience wrapper that picks `num_workers` automatically
/// (`available_parallelism - 2`, min 1).
pub fn compute_vrs_ids_parallel_auto(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
) -> Result<Vec<VrsResult>> {
    let num_workers = std::thread::available_parallelism()
        .map(|n| n.get().saturating_sub(2).max(1))
        .unwrap_or(1);
    compute_vrs_ids_parallel(store, name_to_digest, vcf_path, num_workers)
}

// ── Blockwise parallel API (workers parse) ──────────────────────────────
//
// The pipeline variant above keeps all VCF parsing (field split, POS parse,
// HashMap lookup, ParsedVariant construction with its 5 String allocations)
// on the reader thread. That thread becomes the new bottleneck once
// normalize+digest is offloaded, capping speedup at ~1.3× in practice.
//
// The blockwise variant below batches *raw* VCF lines on the reader and
// hands the parse+normalize+digest pipeline to the workers in its entirety.
// The reader's per-line work collapses to: read_line into a String, skip
// headers, push into a batch Vec. That's ~1 allocation per line, and every
// hot-path operation that was serialized (field split, parse::<u64>,
// HashMap lookup, multi-alt expansion, normalize, digest) now runs across
// the worker pool.

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

/// Compute VRS IDs by batching raw VCF lines on the reader and running the
/// full parse+normalize+digest pipeline on worker threads.
///
/// In contrast to `compute_vrs_ids_parallel`, the reader thread here only
/// reads lines and enqueues them. All VCF field parsing, chromosome
/// lookups, normalization, and digest computation happen on workers. This
/// removes the serial parse/alloc bottleneck observed in the pipeline
/// variant and scales further with core count.
///
/// Preconditions: identical to `compute_vrs_ids_parallel` — every
/// chromosome referenced by the VCF must already be decoded in `store`.
pub fn compute_vrs_ids_parallel_blockwise(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    num_workers: usize,
) -> Result<Vec<VrsResult>> {
    let num_workers = num_workers.max(1);

    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(name, digest)| (name.clone(), format!("SQ.{}", digest)))
        .collect();

    std::thread::scope(|s| -> Result<Vec<VrsResult>> {
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
                                chrom: chrom.to_string(),
                                pos,
                                ref_allele: ref_allele.to_string(),
                                alt_allele: alt.to_string(),
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

        // Collect all result batches, then reorder by batch_id to restore
        // VCF line order. Memory usage is bounded by total results — same
        // as the pipeline variant.
        let mut batches: Vec<LineResultBatch> = Vec::new();
        while let Ok(rb) = result_rx.recv() {
            batches.push(rb);
        }
        batches.sort_by_key(|b| b.batch_id);

        let mut first_err: Option<anyhow::Error> = None;
        let total: usize = batches.iter().map(|b| b.results.len()).sum();
        let mut out: Vec<VrsResult> = Vec::with_capacity(total);
        for batch in batches {
            for r in batch.results {
                match r {
                    Ok(v) => out.push(v),
                    Err(e) => {
                        if first_err.is_none() {
                            first_err = Some(e);
                        }
                    }
                }
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
        Ok(out)
    })
}

/// Convenience wrapper for the blockwise variant with auto-picked worker
/// count (`available_parallelism - 2`, min 1).
pub fn compute_vrs_ids_parallel_blockwise_auto(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
) -> Result<Vec<VrsResult>> {
    let num_workers = std::thread::available_parallelism()
        .map(|n| n.get().saturating_sub(2).max(1))
        .unwrap_or(1);
    compute_vrs_ids_parallel_blockwise(store, name_to_digest, vcf_path, num_workers)
}

// ── BGZF-block parallel API ─────────────────────────────────────────────
//
// The blockwise variant above still routes every VCF byte through one
// reader thread: it decompresses (via noodles-bgzf MT) and line-batches
// sequentially before dispatching. That reader is the new ceiling
// (~3.7× on 14 cores).
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
    if header[0] != 0x1f || header[1] != 0x8b || header[2] != 0x08 {
        return Err(anyhow::anyhow!(
            "Not BGZF: bad gzip magic at block start"
        ));
    }
    if (header[3] & 0x04) == 0 || header[12] != b'B' || header[13] != b'C' {
        return Err(anyhow::anyhow!(
            "Not BGZF: missing BC extra-field marker"
        ));
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
            chrom: chrom.to_string(),
            pos,
            ref_allele: ref_str.to_string(),
            alt_allele: alt.to_string(),
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

/// Compute VRS IDs by distributing BGZF blocks across workers. Each
/// worker handles decompression + parsing + normalize + digest for its
/// own blocks; the reader thread only shuffles raw bytes. Requires BGZF
/// input (plain gzip is rejected — the blockwise variant above is the
/// right choice there).
pub fn compute_vrs_ids_parallel_bgzf(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    num_workers: usize,
) -> Result<Vec<VrsResult>> {
    let mut file = File::open(vcf_path).context(format!("Failed to open VCF: {}", vcf_path))?;
    if !is_bgzf(&mut file)? {
        return Err(anyhow::anyhow!(
            "compute_vrs_ids_parallel_bgzf requires BGZF input; use compute_vrs_ids_parallel_blockwise for plain gzip"
        ));
    }

    let num_workers = num_workers.max(1);
    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(name, digest)| (name.clone(), format!("SQ.{}", digest)))
        .collect();

    std::thread::scope(|s| -> Result<Vec<VrsResult>> {
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

        // Collector: sort blocks, stitch cross-block boundary lines,
        // emit results in VCF order.
        let mut blocks: Vec<BgzfBlockResult> = Vec::new();
        while let Ok(br) = result_rx.recv() {
            blocks.push(br);
        }
        blocks.sort_by_key(|b| b.batch_id);

        let mut first_err: Option<anyhow::Error> = None;
        let mut all_results: Vec<VrsResult> = Vec::new();
        let mut prev_tail: Vec<u8> = Vec::new();
        let mut stitch_writer = DigestWriter::new();

        for (i, block) in blocks.into_iter().enumerate() {
            if !block.had_newline {
                // Block sits entirely inside a single VCF line that spans
                // multiple BGZF blocks. Accumulate its bytes into
                // prev_tail without emitting — the line will complete
                // (and be processed) when some later block produces its
                // first `\n`. Also surface any worker errors (e.g.
                // decompression failure); don't emit bogus results.
                for r in block.results {
                    if let Err(e) = r {
                        if first_err.is_none() {
                            first_err = Some(e);
                        }
                    }
                }
                prev_tail.extend_from_slice(&block.tail_fragment);
                continue;
            }

            // Block contains ≥1 newline. Stitch prev.tail + this.head into
            // the line that straddles the block boundary and process it.
            // For i == 0 with an empty prev_tail the "stitched" bytes are
            // just the first line of the file (typically a header).
            let mut stitched = std::mem::take(&mut prev_tail);
            stitched.extend_from_slice(&block.head_fragment);
            if !stitched.is_empty() || i > 0 {
                let mut tmp: Vec<Result<VrsResult>> = Vec::new();
                process_vcf_line_bytes(
                    &stitched,
                    store,
                    name_to_digest,
                    &chrom_accessions,
                    &mut stitch_writer,
                    &mut tmp,
                );
                for r in tmp {
                    match r {
                        Ok(v) => all_results.push(v),
                        Err(e) => {
                            if first_err.is_none() {
                                first_err = Some(e);
                            }
                        }
                    }
                }
            }
            for r in block.results {
                match r {
                    Ok(v) => all_results.push(v),
                    Err(e) => {
                        if first_err.is_none() {
                            first_err = Some(e);
                        }
                    }
                }
            }
            prev_tail = block.tail_fragment;
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
                    Ok(v) => all_results.push(v),
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
        Ok(all_results)
    })
}

/// Convenience wrapper for the BGZF-block variant with auto-picked worker
/// count (`available_parallelism - 2`, min 1).
pub fn compute_vrs_ids_parallel_bgzf_auto(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
) -> Result<Vec<VrsResult>> {
    let num_workers = std::thread::available_parallelism()
        .map(|n| n.get().saturating_sub(2).max(1))
        .unwrap_or(1);
    compute_vrs_ids_parallel_bgzf(store, name_to_digest, vcf_path, num_workers)
}
