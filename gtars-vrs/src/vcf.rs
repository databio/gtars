//! VCF parsing and VRS ID computation.
//!
//! Reads a VCF file (plain text or gzipped/bgzf), normalizes each variant,
//! and computes GA4GH VRS Allele identifiers using a fast zero-allocation
//! digest path (no serde_json per variant).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::mpsc;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use gtars_refget::digest::alphabet::lookup_alphabet;
use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore, StorageMode};

use crate::digest::DigestWriter;
use crate::normalize::{EncodedSeq, RefView, normalize_ref};

/// Build a mode-aware [`RefView`] over a resident sequence's bytes: decoded bytes
/// for a Raw store, or an on-the-fly [`EncodedSeq`] for an Encoded store. No
/// decoded cache / mmap required — the reference is read through the base accessor.
fn ref_view_for<'a>(store: &'a ReadonlyRefgetStore, raw_digest: &str) -> Result<RefView<'a>> {
    let rec = store
        .get_sequence(raw_digest)
        .context("sequence not found in store")?;
    let meta = rec.metadata();
    let bytes = rec.sequence().context("sequence not resident in store")?;
    Ok(match store.storage_mode() {
        StorageMode::Raw => RefView::Decoded(bytes),
        StorageMode::Encoded => {
            let alphabet = lookup_alphabet(&meta.alphabet);
            RefView::Encoded(EncodedSeq {
                bytes,
                length: meta.length,
                bits_per_symbol: alphabet.bits_per_symbol,
                decoding_array: alphabet.decoding_array,
            })
        }
    })
}

/// Result of computing a VRS identifier for a single VCF variant.
#[derive(Debug, Clone)]
pub struct VrsResult {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub vrs_id: String,
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
                // Ensure the encoded bytes are resident (no decode); in-memory is a no-op.
                let resident = store
                    .get_sequence(digest.as_str())
                    .map(|r| r.sequence().is_some())
                    .unwrap_or(false);
                if !resident {
                    store.load_sequence(digest.as_str())?;
                }
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

        for alt in alt_field.split(',') {
            if alt.starts_with('<') || alt == "*" || alt == "." {
                continue;
            }

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
        let view = ref_view_for(store, raw_digest.as_str())
            .context(format!("Chromosome {} not available in store", chrom))?;

        for alt in alt_field.split(',') {
            if alt.starts_with('<') || alt == "*" || alt == "." {
                continue;
            }

            let norm = normalize_ref(&view, pos, ref_allele.as_bytes(), alt.as_bytes())
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

/// Immutable per-sequence view shared (by `&`) with every worker. Holds the
/// borrowed encoded bytes plus the alphabet bits needed to decode bases on the
/// fly via [`EncodedSeq`]; nothing is decoded up front.
struct SeqView<'a> {
    /// `SQ.<digest>` accession used in the VRS identifier.
    accession: String,
    /// Bit-packed encoded bytes, borrowed from the readonly store.
    encoded: &'a [u8],
    /// Number of bases.
    length: usize,
    /// Bits per symbol (e.g. 2 for DNA_2BIT).
    bits_per_symbol: usize,
    /// Decoding table (encoded code -> symbol byte); `'static` from the alphabet.
    decoding_array: &'static [u8; 256],
}

/// Compute VRS Allele identifiers for every variant in a VCF using an in-process
/// pool of `threads` workers, reading the reference directly from the immutable
/// 2-bit-encoded readonly store (decode-on-the-fly, no decoded cache, no mmap).
///
/// One reader thread parses the VCF (plain or gzip/bgzf via [`open_vcf`]) and
/// dispatches per-record [`WorkItem`]s round-robin to the workers. Each worker
/// owns its own [`DigestWriter`] (which is not `Sync`) and builds an
/// [`EncodedSeq`] per variant from the shared `&ReadonlyRefgetStore`. Results are
/// tagged with the record's input index and reordered by a collector so the
/// callback observes records in exact input order. Within a record, alts are
/// emitted in VCF order.
///
/// The caller must make all referenced sequences resident first (e.g. via
/// `RefgetStore::load_sequence` for each needed digest, then `into_readonly`).
/// Workers only read; the `EncodedSeq` borrows are kept valid by scoped threads.
///
/// Output is identical (same VRS ids, same order) to
/// [`compute_vrs_ids_streaming_readonly`]. Returns the number of results emitted.
pub fn compute_vrs_ids_parallel_encoded(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    threads: usize,
    mut on_result_ordered: impl FnMut(VrsResult),
) -> Result<usize> {
    let n_workers = threads.max(1);

    // Build the shared per-sequence views once, borrowing encoded bytes from the
    // store. We index records by sequence position to avoid per-record hashing in
    // the workers and to give each SeqView a stable index.
    let mut seq_index_of_digest: HashMap<String, usize> = HashMap::new();
    let mut seqs: Vec<SeqView> = Vec::new();
    // Map chrom name -> seq_index, lazily filled as the reader sees new chroms.
    let mut chrom_to_seq_index: HashMap<String, usize> = HashMap::new();
    for (name, digest) in name_to_digest.iter() {
        let idx = if let Some(&i) = seq_index_of_digest.get(digest) {
            i
        } else {
            let rec = store
                .get_sequence(digest.as_str())
                .with_context(|| format!("get_sequence {digest}"))?;
            let meta = rec.metadata();
            let alphabet = lookup_alphabet(&meta.alphabet);
            let encoded = rec.sequence().with_context(|| {
                format!("sequence {digest} not resident/encoded in readonly store")
            })?;
            let i = seqs.len();
            seqs.push(SeqView {
                accession: format!("SQ.{}", meta.sha512t24u),
                encoded,
                length: meta.length,
                bits_per_symbol: alphabet.bits_per_symbol,
                decoding_array: alphabet.decoding_array,
            });
            seq_index_of_digest.insert(digest.clone(), i);
            i
        };
        chrom_to_seq_index.insert(name.clone(), idx);
    }

    let seqs = &seqs; // share by reference with the scoped workers

    let mut count = 0usize;

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
                    let view = &seqs[item.seq_index];
                    let es = EncodedSeq {
                        bytes: view.encoded,
                        length: view.length,
                        bits_per_symbol: view.bits_per_symbol,
                        decoding_array: view.decoding_array,
                    };
                    let mut results = Vec::with_capacity(item.alts.len());
                    for alt in &item.alts {
                        // A variant that fails to normalize (e.g. ref allele past
                        // sequence end) or whose normalized allele is not UTF-8 is
                        // skipped. On well-formed input every variant normalizes,
                        // so output matches the serial path exactly.
                        let norm = match normalize_ref(
                            &es,
                            item.pos,
                            item.ref_allele.as_bytes(),
                            alt.as_bytes(),
                        ) {
                            Ok(n) => n,
                            Err(_) => continue,
                        };
                        let norm_seq = match std::str::from_utf8(&norm.allele) {
                            Ok(s) => s,
                            Err(_) => continue,
                        };
                        let vrs_id = writer.allele_identifier_literal(
                            &view.accession,
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
                let line = line_buf.trim_end_matches('\n').trim_end_matches('\r');
                if line.starts_with('#') || line.is_empty() {
                    continue;
                }
                let fields: Vec<&str> = line.splitn(10, '\t').collect();
                if fields.len() < 5 {
                    continue;
                }
                let chrom = fields[0];
                let seq_index = match chrom_to_seq_index.get(chrom) {
                    Some(&i) => i,
                    None => continue, // unknown chrom: skip (matches serial)
                };
                let pos: u64 = fields[1]
                    .parse::<u64>()
                    .context("Invalid POS field")?
                    .saturating_sub(1);
                let ref_allele = fields[3];
                let alt_field = fields[4];

                let alts: Vec<String> = alt_field
                    .split(',')
                    .filter(|alt| !(alt.starts_with('<') || *alt == "*" || *alt == "."))
                    .map(|s| s.to_string())
                    .collect();
                if alts.is_empty() {
                    continue;
                }

                let item = WorkItem {
                    record_index,
                    seq_index,
                    pos,
                    chrom: chrom.to_string(),
                    ref_allele: ref_allele.to_string(),
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
        Ok(())
    })?;

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
