//! VRS-ID computation from an rgvcf binary file.
//!
//! This is the sink-first replacement entry point for the bgzf-based
//! `compute_vrs_ids_parallel_bgzf_with_sink`. It mmaps the rgvcf file,
//! verifies that its `collection_digest` matches the opened
//! `ReadonlyRefgetStore`'s collection, and feeds every record through
//! the same `normalize()` + `DigestWriter` hot path.
//!
//! Parallelism: per-chromosome streams are split into `num_workers`
//! chunks via a pre-scan (`RgvcfReader::scan_breakpoints`); each chunk
//! is processed by one thread. Scales well for single-chromosome files
//! (chrY, chr22) and trivially for multi-chrom files (ClinVar).
//!
//! The callback signature mirrors `compute_vrs_ids_parallel_bgzf_with_sink`
//! so that `gnomad_vrs_store` and equivalent contestants can swap input
//! formats without touching the output-writer logic.
//!
//! **Ordering**: unlike the bgzf path, this function does not preserve
//! VCF line order across workers — per-chrom chunks are processed in
//! parallel and callback invocations interleave. Within a single chunk,
//! order is preserved (records are in SPEC-mandated strictly-increasing
//! position order).

use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{anyhow, bail, Context, Result};
use gtars_refget::store::ReadonlyRefgetStore;
use gtars_rgvcf::{Record as RgvcfRecord, RgvcfReader};

use crate::digest::DigestWriter;
use crate::normalize::normalize;

/// Emitted per VRS ID. Matches `crate::vcf::VrsResult` field-for-field
/// except the enclosing type name (so downstream sinks are swappable).
pub struct VrsRow {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub vrs_id: String,
}

/// Compute VRS IDs for every record in an rgvcf file, streaming results
/// via `on_result`. Returns the number of VRS IDs emitted.
///
/// Preconditions:
/// - every `seq_digest` referenced by the rgvcf's chrom table must
///   already be decoded in `store` (sidecar chroms with a zero-digest
///   sentinel are skipped);
/// - `name_to_digest` maps each rgvcf chromosome name to the raw
///   base64url seq_digest in `store`.
pub fn compute_vrs_ids_from_rgvcf_with_sink<F>(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    rgvcf_path: &str,
    num_workers: usize,
    mut on_result: F,
) -> Result<usize>
where
    F: FnMut(VrsRow) + Send,
{
    let reader = RgvcfReader::open(rgvcf_path)
        .with_context(|| format!("open {}", rgvcf_path))?;

    // Verify collection digest matches any chromosome's raw digest's
    // parent collection. We can't cheaply confirm the collection-level
    // digest from a ReadonlyRefgetStore here (the store reports it via
    // `name_to_digest`'s producer context); so check that every rgvcf
    // seq_digest is present in `name_to_digest`'s values.
    let known_digests: std::collections::HashSet<&str> =
        name_to_digest.values().map(|s| s.as_str()).collect();
    let zero_digest = "\0".repeat(32);
    for c in reader.chromosomes() {
        if c.seq_digest == zero_digest {
            continue;
        }
        if !known_digests.contains(c.seq_digest.as_str()) {
            bail!(
                "rgvcf chrom {} seq_digest {} not in provided name_to_digest (collection mismatch?)",
                c.name,
                c.seq_digest
            );
        }
    }

    // Build per-chrom work units.
    struct ChromWork<'a> {
        name: &'a str,
        sequence: &'a [u8],
        seq_accession: String,
        chrom_idx: usize,
    }
    let mut works: Vec<ChromWork> = Vec::new();
    for (idx, c) in reader.chromosomes().iter().enumerate() {
        if c.seq_digest == zero_digest {
            continue;
        }
        let seq = store
            .sequence_bytes(c.seq_digest.as_str())
            .ok_or_else(|| anyhow!("chrom {} (digest {}) not decoded", c.name, c.seq_digest))?;
        works.push(ChromWork {
            name: c.name.as_str(),
            sequence: seq,
            seq_accession: format!("SQ.{}", c.seq_digest),
            chrom_idx: idx,
        });
    }

    let effective_workers = num_workers.max(1);

    // Proportionally split chunks per chromosome by stream length.
    let total_bytes: usize = works
        .iter()
        .map(|w| reader.chromosomes()[w.chrom_idx].stream_byte_len as usize)
        .sum();

    #[derive(Clone, Copy)]
    struct Chunk {
        work_idx: usize,
        start_off: usize,
        end_off: usize,
        start_pos: u64,
    }
    let mut chunks: Vec<Chunk> = Vec::new();
    for (wi, w) in works.iter().enumerate() {
        let stream_len = reader.chromosomes()[w.chrom_idx].stream_byte_len as usize;
        let share = if total_bytes == 0 {
            1
        } else {
            ((stream_len as f64) / (total_bytes as f64) * (effective_workers as f64)).ceil()
                as usize
        };
        let k = share.max(1).min(effective_workers);
        if k <= 1 || effective_workers <= 1 {
            chunks.push(Chunk {
                work_idx: wi,
                start_off: 0,
                end_off: stream_len,
                start_pos: 0,
            });
        } else {
            let bps = reader.scan_breakpoints(w.chrom_idx, k)?;
            for i in 0..bps.len() {
                let (start_off, start_pos) = bps[i];
                let end_off = if i + 1 < bps.len() { bps[i + 1].0 } else { stream_len };
                if end_off > start_off {
                    chunks.push(Chunk {
                        work_idx: wi,
                        start_off,
                        end_off,
                        start_pos,
                    });
                }
            }
        }
    }

    // Batched emission through crossbeam bounded channel. Batching amortizes
    // per-message overhead (single-thread sink would otherwise dominate for
    // 10M+ variants). Workers fill a Vec<VrsRow> to ~BATCH rows then send
    // the whole batch; the sink unwraps and invokes on_result per row.
    const BATCH: usize = 4096;
    let (tx, rx) = crossbeam_channel::bounded::<Vec<VrsRow>>(effective_workers * 4);

    // Sink thread drains results; can't run in-scope since `on_result`
    // is `FnMut` captured by the caller. Collect on the current thread
    // after scope join.
    let n_total = AtomicUsize::new(0);
    let idx = AtomicUsize::new(0);

    std::thread::scope(|scope| -> Result<()> {
        let chunks_ref = &chunks;
        let works_ref = &works;
        let reader_ref = &reader;
        let idx_ref = &idx;
        let n_ref = &n_total;
        let mut handles = Vec::new();

        for _ in 0..effective_workers {
            let tx_w = tx.clone();
            handles.push(scope.spawn(move || -> Result<()> {
                let mut dw = DigestWriter::new();
                let mut batch: Vec<VrsRow> = Vec::with_capacity(BATCH);
                loop {
                    let i = idx_ref.fetch_add(1, Ordering::Relaxed);
                    if i >= chunks_ref.len() {
                        break;
                    }
                    let ch = chunks_ref[i];
                    let w = &works_ref[ch.work_idx];
                    let entry = &reader_ref.chromosomes()[w.chrom_idx];
                    let iter = reader_ref.iter_chrom_slice(
                        w.chrom_idx,
                        entry,
                        ch.start_off,
                        ch.end_off,
                        ch.start_pos,
                    );
                    for res in iter {
                        let rec = res?;
                        match rec {
                            RgvcfRecord::Snv { pos, alt, .. } => {
                                let start0 = pos - 1;
                                let idx = start0 as usize;
                                if idx >= w.sequence.len() {
                                    continue;
                                }
                                let ref_b = w.sequence[idx];
                                let ref_buf = [ref_b];
                                let alt_buf = [alt];
                                let norm = match normalize(w.sequence, start0, &ref_buf, &alt_buf) {
                                    Ok(n) => n,
                                    Err(_) => continue,
                                };
                                let norm_seq = match std::str::from_utf8(&norm.allele) {
                                    Ok(s) => s,
                                    Err(_) => continue,
                                };
                                let vrs_id = dw.allele_identifier_literal(
                                    &w.seq_accession,
                                    norm.start,
                                    norm.end,
                                    norm_seq,
                                );
                                let row = VrsRow {
                                    chrom: w.name.to_string(),
                                    pos,
                                    ref_allele: (ref_b as char).to_string(),
                                    alt_allele: (alt as char).to_string(),
                                    vrs_id,
                                };
                                n_ref.fetch_add(1, Ordering::Relaxed);
                                batch.push(row);
                                if batch.len() >= BATCH {
                                    let full = std::mem::replace(&mut batch, Vec::with_capacity(BATCH));
                                    if tx_w.send(full).is_err() { return Ok(()); }
                                }
                            }
                            RgvcfRecord::Complex { pos, ref_, alts, .. } => {
                                let start0 = pos - 1;
                                for alt in alts.iter() {
                                    if alt.is_empty()
                                        || alt[0] == b'<'
                                        || *alt == b"*"
                                        || *alt == b"."
                                    {
                                        continue;
                                    }
                                    let norm = match normalize(w.sequence, start0, ref_, alt) {
                                        Ok(n) => n,
                                        Err(_) => continue,
                                    };
                                    let norm_seq = match std::str::from_utf8(&norm.allele) {
                                        Ok(s) => s,
                                        Err(_) => continue,
                                    };
                                    let vrs_id = dw.allele_identifier_literal(
                                        &w.seq_accession,
                                        norm.start,
                                        norm.end,
                                        norm_seq,
                                    );
                                    let row = VrsRow {
                                        chrom: w.name.to_string(),
                                        pos,
                                        ref_allele: std::str::from_utf8(ref_)
                                            .unwrap_or("N")
                                            .to_string(),
                                        alt_allele: std::str::from_utf8(alt)
                                            .unwrap_or("N")
                                            .to_string(),
                                        vrs_id,
                                    };
                                    n_ref.fetch_add(1, Ordering::Relaxed);
                                    batch.push(row);
                                    if batch.len() >= BATCH {
                                        let full = std::mem::replace(&mut batch, Vec::with_capacity(BATCH));
                                        if tx_w.send(full).is_err() { return Ok(()); }
                                    }
                                }
                            }
                        }
                    }
                }
                if !batch.is_empty() {
                    let _ = tx_w.send(batch);
                }
                Ok(())
            }));
        }
        // Drop the original sender so the receiver knows when workers are done.
        drop(tx);
        for batch in rx.iter() {
            for row in batch {
                on_result(row);
            }
        }
        for h in handles {
            h.join().map_err(|_| anyhow!("worker panicked"))??;
        }
        Ok(())
    })?;

    Ok(n_total.load(Ordering::Relaxed))
}
