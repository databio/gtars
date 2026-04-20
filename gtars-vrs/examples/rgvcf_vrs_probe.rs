//! Plan 9 Phase A throwaway probe: rgvcf -> VRS-ID throughput.
//!
//! Mmaps an rgvcf file, parses header/chrom-table/sidecar, iterates records,
//! and pipes each (chrom, pos, ref, alt) through `normalize` +
//! `DigestWriter::allele_identifier_literal` — the same hot path the
//! `rust_native` contestant uses. Emits the same `variants.tsv` and
//! single-row results TSV as `gnomad_vrs_store`, so numbers can be compared
//! directly.
//!
//! Usage:
//!   cargo run --release --example rgvcf_vrs_probe -- \
//!       <store_dir> <rgvcf_path> <results_tsv> <variants_tsv> [num_workers]
//!
//! num_workers:
//!   0  -> single-threaded (sequential iteration; per-variant TSV written
//!         in stream order)
//!   >0 -> crossbeam-scoped pool. Per-chromosome streams are split into
//!         num_workers chunks via a pre-scan. Per-worker TSV files are
//!         concatenated at the end.

#[cfg(target_os = "linux")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;

use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore};
use gtars_vrs::digest::DigestWriter;
use gtars_vrs::normalize::normalize;

const MAGIC_HEADER: &[u8; 4] = b"RGVC";
const MAGIC_FOOTER: &[u8; 4] = b"CFGR";
const DIGEST_LEN: usize = 32;

struct ChromEntry {
    name: String,
    seq_digest: String,
    #[allow(dead_code)]
    variant_count: u32,
    stream_offset: u64,
    stream_len: u64,
}

struct Rgvcf<'a> {
    mmap: &'a [u8],
    collection_digest: String,
    chroms: Vec<ChromEntry>,
    sidecar_entry_offsets: Vec<u64>,
}

fn read_u8(buf: &[u8], off: &mut usize) -> Result<u8> {
    if *off >= buf.len() {
        bail!("unexpected EOF reading u8");
    }
    let v = buf[*off];
    *off += 1;
    Ok(v)
}
fn read_u16_le(buf: &[u8], off: &mut usize) -> Result<u16> {
    if *off + 2 > buf.len() { bail!("EOF u16"); }
    let v = u16::from_le_bytes(buf[*off..*off + 2].try_into().unwrap());
    *off += 2; Ok(v)
}
fn read_u32_le(buf: &[u8], off: &mut usize) -> Result<u32> {
    if *off + 4 > buf.len() { bail!("EOF u32"); }
    let v = u32::from_le_bytes(buf[*off..*off + 4].try_into().unwrap());
    *off += 4; Ok(v)
}
fn read_u64_le(buf: &[u8], off: &mut usize) -> Result<u64> {
    if *off + 8 > buf.len() { bail!("EOF u64"); }
    let v = u64::from_le_bytes(buf[*off..*off + 8].try_into().unwrap());
    *off += 8; Ok(v)
}

#[inline]
fn read_leb128_u64(buf: &[u8], off: &mut usize) -> Result<u64> {
    let mut result: u64 = 0;
    let mut shift: u32 = 0;
    loop {
        if *off >= buf.len() { bail!("EOF LEB128"); }
        let b = buf[*off];
        *off += 1;
        result |= ((b & 0x7f) as u64) << shift;
        if (b & 0x80) == 0 { return Ok(result); }
        shift += 7;
        if shift >= 64 { bail!("LEB128 overflow"); }
    }
}

fn parse_rgvcf(mmap: &[u8]) -> Result<Rgvcf<'_>> {
    if mmap.len() < 12 || &mmap[0..4] != MAGIC_HEADER {
        bail!("bad header magic");
    }
    if mmap[4] != 1 { bail!("bad version byte: {}", mmap[4]); }
    let mut off = 8usize;
    let digest_len = read_u8(mmap, &mut off)? as usize;
    if digest_len != DIGEST_LEN { bail!("digest len != 32"); }
    let collection_digest =
        String::from_utf8(mmap[off..off + DIGEST_LEN].to_vec()).context("coll digest utf-8")?;
    off += DIGEST_LEN;

    let num_chroms = read_u16_le(mmap, &mut off)? as usize;
    let mut chroms = Vec::with_capacity(num_chroms);
    for _ in 0..num_chroms {
        let name_len = read_u8(mmap, &mut off)? as usize;
        let name = String::from_utf8(mmap[off..off + name_len].to_vec())?;
        off += name_len;
        let seq_digest =
            String::from_utf8(mmap[off..off + DIGEST_LEN].to_vec())?;
        off += DIGEST_LEN;
        let variant_count = read_u32_le(mmap, &mut off)?;
        let stream_offset = read_u64_le(mmap, &mut off)?;
        let stream_len = read_u64_le(mmap, &mut off)?;
        chroms.push(ChromEntry { name, seq_digest, variant_count, stream_offset, stream_len });
    }

    let fsz = mmap.len();
    if fsz < 12 || &mmap[fsz - 4..fsz] != MAGIC_FOOTER {
        bail!("bad footer magic");
    }
    let sidecar_byte_offset =
        u64::from_le_bytes(mmap[fsz - 12..fsz - 4].try_into().unwrap()) as usize;

    let mut scur = sidecar_byte_offset;
    let num_sidecar = read_u32_le(mmap, &mut scur)? as usize;
    let mut sidecar_entry_offsets = Vec::with_capacity(num_sidecar);
    for _ in 0..num_sidecar {
        sidecar_entry_offsets.push(read_u64_le(mmap, &mut scur)?);
    }

    Ok(Rgvcf { mmap, collection_digest, chroms, sidecar_entry_offsets })
}

struct SidecarEntry<'a> {
    pos: u64,
    ref_: &'a [u8],
    alts: Vec<&'a [u8]>,
}

fn parse_sidecar_entry<'a>(mmap: &'a [u8], off: usize) -> Result<SidecarEntry<'a>> {
    let mut cur = off;
    let _chrom_index = read_u16_le(mmap, &mut cur)?;
    let pos = read_u64_le(mmap, &mut cur)?;
    let ref_len = read_u16_le(mmap, &mut cur)? as usize;
    let ref_ = &mmap[cur..cur + ref_len];
    cur += ref_len;
    let num_alts = read_u8(mmap, &mut cur)? as usize;
    let mut alts = Vec::with_capacity(num_alts);
    for _ in 0..num_alts {
        let alt_len = read_u16_le(mmap, &mut cur)? as usize;
        alts.push(&mmap[cur..cur + alt_len]);
        cur += alt_len;
    }
    Ok(SidecarEntry { pos, ref_, alts })
}

/// Scan a chrom stream once, emitting (record_index, byte_offset, abs_pos)
/// for each record. Used to carve chunks for parallel workers.
fn scan_chrom_breakpoints(stream: &[u8], num_chunks: usize) -> Result<Vec<(usize, u64)>> {
    // First pass: count records + store per-record (offset, pos) for the
    // requested split points.
    // Two-pass: pass 1 counts records, pass 2 picks the records at
    // i * n_records / num_chunks.
    let mut n = 0usize;
    let mut cursor = 0usize;
    while cursor < stream.len() {
        // Skip LEB128.
        while cursor < stream.len() && stream[cursor] & 0x80 != 0 { cursor += 1; }
        cursor += 1; // last LEB128 byte
        if cursor > stream.len() { bail!("bad stream"); }
        let code = stream[cursor] & 0x07;
        cursor += 1;
        if code == 4 { cursor += 4; }
        n += 1;
    }
    if n == 0 {
        return Ok(vec![(0, 0)]);
    }
    let mut targets: Vec<usize> = (0..num_chunks).map(|i| i * n / num_chunks).collect();
    targets.dedup();

    let mut bps: Vec<(usize, u64)> = Vec::with_capacity(targets.len());
    let mut ti = 0usize;
    let mut cursor = 0usize;
    let mut prev_pos: u64 = 0;
    let mut rec_idx = 0usize;
    while cursor < stream.len() && ti < targets.len() {
        if rec_idx == targets[ti] {
            bps.push((cursor, prev_pos));
            ti += 1;
        }
        // Read LEB128
        let mut val: u64 = 0;
        let mut shift = 0;
        loop {
            let b = stream[cursor];
            cursor += 1;
            val |= ((b & 0x7f) as u64) << shift;
            if b & 0x80 == 0 { break; }
            shift += 7;
        }
        let code = stream[cursor] & 0x07;
        cursor += 1;
        if code == 4 { cursor += 4; }
        prev_pos += val;
        rec_idx += 1;
    }
    Ok(bps)
}

/// Process a chunk of a chromosome stream, writing per-variant TSV rows
/// into `vw`. Returns number of variants emitted.
#[allow(clippy::too_many_arguments)]
fn process_chunk<W: Write>(
    stream: &[u8],
    start_off: usize,
    end_off: usize,       // exclusive; == stream.len() for last chunk
    start_pos: u64,       // absolute pos BEFORE the first record in this chunk
    sequence: &[u8],
    seq_accession: &str,
    sidecar_entry_offsets: &[u64],
    mmap: &[u8],
    chrom_name: &str,
    vw: &mut W,
) -> Result<usize> {
    let mut n = 0usize;
    let mut cursor = start_off;
    let mut prev_pos = start_pos;
    let mut digest_writer = DigestWriter::new();
    while cursor < end_off {
        let pos_delta = read_leb128_u64(stream, &mut cursor)?;
        let pos1 = prev_pos + pos_delta;
        prev_pos = pos1;
        let code = read_u8(stream, &mut cursor)?;
        if code & 0xF8 != 0 { bail!("invalid code upper bits"); }
        let code_low = code & 0x07;
        if code_low < 4 {
            let start0 = pos1 - 1;
            let idx = start0 as usize;
            if idx >= sequence.len() { bail!("pos {} OOB for {}", pos1, chrom_name); }
            let ref_b = sequence[idx];
            let alt_b = b"ACGT"[code_low as usize];
            let ref_buf = [ref_b];
            let alt_buf = [alt_b];
            let norm = match normalize(sequence, start0, &ref_buf, &alt_buf) {
                Ok(n) => n, Err(_) => continue,
            };
            let norm_seq = std::str::from_utf8(&norm.allele)?;
            let vrs_id = digest_writer.allele_identifier_literal(
                seq_accession, norm.start, norm.end, norm_seq,
            );
            writeln!(vw, "{}\t{}\t{}\t{}\t{}", chrom_name, pos1, ref_b as char, alt_b as char, vrs_id)?;
            n += 1;
        } else if code_low == 4 {
            let sidecar_index = read_u32_le(stream, &mut cursor)? as usize;
            let entry_off = sidecar_entry_offsets[sidecar_index] as usize;
            let entry = parse_sidecar_entry(mmap, entry_off)?;
            let start0 = entry.pos - 1;
            for alt in &entry.alts {
                if alt.is_empty() || alt[0] == b'<' || alt == b"*" || alt == b"." { continue; }
                let norm = match normalize(sequence, start0, entry.ref_, alt) {
                    Ok(n) => n, Err(_) => continue,
                };
                let norm_seq = std::str::from_utf8(&norm.allele)?;
                let vrs_id = digest_writer.allele_identifier_literal(
                    seq_accession, norm.start, norm.end, norm_seq,
                );
                writeln!(vw, "{}\t{}\t{}\t{}\t{}",
                    chrom_name, entry.pos,
                    std::str::from_utf8(entry.ref_)?,
                    std::str::from_utf8(alt)?,
                    vrs_id)?;
                n += 1;
            }
        } else {
            bail!("reserved code {}", code_low);
        }
    }
    Ok(n)
}

fn max_rss_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status").ok().and_then(|s| {
        s.lines().find(|l| l.starts_with("VmHWM:"))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|n| n.parse::<u64>().ok())
            .map(|kb| (kb as f64 / 1024.0 * 10.0).round() / 10.0)
    }).unwrap_or(0.0)
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if !(5..=6).contains(&args.len()) {
        eprintln!("Usage: {} <store_dir> <rgvcf_path> <results_tsv> <variants_tsv> [num_workers]", args[0]);
        std::process::exit(2);
    }
    let store_dir = &args[1];
    let rgvcf_path = &args[2];
    let results_tsv = &args[3];
    let variants_tsv = &args[4];
    let num_workers_arg: usize = if args.len() == 6 { args[5].parse().unwrap_or(0) } else { 0 };
    let effective_workers = if num_workers_arg == 0 {
        std::thread::available_parallelism().map(|n| n.get().saturating_sub(2).max(1)).unwrap_or(1)
    } else { num_workers_arg };

    let t0 = Instant::now();
    let mut store = RefgetStore::open_local(store_dir).context("open store")?;
    let t_open = t0.elapsed();

    let paged = store.list_collections(0, 10, &[]).context("list collections")?;
    let collections = &paged.results;
    if collections.is_empty() { bail!("No collections"); }
    let collection_digest = collections[0].digest.clone();
    store.load_collection(collection_digest.as_str()).context("load collection")?;

    let t1 = Instant::now();
    let file = File::open(rgvcf_path).with_context(|| format!("open {}", rgvcf_path))?;
    let mmap = unsafe { Mmap::map(&file)? };
    let rgvcf = parse_rgvcf(&mmap[..])?;

    let store_coll_core = collection_digest.trim_start_matches("SQ.");
    if rgvcf.collection_digest != store_coll_core {
        bail!("collection digest mismatch: rgvcf={} store={}", rgvcf.collection_digest, store_coll_core);
    }

    let zero_digest = "\0".repeat(DIGEST_LEN);
    for c in &rgvcf.chroms {
        if c.seq_digest == zero_digest {
            eprintln!("[rgvcf_probe] WARN chrom {} zero digest; skip", c.name);
            continue;
        }
        store.load_sequence(c.seq_digest.as_str())?;
        store.ensure_decoded(c.seq_digest.as_str())?;
    }
    let t_preload = t1.elapsed();
    let readonly: ReadonlyRefgetStore = store.into_readonly();

    // Collect (chrom, sequence bytes, stream slice, seq_accession) for workers.
    // Decide chunking.
    let t2 = Instant::now();

    // Prepare chrom-level work units.
    struct ChromWork<'a> {
        name: &'a str,
        stream: &'a [u8],
        sequence: &'a [u8],
        seq_accession: String,
    }
    let mut works: Vec<ChromWork> = Vec::new();
    for c in &rgvcf.chroms {
        if c.seq_digest == zero_digest { continue; }
        let raw_digest = c.seq_digest.as_str();
        let sequence = match readonly.sequence_bytes(raw_digest) {
            Some(b) => b,
            None => bail!("chrom {} not decoded", c.name),
        };
        let s_off = c.stream_offset as usize;
        let s_end = s_off + c.stream_len as usize;
        works.push(ChromWork {
            name: c.name.as_str(),
            stream: &rgvcf.mmap[s_off..s_end],
            sequence,
            seq_accession: format!("SQ.{}", raw_digest),
        });
    }

    // Build chunks: per chrom, decide how many sub-chunks based on stream
    // size and total workers.
    // Simple heuristic: if only one chrom with a big stream, split into
    // effective_workers chunks. Otherwise, chunks per chrom proportional
    // to stream size.
    let total_bytes: usize = works.iter().map(|w| w.stream.len()).sum();
    #[allow(clippy::type_complexity)]
    let mut chunks: Vec<(usize /* work idx */, usize /* start_off */, usize /* end_off */, u64 /* start_pos */)> = Vec::new();
    for (wi, w) in works.iter().enumerate() {
        let chunks_for_this = if total_bytes == 0 {
            1
        } else {
            let share = ((w.stream.len() as f64) / (total_bytes as f64) * (effective_workers as f64)).ceil() as usize;
            share.max(1).min(effective_workers)
        };
        if chunks_for_this == 1 || num_workers_arg == 0 {
            chunks.push((wi, 0, w.stream.len(), 0));
        } else {
            let bps = scan_chrom_breakpoints(w.stream, chunks_for_this)?;
            // bps are (byte_offset, abs_pos_BEFORE_this_record)
            for i in 0..bps.len() {
                let (start_off, start_pos) = bps[i];
                let end_off = if i + 1 < bps.len() { bps[i + 1].0 } else { w.stream.len() };
                if end_off > start_off {
                    chunks.push((wi, start_off, end_off, start_pos));
                }
            }
        }
    }

    // Each chunk writes its own TSV to a tmp path; we concat at the end.
    let tmpdir = std::env::temp_dir();
    let run_id = std::process::id();
    let chunk_paths: Vec<PathBuf> = (0..chunks.len())
        .map(|i| tmpdir.join(format!("rgvcf_probe_{}_{}.tsv", run_id, i)))
        .collect();

    // Parallel processing using scoped threads.
    let n_total = std::sync::atomic::AtomicUsize::new(0);
    let idx = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| -> Result<()> {
        let sidecar_offsets = &rgvcf.sidecar_entry_offsets;
        let mmap_slice = &rgvcf.mmap[..];
        let works_ref = &works;
        let chunks_ref = &chunks;
        let chunk_paths_ref = &chunk_paths;
        let n_total_ref = &n_total;

        // Simple work-stealing via an atomic index.
        let idx_ref = &idx;
        let mut handles = Vec::new();
        let nthreads = if num_workers_arg == 0 { 1 } else { effective_workers };
        for _ in 0..nthreads {
            let handle = scope.spawn(move || -> Result<()> {
                loop {
                    let i = idx_ref.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if i >= chunks_ref.len() { break; }
                    let (wi, start_off, end_off, start_pos) = chunks_ref[i];
                    let w = &works_ref[wi];
                    let f = File::create(&chunk_paths_ref[i])?;
                    let mut bw = BufWriter::with_capacity(1 << 20, f);
                    let n = process_chunk(
                        w.stream, start_off, end_off, start_pos,
                        w.sequence, &w.seq_accession,
                        sidecar_offsets, mmap_slice,
                        w.name, &mut bw,
                    )?;
                    bw.flush()?;
                    n_total_ref.fetch_add(n, std::sync::atomic::Ordering::Relaxed);
                }
                Ok(())
            });
            handles.push(handle);
        }
        for h in handles {
            h.join().map_err(|_| anyhow!("worker panicked"))??;
        }
        Ok(())
    })?;
    let t_compute = t2.elapsed();
    let n = n_total.load(std::sync::atomic::Ordering::Relaxed);

    // Concat chunk TSVs into variants_tsv.
    let vf = File::create(variants_tsv)?;
    let mut vw = BufWriter::with_capacity(1 << 20, vf);
    writeln!(vw, "chrom\tpos\tref\talt\tvrs_id")?;
    for p in &chunk_paths {
        let mut src = File::open(p)?;
        std::io::copy(&mut src, &mut vw)?;
        let _ = std::fs::remove_file(p);
    }
    vw.flush()?;

    let open_s = t_open.as_secs_f64();
    let preload_s = t_preload.as_secs_f64();
    let compute_s = t_compute.as_secs_f64();
    let total_s = open_s + preload_s + compute_s;
    let vps = if compute_s > 0.0 { n as f64 / compute_s } else { 0.0 };
    let rss = max_rss_mb();

    let vcf_basename = Path::new(rgvcf_path)
        .file_name().map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| rgvcf_path.to_string());
    let nthreads_used = if num_workers_arg == 0 { 1 } else { effective_workers };

    let rf = File::create(results_tsv)?;
    let mut rw = BufWriter::new(rf);
    writeln!(rw, "contestant\tvcf\tn_variants\tnum_workers\topen_s\tpreload_s\tcompute_s\ttotal_s\tvariants_per_sec\tmax_rss_mb")?;
    writeln!(rw, "rust_rgvcf_probe\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.0}\t{:.1}",
        vcf_basename, n, nthreads_used, open_s, preload_s, compute_s, total_s, vps, rss)?;
    rw.flush()?;

    eprintln!(
        "[rgvcf_probe] n={} workers={} open={:.2}s preload={:.2}s compute={:.2}s total={:.2}s  {:.0} v/s  rss={:.1} MB",
        n, nthreads_used, open_s, preload_s, compute_s, total_s, vps, rss
    );
    Ok(())
}
