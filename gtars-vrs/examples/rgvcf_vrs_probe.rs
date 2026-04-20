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
//! num_workers is currently ignored; this probe is single-threaded. Phase A.2
//! goal is "how fast is rgvcf single-threaded compared to bgzf single-thread
//! (rust_serial)?" If single-thread beats bgzf serial by a large margin, we
//! promote; Phase A.3 adds per-chromosome parallelism.

#[cfg(target_os = "linux")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;

use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore};
use gtars_vrs::digest::DigestWriter;
use gtars_vrs::normalize::normalize;
use gtars_vrs::vcf::build_name_to_digest_readonly;

const MAGIC_HEADER: &[u8; 4] = b"RGVC";
const MAGIC_FOOTER: &[u8; 4] = b"CFGR";
const DIGEST_LEN: usize = 32;

struct ChromEntry {
    name: String,
    seq_digest: String, // 32 ASCII base64url (no SQ. prefix)
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
    if *off + 2 > buf.len() {
        bail!("unexpected EOF reading u16");
    }
    let v = u16::from_le_bytes(buf[*off..*off + 2].try_into().unwrap());
    *off += 2;
    Ok(v)
}
fn read_u32_le(buf: &[u8], off: &mut usize) -> Result<u32> {
    if *off + 4 > buf.len() {
        bail!("unexpected EOF reading u32");
    }
    let v = u32::from_le_bytes(buf[*off..*off + 4].try_into().unwrap());
    *off += 4;
    Ok(v)
}
fn read_u64_le(buf: &[u8], off: &mut usize) -> Result<u64> {
    if *off + 8 > buf.len() {
        bail!("unexpected EOF reading u64");
    }
    let v = u64::from_le_bytes(buf[*off..*off + 8].try_into().unwrap());
    *off += 8;
    Ok(v)
}

#[inline]
fn read_leb128_u64(buf: &[u8], off: &mut usize) -> Result<u64> {
    let mut result: u64 = 0;
    let mut shift: u32 = 0;
    loop {
        if *off >= buf.len() {
            bail!("unexpected EOF reading LEB128");
        }
        let b = buf[*off];
        *off += 1;
        result |= ((b & 0x7f) as u64) << shift;
        if (b & 0x80) == 0 {
            return Ok(result);
        }
        shift += 7;
        if shift >= 64 {
            bail!("LEB128 overflow");
        }
    }
}

fn parse_rgvcf(mmap: &[u8]) -> Result<Rgvcf<'_>> {
    // Header
    if mmap.len() < 12 || &mmap[0..4] != MAGIC_HEADER {
        bail!("bad header magic");
    }
    if mmap[4] != 1 {
        bail!("bad version byte: {}", mmap[4]);
    }
    // flags byte (mmap[5]) and reserved (mmap[6..8]) ignored by probe.
    let mut off = 8usize;
    let digest_len = read_u8(mmap, &mut off)? as usize;
    if digest_len != DIGEST_LEN {
        bail!("digest len != 32");
    }
    let collection_digest =
        String::from_utf8(mmap[off..off + DIGEST_LEN].to_vec()).context("collection digest utf-8")?;
    off += DIGEST_LEN;

    // Chrom table
    let num_chroms = read_u16_le(mmap, &mut off)? as usize;
    let mut chroms = Vec::with_capacity(num_chroms);
    for _ in 0..num_chroms {
        let name_len = read_u8(mmap, &mut off)? as usize;
        let name = String::from_utf8(mmap[off..off + name_len].to_vec())
            .context("chrom name utf-8")?;
        off += name_len;
        let seq_digest = String::from_utf8(mmap[off..off + DIGEST_LEN].to_vec())
            .context("seq digest utf-8")?;
        off += DIGEST_LEN;
        let variant_count = read_u32_le(mmap, &mut off)?;
        let stream_offset = read_u64_le(mmap, &mut off)?;
        let stream_len = read_u64_le(mmap, &mut off)?;
        chroms.push(ChromEntry {
            name,
            seq_digest,
            variant_count,
            stream_offset,
            stream_len,
        });
    }

    // Footer: last 12 bytes.
    let fsz = mmap.len();
    if fsz < 12 || &mmap[fsz - 4..fsz] != MAGIC_FOOTER {
        bail!("bad footer magic");
    }
    let sidecar_byte_offset =
        u64::from_le_bytes(mmap[fsz - 12..fsz - 4].try_into().unwrap()) as usize;

    // Sidecar header: u32 num_entries + num_entries * u64 offsets.
    let mut scur = sidecar_byte_offset;
    let num_sidecar = read_u32_le(mmap, &mut scur)? as usize;
    let mut sidecar_entry_offsets = Vec::with_capacity(num_sidecar);
    for _ in 0..num_sidecar {
        sidecar_entry_offsets.push(read_u64_le(mmap, &mut scur)?);
    }

    Ok(Rgvcf {
        mmap,
        collection_digest,
        chroms,
        sidecar_entry_offsets,
    })
}

struct SidecarEntry<'a> {
    chrom_index: u16,
    pos: u64,
    ref_: &'a [u8],
    alts: Vec<&'a [u8]>,
}

fn parse_sidecar_entry<'a>(mmap: &'a [u8], off: usize) -> Result<SidecarEntry<'a>> {
    let mut cur = off;
    let chrom_index = read_u16_le(mmap, &mut cur)?;
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
    Ok(SidecarEntry {
        chrom_index,
        pos,
        ref_,
        alts,
    })
}

fn max_rss_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|n| n.parse::<u64>().ok())
                .map(|kb| (kb as f64 / 1024.0 * 10.0).round() / 10.0)
        })
        .unwrap_or(0.0)
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if !(5..=6).contains(&args.len()) {
        eprintln!(
            "Usage: {} <store_dir> <rgvcf_path> <results_tsv> <variants_tsv> [num_workers]",
            args[0]
        );
        std::process::exit(2);
    }
    let store_dir = &args[1];
    let rgvcf_path = &args[2];
    let results_tsv = &args[3];
    let variants_tsv = &args[4];

    // --- Phase 1: open store ---
    let t0 = Instant::now();
    let mut store = RefgetStore::open_local(store_dir).context("open store")?;
    let t_open = t0.elapsed();

    // Pick first collection (mirror rust_native).
    let paged = store.list_collections(0, 10, &[]).context("list collections")?;
    let collections = &paged.results;
    if collections.is_empty() {
        bail!("No collections in store");
    }
    let collection_digest = collections[0].digest.clone();
    store
        .load_collection(collection_digest.as_str())
        .context("load_collection")?;

    // --- Phase 2: mmap rgvcf + parse header + decode needed chroms ---
    let t1 = Instant::now();
    let file = File::open(rgvcf_path).with_context(|| format!("open {}", rgvcf_path))?;
    let mmap = unsafe { Mmap::map(&file).context("mmap rgvcf")? };
    let rgvcf = parse_rgvcf(&mmap[..])?;

    // Verify collection digest matches.
    // rgvcf stores the raw 32-byte base64url digest (no "SQ." prefix).
    let store_coll_core = collection_digest.trim_start_matches("SQ.");
    if rgvcf.collection_digest != store_coll_core {
        bail!(
            "rgvcf collection digest {} != store collection {}",
            rgvcf.collection_digest,
            store_coll_core
        );
    }

    // Decode every chrom in the rgvcf's chrom table. Store APIs use the
    // raw (no "SQ.") base64url digest.
    for c in &rgvcf.chroms {
        store
            .load_sequence(c.seq_digest.as_str())
            .with_context(|| format!("load {}", c.name))?;
        store
            .ensure_decoded(c.seq_digest.as_str())
            .with_context(|| format!("decode {}", c.name))?;
    }
    let t_preload = t1.elapsed();

    // --- Phase 3: convert to readonly, build helpers, iterate records ---
    let readonly = store.into_readonly();
    let name_to_digest = build_name_to_digest_readonly(&readonly, &collection_digest)?;
    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(n, d)| (n.clone(), format!("SQ.{}", d)))
        .collect();

    let vf = File::create(variants_tsv).context("create variants tsv")?;
    let mut vw = BufWriter::with_capacity(1 << 20, vf);
    writeln!(vw, "chrom\tpos\tref\talt\tvrs_id")?;

    let mut digest_writer = DigestWriter::new();
    let mut n: usize = 0;

    let t2 = Instant::now();
    for chrom in &rgvcf.chroms {
        // Grab sequence bytes for this chrom from the readonly store.
        let raw_digest = match name_to_digest.get(&chrom.name) {
            Some(d) => d.as_str(),
            None => {
                eprintln!("[rgvcf_probe] WARN: chrom {} not in name_to_digest; skipping", chrom.name);
                continue;
            }
        };
        let sequence = match readonly.sequence_bytes(raw_digest) {
            Some(b) => b,
            None => {
                bail!("chrom {} not decoded after preload", chrom.name);
            }
        };
        let seq_accession = chrom_accessions
            .get(&chrom.name)
            .ok_or_else(|| anyhow!("no accession for chrom {}", chrom.name))?
            .as_str();

        let s_off = chrom.stream_offset as usize;
        let s_end = s_off + chrom.stream_len as usize;
        let stream = &rgvcf.mmap[s_off..s_end];

        let mut cursor = 0usize;
        let mut prev_pos: u64 = 0;
        while cursor < stream.len() {
            let pos_delta = read_leb128_u64(stream, &mut cursor)?;
            let pos1 = prev_pos + pos_delta;
            prev_pos = pos1;
            let code = read_u8(stream, &mut cursor)?;
            if code & 0xF8 != 0 {
                bail!(
                    "invalid code upper bits at {} offset {}: {}",
                    chrom.name,
                    cursor - 1,
                    code
                );
            }
            let code_low = code & 0x07;
            if code_low < 4 {
                // Inline SNV. REF = sequence[pos1-1], ALT = ACGT[code_low].
                let start0: u64 = pos1 - 1; // 0-based VCF-style pos
                let idx = start0 as usize;
                if idx >= sequence.len() {
                    bail!(
                        "pos {} out of bounds for {} (len {})",
                        pos1,
                        chrom.name,
                        sequence.len()
                    );
                }
                let ref_b = sequence[idx];
                let alt_b = b"ACGT"[code_low as usize];
                let ref_buf = [ref_b];
                let alt_buf = [alt_b];

                let norm = match normalize(sequence, start0, &ref_buf, &alt_buf) {
                    Ok(n) => n,
                    Err(e) => {
                        eprintln!("[rgvcf_probe] normalize failed at {}:{}: {}", chrom.name, pos1, e);
                        continue;
                    }
                };
                let norm_seq = std::str::from_utf8(&norm.allele)?;
                let vrs_id = digest_writer.allele_identifier_literal(
                    seq_accession,
                    norm.start,
                    norm.end,
                    norm_seq,
                );
                writeln!(
                    vw,
                    "{}\t{}\t{}\t{}\t{}",
                    chrom.name,
                    pos1,
                    ref_b as char,
                    alt_b as char,
                    vrs_id
                )?;
                n += 1;
            } else if code_low == 4 {
                // Escape: next 4 bytes sidecar_index u32 LE.
                let sidecar_index = read_u32_le(stream, &mut cursor)? as usize;
                let entry_off = rgvcf.sidecar_entry_offsets[sidecar_index] as usize;
                let entry = parse_sidecar_entry(&rgvcf.mmap[..], entry_off)?;
                let _ = entry.chrom_index;
                let start0 = entry.pos - 1;
                for alt in &entry.alts {
                    // Skip symbolic / spanning-del per bgzf path semantics.
                    if alt.is_empty() || alt[0] == b'<' || alt == b"*" || alt == b"." {
                        continue;
                    }
                    let norm = match normalize(sequence, start0, entry.ref_, alt) {
                        Ok(n) => n,
                        Err(e) => {
                            eprintln!(
                                "[rgvcf_probe] normalize failed at {}:{}: {}",
                                chrom.name, entry.pos, e
                            );
                            continue;
                        }
                    };
                    let norm_seq = std::str::from_utf8(&norm.allele)?;
                    let vrs_id = digest_writer.allele_identifier_literal(
                        seq_accession,
                        norm.start,
                        norm.end,
                        norm_seq,
                    );
                    writeln!(
                        vw,
                        "{}\t{}\t{}\t{}\t{}",
                        chrom.name,
                        entry.pos,
                        std::str::from_utf8(entry.ref_)?,
                        std::str::from_utf8(alt)?,
                        vrs_id
                    )?;
                    n += 1;
                }
            } else {
                bail!("reserved code {} at {} offset {}", code_low, chrom.name, cursor);
            }
        }
    }
    let t_compute = t2.elapsed();
    vw.flush()?;

    let open_s = t_open.as_secs_f64();
    let preload_s = t_preload.as_secs_f64();
    let compute_s = t_compute.as_secs_f64();
    let total_s = open_s + preload_s + compute_s;
    let vps = if compute_s > 0.0 { n as f64 / compute_s } else { 0.0 };
    let rss = max_rss_mb();

    let vcf_basename = Path::new(rgvcf_path)
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| rgvcf_path.to_string());

    let rf = File::create(results_tsv).context("create results tsv")?;
    let mut rw = BufWriter::new(rf);
    writeln!(
        rw,
        "contestant\tvcf\tn_variants\tnum_workers\topen_s\tpreload_s\tcompute_s\ttotal_s\tvariants_per_sec\tmax_rss_mb"
    )?;
    writeln!(
        rw,
        "rust_rgvcf_probe\t{}\t{}\t1\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.0}\t{:.1}",
        vcf_basename, n, open_s, preload_s, compute_s, total_s, vps, rss
    )?;
    rw.flush()?;

    eprintln!(
        "[rgvcf_probe] n={} open={:.2}s preload={:.2}s compute={:.2}s total={:.2}s  {:.0} v/s  rss={:.1} MB",
        n, open_s, preload_s, compute_s, total_s, vps, rss
    );
    Ok(())
}
