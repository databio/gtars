//! Benchmark: encoded-on-the-fly vs decoded-in-memory reference access, multithreaded.
//!
//! Answers: when computing VRS IDs across N in-process threads (shared store via
//! `&`, the `Arc`-equivalent), is it just as fast to keep the reference 2-bit
//! *encoded* and decode single bases on the fly as it is to hold a fully *decoded*
//! copy in RAM? If yes, the decoded cache + file-backed mmap variant + `unsafe`
//! can be removed in favor of the simpler, smaller encoded store.
//!
//! Both arms use identical `std::thread::scope` threading and per-thread
//! `DigestWriter`. The ONLY difference is the reference access path:
//!   - `encoded`: `EncodedSeq` decodes each base on demand from the 2-bit buffer.
//!   - `decoded`: a `Vec<u8>` decoded once up front, indexed directly.
//!
//! This is an in-process comparison (arm 1 vs arm 2). The multi-process mmap arm
//! (arm 3) only differs cross-process and is measured separately (PSS).
//!
//! Usage:
//!   cargo run --release --example bench_encoded -- <store_path> <vcf_path> [mode] [threads_csv]
//!     mode:        encoded | decoded | both | verify   (default: both)
//!     threads_csv: e.g. 1,2,4,8                          (default: 1,<num_cpus>)

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use flate2::read::MultiGzDecoder;
use gtars_refget::digest::alphabet::lookup_alphabet;
use gtars_refget::digest::decode_substring_from_bytes;
use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore};
use gtars_vrs::digest::DigestWriter;
use gtars_vrs::normalize::{normalize_ref, EncodedSeq};
use gtars_vrs::vcf::parse_vcf_record;

/// One reference sequence used by the VCF, with everything needed for both arms.
struct SeqEntry<'a> {
    /// `SQ.<digest>` refget accession (used in the VRS digest).
    accession: String,
    /// Encoded (bit-packed) bytes, borrowed from the store.
    encoded: &'a [u8],
    /// Number of bases.
    length: usize,
    /// Bits per symbol (e.g. 2 for DNA_2BIT).
    bits_per_symbol: usize,
    /// Decoding table (encoded code -> symbol byte).
    decoding_array: &'static [u8; 256],
}

/// A single VCF variant, with ref/alt bytes stored in a shared arena.
#[derive(Clone, Copy)]
struct Variant {
    entry: u32,
    pos: u64,
    ref_off: u32,
    ref_len: u32,
    alt_off: u32,
    alt_len: u32,
}

fn peak_rss_mb() -> f64 {
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            let kb: f64 = rest.split_whitespace().next().and_then(|s| s.parse().ok()).unwrap_or(0.0);
            return kb / 1024.0;
        }
    }
    0.0
}

fn open_vcf(path: &str) -> Box<dyn BufRead> {
    let file = File::open(path).unwrap_or_else(|e| panic!("open {path}: {e}"));
    let cap = 1 << 20;
    if path.ends_with(".gz") || path.ends_with(".bgz") {
        Box::new(BufReader::with_capacity(cap, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(cap, file))
    }
}

fn read_line(reader: &mut dyn BufRead, buf: &mut String) -> bool {
    buf.clear();
    match reader.read_line(buf) {
        Ok(0) => false,
        Ok(_) => true,
        // BGZF EOF empty block can surface as InvalidInput; treat as EOF.
        Err(e) if e.kind() == std::io::ErrorKind::InvalidInput => false,
        Err(e) => panic!("read vcf line: {e}"),
    }
}

/// Run one timed parallel pass. Returns (count, checksum). The checksum is a
/// wrapping byte-sum of every VRS id; it MUST match between encoded and decoded
/// arms (strong correctness cross-check) and prevents dead-code elimination.
fn run_pass(
    variants: &[Variant],
    arena: &[u8],
    entries: &[SeqEntry<'_>],
    decoded: &[Vec<u8>],
    encoded_mode: bool,
    threads: usize,
) -> (u64, u64) {
    let chunk = variants.len().div_ceil(threads.max(1));
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for t in 0..threads {
            let start = t * chunk;
            if start >= variants.len() {
                break;
            }
            let end = (start + chunk).min(variants.len());
            let slice = &variants[start..end];
            let h = scope.spawn(move || {
                let mut writer = DigestWriter::new();
                let mut count: u64 = 0;
                let mut checksum: u64 = 0;
                for v in slice {
                    let e = &entries[v.entry as usize];
                    let ref_allele = &arena[v.ref_off as usize..(v.ref_off + v.ref_len) as usize];
                    let alt_allele = &arena[v.alt_off as usize..(v.alt_off + v.alt_len) as usize];

                    let norm = if encoded_mode {
                        let es = EncodedSeq {
                            bytes: e.encoded,
                            length: e.length,
                            bits_per_symbol: e.bits_per_symbol,
                            decoding_array: e.decoding_array,
                        };
                        normalize_ref(&es, v.pos, ref_allele, alt_allele)
                    } else {
                        let dec = decoded[v.entry as usize].as_slice();
                        normalize_ref(dec, v.pos, ref_allele, alt_allele)
                    };
                    let norm = match norm {
                        Ok(n) => n,
                        Err(_) => continue,
                    };
                    let norm_seq = match std::str::from_utf8(&norm.allele) {
                        Ok(s) => s,
                        Err(_) => continue,
                    };
                    let vrs = writer.allele_identifier_literal(&e.accession, norm.start, norm.end, norm_seq);
                    for &b in vrs.as_bytes() {
                        checksum = checksum.wrapping_add(b as u64);
                    }
                    count += 1;
                }
                (count, checksum)
            });
            handles.push(h);
        }
        let mut count = 0u64;
        let mut checksum = 0u64;
        for h in handles {
            let (c, k) = h.join().expect("worker panicked");
            count += c;
            checksum = checksum.wrapping_add(k);
        }
        (count, checksum)
    })
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <store_path> <vcf_path> [encoded|decoded|both|verify] [threads_csv]", args[0]);
        std::process::exit(1);
    }
    let store_path = &args[1];
    let vcf_path = &args[2];
    let mode = args.get(3).map(|s| s.as_str()).unwrap_or("both");
    let max_threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8);
    let thread_counts: Vec<usize> = match args.get(4) {
        Some(csv) => csv.split(',').filter_map(|s| s.trim().parse().ok()).collect(),
        None => {
            let mut v = vec![1usize];
            if max_threads > 1 {
                v.push(max_threads);
            }
            v
        }
    };

    // ── Setup (untimed) ───────────────────────────────────────────────────
    eprintln!("Opening store: {store_path}");
    let t0 = Instant::now();
    let mut store = RefgetStore::open_local(store_path).expect("open_local");
    store.load_all_collections().expect("load_all_collections");
    let paged = store.list_collections(0, 10, &[]).expect("list_collections");
    let coll_digest = paged.results.first().expect("no collections").digest.clone();
    let collection = store.get_collection(&coll_digest).expect("get_collection");
    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    for rec in &collection.sequences {
        let m = rec.metadata();
        name_to_digest.insert(m.name.clone(), m.sha512t24u.clone());
    }
    eprintln!(
        "Store opened in {:.1}s ({} sequences in collection {})",
        t0.elapsed().as_secs_f64(),
        collection.sequences.len(),
        coll_digest
    );

    // Parse VCF -> variants + arena, interning digests as entry indices.
    eprintln!("Parsing VCF: {vcf_path}");
    let tp = Instant::now();
    let mut digest_index: HashMap<String, u32> = HashMap::new();
    let mut digests: Vec<String> = Vec::new();
    let mut variants: Vec<Variant> = Vec::new();
    let mut arena: Vec<u8> = Vec::new();
    let mut line = String::new();
    let mut reader = open_vcf(vcf_path);
    let mut skipped_unknown_chrom = 0u64;
    while read_line(&mut *reader, &mut line) {
        let Some(rec) = parse_vcf_record(&line) else {
            continue;
        };
        if rec.ref_allele.is_empty() || rec.alts.is_empty() {
            continue;
        }
        let raw_digest = match name_to_digest.get(rec.chrom) {
            Some(d) => d,
            None => {
                skipped_unknown_chrom += 1;
                continue;
            }
        };
        let pos = rec.pos;
        let entry = *digest_index.entry(raw_digest.clone()).or_insert_with(|| {
            digests.push(raw_digest.clone());
            (digests.len() - 1) as u32
        });
        let ref_off = arena.len() as u32;
        arena.extend_from_slice(rec.ref_allele.as_bytes());
        let ref_len = rec.ref_allele.len() as u32;
        for alt in rec.real_alts() {
            let alt_off = arena.len() as u32;
            arena.extend_from_slice(alt.as_bytes());
            variants.push(Variant {
                entry,
                pos,
                ref_off,
                ref_len,
                alt_off,
                alt_len: alt.len() as u32,
            });
        }
    }
    eprintln!(
        "Parsed {} variants over {} sequences in {:.1}s ({} rows skipped: unknown chrom)",
        variants.len(),
        digests.len(),
        tp.elapsed().as_secs_f64(),
        skipped_unknown_chrom
    );

    // Make encoded bytes resident for the needed sequences, then freeze to readonly.
    for d in &digests {
        store.load_sequence(d).unwrap_or_else(|e| panic!("load_sequence {d}: {e}"));
    }
    let store: ReadonlyRefgetStore = store.into_readonly();

    // Build per-sequence views (encoded borrows from store; decoded built once).
    let mut entries: Vec<SeqEntry> = Vec::with_capacity(digests.len());
    let mut decoded: Vec<Vec<u8>> = Vec::with_capacity(digests.len());
    let mut encoded_resident: u64 = 0;
    let mut decoded_resident: u64 = 0;
    for d in &digests {
        let rec = store.get_sequence(d.as_str()).unwrap_or_else(|e| panic!("get_sequence {d}: {e}"));
        let meta = rec.metadata();
        let alphabet = lookup_alphabet(&meta.alphabet);
        let enc = rec.sequence().unwrap_or_else(|| panic!("sequence {d} not resident/encoded"));
        encoded_resident += enc.len() as u64;
        decoded_resident += meta.length as u64;
        // Decode once for the decoded arm (and for verify).
        let dec = decode_substring_from_bytes(enc, 0, meta.length, alphabet);
        decoded.push(dec);
        entries.push(SeqEntry {
            accession: format!("SQ.{}", meta.sha512t24u),
            encoded: enc,
            length: meta.length,
            bits_per_symbol: alphabet.bits_per_symbol,
            decoding_array: alphabet.decoding_array,
        });
    }
    eprintln!(
        "Reference resident bytes: encoded = {:.1} MB, decoded = {:.1} MB ({:.2}x)",
        encoded_resident as f64 / 1e6,
        decoded_resident as f64 / 1e6,
        decoded_resident as f64 / encoded_resident.max(1) as f64
    );

    // ── Verify: encoded and decoded must agree (sample) ───────────────────
    if mode == "verify" || mode == "both" {
        let sample = variants.len().min(200_000);
        let (_, ck_enc) = run_pass(&variants[..sample], &arena, &entries, &decoded, true, 1);
        let (_, ck_dec) = run_pass(&variants[..sample], &arena, &entries, &decoded, false, 1);
        if ck_enc == ck_dec {
            eprintln!("VERIFY ✓ encoded == decoded VRS-id checksum over {sample} variants ({ck_enc})");
        } else {
            eprintln!("VERIFY ✗ MISMATCH encoded={ck_enc} decoded={ck_dec} over {sample} variants");
            std::process::exit(2);
        }
        if mode == "verify" {
            return;
        }
    }

    // ── Timed parallel passes ─────────────────────────────────────────────
    println!("mode\tthreads\tn_variants\tcompute_s\tvariants_per_sec\tchecksum");
    let modes: Vec<bool> = match mode {
        "encoded" => vec![true],
        "decoded" => vec![false],
        _ => vec![true, false],
    };
    for &enc_mode in &modes {
        let label = if enc_mode { "encoded" } else { "decoded" };
        for &t in &thread_counts {
            // warm pass discarded to stabilize caches, then a timed pass.
            let _ = run_pass(&variants, &arena, &entries, &decoded, enc_mode, t);
            let start = Instant::now();
            let (count, checksum) = run_pass(&variants, &arena, &entries, &decoded, enc_mode, t);
            let secs = start.elapsed().as_secs_f64();
            println!(
                "{label}\t{t}\t{count}\t{secs:.3}\t{:.0}\t{checksum}",
                count as f64 / secs
            );
        }
    }
    eprintln!("Peak RSS: {:.0} MB", peak_rss_mb());
}
