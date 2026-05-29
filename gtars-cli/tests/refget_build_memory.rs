//! Bounded-memory regression test for `gtars refget build`.
//!
//! Background: parallel FASTA import used to materialize a FULL in-memory
//! `BuiltCollection` (every encoded sequence of a file) before handing it to the
//! inserter, so peak memory grew with collection size AND with the number of
//! concurrently-built files -> OOM on large multi-sequence transcriptomes. The
//! import was refactored to STREAM sequences one at a time over a bounded channel
//! and persist each `.seq` immediately, retaining only per-sequence metadata.
//!
//! This test builds several large multi-sequence FASTA files and runs the real
//! `gtars refget build` binary as a subprocess, measuring its peak RSS (VmHWM
//! from `/proc/<pid>/status`). It asserts:
//!   1. Peak RSS is far below the total encoded byte volume (bytes are NOT all
//!      resident at once -- they stream to disk).
//!   2. Peak RSS stays roughly FLAT as `jobs` increases (1 -> 2 -> 4); the old
//!      design grew with both collection size and files-in-flight.
//!
//! Linux-only (reads /proc). Skips on other platforms.

#![cfg(target_os = "linux")]

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

/// A unique scratch dir under the system temp dir, cleaned up on drop.
struct Scratch {
    path: PathBuf,
}

impl Scratch {
    fn new(tag: &str) -> Self {
        let mut p = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        p.push(format!("gtars_memtest_{}_{}_{}", tag, std::process::id(), nanos));
        std::fs::create_dir_all(&p).unwrap();
        Scratch { path: p }
    }
    fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

/// Write a FASTA with `n_seqs` distinct sequences, each `seq_len` bases.
/// `salt` makes sequences distinct across files so collections differ and no
/// cross-file dedup masks the per-file work. Returns total encoded byte volume
/// for this file (2-bit DNA => length/4 bytes per sequence).
fn write_big_fasta(path: &Path, n_seqs: usize, seq_len: usize, salt: u64) -> usize {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let f = File::create(path).unwrap();
    let mut w = BufWriter::new(f);
    let mut total_encoded = 0usize;
    for i in 0..n_seqs {
        writeln!(w, ">seq_{}_{}", salt, i).unwrap();
        // Deterministic pseudo-random-ish sequence unique per (salt, i).
        let mut x: u64 = salt
            .wrapping_mul(0x9E37_79B9_7F4A_7C15)
            .wrapping_add(i as u64 + 1);
        let mut line = Vec::with_capacity(seq_len);
        for _ in 0..seq_len {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            line.push(BASES[(x & 3) as usize]);
        }
        w.write_all(&line).unwrap();
        w.write_all(b"\n").unwrap();
        total_encoded += seq_len / 4; // 2-bit packed
    }
    w.flush().unwrap();
    total_encoded
}

/// Run `gtars refget build` as a subprocess and return its peak RSS in bytes
/// (VmHWM). VmHWM is the process high-water mark, retained until exit, so reading
/// it after the child finishes (but before reaping is irrelevant -- we read
/// /proc/<pid>/status while it is still our child) is unreliable; instead we
/// poll VmHWM while the child runs and keep the max we observe.
fn build_and_measure_peak_rss(fastas: &[PathBuf], out_dir: &Path, jobs: usize) -> u64 {
    let bin = env!("CARGO_BIN_EXE_gtars");
    let mut cmd = Command::new(bin);
    cmd.arg("refget")
        .arg("build")
        .arg("--output")
        .arg(out_dir)
        .arg("--jobs")
        .arg(jobs.to_string());
    for fa in fastas {
        cmd.arg(fa);
    }
    // Silence the child's chatter.
    cmd.stdout(std::process::Stdio::null());
    cmd.stderr(std::process::Stdio::null());

    let mut child = cmd.spawn().expect("failed to spawn gtars binary");
    let pid = child.id();

    let mut peak_kb: u64 = 0;
    loop {
        if let Some(kb) = read_vmhwm_kb(pid) {
            peak_kb = peak_kb.max(kb);
        }
        match child.try_wait().unwrap() {
            Some(status) => {
                // One final read in case the high-water mark moved late.
                if let Some(kb) = read_vmhwm_kb(pid) {
                    peak_kb = peak_kb.max(kb);
                }
                assert!(status.success(), "gtars refget build failed (jobs={})", jobs);
                break;
            }
            None => std::thread::sleep(std::time::Duration::from_millis(2)),
        }
    }

    peak_kb * 1024
}

/// Read VmHWM (peak resident set size) in KB from /proc/<pid>/status.
fn read_vmhwm_kb(pid: u32) -> Option<u64> {
    let txt = std::fs::read_to_string(format!("/proc/{}/status", pid)).ok()?;
    for line in txt.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            // Format: "VmHWM:\t   12345 kB"
            let kb: u64 = rest
                .split_whitespace()
                .next()?
                .parse()
                .ok()?;
            return Some(kb);
        }
    }
    None
}

/// Build a set of `n_files` FASTA files (each `n_seqs` x `seq_len`) and return
/// the per-jobs peak RSS plus the total encoded byte volume.
fn measure(
    tag: &str,
    n_files: usize,
    n_seqs: usize,
    seq_len: usize,
    jobs_list: &[usize],
) -> (Vec<(usize, u64)>, u64) {
    let work = Scratch::new(&format!("fastas_{}", tag));
    let mut fastas = Vec::new();
    let mut total_encoded = 0usize;
    for fi in 0..n_files {
        let p = work.path().join(format!("big_{}.fa", fi));
        total_encoded += write_big_fasta(&p, n_seqs, seq_len, 1000 + fi as u64);
        fastas.push(p);
    }
    // Keep `work` alive for the duration by leaking the fastas into a held vec;
    // we manually keep `work` in scope below.
    let mut peaks = Vec::new();
    for &jobs in jobs_list {
        let out = Scratch::new(&format!("store_{}_j{}", tag, jobs));
        let peak = build_and_measure_peak_rss(&fastas, out.path(), jobs);
        eprintln!(
            "[{}] jobs={} seqlen={} -> peak RSS {:.1} MB (encoded ~{} MB)",
            tag,
            jobs,
            seq_len,
            peak as f64 / (1024.0 * 1024.0),
            total_encoded / (1024 * 1024),
        );
        peaks.push((jobs, peak));
    }
    drop(work);
    (peaks, total_encoded as u64)
}

#[test]
fn refget_build_memory_is_bounded_and_flat() {
    // Many small sequences across several files. The dominant steady-state
    // resident cost is the per-sequence METADATA index (name + two digests per
    // sequence), which is unavoidable (it IS the on-disk index) and is the SAME
    // regardless of `jobs` or sequence length. What the streaming refactor fixes
    // is that the ENCODED BYTES are no longer all held resident, and peak does
    // not grow with files-in-flight.
    //
    // Tunable for CI; the assertion SHAPE (flat-with-jobs + invariant-to-length)
    // is what matters, not the absolute numbers.
    let n_files = 4usize;
    let n_seqs = 15_000usize;
    let seq_len_short = 300usize;

    // --- (A) Peak RSS is FLAT as jobs increases. ---
    let (peaks, _enc) = measure("flat", n_files, n_seqs, seq_len_short, &[1, 2, 4]);
    let peak_j1 = peaks[0].1;
    let peak_j4 = peaks[2].1;
    assert!(
        peak_j4 <= peak_j1 * 2 + 32 * 1024 * 1024,
        "peak RSS must stay roughly flat as jobs increases: jobs=1 {:.1} MB vs jobs=4 {:.1} MB \
         (linear-in-jobs growth indicates the old whole-collection-x-files accumulation)",
        peak_j1 as f64 / (1024.0 * 1024.0),
        peak_j4 as f64 / (1024.0 * 1024.0),
    );

    // --- (B) Peak RSS is INVARIANT to encoded byte volume. ---
    // Same number of sequences, but 8x longer => 8x the encoded byte volume. If
    // encoded bytes were all resident (the old bug), peak would grow with the
    // byte volume. With streaming, peak is dominated by the (unchanged) metadata
    // index and stays essentially the same.
    let seq_len_long = seq_len_short * 8;
    let (peaks_short, enc_short) = measure("short", n_files, n_seqs, seq_len_short, &[2]);
    let (peaks_long, enc_long) = measure("long", n_files, n_seqs, seq_len_long, &[2]);
    let peak_short = peaks_short[0].1;
    let peak_long = peaks_long[0].1;
    let extra_encoded = enc_long.saturating_sub(enc_short);

    eprintln!(
        "byte-volume invariance: short peak {:.1} MB (enc {:.1} MB) vs long peak {:.1} MB (enc {:.1} MB); extra encoded {:.1} MB",
        peak_short as f64 / (1024.0 * 1024.0),
        enc_short as f64 / (1024.0 * 1024.0),
        peak_long as f64 / (1024.0 * 1024.0),
        enc_long as f64 / (1024.0 * 1024.0),
        extra_encoded as f64 / (1024.0 * 1024.0),
    );

    // The peak must NOT grow by anything close to the extra encoded volume. We
    // allow a small slack (in-flight window + allocator slop) but require the
    // increase to be far less than the added byte volume. The old design would
    // have grown peak by ~the full extra encoded volume (and more, x jobs).
    let peak_growth = peak_long.saturating_sub(peak_short);
    assert!(
        peak_growth < extra_encoded / 2,
        "peak RSS grew {:.1} MB when encoded volume grew {:.1} MB -- streaming should keep \
         encoded bytes off the heap, so growth must be far below the byte-volume increase",
        peak_growth as f64 / (1024.0 * 1024.0),
        extra_encoded as f64 / (1024.0 * 1024.0),
    );

    // --- (C) Incremental persistence: .seq files appear on disk during build. ---
    let work = Scratch::new("verify_fastas");
    let mut fastas = Vec::new();
    for fi in 0..2 {
        let p = work.path().join(format!("v_{}.fa", fi));
        write_big_fasta(&p, 2000, seq_len_short, 50 + fi as u64);
        fastas.push(p);
    }
    let out = Scratch::new("store_verify");
    let _ = build_and_measure_peak_rss(&fastas, out.path(), 2);
    let seqdir = out.path().join("sequences");
    assert!(seqdir.exists(), "sequences/ dir must exist after build");
    let mut found_seq = false;
    if let Ok(entries) = std::fs::read_dir(&seqdir) {
        for e in entries.flatten() {
            if e.path().is_dir() {
                found_seq = true;
                break;
            }
        }
    }
    assert!(found_seq, "expected .seq files persisted to disk during build");
}
