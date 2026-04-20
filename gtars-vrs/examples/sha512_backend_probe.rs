//! SHA-512 backend probe (plan 7: sha_ni investigation).
//!
//! Runs SHA-512 on a fixed 4 KiB buffer in a tight loop, reports ns/byte and
//! GB/s. Also prints runtime CPU feature detection for `avx2`, `avx512f`, and
//! `sha` so we can tell which backend `sha2 0.10.x` is dispatching to.
//!
//! This file is intended as a one-off diagnostic for plan `sha_ni_plan_v1.md`.
//! It stays in tree as a historical probe but is not wired into the CLI.

use sha2::{Digest, Sha512};
use std::time::Instant;

fn main() {
    #[cfg(target_arch = "x86_64")]
    {
        println!("cpu features:");
        println!("  avx2     = {}", is_x86_feature_detected!("avx2"));
        println!("  avx512f  = {}", is_x86_feature_detected!("avx512f"));
        println!("  sha      = {}", is_x86_feature_detected!("sha"));
    }

    let buf = vec![0x5au8; 4096];
    // Warmup
    for _ in 0..1000 {
        let _ = Sha512::digest(&buf);
    }

    // Scale iterations to run for at least ~500 ms.
    let iters: u64 = 200_000;
    let start = Instant::now();
    let mut acc: u64 = 0;
    for _ in 0..iters {
        let d = Sha512::digest(&buf);
        acc = acc.wrapping_add(d[0] as u64);
    }
    let elapsed = start.elapsed();
    std::hint::black_box(acc);

    let bytes = (iters as u128) * (buf.len() as u128);
    let ns = elapsed.as_nanos();
    let ns_per_byte = ns as f64 / bytes as f64;
    let gb_per_s = (bytes as f64) / 1e9 / elapsed.as_secs_f64();

    println!(
        "sha512: iters={} buflen={} total_bytes={} elapsed={:.3?}",
        iters,
        buf.len(),
        bytes,
        elapsed
    );
    println!("  ns/byte = {:.3}", ns_per_byte);
    println!("  GB/s    = {:.3}", gb_per_s);
}
