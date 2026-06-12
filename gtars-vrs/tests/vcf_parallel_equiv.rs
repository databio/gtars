//! VCF parallel-vs-serial VRS equivalence — single-reader + BGZF-block paths.
//!
//! Merged from the former `test_parallel_encoded` (single-reader parallel) and
//! `test_parallel_bgzf_encoded` (BGZF-block-parallel) binaries. Both prove their
//! parallel path produces exactly the same VRS ids in exactly the same order as
//! the serial `compute_vrs_ids_streaming_readonly` reference path.
//!
//! The shared disk-backed encoded-store fixture (`EncodedVcfFixture`), the
//! reopen helpers, the serial runner, and the BGZF writer helpers are
//! deduplicated into `tests/common/mod.rs`. The fixture writes the SAME VCF text
//! both as a plain `.vcf` (single-reader + serial reference) and as a multi-block
//! BGZF `.vcf.bgz` (BGZF path); the BGZF file's longer INFO padding only affects
//! block-straddling, not VRS output.

use std::collections::HashMap;
use std::io::Write;

use gtars_refget::store::{FastaImportOptions, ReadonlyRefgetStore, RefgetStore};
use gtars_vrs::vcf::{
    VrsResult, compute_vrs_ids_parallel_bgzf_encoded_with_sink, compute_vrs_ids_parallel_encoded,
    compute_vrs_ids_streaming_readonly,
};
use tempfile::tempdir;

mod common;
use common::{build_encoded_vcf_fixture, write_bgzf_file};

fn serial_results(
    store: &ReadonlyRefgetStore,
    n2d: &HashMap<String, String>,
    vcf: &str,
) -> Vec<VrsResult> {
    let mut out = Vec::new();
    let n = compute_vrs_ids_streaming_readonly(store, n2d, vcf, |r| out.push(r)).unwrap();
    assert_eq!(n, out.len());
    out
}

fn parallel_results(
    store: &ReadonlyRefgetStore,
    n2d: &HashMap<String, String>,
    vcf: &str,
    threads: usize,
) -> Vec<VrsResult> {
    let mut out = Vec::new();
    let n = compute_vrs_ids_parallel_encoded(store, n2d, vcf, threads, |r| out.push(r)).unwrap();
    assert_eq!(n, out.len());
    out
}

fn bgzf_results(
    store: &ReadonlyRefgetStore,
    n2d: &HashMap<String, String>,
    vcf: &str,
    workers: usize,
) -> Vec<VrsResult> {
    let mut out = Vec::new();
    let n =
        compute_vrs_ids_parallel_bgzf_encoded_with_sink(store, n2d, vcf, workers, |r| out.push(r))
            .unwrap();
    assert_eq!(n, out.len());
    out
}

// ============================================================================
// Single-reader parallel path (plain VCF)
// ============================================================================

#[test]
fn parallel_encoded_matches_serial() {
    let fx = build_encoded_vcf_fixture();
    let n2d = &fx.name_to_digest;
    let vcf = fx.vcf.as_str();

    let serial_store = fx.open_serial();
    let serial = serial_results(&serial_store, n2d, vcf);
    assert!(!serial.is_empty(), "fixture produced no results");

    let encoded_store = fx.open_encoded();
    for &threads in &[1usize, 2, 4, 8] {
        let par = parallel_results(&encoded_store, n2d, vcf, threads);
        assert_eq!(
            par.len(),
            serial.len(),
            "count mismatch at {threads} thread(s): serial={} parallel={}",
            serial.len(),
            par.len()
        );
        for (i, (s, p)) in serial.iter().zip(par.iter()).enumerate() {
            assert_eq!(
                s.vrs_id, p.vrs_id,
                "VRS id mismatch at record {i} with {threads} thread(s): serial={} parallel={}",
                s.vrs_id, p.vrs_id
            );
            assert_eq!(
                s.chrom, p.chrom,
                "chrom order mismatch at {i} ({threads} threads)"
            );
            assert_eq!(
                s.pos, p.pos,
                "pos order mismatch at {i} ({threads} threads)"
            );
            assert_eq!(
                s.ref_allele, p.ref_allele,
                "ref mismatch at {i} ({threads} threads)"
            );
            assert_eq!(
                s.alt_allele, p.alt_allele,
                "alt mismatch at {i} ({threads} threads)"
            );
        }
    }
}

#[test]
fn vrs_ids_from_raw_mode_store_match_serial() {
    // A store opened with encoding disabled (Raw, 1 byte/base) must yield the
    // same VRS ids from the parallel path as from the serial reference path.
    // This locks in the mode-aware reference-view selection in the parallel path:
    // before the fix, the worker always built a 2-bit `EncodedSeq` regardless of
    // store mode, so it would 2-bit-decode the already-decoded ASCII bytes of a
    // Raw store and emit garbage bases / wrong VRS ids, diverging from serial.
    let fx = build_encoded_vcf_fixture();
    let n2d = &fx.name_to_digest;
    let vcf = fx.vcf.as_str();

    let serial_store = fx.open_raw();
    let serial = serial_results(&serial_store, n2d, vcf);
    assert!(!serial.is_empty(), "fixture produced no results");

    let raw_store = fx.open_raw();
    for &threads in &[1usize, 2, 4, 8] {
        let par = parallel_results(&raw_store, n2d, vcf, threads);
        assert_eq!(
            par.len(),
            serial.len(),
            "count mismatch at {threads} thread(s)"
        );
        for (i, (s, p)) in serial.iter().zip(par.iter()).enumerate() {
            assert_eq!(
                s.vrs_id, p.vrs_id,
                "VRS id mismatch at record {i} with {threads} thread(s): serial={} parallel={}",
                s.vrs_id, p.vrs_id
            );
        }
    }
}

#[test]
fn parallel_encoded_is_deterministic() {
    // Same input, repeated runs at a high thread count, must be byte-identical.
    let fx = build_encoded_vcf_fixture();
    let store = fx.open_encoded();
    let a = parallel_results(&store, &fx.name_to_digest, &fx.vcf, 8);
    let b = parallel_results(&store, &fx.name_to_digest, &fx.vcf, 8);
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b.iter()) {
        assert_eq!(x.vrs_id, y.vrs_id);
        assert_eq!(x.pos, y.pos);
        assert_eq!(x.alt_allele, y.alt_allele);
    }
}

// FIX 3: A VCF record whose REF allele extends past the sequence end must
// cause BOTH serial and parallel paths to return Err (not silently skip).
#[test]
fn fix3_ref_past_end_returns_err_both_paths() {
    let dir = tempdir().unwrap();

    // A single short sequence "chr1" of length 10.
    let fasta_path = dir.path().join("short.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">chr1\nACGTACGTAC").unwrap(); // length 10
    }

    // A VCF record whose REF (ACGTACGTAC, len=10) starting at 1-based pos=5
    // (0-based=4) would need bytes [4..14], extending past the end (len=10).
    let vcf_path = dir.path().join("bad.vcf");
    {
        let mut f = std::fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        // POS=5 (0-based=4), REF=ACGTACGTAC (len=10) → end=14 > seq_len=10 → error
        writeln!(f, "chr1\t5\t.\tACGTACGTAC\tA\t.\tPASS\t.").unwrap();
    }

    let store_dir = dir.path().join("store");
    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
        .unwrap();
    store.load_all_sequences().unwrap();

    let paged = store.list_collections(0, 10, &[]).unwrap();
    let coll_digest = paged.results[0].digest.clone();
    let collection = store.get_collection(&coll_digest).unwrap();
    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    for rec in &collection.sequences {
        let m = rec.metadata();
        name_to_digest.insert(m.name.clone(), m.sha512t24u.clone());
    }
    let store_dir_str = store_dir.to_str().unwrap().to_string();
    drop(store);

    // Serial path must return Err.
    let serial_store = {
        let mut s = RefgetStore::open_local(&store_dir_str).unwrap();
        s.load_all_collections().unwrap();
        for d in name_to_digest.values() {
            s.load_sequence(d.as_str()).unwrap();
        }
        s.into_readonly()
    };
    let serial_result = compute_vrs_ids_streaming_readonly(
        &serial_store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        |_| {},
    );
    assert!(
        serial_result.is_err(),
        "serial path must return Err for ref-past-end variant"
    );

    // Parallel path must also return Err.
    let parallel_store = {
        let mut s = RefgetStore::open_local(&store_dir_str).unwrap();
        s.load_all_collections().unwrap();
        for d in name_to_digest.values() {
            s.load_sequence(d.as_str()).unwrap();
        }
        s.into_readonly()
    };
    let parallel_result = compute_vrs_ids_parallel_encoded(
        &parallel_store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        2,
        |_| {},
    );
    assert!(
        parallel_result.is_err(),
        "parallel path must return Err for ref-past-end variant (Fix 3)"
    );
}

// ============================================================================
// BGZF-block-parallel path (multi-block BGZF VCF)
// ============================================================================

#[test]
fn bgzf_parallel_matches_serial() {
    let fx = build_encoded_vcf_fixture();
    let n2d = &fx.name_to_digest;
    let vcf = fx.vcf_bgz.as_str();

    let store = fx.open_readonly();
    let serial = serial_results(&store, n2d, vcf);
    assert!(!serial.is_empty(), "fixture produced no results");

    for &workers in &[1usize, 2, 4, 8] {
        let par = bgzf_results(&store, n2d, vcf, workers);
        assert_eq!(
            par.len(),
            serial.len(),
            "count mismatch at {workers} worker(s): serial={} parallel={}",
            serial.len(),
            par.len()
        );
        for (i, (s, p)) in serial.iter().zip(par.iter()).enumerate() {
            assert_eq!(
                s.vrs_id, p.vrs_id,
                "VRS id mismatch at record {i} with {workers} worker(s)"
            );
            assert_eq!(
                s.chrom, p.chrom,
                "chrom order mismatch at {i} ({workers} workers)"
            );
            assert_eq!(
                s.pos, p.pos,
                "pos order mismatch at {i} ({workers} workers)"
            );
            assert_eq!(
                s.ref_allele, p.ref_allele,
                "ref mismatch at {i} ({workers} workers)"
            );
            assert_eq!(
                s.alt_allele, p.alt_allele,
                "alt mismatch at {i} ({workers} workers)"
            );
        }
    }
}

#[test]
fn bgzf_parallel_is_deterministic() {
    let fx = build_encoded_vcf_fixture();
    let store = fx.open_readonly();
    let a = bgzf_results(&store, &fx.name_to_digest, &fx.vcf_bgz, 8);
    let b = bgzf_results(&store, &fx.name_to_digest, &fx.vcf_bgz, 8);
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b.iter()) {
        assert_eq!(x.vrs_id, y.vrs_id);
        assert_eq!(x.pos, y.pos);
        assert_eq!(x.alt_allele, y.alt_allele);
    }
}

#[test]
fn bgzf_parallel_rejects_plain_gzip() {
    // A non-BGZF (plain) file must be rejected with a clear Err, not silently
    // mis-read. We point the BGZF path at the uncompressed .vcf, which is not a
    // BGZF stream.
    let fx = build_encoded_vcf_fixture();
    let store = fx.open_readonly();
    let res = compute_vrs_ids_parallel_bgzf_encoded_with_sink(
        &store,
        &fx.name_to_digest,
        &fx.vcf,
        4,
        |_| {},
    );
    assert!(
        res.is_err(),
        "block-parallel path must reject non-BGZF input"
    );
}

#[test]
fn bgzf_parallel_ref_past_end_returns_err() {
    // A VCF record whose REF allele extends past the sequence end must cause the
    // block-parallel path to return Err (matching the serial path), not skip.
    let dir = tempdir().unwrap();

    let fasta_path = dir.path().join("short.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">chr1\nACGTACGTAC").unwrap(); // length 10
    }

    let mut vcf = String::new();
    vcf.push_str("##fileformat=VCFv4.2\n");
    vcf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    // POS=5 (0-based=4), REF=ACGTACGTAC (len=10) -> end=14 > seq_len=10 -> error
    vcf.push_str("chr1\t5\t.\tACGTACGTAC\tA\t.\tPASS\t.\n");

    let vcf_bgz = dir.path().join("bad.vcf.bgz");
    write_bgzf_file(&vcf_bgz, vcf.as_bytes(), 250);

    let store_dir = dir.path().join("store");
    let mut store = RefgetStore::on_disk(&store_dir).unwrap();
    store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
        .unwrap();
    store.load_all_sequences().unwrap();
    let paged = store.list_collections(0, 10, &[]).unwrap();
    let coll_digest = paged.results[0].digest.clone();
    let collection = store.get_collection(&coll_digest).unwrap();
    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    for rec in &collection.sequences {
        let m = rec.metadata();
        name_to_digest.insert(m.name.clone(), m.sha512t24u.clone());
    }
    drop(store);

    let mut s = RefgetStore::open_local(store_dir.to_str().unwrap()).unwrap();
    s.load_all_collections().unwrap();
    for d in name_to_digest.values() {
        s.load_sequence(d.as_str()).unwrap();
    }
    let store = s.into_readonly();

    let res = compute_vrs_ids_parallel_bgzf_encoded_with_sink(
        &store,
        &name_to_digest,
        vcf_bgz.to_str().unwrap(),
        2,
        |_| {},
    );
    assert!(
        res.is_err(),
        "block-parallel path must return Err for ref-past-end variant"
    );
}
