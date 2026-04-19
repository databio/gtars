//! Tests for the parallel VCF → VRS pipeline.
//!
//! Verifies that `compute_vrs_ids_parallel_blockwise_with_sink` and
//! `compute_vrs_ids_parallel_bgzf_with_sink` produce the same results as
//! `compute_vrs_ids_from_vcf_readonly` across a variety of VCF shapes.

use std::io::Write;

use gtars_refget::store::{FastaImportOptions, RefgetStore};
use gtars_vrs::vcf::{
    build_name_to_digest_readonly, compute_vrs_ids_from_vcf_readonly,
    compute_vrs_ids_parallel_bgzf_with_sink, compute_vrs_ids_parallel_blockwise_with_sink,
    decode_vcf_chroms, VrsResult,
};
use tempfile::tempdir;

const CHR1: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACAAAAAAAAAAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCTTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
const CHR2: &str = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG";

/// Build a preloaded readonly store + name_to_digest map from two canned chroms.
fn build_readonly_fixture(dir: &std::path::Path) -> (
    gtars_refget::store::ReadonlyRefgetStore,
    std::collections::HashMap<String, String>,
) {
    let fasta_path = dir.join("test.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        write!(f, ">chr1\n{}\n>chr2\n{}\n", CHR1, CHR2).unwrap();
    }

    let mut store = RefgetStore::in_memory();
    store.disable_encoding();
    store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
        .unwrap();

    let paged = store.list_collections(0, 10, &[]).unwrap();
    let digest = paged.results[0].digest.clone();

    store.load_all_sequences().unwrap();
    // Populate the decoded cache so sequence_bytes(&self) returns Some().
    let collection = store.get_collection(&digest).unwrap().clone();
    for seq in &collection.sequences {
        store.ensure_decoded(seq.metadata().sha512t24u.as_str()).unwrap();
    }
    let readonly = store.into_readonly();
    let name_to_digest = build_name_to_digest_readonly(&readonly, &digest).unwrap();
    (readonly, name_to_digest)
}

fn write_vcf(path: &std::path::Path, body: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    f.write_all(body.as_bytes()).unwrap();
}

/// Write a BGZF-encoded VCF using noodles-bgzf. The test exercises
/// boundary stitching by forcing small blocks (the caller picks the
/// content; noodles chooses block boundaries naturally).
fn write_vcf_bgzf(path: &std::path::Path, body: &str) {
    use noodles_bgzf::io::writer::Writer as BgzfWriter;
    let file = std::fs::File::create(path).unwrap();
    let mut w = BgzfWriter::new(file);
    write!(w, "##fileformat=VCFv4.2\n").unwrap();
    write!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n").unwrap();
    w.write_all(body.as_bytes()).unwrap();
    w.flush().unwrap();
}

// ── Blockwise variant tests ────────────────────────────────────────────

#[test]
fn test_blockwise_matches_sequential() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());

    let vcf_path = dir.path().join("mixed_bw.vcf");
    write_vcf(
        &vcf_path,
        "\
chr1\t5\t.\tA\tT\t.\tPASS\t.\n\
chr1\t10\t.\tG\tC\t.\tPASS\t.\n\
chr1\t51\t.\tA\tAA\t.\tPASS\t.\n\
chr1\t51\t.\tAA\tA\t.\tPASS\t.\n\
chr2\t5\t.\tG\tA,T\t.\tPASS\t.\n\
chr1\t20\t.\tA\t<DEL>\t.\tPASS\t.\n",
    );

    let serial =
        compute_vrs_ids_from_vcf_readonly(&store, &name_to_digest, vcf_path.to_str().unwrap())
            .unwrap();
    let mut blockwise: Vec<VrsResult> = Vec::new();
    compute_vrs_ids_parallel_blockwise_with_sink(
        &store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        4,
        |r| blockwise.push(r),
    )
    .unwrap();

    assert_eq!(serial.len(), blockwise.len());
    for (a, b) in serial.iter().zip(blockwise.iter()) {
        assert_eq!(a.vrs_id, b.vrs_id);
        assert_eq!(a.chrom, b.chrom);
        assert_eq!(a.pos, b.pos);
    }
}

#[test]
fn test_blockwise_single_worker() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());
    let vcf_path = dir.path().join("single_bw.vcf");
    write_vcf(
        &vcf_path,
        "chr1\t5\t.\tA\tT\t.\tPASS\t.\nchr2\t5\t.\tG\tA\t.\tPASS\t.\n",
    );
    let mut out: Vec<VrsResult> = Vec::new();
    compute_vrs_ids_parallel_blockwise_with_sink(
        &store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        1,
        |r| out.push(r),
    )
    .unwrap();
    assert_eq!(out.len(), 2);
    for r in &out {
        assert!(r.vrs_id.starts_with("ga4gh:VA."));
    }
}

#[test]
fn test_blockwise_empty_vcf() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());
    let vcf_path = dir.path().join("empty_bw.vcf");
    write_vcf(&vcf_path, "");
    let mut out: Vec<VrsResult> = Vec::new();
    compute_vrs_ids_parallel_blockwise_with_sink(
        &store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        4,
        |r| out.push(r),
    )
    .unwrap();
    assert!(out.is_empty());
}

#[test]
fn test_blockwise_skips_symbolic_alleles() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());
    let vcf_path = dir.path().join("sym_bw.vcf");
    write_vcf(
        &vcf_path,
        "\
chr1\t10\t.\tA\t<DEL>\t.\tPASS\t.\n\
chr1\t10\t.\tA\t*\t.\tPASS\t.\n\
chr1\t10\t.\tA\t.\t.\tPASS\t.\n\
chr1\t11\t.\tC\tG\t.\tPASS\t.\n",
    );
    let mut out: Vec<VrsResult> = Vec::new();
    compute_vrs_ids_parallel_blockwise_with_sink(
        &store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        4,
        |r| out.push(r),
    )
    .unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].alt_allele, "G");
}

// ── BGZF-block parallel tests ──────────────────────────────────────────

#[test]
fn test_bgzf_matches_sequential_small() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());
    let plain = dir.path().join("mixed_plain.vcf");
    let bgz = dir.path().join("mixed.vcf.gz");
    let body = "\
chr1\t5\t.\tA\tT\t.\tPASS\t.\n\
chr1\t10\t.\tG\tC\t.\tPASS\t.\n\
chr1\t51\t.\tA\tAA\t.\tPASS\t.\n\
chr1\t51\t.\tAA\tA\t.\tPASS\t.\n\
chr2\t5\t.\tG\tA,T\t.\tPASS\t.\n";
    write_vcf(&plain, body);
    write_vcf_bgzf(&bgz, body);

    let serial =
        compute_vrs_ids_from_vcf_readonly(&store, &name_to_digest, plain.to_str().unwrap())
            .unwrap();
    let mut bgzf_parallel: Vec<VrsResult> = Vec::new();
    compute_vrs_ids_parallel_bgzf_with_sink(
        &store,
        &name_to_digest,
        bgz.to_str().unwrap(),
        4,
        |r| bgzf_parallel.push(r),
    )
    .unwrap();

    assert_eq!(serial.len(), bgzf_parallel.len());
    for (a, b) in serial.iter().zip(bgzf_parallel.iter()) {
        assert_eq!(a.vrs_id, b.vrs_id);
        assert_eq!(a.chrom, b.chrom);
        assert_eq!(a.pos, b.pos);
    }
}

#[test]
fn test_bgzf_rejects_plain_gzip() {
    use flate2::write::GzEncoder;
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());
    let plain_gz = dir.path().join("plain.vcf.gz");
    {
        let f = std::fs::File::create(&plain_gz).unwrap();
        let mut w = GzEncoder::new(f, flate2::Compression::default());
        writeln!(w, "##fileformat=VCFv4.2").unwrap();
        writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(w, "chr1\t5\t.\tA\tT\t.\tPASS\t.").unwrap();
        w.finish().unwrap();
    }
    let err = compute_vrs_ids_parallel_bgzf_with_sink(
        &store,
        &name_to_digest,
        plain_gz.to_str().unwrap(),
        4,
        |_r| {},
    );
    assert!(err.is_err(), "plain gzip must be rejected");
}

#[test]
fn test_bgzf_spans_multiple_blocks() {
    // Generate enough content that noodles will emit multiple BGZF blocks
    // (BGZF default block size ≈ 64 KB). 3000 lines at ~30 bytes ≈ 90 KB.
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());
    let bgz = dir.path().join("many.vcf.gz");
    let plain = dir.path().join("many.vcf");
    let mut body = String::new();
    for i in 1..=3000 {
        let pos = ((i - 1) % 150) + 1;
        let alt = if i % 2 == 0 { "T" } else { "G" };
        body.push_str(&format!("chr1\t{}\t.\tN\t{}\t.\tPASS\t.\n", pos, alt));
    }
    write_vcf(&plain, &body);
    write_vcf_bgzf(&bgz, &body);

    let serial =
        compute_vrs_ids_from_vcf_readonly(&store, &name_to_digest, plain.to_str().unwrap())
            .unwrap();
    let mut bgzf_parallel: Vec<VrsResult> = Vec::new();
    compute_vrs_ids_parallel_bgzf_with_sink(
        &store,
        &name_to_digest,
        bgz.to_str().unwrap(),
        8,
        |r| bgzf_parallel.push(r),
    )
    .unwrap();

    assert_eq!(
        serial.len(),
        bgzf_parallel.len(),
        "BGZF parallel count mismatch (likely lost a boundary line)"
    );
    for (a, b) in serial.iter().zip(bgzf_parallel.iter()) {
        assert_eq!(a.vrs_id, b.vrs_id, "order/id mismatch at pos {}", a.pos);
    }
}

#[test]
fn test_bgzf_line_spans_multiple_blocks() {
    // Regression: a BGZF *block* (≤64 KB decompressed) can sit entirely
    // inside one VCF *line* when INFO payloads are large. Previously,
    // the block's whole content was mistakenly treated as a
    // "head_fragment" and emitted as a bogus complete line, dropping
    // the real line's continuation. Correct behavior: a block with no
    // newline contributes nothing to output and accumulates into the
    // running tail until a later block's first `\n` completes the line.
    //
    // We pad a single data line with a 200 KB INFO field to guarantee
    // that at least one full BGZF block (≤64 KB) decompresses with no
    // newline inside. Expected: 2 complete lines → 2 results.

    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());
    let bgz = dir.path().join("long_info.vcf.gz");
    let plain = dir.path().join("long_info.vcf");

    let big = "X".repeat(200_000);
    let body = format!(
        "chr1\t5\t.\tA\tT\t.\tPASS\tBIG={}\n\
chr1\t10\t.\tG\tC\t.\tPASS\t.\n",
        big
    );
    write_vcf(&plain, &body);
    write_vcf_bgzf(&bgz, &body);

    let serial =
        compute_vrs_ids_from_vcf_readonly(&store, &name_to_digest, plain.to_str().unwrap())
            .unwrap();
    let mut bgzf_parallel: Vec<VrsResult> = Vec::new();
    compute_vrs_ids_parallel_bgzf_with_sink(
        &store,
        &name_to_digest,
        bgz.to_str().unwrap(),
        4,
        |r| bgzf_parallel.push(r),
    )
    .unwrap();

    assert_eq!(
        serial.len(),
        bgzf_parallel.len(),
        "BGZF parallel lost or duplicated lines when a line spans multiple blocks"
    );
    assert_eq!(bgzf_parallel.len(), 2);
    for (a, b) in serial.iter().zip(bgzf_parallel.iter()) {
        assert_eq!(a.chrom, b.chrom);
        assert_eq!(a.pos, b.pos);
        assert_eq!(a.ref_allele, b.ref_allele);
        assert_eq!(a.alt_allele, b.alt_allele);
        assert_eq!(a.vrs_id, b.vrs_id);
    }
}

#[test]
fn test_bgzf_unterminated_final_line() {
    // File ends without a trailing `\n`: the last block's tail_fragment
    // must still be processed after the main loop.
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());
    let bgz = dir.path().join("noeol.vcf.gz");
    let plain = dir.path().join("noeol.vcf");
    // NB: write_vcf/write_vcf_bgzf already emit a terminating \n on
    // header rows; the data line we pass here has no trailing \n.
    let body = "chr1\t5\t.\tA\tT\t.\tPASS\t.";
    write_vcf(&plain, body);
    write_vcf_bgzf(&bgz, body);
    let serial =
        compute_vrs_ids_from_vcf_readonly(&store, &name_to_digest, plain.to_str().unwrap())
            .unwrap();
    let mut bgzf_parallel: Vec<VrsResult> = Vec::new();
    compute_vrs_ids_parallel_bgzf_with_sink(
        &store,
        &name_to_digest,
        bgz.to_str().unwrap(),
        4,
        |r| bgzf_parallel.push(r),
    )
    .unwrap();
    assert_eq!(serial.len(), 1);
    assert_eq!(bgzf_parallel.len(), 1);
    assert_eq!(serial[0].vrs_id, bgzf_parallel[0].vrs_id);
}

#[test]
fn test_bgzf_multiallelic() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());
    let bgz = dir.path().join("multi.vcf.gz");
    write_vcf_bgzf(&bgz, "chr2\t5\t.\tG\tA,T,C\t.\tPASS\t.\n");
    let mut out: Vec<VrsResult> = Vec::new();
    compute_vrs_ids_parallel_bgzf_with_sink(
        &store,
        &name_to_digest,
        bgz.to_str().unwrap(),
        4,
        |r| out.push(r),
    )
    .unwrap();
    assert_eq!(out.len(), 3);
    assert_eq!(out[0].alt_allele, "A");
    assert_eq!(out[1].alt_allele, "T");
    assert_eq!(out[2].alt_allele, "C");
}

// ── decode_vcf_chroms tests ─────────────────────────────────────────────

/// Build a RefgetStore (not yet loaded/decoded) with two chroms and
/// return the store + its collection digest plus the per-chrom digests.
/// In contrast to `build_readonly_fixture`, this one does NOT preload or
/// decode — the caller is expected to exercise `decode_vcf_chroms`.
fn build_unloaded_store(dir: &std::path::Path) -> (RefgetStore, String, String, String) {
    let fasta_path = dir.join("test_decode.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        write!(f, ">chr1\n{}\n>chr2\n{}\n", CHR1, CHR2).unwrap();
    }

    // Use a disk-backed store so we can write then reopen, which puts
    // sequences in the Stub state (unloaded + undecoded) and gives us a
    // realistic "slow path" starting point.
    let store_dir = dir.join("store_unloaded");
    {
        let mut store = RefgetStore::on_disk(&store_dir).unwrap();
        store
            .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
            .unwrap();
        // Drop store to flush any pending writes.
    }

    let mut store = RefgetStore::open_local(&store_dir).unwrap();
    let paged = store.list_collections(0, 10, &[]).unwrap();
    let coll_digest = paged.results[0].digest.clone();

    let coll = store.get_collection(&coll_digest).unwrap().clone();
    let chr1_digest = coll
        .sequences
        .iter()
        .find(|sr| sr.metadata().name == "chr1")
        .unwrap()
        .metadata()
        .sha512t24u
        .clone();
    let chr2_digest = coll
        .sequences
        .iter()
        .find(|sr| sr.metadata().name == "chr2")
        .unwrap()
        .metadata()
        .sha512t24u
        .clone();
    (store, coll_digest, chr1_digest, chr2_digest)
}

#[test]
fn test_decode_vcf_chroms_slow_path_only_referenced() {
    let dir = tempdir().unwrap();
    let (mut store, coll_digest, chr1_digest, chr2_digest) = build_unloaded_store(dir.path());

    // Confirm neither sequence is decoded yet.
    assert!(!store.is_sequence_decoded(chr1_digest.as_str()));
    assert!(!store.is_sequence_decoded(chr2_digest.as_str()));

    // VCF references only chr1.
    let vcf_path = dir.path().join("chr1_only.vcf");
    write_vcf(&vcf_path, "chr1\t10\t.\tA\tG\t.\tPASS\t.\n");

    decode_vcf_chroms(&mut store, &coll_digest, vcf_path.to_str().unwrap()).unwrap();

    // chr1 is decoded (referenced); chr2 is not (not referenced).
    assert!(
        store.is_sequence_decoded(chr1_digest.as_str()),
        "chr1 should be decoded after decode_vcf_chroms"
    );
    assert!(
        !store.is_sequence_decoded(chr2_digest.as_str()),
        "chr2 should NOT be decoded (not referenced by VCF)"
    );
}

#[test]
fn test_decode_vcf_chroms_unknown_chrom_silently_skipped() {
    let dir = tempdir().unwrap();
    let (mut store, coll_digest, chr1_digest, chr2_digest) = build_unloaded_store(dir.path());

    // VCF references a chrom not in the collection.
    let vcf_path = dir.path().join("unknown.vcf");
    write_vcf(&vcf_path, "chrX\t10\t.\tA\tG\t.\tPASS\t.\n");

    // Should succeed — unknown chroms are silently skipped.
    decode_vcf_chroms(&mut store, &coll_digest, vcf_path.to_str().unwrap()).unwrap();

    // Nothing got decoded.
    assert!(!store.is_sequence_decoded(chr1_digest.as_str()));
    assert!(!store.is_sequence_decoded(chr2_digest.as_str()));
}

#[test]
fn test_decode_vcf_chroms_fast_path_all_decoded() {
    let dir = tempdir().unwrap();
    let (mut store, coll_digest, chr1_digest, chr2_digest) = build_unloaded_store(dir.path());

    // Pre-decode everything.
    store.load_all_sequences().unwrap();
    store.ensure_decoded(chr1_digest.as_str()).unwrap();
    store.ensure_decoded(chr2_digest.as_str()).unwrap();
    assert!(store.is_sequence_decoded(chr1_digest.as_str()));
    assert!(store.is_sequence_decoded(chr2_digest.as_str()));

    // Point at a nonexistent VCF path; the fast path must avoid opening it.
    let nonexistent = dir.path().join("does_not_exist.vcf");
    decode_vcf_chroms(&mut store, &coll_digest, nonexistent.to_str().unwrap()).unwrap();

    // Still decoded.
    assert!(store.is_sequence_decoded(chr1_digest.as_str()));
    assert!(store.is_sequence_decoded(chr2_digest.as_str()));
}

#[test]
fn test_decode_vcf_chroms_fast_path_all_loaded() {
    let dir = tempdir().unwrap();
    let (mut store, coll_digest, chr1_digest, chr2_digest) = build_unloaded_store(dir.path());

    // Load (but don't decode) all sequences.
    store.load_all_sequences().unwrap();
    assert!(store.is_sequence_loaded(chr1_digest.as_str()));
    assert!(store.is_sequence_loaded(chr2_digest.as_str()));
    assert!(!store.is_sequence_decoded(chr1_digest.as_str()));
    assert!(!store.is_sequence_decoded(chr2_digest.as_str()));

    // Point at a nonexistent VCF path; fast-path-2 must avoid opening it.
    let nonexistent = dir.path().join("does_not_exist.vcf");
    decode_vcf_chroms(&mut store, &coll_digest, nonexistent.to_str().unwrap()).unwrap();

    // Both sequences are now decoded (fast path 2 decodes everything).
    assert!(store.is_sequence_decoded(chr1_digest.as_str()));
    assert!(store.is_sequence_decoded(chr2_digest.as_str()));
}

#[test]
fn test_blockwise_preserves_order_across_batches() {
    // Write enough lines to span multiple 1024-line batches, then assert
    // the VRS IDs come back in VCF order (batch reordering matters).
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());
    let vcf_path = dir.path().join("many_bw.vcf");
    let mut body = String::new();
    // chr1 is 193 bp in the fixture; use positions 1..=150 with simple SNVs.
    for i in 1..=3000 {
        // Cycle positions through 1..150 to avoid ref-base mismatches.
        let pos = ((i - 1) % 150) + 1;
        // Toggle alt so multi-allelic-free lines produce unique results.
        let alt = if i % 2 == 0 { "T" } else { "G" };
        body.push_str(&format!("chr1\t{}\t.\tN\t{}\t.\tPASS\t.\n", pos, alt));
    }
    write_vcf(&vcf_path, &body);

    let serial =
        compute_vrs_ids_from_vcf_readonly(&store, &name_to_digest, vcf_path.to_str().unwrap())
            .unwrap();
    let mut blockwise: Vec<VrsResult> = Vec::new();
    compute_vrs_ids_parallel_blockwise_with_sink(
        &store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        8,
        |r| blockwise.push(r),
    )
    .unwrap();

    // Skip ref-mismatch failures deterministically by comparing in lockstep.
    assert_eq!(serial.len(), blockwise.len(), "count mismatch");
    for (a, b) in serial.iter().zip(blockwise.iter()) {
        assert_eq!(a.vrs_id, b.vrs_id, "order mismatch at pos {}", a.pos);
    }
}
