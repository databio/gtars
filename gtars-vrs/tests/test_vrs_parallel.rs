//! Tests for the parallel VCF → VRS pipeline.
//!
//! Verifies that `compute_vrs_ids_parallel` produces the same results as
//! `compute_vrs_ids_from_vcf_readonly` (same IDs, same order) across a
//! variety of VCF shapes.

use std::io::Write;

use gtars_refget::store::{FastaImportOptions, RefgetStore};
use gtars_vrs::vcf::{
    build_name_to_digest_readonly, compute_vrs_ids_from_vcf_readonly,
    compute_vrs_ids_parallel, compute_vrs_ids_parallel_bgzf,
    compute_vrs_ids_parallel_blockwise,
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

#[test]
fn test_parallel_matches_sequential() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());

    let vcf_path = dir.path().join("mixed.vcf");
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
    let parallel =
        compute_vrs_ids_parallel(&store, &name_to_digest, vcf_path.to_str().unwrap(), 4).unwrap();

    assert_eq!(serial.len(), parallel.len());
    for (a, b) in serial.iter().zip(parallel.iter()) {
        assert_eq!(a.chrom, b.chrom);
        assert_eq!(a.pos, b.pos);
        assert_eq!(a.ref_allele, b.ref_allele);
        assert_eq!(a.alt_allele, b.alt_allele);
        assert_eq!(a.vrs_id, b.vrs_id);
    }
    assert_eq!(parallel.len(), 6, "symbolic <DEL> should be skipped");
}

#[test]
fn test_parallel_single_worker() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());

    let vcf_path = dir.path().join("single.vcf");
    write_vcf(
        &vcf_path,
        "\
chr1\t5\t.\tA\tT\t.\tPASS\t.\n\
chr2\t5\t.\tG\tA\t.\tPASS\t.\n",
    );

    let out =
        compute_vrs_ids_parallel(&store, &name_to_digest, vcf_path.to_str().unwrap(), 1).unwrap();
    assert_eq!(out.len(), 2);
    for r in &out {
        assert!(r.vrs_id.starts_with("ga4gh:VA."));
        assert_eq!(r.vrs_id.len(), 9 + 32);
    }
}

#[test]
fn test_parallel_empty_vcf() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());

    let vcf_path = dir.path().join("empty.vcf");
    write_vcf(&vcf_path, "");

    let out =
        compute_vrs_ids_parallel(&store, &name_to_digest, vcf_path.to_str().unwrap(), 4).unwrap();
    assert!(out.is_empty());
}

#[test]
fn test_parallel_skips_symbolic_alleles() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());

    let vcf_path = dir.path().join("sym.vcf");
    write_vcf(
        &vcf_path,
        "\
chr1\t10\t.\tA\t<DEL>\t.\tPASS\t.\n\
chr1\t10\t.\tA\t*\t.\tPASS\t.\n\
chr1\t10\t.\tA\t.\t.\tPASS\t.\n\
chr1\t11\t.\tC\tG\t.\tPASS\t.\n",
    );

    let out =
        compute_vrs_ids_parallel(&store, &name_to_digest, vcf_path.to_str().unwrap(), 4).unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].alt_allele, "G");
}

#[test]
fn test_parallel_multiallelic_preserves_order() {
    let dir = tempdir().unwrap();
    let (store, name_to_digest) = build_readonly_fixture(dir.path());

    let vcf_path = dir.path().join("multi.vcf");
    write_vcf(&vcf_path, "chr2\t5\t.\tG\tA,T,C\t.\tPASS\t.\n");

    let out =
        compute_vrs_ids_parallel(&store, &name_to_digest, vcf_path.to_str().unwrap(), 8).unwrap();
    assert_eq!(out.len(), 3);
    assert_eq!(out[0].alt_allele, "A");
    assert_eq!(out[1].alt_allele, "T");
    assert_eq!(out[2].alt_allele, "C");
    // Distinct IDs
    assert_ne!(out[0].vrs_id, out[1].vrs_id);
    assert_ne!(out[1].vrs_id, out[2].vrs_id);
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
    let blockwise = compute_vrs_ids_parallel_blockwise(
        &store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        4,
    )
    .unwrap();
    let pipeline =
        compute_vrs_ids_parallel(&store, &name_to_digest, vcf_path.to_str().unwrap(), 4).unwrap();

    assert_eq!(serial.len(), blockwise.len());
    assert_eq!(serial.len(), pipeline.len());
    for ((a, b), c) in serial.iter().zip(blockwise.iter()).zip(pipeline.iter()) {
        assert_eq!(a.vrs_id, b.vrs_id);
        assert_eq!(a.vrs_id, c.vrs_id);
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
    let out = compute_vrs_ids_parallel_blockwise(
        &store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        1,
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
    let out = compute_vrs_ids_parallel_blockwise(
        &store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        4,
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
    let out = compute_vrs_ids_parallel_blockwise(
        &store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        4,
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
    let bgzf_parallel =
        compute_vrs_ids_parallel_bgzf(&store, &name_to_digest, bgz.to_str().unwrap(), 4).unwrap();

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
    let err = compute_vrs_ids_parallel_bgzf(
        &store,
        &name_to_digest,
        plain_gz.to_str().unwrap(),
        4,
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
    let bgzf_parallel =
        compute_vrs_ids_parallel_bgzf(&store, &name_to_digest, bgz.to_str().unwrap(), 8).unwrap();

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
    let bgzf_parallel =
        compute_vrs_ids_parallel_bgzf(&store, &name_to_digest, bgz.to_str().unwrap(), 4).unwrap();

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
    let bgzf_parallel =
        compute_vrs_ids_parallel_bgzf(&store, &name_to_digest, bgz.to_str().unwrap(), 4).unwrap();
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
    let out =
        compute_vrs_ids_parallel_bgzf(&store, &name_to_digest, bgz.to_str().unwrap(), 4).unwrap();
    assert_eq!(out.len(), 3);
    assert_eq!(out[0].alt_allele, "A");
    assert_eq!(out[1].alt_allele, "T");
    assert_eq!(out[2].alt_allele, "C");
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
    let blockwise = compute_vrs_ids_parallel_blockwise(
        &store,
        &name_to_digest,
        vcf_path.to_str().unwrap(),
        8,
    )
    .unwrap();

    // Skip ref-mismatch failures deterministically by comparing in lockstep.
    assert_eq!(serial.len(), blockwise.len(), "count mismatch");
    for (a, b) in serial.iter().zip(blockwise.iter()) {
        assert_eq!(a.vrs_id, b.vrs_id, "order mismatch at pos {}", a.pos);
    }
}
