//! Equivalence test: the BGZF-block-parallel, encoded-on-the-fly VRS path
//! (`compute_vrs_ids_parallel_bgzf_encoded_with_sink`) must produce exactly the
//! same VRS ids in exactly the same order as the serial
//! `compute_vrs_ids_streaming_readonly` reference path.
//!
//! Unlike the single-reader parallel path, this one reads raw BGZF blocks on
//! the reader thread and decompresses/parses on the workers, then stitches
//! cross-block boundary lines back together. To exercise that stitching we
//! write the fixture VCF as a real multi-block BGZF stream with deliberately
//! small blocks (so VCF lines straddle block boundaries), plus a final BGZF EOF
//! block. Byte-identical output proves the block-parallel path is faithful to
//! the serial reference at 1/2/4/8 workers, that it is deterministic, and that
//! worker errors (ref-past-end) surface as `Err`.

use std::collections::HashMap;
use std::io::Write;

use flate2::Compression;
use flate2::write::DeflateEncoder;
use gtars_refget::store::{FastaImportOptions, ReadonlyRefgetStore, RefgetStore};
use gtars_vrs::vcf::{
    compute_vrs_ids_parallel_bgzf_encoded_with_sink, compute_vrs_ids_streaming_readonly, VrsResult,
};
use tempfile::{tempdir, TempDir};

/// Write `payload` as one BGZF block (gzip member with the BGZF `BC` extra
/// subfield) onto `out`.
fn write_bgzf_block(out: &mut Vec<u8>, payload: &[u8]) {
    // Raw DEFLATE of the payload.
    let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
    enc.write_all(payload).unwrap();
    let deflated = enc.finish().unwrap();

    // Total block size = 12 (fixed header) + 6 (extra: BC subfield) + deflate +
    // 8 (CRC32 + ISIZE). BSIZE stored in the extra field is total - 1.
    let bsize = 12 + 6 + deflated.len() + 8 - 1;
    assert!(bsize <= u16::MAX as usize, "block too large for BGZF");

    // Fixed gzip header with FEXTRA set.
    out.extend_from_slice(&[
        0x1f, 0x8b, // magic
        0x08, // CM = deflate
        0x04, // FLG = FEXTRA
        0x00, 0x00, 0x00, 0x00, // MTIME
        0x00, // XFL
        0xff, // OS = unknown
    ]);
    // XLEN = 6 (one BC subfield).
    out.extend_from_slice(&(6u16).to_le_bytes());
    // Extra subfield: SI1='B', SI2='C', SLEN=2, BSIZE (u16 LE).
    out.extend_from_slice(&[b'B', b'C', 0x02, 0x00]);
    out.extend_from_slice(&(bsize as u16).to_le_bytes());
    // Deflate payload.
    out.extend_from_slice(&deflated);
    // CRC32 + ISIZE (uncompressed length mod 2^32).
    let crc = crc32(payload);
    out.extend_from_slice(&crc.to_le_bytes());
    out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
}

/// The 28-byte BGZF EOF marker block (empty payload).
fn write_bgzf_eof(out: &mut Vec<u8>) {
    write_bgzf_block(out, &[]);
}

/// Minimal CRC32 (IEEE) implementation so the test has no extra deps.
fn crc32(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;
    for &b in data {
        crc ^= b as u32;
        for _ in 0..8 {
            let mask = (crc & 1).wrapping_neg();
            crc = (crc >> 1) ^ (0xEDB8_8320 & mask);
        }
    }
    !crc
}

/// Write `data` to `path` as a multi-block BGZF stream with ~`block_payload`
/// byte payloads, so VCF lines straddle block boundaries, then a BGZF EOF block.
fn write_bgzf_file(path: &std::path::Path, data: &[u8], block_payload: usize) {
    let mut out = Vec::new();
    let mut i = 0;
    while i < data.len() {
        let end = (i + block_payload).min(data.len());
        write_bgzf_block(&mut out, &data[i..end]);
        i = end;
    }
    write_bgzf_eof(&mut out);
    std::fs::write(path, &out).unwrap();
}

/// Holds the fixture alive: the tempdir (store + BGZF VCF on disk), the
/// name->digest map, the BGZF VCF path, and the store directory.
struct Fixture {
    _dir: TempDir,
    store_dir: String,
    name_to_digest: HashMap<String, String>,
    vcf_bgz: String,
}

fn open_readonly(fx: &Fixture) -> ReadonlyRefgetStore {
    let mut store = RefgetStore::open_local(&fx.store_dir).unwrap();
    store.load_all_collections().unwrap();
    for d in fx.name_to_digest.values() {
        store.load_sequence(d.as_str()).unwrap();
    }
    store.into_readonly()
}

/// Build a small disk-backed encoded store plus the name->digest map and a real
/// multi-block BGZF VCF (small blocks so lines span boundaries).
fn build_fixture() -> Fixture {
    let dir = tempdir().unwrap();

    let chr1_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACAAAAAAAAAAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCTTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let chr2_seq = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG";

    let fasta_path = dir.path().join("test.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        write!(f, ">chr1\n{}\n>chr2\n{}\n", chr1_seq, chr2_seq).unwrap();
    }

    // Build the VCF text in memory.
    let mut vcf = String::new();
    vcf.push_str("##fileformat=VCFv4.2\n");
    vcf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for i in 0..200 {
        let pos = 5 + (i % 40);
        // A long-ish INFO field makes some lines straddle small BGZF blocks.
        vcf.push_str(&format!(
            "chr1\t{pos}\t.\tA\tT\t.\tPASS\tAF=0.1;INFO_PAD=XXXXXXXXXXXXXXXXXXXX\n"
        ));
        vcf.push_str(&format!("chr2\t{pos}\t.\tG\tA,T\t.\tPASS\t.\n"));
    }
    vcf.push_str("chr1\t51\t.\tA\tAA\t.\tPASS\t.\n");
    vcf.push_str("chr1\t51\t.\tAA\tA\t.\tPASS\t.\n");
    vcf.push_str("chr1\t10\t.\tG\tC\t.\tPASS\t.\n");
    vcf.push_str("chr1\t20\t.\tA\t<DEL>\t.\tPASS\t.\n");
    vcf.push_str("chrZ\t5\t.\tA\tT\t.\tPASS\t.\n");

    // Plain-text VCF (for the serial reference, which reads via open_vcf).
    let vcf_plain = dir.path().join("test.vcf");
    std::fs::write(&vcf_plain, vcf.as_bytes()).unwrap();

    // Multi-block BGZF VCF with small (250-byte) payloads to force lines to span
    // block boundaries, exercising the head/tail stitching.
    let vcf_bgz = dir.path().join("test.vcf.bgz");
    write_bgzf_file(&vcf_bgz, vcf.as_bytes(), 250);

    // Build a disk-backed encoded store.
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

    Fixture {
        _dir: dir,
        store_dir: store_dir.to_str().unwrap().to_string(),
        name_to_digest,
        vcf_bgz: vcf_bgz.to_str().unwrap().to_string(),
    }
}

fn serial_results(
    store: &ReadonlyRefgetStore,
    n2d: &HashMap<String, String>,
    vcf: &str,
) -> Vec<VrsResult> {
    let mut out = Vec::new();
    // The serial path reads the BGZF file too (via MultiGzDecoder in open_vcf).
    let n = compute_vrs_ids_streaming_readonly(store, n2d, vcf, |r| out.push(r)).unwrap();
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
    let n = compute_vrs_ids_parallel_bgzf_encoded_with_sink(store, n2d, vcf, workers, |r| {
        out.push(r)
    })
    .unwrap();
    assert_eq!(n, out.len());
    out
}

#[test]
fn bgzf_parallel_matches_serial() {
    let fx = build_fixture();
    let n2d = &fx.name_to_digest;
    let vcf = fx.vcf_bgz.as_str();

    let store = open_readonly(&fx);
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
            assert_eq!(s.chrom, p.chrom, "chrom order mismatch at {i} ({workers} workers)");
            assert_eq!(s.pos, p.pos, "pos order mismatch at {i} ({workers} workers)");
            assert_eq!(s.ref_allele, p.ref_allele, "ref mismatch at {i} ({workers} workers)");
            assert_eq!(s.alt_allele, p.alt_allele, "alt mismatch at {i} ({workers} workers)");
        }
    }
}

#[test]
fn bgzf_parallel_is_deterministic() {
    let fx = build_fixture();
    let store = open_readonly(&fx);
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
    let fx = build_fixture();
    let store = open_readonly(&fx);
    let plain_vcf = fx.vcf_bgz.replace("test.vcf.bgz", "test.vcf");
    let res = compute_vrs_ids_parallel_bgzf_encoded_with_sink(
        &store,
        &fx.name_to_digest,
        &plain_vcf,
        4,
        |_| {},
    );
    assert!(res.is_err(), "block-parallel path must reject non-BGZF input");
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
