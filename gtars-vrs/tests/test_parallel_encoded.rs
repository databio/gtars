//! Equivalence test: the in-process parallel, encoded-on-the-fly VRS path must
//! produce exactly the same VRS ids in exactly the same order as the serial
//! `compute_vrs_ids_streaming_readonly` reference path.
//!
//! Builds a tiny disk-backed, 2-bit-encoded refget store from a few short
//! sequences (including A-repeat regions to exercise rolling normalization),
//! writes a small VCF with SNVs, insertions, deletions, multi-allelics, and a
//! symbolic allele, then compares serial vs parallel output at several thread
//! counts.
//!
//! Both paths read from the same on-disk store using the encoded 2-bit store;
//! decoding happens on the fly in both cases. Identical output proves the paths
//! are byte-faithful to each other.

use std::collections::HashMap;
use std::io::Write;

use gtars_refget::store::{FastaImportOptions, ReadonlyRefgetStore, RefgetStore};
use gtars_vrs::vcf::{compute_vrs_ids_parallel_encoded, compute_vrs_ids_streaming_readonly, VrsResult};
use tempfile::{tempdir, TempDir};

/// Holds the fixture alive: the tempdir (store + VCF on disk), the name->digest
/// map, the VCF path, and the store directory so callers can open fresh views.
struct Fixture {
    _dir: TempDir,
    store_dir: String,
    name_to_digest: HashMap<String, String>,
    vcf: String,
}

/// A serial (encoded, decode-on-the-fly) readonly view of the fixture store.
fn open_serial(fx: &Fixture) -> ReadonlyRefgetStore {
    let mut store = RefgetStore::open_local(&fx.store_dir).unwrap();
    store.load_all_collections().unwrap();
    for d in fx.name_to_digest.values() {
        store.load_sequence(d.as_str()).unwrap();
    }
    store.into_readonly()
}

/// A parallel (encoded, decode-on-the-fly) readonly view of the fixture store.
fn open_encoded(fx: &Fixture) -> ReadonlyRefgetStore {
    let mut store = RefgetStore::open_local(&fx.store_dir).unwrap();
    store.load_all_collections().unwrap();
    for d in fx.name_to_digest.values() {
        store.load_sequence(d.as_str()).unwrap();
    }
    store.into_readonly()
}

/// A Raw-mode (decoded, 1 byte/base) readonly view of the fixture store.
///
/// `set_encoding_mode` only re-encodes/decodes sequences that are already
/// resident, so we load every needed sequence first and then `disable_encoding`
/// to convert the in-memory bytes to Raw (1 byte/base, `bytes.len() == length`).
fn open_raw(fx: &Fixture) -> ReadonlyRefgetStore {
    let mut store = RefgetStore::open_local(&fx.store_dir).unwrap();
    store.load_all_collections().unwrap();
    for d in fx.name_to_digest.values() {
        store.load_sequence(d.as_str()).unwrap();
    }
    store.disable_encoding(); // StorageMode::Raw: 1 byte/base, len == length
    store.into_readonly()
}

/// Build a small disk-backed encoded store plus the name->digest map and VCF.
fn build_fixture() -> Fixture {
    let dir = tempdir().unwrap();

    // chr1 has an A-repeat at positions 50..60 to exercise roll left/right.
    let chr1_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACAAAAAAAAAAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCTTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let chr2_seq = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG";

    let fasta_path = dir.path().join("test.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        write!(f, ">chr1\n{}\n>chr2\n{}\n", chr1_seq, chr2_seq).unwrap();
    }

    let vcf_path = dir.path().join("test.vcf");
    {
        let mut f = std::fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        // Many records so multi-thread reordering is actually exercised.
        for i in 0..200 {
            let pos = 5 + (i % 40);
            writeln!(f, "chr1\t{}\t.\tA\tT\t.\tPASS\t.", pos).unwrap();
            writeln!(f, "chr2\t{}\t.\tG\tA,T\t.\tPASS\t.", pos).unwrap();
        }
        // Insertion + deletion in the A-repeat (1-based 51 == 0-based 50).
        writeln!(f, "chr1\t51\t.\tA\tAA\t.\tPASS\t.").unwrap();
        writeln!(f, "chr1\t51\t.\tAA\tA\t.\tPASS\t.").unwrap();
        // SNV.
        writeln!(f, "chr1\t10\t.\tG\tC\t.\tPASS\t.").unwrap();
        // Symbolic allele: must be skipped by both paths.
        writeln!(f, "chr1\t20\t.\tA\t<DEL>\t.\tPASS\t.").unwrap();
        // Unknown chromosome: must be skipped by both paths.
        writeln!(f, "chrZ\t5\t.\tA\tT\t.\tPASS\t.").unwrap();
    }

    // Build a disk-backed encoded store (on_disk() defaults to StorageMode::Encoded).
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
    drop(store); // flush/close so fresh views can reopen the on-disk store

    Fixture {
        _dir: dir,
        store_dir: store_dir.to_str().unwrap().to_string(),
        name_to_digest,
        vcf: vcf_path.to_str().unwrap().to_string(),
    }
}

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
    let n =
        compute_vrs_ids_parallel_encoded(store, n2d, vcf, threads, |r| out.push(r)).unwrap();
    assert_eq!(n, out.len());
    out
}

#[test]
fn parallel_encoded_matches_serial() {
    let fx = build_fixture();
    let n2d = &fx.name_to_digest;
    let vcf = fx.vcf.as_str();

    let serial_store = open_serial(&fx);
    let serial = serial_results(&serial_store, n2d, vcf);
    assert!(!serial.is_empty(), "fixture produced no results");

    let encoded_store = open_encoded(&fx);
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
            assert_eq!(s.chrom, p.chrom, "chrom order mismatch at {i} ({threads} threads)");
            assert_eq!(s.pos, p.pos, "pos order mismatch at {i} ({threads} threads)");
            assert_eq!(s.ref_allele, p.ref_allele, "ref mismatch at {i} ({threads} threads)");
            assert_eq!(s.alt_allele, p.alt_allele, "alt mismatch at {i} ({threads} threads)");
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
    let fx = build_fixture();
    let n2d = &fx.name_to_digest;
    let vcf = fx.vcf.as_str();

    let serial_store = open_raw(&fx);
    let serial = serial_results(&serial_store, n2d, vcf);
    assert!(!serial.is_empty(), "fixture produced no results");

    let raw_store = open_raw(&fx);
    for &threads in &[1usize, 2, 4, 8] {
        let par = parallel_results(&raw_store, n2d, vcf, threads);
        assert_eq!(par.len(), serial.len(), "count mismatch at {threads} thread(s)");
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
    let fx = build_fixture();
    let store = open_encoded(&fx);
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
    use std::io::Write;
    use tempfile::tempdir;
    use gtars_refget::store::{FastaImportOptions, RefgetStore};

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
