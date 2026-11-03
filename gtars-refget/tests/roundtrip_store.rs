//! End-to-end round-trip tests through `RefgetStore`.
//!
//! The rule under test: every store read path, in every storage mode, must
//! return bytes identical to the original (uppercased) FASTA input, and
//! those bytes must hash to the digest recorded in the sequence metadata.
//! Unlike `test_roundtrip_regressions.rs` (targeted, named bugs) and
//! `roundtrip_alphabet_bytes.rs` (byte-level, encoder-only), this file
//! exercises the full store: multi-line and mixed-case FASTA import,
//! `stream_sequence` (full/ranged/via MD5), `get_substring_bytes`, every
//! `StorageMode`, mode switches, on-disk save/reopen, and export.

use std::collections::HashMap;
use std::io::{Read, Write};
use std::path::Path;

use gtars_refget::digest::{SequenceMetadata, md5, sha512t24u};
use gtars_refget::store::{FastaImportOptions, ReadonlyRefgetStore, RefgetStore, StorageMode};

// =========================================================================
// Fixture data
// =========================================================================

/// (name, raw FASTA sequence as written -- some are mixed/lower case) for
/// every fixture record. Each record exercises a distinct alphabet or a
/// known-bad symbol; see the plan for the rationale behind each one.
fn fixture_records() -> Vec<(&'static str, &'static str)> {
    vec![
        ("dna2bit", "ACGTACGTACacgt"),
        ("dna3bit", "ACGTNNNRYRYACGT"),
        (
            "iupac_all",
            "ACGTURYSWKMBDHVNACGTURYSWKMBDHVNACGTURYSWKMBDHVNA",
        ),
        ("iupac_dh", "DDDHHHVVVDHVDHV"),
        ("rna", "ACGUACGUUUUAGCU"),
        ("rna_lower", "acguacguuuu"),
        ("protein_std", "MEFILPQACDEGHKNRSTVWY*"),
        ("selenoprotein", "MPRLLUGSEEAUVLLK*"),
        ("protein_bzoj", "MEBZOJXBZOJK*"),
        ("protein_b_only", "MEFBBBILPQ"),
        ("protein_iupac_letters", "MCVTYHNGTGYC"),
        ("protein_gaps", "MEF-ILP.QX*"),
        ("ascii_misc", "MEF123ILP"),
    ]
}

/// Write the fixture records to a FASTA file, wrapped at 10 chars/line so
/// the multi-line join path in the parser runs.
fn write_fixture_fasta(path: &Path) {
    let mut f = std::fs::File::create(path).expect("create fixture fasta");
    for (name, seq) in fixture_records() {
        writeln!(f, ">{name}").unwrap();
        for chunk in seq.as_bytes().chunks(10) {
            f.write_all(chunk).unwrap();
            f.write_all(b"\n").unwrap();
        }
    }
}

/// Expected bytes for each fixture: the uppercased, line-joined sequence --
/// what the ingest pipeline uppercases every sequence byte to before it
/// digests, guesses or encodes.
fn expected_sequences() -> HashMap<&'static str, Vec<u8>> {
    fixture_records()
        .into_iter()
        .map(|(name, seq)| (name, seq.to_ascii_uppercase().into_bytes()))
        .collect()
}

// =========================================================================
// Assertion helpers
// =========================================================================

/// Compare `got` against `expected` and the digests in `meta`, pushing a
/// readable failure message instead of panicking so one run reports every
/// broken (label, name) combination.
fn assert_lossless(
    label: &str,
    name: &str,
    got: &[u8],
    meta: &SequenceMetadata,
    expected: &[u8],
    failures: &mut Vec<String>,
) {
    if got != expected {
        let diff_at = got
            .iter()
            .zip(expected.iter())
            .position(|(a, b)| a != b)
            .unwrap_or_else(|| got.len().min(expected.len()));
        let ctx = |b: &[u8], i: usize| {
            let start = i.saturating_sub(3);
            let end = (i + 4).min(b.len());
            String::from_utf8_lossy(&b[start..end]).into_owned()
        };
        failures.push(format!(
            "{label}/{name}: mismatch at index {diff_at} (got len={}, expected len={}; \
             got ...{}..., expected ...{}...)",
            got.len(),
            expected.len(),
            ctx(got, diff_at),
            ctx(expected, diff_at),
        ));
        return;
    }
    let got_sha = sha512t24u(got);
    if got_sha != meta.sha512t24u {
        failures.push(format!(
            "{label}/{name}: sha512t24u(got)={got_sha} != metadata.sha512t24u={}",
            meta.sha512t24u
        ));
    }
    let got_md5 = md5(got);
    if got_md5 != meta.md5 {
        failures.push(format!(
            "{label}/{name}: md5(got)={got_md5} != metadata.md5={}",
            meta.md5
        ));
    }
    if got.len() != meta.length {
        failures.push(format!(
            "{label}/{name}: got.len()={} != metadata.length={}",
            got.len(),
            meta.length
        ));
    }
}

/// Run every store read path against one sequence and record failures.
fn read_paths(
    store: &ReadonlyRefgetStore,
    meta: &SequenceMetadata,
    expected: &[u8],
    label: &str,
    failures: &mut Vec<String>,
) {
    let digest = meta.sha512t24u.clone();
    let len = expected.len();

    // Full stream_sequence.
    {
        let mut reader = store
            .stream_sequence(&digest, None, None)
            .unwrap_or_else(|e| {
                panic!("{label}/{}: stream_sequence(full) failed: {e}", meta.name)
            });
        let mut got = Vec::new();
        reader
            .read_to_end(&mut got)
            .unwrap_or_else(|e| panic!("{label}/{}: read_to_end failed: {e}", meta.name));
        assert_lossless(
            &format!("{label} stream_sequence(full)"),
            &meta.name,
            &got,
            meta,
            expected,
            failures,
        );
    }

    // Ranged stream_sequence: (1,len), (3,len-2), every single-base window.
    let mut ranges: Vec<(u64, u64)> = Vec::new();
    if len >= 1 {
        ranges.push((1, len as u64));
    }
    if len >= 5 {
        ranges.push((3, (len - 2) as u64));
    }
    for i in 0..len {
        ranges.push((i as u64, (i + 1) as u64));
    }
    for (s, e) in ranges {
        if s > e || e as usize > len {
            continue;
        }
        let mut reader = store
            .stream_sequence(&digest, Some(s), Some(e))
            .unwrap_or_else(|err| {
                panic!(
                    "{label}/{}: stream_sequence({s},{e}) failed: {err}",
                    meta.name
                )
            });
        let mut got = Vec::new();
        reader
            .read_to_end(&mut got)
            .unwrap_or_else(|err| panic!("{label}/{}: read failed: {err}", meta.name));
        let want = &expected[s as usize..e as usize];
        if got != want {
            failures.push(format!(
                "{label}/{}: stream_sequence({s},{e}) = {:?}, want {:?}",
                meta.name,
                String::from_utf8_lossy(&got),
                String::from_utf8_lossy(want),
            ));
        }
    }

    // stream_sequence via the MD5 digest (lookup fallback path), full read.
    {
        let mut reader = store
            .stream_sequence(&meta.md5, None, None)
            .unwrap_or_else(|e| panic!("{label}/{}: stream_sequence(md5) failed: {e}", meta.name));
        let mut got = Vec::new();
        reader
            .read_to_end(&mut got)
            .unwrap_or_else(|e| panic!("{label}/{}: read_to_end(md5) failed: {e}", meta.name));
        assert_lossless(
            &format!("{label} stream_sequence(md5)"),
            &meta.name,
            &got,
            meta,
            expected,
            failures,
        );
    }

    // get_substring_bytes: (0,len) and (3,len-2).
    {
        let got = store
            .get_substring_bytes(&digest, 0, len)
            .unwrap_or_else(|e| {
                panic!(
                    "{label}/{}: get_substring_bytes(0,{len}) failed: {e}",
                    meta.name
                )
            });
        assert_lossless(
            &format!("{label} get_substring_bytes(full)"),
            &meta.name,
            &got,
            meta,
            expected,
            failures,
        );
    }
    if len >= 5 {
        let (s, e) = (3usize, len - 2);
        let got = store.get_substring_bytes(&digest, s, e).unwrap_or_else(|err| {
            panic!(
                "{label}/{}: get_substring_bytes({s},{e}) failed: {err}",
                meta.name
            )
        });
        let want = &expected[s..e];
        if got != want {
            failures.push(format!(
                "{label}/{}: get_substring_bytes({s},{e}) = {:?}, want {:?}",
                meta.name,
                String::from_utf8_lossy(&got),
                String::from_utf8_lossy(want),
            ));
        }
    }

    // Resident, Encoded-mode-only: the SequenceRecord::decode() path. Zstd
    // is skipped because decode() does not understand zstd frames -- a
    // separate gap, not addressed by this plan.
    if store.storage_mode() == StorageMode::Encoded && store.is_sequence_loaded(&digest) {
        if let Ok(record) = store.get_sequence(&digest) {
            match record.decode() {
                Some(decoded) => assert_lossless(
                    &format!("{label} SequenceRecord::decode()"),
                    &meta.name,
                    decoded.as_bytes(),
                    meta,
                    expected,
                    failures,
                ),
                None => failures.push(format!(
                    "{label}/{}: SequenceRecord::decode() returned None for a resident record",
                    meta.name
                )),
            }
        }
    }
}

/// Run `read_paths` for every sequence in `store`.
fn run_all(store: &ReadonlyRefgetStore, label: &str, failures: &mut Vec<String>) {
    let expected_map = expected_sequences();
    let metas: Vec<SequenceMetadata> = store.sequence_metadata().cloned().collect();
    assert_eq!(
        metas.len(),
        fixture_records().len(),
        "{label}: expected {} sequences in the store, found {}",
        fixture_records().len(),
        metas.len()
    );
    for meta in &metas {
        let expected = expected_map
            .get(meta.name.as_str())
            .unwrap_or_else(|| panic!("{label}: no expected sequence for {}", meta.name));
        read_paths(store, meta, expected, label, failures);
    }
}

fn build_in_memory(mode: StorageMode, fasta_path: &Path) -> RefgetStore {
    let mut store = RefgetStore::in_memory();
    store.set_encoding_mode(mode).unwrap();
    store
        .add_sequence_collection_from_fasta(fasta_path, FastaImportOptions::new())
        .expect("import fasta");
    store
}

// =========================================================================
// Step 9: in-memory tests, one per mode
// =========================================================================

#[test]
fn in_memory_encoded_is_lossless() {
    let dir = tempfile::tempdir().unwrap();
    let fasta_path = dir.path().join("fixtures.fa");
    write_fixture_fasta(&fasta_path);

    let store = build_in_memory(StorageMode::Encoded, &fasta_path);
    let ro = store.into_readonly();

    let mut failures = Vec::new();
    run_all(&ro, "in_memory/encoded", &mut failures);
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn in_memory_raw_is_lossless() {
    let dir = tempfile::tempdir().unwrap();
    let fasta_path = dir.path().join("fixtures.fa");
    write_fixture_fasta(&fasta_path);

    let store = build_in_memory(StorageMode::Raw, &fasta_path);
    let ro = store.into_readonly();

    let mut failures = Vec::new();
    run_all(&ro, "in_memory/raw", &mut failures);
    assert!(
        failures.is_empty(),
        "Raw is the control mode (no bit-packing table involved) and must be \
         lossless. If this fails, the test fixture or harness is wrong, not \
         the encoding tables. {} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn in_memory_zstd_is_lossless() {
    let dir = tempfile::tempdir().unwrap();
    let fasta_path = dir.path().join("fixtures.fa");
    write_fixture_fasta(&fasta_path);

    let store = build_in_memory(StorageMode::Zstd, &fasta_path);
    let ro = store.into_readonly();

    let mut failures = Vec::new();
    run_all(&ro, "in_memory/zstd", &mut failures);
    assert!(
        failures.is_empty(),
        "Zstd stores whole ASCII records (no alphabet bit-packing table \
         involved) and must be lossless. {} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

// =========================================================================
// Step 10: mode-switch tests
// =========================================================================

#[test]
fn mode_switches_are_lossless() {
    let dir = tempfile::tempdir().unwrap();
    let fasta_path = dir.path().join("fixtures.fa");
    write_fixture_fasta(&fasta_path);

    let mut failures = Vec::new();

    // Encoded -> Raw: set_encoding_mode decodes on the way, so this is its
    // own lossy path independent of the read paths tested elsewhere.
    {
        let mut store = build_in_memory(StorageMode::Encoded, &fasta_path);
        store.set_encoding_mode(StorageMode::Raw).unwrap();
        let ro = store.into_readonly();
        run_all(&ro, "mode_switch/encoded_to_raw", &mut failures);
    }

    // Encoded -> Zstd.
    {
        let mut store = build_in_memory(StorageMode::Encoded, &fasta_path);
        store.set_encoding_mode(StorageMode::Zstd).unwrap();
        let ro = store.into_readonly();
        run_all(&ro, "mode_switch/encoded_to_zstd", &mut failures);
    }

    // Raw -> Encoded -> Raw: checks that encode_sequence and decode are
    // exact inverses on real data.
    {
        let mut store = build_in_memory(StorageMode::Raw, &fasta_path);
        store.set_encoding_mode(StorageMode::Encoded).unwrap();
        store.set_encoding_mode(StorageMode::Raw).unwrap();
        let ro = store.into_readonly();
        run_all(&ro, "mode_switch/raw_encoded_raw", &mut failures);
    }

    // Zstd -> Encoded -> Zstd.
    {
        let mut store = build_in_memory(StorageMode::Zstd, &fasta_path);
        store.set_encoding_mode(StorageMode::Encoded).unwrap();
        store.set_encoding_mode(StorageMode::Zstd).unwrap();
        let ro = store.into_readonly();
        run_all(&ro, "mode_switch/zstd_encoded_zstd", &mut failures);
    }

    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

// =========================================================================
// Step 11: on-disk tests, one per mode
// =========================================================================

fn check_on_disk_reopen(store_dir: &Path, mode: StorageMode, label: &str, failures: &mut Vec<String>) {
    // Reopen fresh; the mode must survive the round trip through
    // rgstore.json.
    let ro = RefgetStore::open_local(store_dir)
        .expect("open_local")
        .into_readonly();
    assert_eq!(
        ro.storage_mode(),
        mode,
        "{label}: storage_mode did not survive reopen"
    );

    let expected_map = expected_sequences();
    let metas: Vec<SequenceMetadata> = ro.sequence_metadata().cloned().collect();
    for meta in &metas {
        assert!(
            !ro.is_sequence_loaded(&meta.sha512t24u),
            "{label}/{}: expected a stub (disk-backed) record right after reopen, \
             so the on-disk read paths are actually exercised",
            meta.name
        );
    }
    for meta in &metas {
        let expected = expected_map.get(meta.name.as_str()).unwrap();
        read_paths(&ro, meta, expected, &format!("{label}/stub"), failures);
    }

    // Load everything resident and re-run, to cover the resident-after-load
    // branch (and, for Encoded mode, the SequenceRecord::decode() path).
    let mut store2 = RefgetStore::open_local(store_dir).expect("open_local #2");
    store2.load_all_sequences().expect("load_all_sequences");
    let ro2 = store2.into_readonly();
    for meta in &metas {
        let expected = expected_map.get(meta.name.as_str()).unwrap();
        read_paths(&ro2, meta, expected, &format!("{label}/resident"), failures);
    }
}

fn on_disk_mode_test(mode: StorageMode, label: &str) {
    let src_dir = tempfile::tempdir().unwrap();
    let fasta_path = src_dir.path().join("fixtures.fa");
    write_fixture_fasta(&fasta_path);

    let mut failures = Vec::new();

    // Variant A: build in memory, then enable_persistence (persist-after).
    {
        let store_dir = tempfile::tempdir().unwrap();
        let mut store = build_in_memory(mode, &fasta_path);
        store
            .enable_persistence(store_dir.path())
            .expect("enable_persistence");
        check_on_disk_reopen(
            store_dir.path(),
            mode,
            &format!("{label}/persist_after"),
            &mut failures,
        );
    }

    // Variant B: RefgetStore::on_disk + set_encoding_mode + import, which
    // writes through the deferred-write seam instead of persisting after
    // the fact.
    {
        let store_dir = tempfile::tempdir().unwrap();
        let mut store = RefgetStore::on_disk(store_dir.path()).expect("on_disk");
        store.set_encoding_mode(mode).unwrap();
        store
            .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
            .expect("import fasta (on_disk)");
        check_on_disk_reopen(
            store_dir.path(),
            mode,
            &format!("{label}/on_disk_deferred"),
            &mut failures,
        );
    }

    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn on_disk_encoded_is_lossless() {
    on_disk_mode_test(StorageMode::Encoded, "on_disk/encoded");
}

#[test]
fn on_disk_raw_is_lossless() {
    on_disk_mode_test(StorageMode::Raw, "on_disk/raw");
}

#[test]
fn on_disk_zstd_is_lossless() {
    on_disk_mode_test(StorageMode::Zstd, "on_disk/zstd");
}

// =========================================================================
// Step 12: export round-trip
// =========================================================================

#[test]
fn export_fasta_is_lossless() {
    let mut failures = Vec::new();

    for mode in [StorageMode::Encoded, StorageMode::Zstd] {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("fixtures.fa");
        write_fixture_fasta(&fasta_path);

        let mut store = build_in_memory(mode, &fasta_path);
        let collection_digest = store.list_collections(0, usize::MAX, &[]).unwrap().results[0]
            .digest
            .clone();

        // Line width of 7 so records wrap at odd lengths.
        let out_path = dir.path().join("exported.fa");
        store
            .export_fasta(&collection_digest, &out_path, None, Some(7))
            .expect("export_fasta");

        // Import the exported FASTA into a fresh Raw store and check that
        // the collection- and sequence-level digests match exactly. This is
        // the full seqcol-level check that export did not lose data.
        let mut reimported = RefgetStore::in_memory();
        reimported.set_encoding_mode(StorageMode::Raw).unwrap();
        reimported
            .add_sequence_collection_from_fasta(&out_path, FastaImportOptions::new())
            .expect("reimport exported fasta");
        let reimported_digest = reimported.list_collections(0, usize::MAX, &[]).unwrap().results[0]
            .digest
            .clone();

        if reimported_digest != collection_digest {
            failures.push(format!(
                "export/{mode:?}: collection digest changed across export/reimport: {collection_digest} -> {reimported_digest}"
            ));
        }

        let orig_metas: Vec<SequenceMetadata> = store.sequence_metadata().cloned().collect();
        let new_metas: HashMap<String, SequenceMetadata> = reimported
            .sequence_metadata()
            .cloned()
            .map(|m| (m.name.clone(), m))
            .collect();
        for m in &orig_metas {
            match new_metas.get(&m.name) {
                Some(nm) if nm.sha512t24u == m.sha512t24u => {}
                Some(nm) => failures.push(format!(
                    "export/{mode:?}/{}: sha512t24u changed across export/reimport: {} -> {}",
                    m.name, m.sha512t24u, nm.sha512t24u
                )),
                None => failures.push(format!(
                    "export/{mode:?}/{}: missing from the reimported collection",
                    m.name
                )),
            }
        }
    }

    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
