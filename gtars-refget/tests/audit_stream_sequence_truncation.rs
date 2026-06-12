//! Audit test for FINDING #7 (MEDIUM): silent truncation in the local-file
//! branch of `ReadonlyRefgetStore::stream_sequence`.
//!
//! In `gtars-refget/src/store/readonly.rs`, `stream_sequence(digest, start, end)`
//! has a LOCAL-FILE branch that opens the `.seq` file, does
//! `seek(SeekFrom::Start(byte_start))` then `BufReader::new(file).take(byte_len)`.
//!
//! Unlike:
//!   - the REMOTE path, which validates the response is HTTP 206 Partial Content
//!     and refuses to stream otherwise (avoiding silent corruption), and
//!   - `posread::read_exact_window` (used by the partial-read `get_substring`
//!     path), which errors on a short read,
//!
//! the `.take(byte_len)` adapter SILENTLY yields FEWER bytes (no error) when the
//! `.seq` file on disk is shorter than expected (truncated / corrupt). For a
//! Raw-mode sequence (1 ASCII byte == 1 base) this returns a silently truncated
//! sequence rather than erroring.
//!
//! This test demonstrates the discrepancy:
//!   1. Build a Raw-mode store, import a single chromosome of known sequence,
//!      persist it to disk.
//!   2. Reopen as a readonly (disk-backed, stub) store and confirm
//!      `stream_sequence` returns the correct full sequence.
//!   3. Truncate the on-disk `.seq` file so it is shorter than the sequence.
//!   4. Reopen a FRESH readonly store (bypass any in-memory cache) and call
//!      `stream_sequence` for a range that extends into the truncated region.
//!   5. Assert the EXPECTED-CORRECT behavior: the read must either return an
//!      error OR yield the full requested byte count. If it instead returns
//!      fewer bytes with `Ok` (silent truncation), the assertion FAILS,
//!      confirming the bug.
//!
//! Requires the `filesystem` feature (disk-backed store / local-file streaming).

use std::fs;
use std::io::Read;

use gtars_refget::store::{FastaImportOptions, RefgetStore, StorageMode};

/// Build a Raw-mode, disk-backed store from a single-chrom FASTA, persist it,
/// and return (sequence digest, sequence length in bases).
fn build_raw_disk_store(
    store_dir: &std::path::Path,
    fasta_content: &str,
) -> (String, usize) {
    // Write the FASTA to a temp file.
    let fasta_path = store_dir.join("input.fa");
    fs::write(&fasta_path, fasta_content).expect("write FASTA");

    // Build in-memory in RAW mode (so on-disk bytes are 1 ASCII byte per base,
    // making the truncation silent rather than surfaced by a decoder EOF).
    let mut store = RefgetStore::in_memory();
    store.disable_encoding();
    assert_eq!(store.storage_mode(), StorageMode::Raw);

    store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
        .expect("import FASTA");

    // Persist to disk: writes `.seq` files (Raw bytes) and index files.
    store
        .enable_persistence(store_dir)
        .expect("enable persistence (persist to disk)");

    // Grab the single sequence's digest + length from the readonly view.
    let ro = store.into_readonly();
    let metas: Vec<_> = ro.sequence_metadata().cloned().collect();
    assert_eq!(metas.len(), 1, "expected exactly one sequence");
    let digest = metas[0].sha512t24u.clone();
    let length = metas[0].length;

    (digest, length)
}

/// Read a `stream_sequence` reader to completion into a String.
fn drain_to_string(mut reader: Box<dyn Read + Send>) -> std::io::Result<String> {
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;
    Ok(String::from_utf8_lossy(&buf).into_owned())
}

#[test]
fn audit_stream_sequence_silent_truncation_local_file() {
    // A known single chromosome. 64 bases, plain ACGT so it is valid in Raw mode.
    let full_seq: String = "ACGT".repeat(16); // 64 bases
    assert_eq!(full_seq.len(), 64);
    let fasta_content = format!(">chr1\n{}\n", full_seq);

    let store_tmp = tempfile::tempdir().expect("tempdir");
    let store_dir = store_tmp.path();

    let (digest, length) = build_raw_disk_store(store_dir, &fasta_content);
    assert_eq!(length, 64);

    // ----- Sanity check: reopen readonly from disk, stream the full sequence. -----
    {
        let ro = RefgetStore::open_local(store_dir)
            .expect("open_local")
            .into_readonly();
        // Sequence should be a stub (data on disk), so stream_sequence exercises
        // the disk-backed local-file branch rather than the in-memory ArcSlice.
        assert!(
            !ro.is_sequence_loaded(&digest),
            "expected sequence to be a stub (disk-backed), not resident in memory; \
             the disk-streaming branch would not be exercised otherwise"
        );

        let reader = ro
            .stream_sequence(&digest, None, None)
            .expect("stream full sequence");
        let got = drain_to_string(reader).expect("drain full");
        assert_eq!(
            got, full_seq,
            "sanity: full streamed sequence must match the imported sequence"
        );
    }

    // ----- Locate the on-disk .seq file and TRUNCATE it. -----
    let seq_path = {
        let ro = RefgetStore::open_local(store_dir)
            .expect("open_local for path")
            .into_readonly();
        ro.sequence_file_path(&digest)
            .expect("sequence_file_path should resolve for a disk-backed store")
    };
    assert!(
        seq_path.exists(),
        ".seq file should exist at {}",
        seq_path.display()
    );

    let original_len = fs::metadata(&seq_path).unwrap().len();
    // In Raw mode the .seq file is exactly `length` bytes (1 byte/base).
    assert_eq!(
        original_len, 64,
        "Raw .seq file should be 64 bytes (1 byte/base) before truncation"
    );

    // Truncate the file to 32 bytes (i.e. only the first 32 bases survive on disk).
    let truncated_len: u64 = 32;
    let f = fs::OpenOptions::new()
        .write(true)
        .open(&seq_path)
        .expect("open .seq for truncate");
    f.set_len(truncated_len).expect("truncate .seq file");
    drop(f);
    assert_eq!(
        fs::metadata(&seq_path).unwrap().len(),
        truncated_len,
        ".seq file should now be truncated to 32 bytes"
    );

    // ----- Reopen a FRESH readonly store (bypass any in-memory cache) and stream
    //       a range that extends past the truncation point. -----
    let ro = RefgetStore::open_local(store_dir)
        .expect("open_local after truncation")
        .into_readonly();

    // Request bases [16, 48): byte range [16, 48) in Raw mode. The file only
    // has bytes [0, 32) on disk now, so bytes [32, 48) are missing.
    let req_start: u64 = 16;
    let req_end: u64 = 48;
    let requested_bases = (req_end - req_start) as usize; // 32 bases requested

    let stream_result = ro.stream_sequence(&digest, Some(req_start), Some(req_end));

    match stream_result {
        Err(_e) => {
            // ACCEPTABLE / CORRECT: streaming refused to open a range it cannot
            // satisfy. (This is what the remote path does via the 206 check.)
            // Bug is NOT reproduced in this configuration.
        }
        Ok(reader) => {
            // The reader opened. Drain it to completion and count the bytes.
            let drained = drain_to_string(reader);
            match drained {
                Err(_e) => {
                    // ACCEPTABLE / CORRECT: the read surfaced an EOF/short-read
                    // error. Bug is NOT reproduced.
                }
                Ok(got) => {
                    let returned_bases = got.len();
                    // EXPECTED-CORRECT behavior: a successful read of a 32-base
                    // range must yield exactly 32 bases. If we got FEWER bases
                    // with Ok, that is the SILENT TRUNCATION bug.
                    assert_eq!(
                        returned_bases, requested_bases,
                        "FINDING #7 CONFIRMED: stream_sequence silently returned {} bases \
                         (Ok) for a {}-base request against a truncated .seq file, instead \
                         of erroring or returning the full requested count. The local-file \
                         branch's `.take(byte_len)` yields fewer bytes without error when \
                         the .seq file is shorter than expected.",
                        returned_bases, requested_bases
                    );
                }
            }
        }
    }
}
