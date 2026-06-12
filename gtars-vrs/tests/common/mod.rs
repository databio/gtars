//! Shared test helpers for the gtars-vrs integration test suite.
//!
//! Cargo treats `tests/common/mod.rs` as a non-test support module: because it
//! lives in a subdirectory (and is named `mod.rs`), it is NOT compiled as its
//! own test binary and does not trigger the "unused integration test" warning.
//! Files that need these helpers declare `mod common;` and call into it.
//!
//! Organization:
//!   - Pure VRS digest helpers (`vcf_equiv_vrs_id_*`) are always available.
//!   - Filesystem-only helpers (synthetic harness + parallel/VCF fixtures) are
//!     gated behind `#[cfg(feature = "filesystem")]` so this module still
//!     compiles in the `transcripts`-only (WASM-ish) configuration.

// Each test binary that does `mod common;` only uses a subset of these helpers,
// so per-binary `dead_code` warnings are expected and benign for a shared
// support module. This is the standard idiom for `tests/common/mod.rs`.
#![allow(dead_code)]

// ============================================================================
// Pure VRS digest helpers (always available)
// ============================================================================

use gtars_vrs::digest::DigestWriter;
use gtars_vrs::normalize::normalize;

/// Compute the equivalent VCF-derived VRS ID for comparison, taking a
/// pre-formatted `SQ.<digest>` sequence reference (the `hgvs_bridge` form).
#[allow(dead_code)]
pub fn vcf_equiv_vrs_id_with_sq(
    seq: &[u8],
    sq: &str,
    pos_ib: u64,
    refb: &[u8],
    alt: &[u8],
) -> String {
    let norm = normalize(seq, pos_ib, refb, alt).unwrap();
    let mut writer = DigestWriter::new();
    writer.allele_identifier_literal(
        sq,
        norm.start,
        norm.end,
        std::str::from_utf8(&norm.allele).unwrap(),
    )
}

/// Compute the canonical VRS id via the normalize+digest path, deriving the
/// `SQ.<digest>` from the given bases (the readonly/WASM form). The bridge's
/// SequenceReference accession is the raw digest of the bases.
#[allow(dead_code)]
pub fn vcf_equiv_vrs_id_from_bases(
    name: &str,
    bases: &[u8],
    pos_ib: u64,
    refb: &[u8],
    alt: &[u8],
) -> String {
    let raw = gtars_refget::digest::digest_sequence(name, bases)
        .metadata()
        .sha512t24u
        .clone();
    let sq = format!("SQ.{raw}");
    let norm = normalize(bases, pos_ib, refb, alt).unwrap();
    let mut writer = DigestWriter::new();
    writer.allele_identifier_literal(
        &sq,
        norm.start,
        norm.end,
        std::str::from_utf8(&norm.allele).unwrap(),
    )
}

// ============================================================================
// Synthetic HGVS harness helpers (filesystem-only)
// ============================================================================

#[cfg(feature = "filesystem")]
mod synthetic {
    use std::path::PathBuf;
    use std::sync::Arc;

    use gtars_refget::store::RefgetStore;
    use gtars_refget::transcripts::{TxStore, TxStoreBuilder};
    use gtars_vrs::TxProvider;

    pub const CHROM_NAME: &str = "chr_synth";

    pub fn synthetic_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/hgvs/synthetic")
    }

    pub fn fixtures_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
    }

    /// Build the synthetic RefgetStore + readonly TxStore + provider used by
    /// the layer 1-3 synthetic tests. Byte-identical setup across both former
    /// synthetic test files.
    pub fn synthetic_build_fixture() -> (
        RefgetStore,
        Arc<gtars_refget::ReadonlyTxStore>,
        TxProvider,
        String,
    ) {
        let dir = synthetic_dir();
        let fasta = dir.join("synthetic.fa");

        let mut store = RefgetStore::in_memory();
        store.disable_encoding();
        store.set_quiet(true);
        store
            .add_sequence_collection_from_fasta(
                &fasta,
                gtars_refget::store::FastaImportOptions::new(),
            )
            .expect("import synthetic.fa");
        store.load_all_sequences().expect("load sequences");

        let collection_digest = store
            .iter_collections()
            .next()
            .map(|c| c.metadata.digest.clone())
            .expect("collection");

        let coll = store.get_collection(&collection_digest).unwrap();
        let mut chrom_digest_b64 = String::new();
        for r in &coll.sequences {
            let m = r.metadata();
            if m.name == CHROM_NAME {
                chrom_digest_b64 = m.sha512t24u.clone();
            }
        }
        assert!(!chrom_digest_b64.is_empty());

        let digest_bytes = base64_url::decode(&chrom_digest_b64).expect("decode digest");
        assert_eq!(digest_bytes.len(), 24);
        let mut digest_arr = [0u8; 24];
        digest_arr.copy_from_slice(&digest_bytes);

        let cdot_path = dir.join("synthetic_transcripts.json");
        let tmpdir = tempfile::tempdir().unwrap();
        let bin_path = tmpdir.path().join("synth.reftx");

        let mut builder = TxStoreBuilder::new();
        builder.add_chrom_mapping(CHROM_NAME, digest_arr);
        let n = builder.ingest_cdot(&cdot_path).expect("ingest cdot");
        assert!(n >= 5, "expected >=5 transcripts, got {}", n);
        builder.build(&bin_path).expect("build reftx");

        let txstore = Arc::new(TxStore::open(&bin_path).unwrap().into_readonly());
        // Leak the tempdir to keep the mmap'd binary file alive for the test process.
        std::mem::forget(tmpdir);

        let provider = TxProvider::new(Arc::clone(&txstore));
        (store, txstore, provider, collection_digest)
    }

    pub fn resolve_seq_bytes(store: &mut RefgetStore, coll: &str, name: &str) -> Option<Vec<u8>> {
        let c = store.get_collection(coll).ok()?;
        let sha = c
            .sequences
            .iter()
            .find(|r| r.metadata().name == name)
            .map(|r| r.metadata().sha512t24u.clone())?;
        // Mode-agnostic decode straight from the resident record (works for Raw
        // and Encoded stores alike; no decoded cache / mmap required).
        store.load_sequence(&sha).ok()?;
        store
            .get_sequence(&sha)
            .ok()?
            .decode()
            .map(|s| s.into_bytes())
    }

    pub fn resolve_raw_digest(store: &mut RefgetStore, coll: &str, name: &str) -> Option<String> {
        let c = store.get_collection(coll).ok()?;
        c.sequences
            .iter()
            .find(|r| r.metadata().name == name)
            .map(|r| r.metadata().sha512t24u.clone())
    }
}

#[cfg(feature = "filesystem")]
#[allow(unused_imports)]
pub use synthetic::{
    CHROM_NAME, fixtures_dir, resolve_raw_digest, resolve_seq_bytes, synthetic_build_fixture,
    synthetic_dir,
};

// ============================================================================
// Parallel/VCF equivalence fixtures (filesystem-only)
// ============================================================================

#[cfg(feature = "filesystem")]
mod parallel {
    use std::collections::HashMap;
    use std::io::Write;

    use flate2::Compression;
    use flate2::write::DeflateEncoder;
    use gtars_refget::store::{FastaImportOptions, ReadonlyRefgetStore, RefgetStore};
    use tempfile::{TempDir, tempdir};

    // ── BGZF writer helpers ─────────────────────────────────────────────────

    /// Write `payload` as one BGZF block (gzip member with the BGZF `BC` extra
    /// subfield) onto `out`.
    pub fn write_bgzf_block(out: &mut Vec<u8>, payload: &[u8]) {
        // Raw DEFLATE of the payload.
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
        enc.write_all(payload).unwrap();
        let deflated = enc.finish().unwrap();

        // Total block size = 12 (fixed header) + 6 (extra: BC subfield) +
        // deflate + 8 (CRC32 + ISIZE). BSIZE stored in the extra field is
        // total - 1.
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
    pub fn write_bgzf_eof(out: &mut Vec<u8>) {
        write_bgzf_block(out, &[]);
    }

    /// Minimal CRC32 (IEEE) implementation so the test has no extra deps.
    pub fn crc32(data: &[u8]) -> u32 {
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
    /// byte payloads, so VCF lines straddle block boundaries, then a BGZF EOF
    /// block.
    pub fn write_bgzf_file(path: &std::path::Path, data: &[u8], block_payload: usize) {
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

    /// Disk-backed encoded refget store fixture plus a plain `.vcf` and a
    /// multi-block BGZF `.vcf.bgz` built from the SAME VCF text. Both the
    /// single-reader parallel suite and the BGZF-block-parallel suite open this.
    pub struct EncodedVcfFixture {
        _dir: TempDir,
        pub store_dir: String,
        pub name_to_digest: HashMap<String, String>,
        pub vcf: String,
        pub vcf_bgz: String,
    }

    impl EncodedVcfFixture {
        /// A serial / single-reader-parallel (encoded, decode-on-the-fly)
        /// readonly view of the fixture store.
        pub fn open_serial(&self) -> ReadonlyRefgetStore {
            let mut store = RefgetStore::open_local(&self.store_dir).unwrap();
            store.load_all_collections().unwrap();
            for d in self.name_to_digest.values() {
                store.load_sequence(d.as_str()).unwrap();
            }
            store.into_readonly()
        }

        /// An encoded (decode-on-the-fly) readonly view of the fixture store.
        pub fn open_encoded(&self) -> ReadonlyRefgetStore {
            self.open_serial()
        }

        /// Alias used by the BGZF suite.
        pub fn open_readonly(&self) -> ReadonlyRefgetStore {
            self.open_serial()
        }

        /// A Raw-mode (decoded, 1 byte/base) readonly view of the fixture store.
        ///
        /// `set_encoding_mode` only re-encodes/decodes sequences that are
        /// already resident, so we load every needed sequence first and then
        /// `disable_encoding` to convert the in-memory bytes to Raw (1
        /// byte/base, `bytes.len() == length`).
        pub fn open_raw(&self) -> ReadonlyRefgetStore {
            let mut store = RefgetStore::open_local(&self.store_dir).unwrap();
            store.load_all_collections().unwrap();
            for d in self.name_to_digest.values() {
                store.load_sequence(d.as_str()).unwrap();
            }
            store.disable_encoding(); // StorageMode::Raw: 1 byte/base, len == length
            store.into_readonly()
        }
    }

    /// Build the shared encoded-store fixture. The VCF record set contains SNVs,
    /// insertions, deletions, multi-allelics, a symbolic allele, and an
    /// unknown-chromosome record (exercising skip behavior). The same VCF text
    /// is written both as a plain `.vcf` and as a multi-block BGZF `.vcf.bgz`
    /// with small (250-byte) payloads so lines straddle block boundaries.
    pub fn build_encoded_vcf_fixture() -> EncodedVcfFixture {
        let dir = tempdir().unwrap();

        // chr1 has an A-repeat at positions 50..60 to exercise roll left/right.
        let chr1_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACAAAAAAAAAAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCTTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        let chr2_seq = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG";

        let fasta_path = dir.path().join("test.fa");
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            write!(f, ">chr1\n{}\n>chr2\n{}\n", chr1_seq, chr2_seq).unwrap();
        }

        // Build the VCF text in memory. A long-ish INFO field on the chr1/A>T
        // records makes some lines straddle small BGZF blocks (this only affects
        // block-straddling, not VRS output).
        let chr1_bytes = chr1_seq.as_bytes();
        let chr2_bytes = chr2_seq.as_bytes();
        let mut vcf = String::new();
        vcf.push_str("##fileformat=VCFv4.2\n");
        vcf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for i in 0..200 {
            let pos = 5 + (i % 40);
            // Derive REF from the actual reference base at this position (VCF POS
            // is 1-based, so the reference base is seq[pos-1]). `normalize_ref`
            // validates REF against the reference, so a hardcoded REF would only
            // match at some positions. Pick an ALT that differs from REF.
            let r1 = chr1_bytes[pos - 1] as char;
            let a1 = if r1 == 'A' { 'C' } else { 'A' };
            vcf.push_str(&format!(
                "chr1\t{pos}\t.\t{r1}\t{a1}\t.\tPASS\tAF=0.1;INFO_PAD=XXXXXXXXXXXXXXXXXXXX\n"
            ));
            // chr2 records stay multi-allelic; both ALTs must differ from REF.
            let r2 = chr2_bytes[pos - 1] as char;
            let (a2a, a2b) = match r2 {
                'A' => ('C', 'G'),
                'C' => ('A', 'G'),
                'G' => ('A', 'C'),
                _ => ('A', 'C'), // r2 == 'T'
            };
            vcf.push_str(&format!("chr2\t{pos}\t.\t{r2}\t{a2a},{a2b}\t.\tPASS\t.\n"));
        }
        // Insertion + deletion in the A-repeat (1-based 51 == 0-based 50).
        vcf.push_str("chr1\t51\t.\tA\tAA\t.\tPASS\t.\n");
        vcf.push_str("chr1\t51\t.\tAA\tA\t.\tPASS\t.\n");
        // SNV. REF must match chr1[9] (0-based) == 'C'.
        vcf.push_str("chr1\t10\t.\tC\tA\t.\tPASS\t.\n");
        // Symbolic allele: must be skipped by both paths.
        vcf.push_str("chr1\t20\t.\tA\t<DEL>\t.\tPASS\t.\n");
        // Unknown chromosome: must be skipped by both paths.
        vcf.push_str("chrZ\t5\t.\tA\tT\t.\tPASS\t.\n");

        // Plain-text VCF (for the serial reference and single-reader path).
        let vcf_plain = dir.path().join("test.vcf");
        std::fs::write(&vcf_plain, vcf.as_bytes()).unwrap();

        // Multi-block BGZF VCF with small (250-byte) payloads to force lines to
        // span block boundaries, exercising the head/tail stitching.
        let vcf_bgz = dir.path().join("test.vcf.bgz");
        write_bgzf_file(&vcf_bgz, vcf.as_bytes(), 250);

        // Build a disk-backed encoded store (on_disk() defaults to
        // StorageMode::Encoded).
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

        EncodedVcfFixture {
            _dir: dir,
            store_dir: store_dir.to_str().unwrap().to_string(),
            name_to_digest,
            vcf: vcf_plain.to_str().unwrap().to_string(),
            vcf_bgz: vcf_bgz.to_str().unwrap().to_string(),
        }
    }
}

#[cfg(feature = "filesystem")]
#[allow(unused_imports)]
pub use parallel::{
    EncodedVcfFixture, build_encoded_vcf_fixture, crc32, write_bgzf_block, write_bgzf_eof,
    write_bgzf_file,
};
