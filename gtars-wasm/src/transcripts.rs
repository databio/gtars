//! WASM bindings for a persistent, reusable transcript store.
//!
//! Holds a transcript annotation store **once** so a whole batch of `c.`/`n.`
//! HGVS variants (or a transcript-coordinate VCF) can reuse it, mirroring the
//! [`RefgetStore`](crate::refget::RefgetStore) holder pattern. JS fetches the
//! `.reftx` store bytes asynchronously (there is no filesystem in the browser),
//! then hands the resident bytes in synchronously — exactly the contract used by
//! `hgvs.rs` for bases and `RefgetStore::add_encoded_sequence` for the genome
//! blob.
//!
//! The store is built via the WASM-safe in-memory backend
//! ([`gtars_refget::transcripts::ReadonlyTxStore::from_bytes`], which wraps the
//! bytes in `TxBytes::InMemory` — no `memmap2`, no `std::fs`) and wrapped in the
//! gtars-vrs [`TxProvider`], the same `TranscriptProvider` used on native. The
//! provider lives for the lifetime of the holder.
//!
//! # Example (JavaScript)
//! ```javascript
//! import init, { TranscriptStore, hgvs_to_vrs_id_with_transcripts } from "gtars-js";
//! await init();
//!
//! // JS owns the async fetch of the .reftx blob.
//! const reftxBytes = new Uint8Array(await (await fetch(reftxUrl)).arrayBuffer());
//! const tx = new TranscriptStore(reftxBytes);
//!
//! // The referenced chromosome's bases must be resident in the RefgetStore /
//! // passed to the standalone entry for the c./n. variant to resolve.
//! const id = hgvs_to_vrs_id_with_transcripts(
//!     "NM_004333.6:c.1799T>A", "chr7", chr7Bases, tx);
//! // -> "ga4gh:VA.<digest>"
//! ```

use std::sync::Arc;

use gtars_refget::transcripts::ReadonlyTxStore;
use gtars_vrs::provider::{TranscriptProvider, TxProvider};
use wasm_bindgen::prelude::*;

/// A persistent, reusable transcript store for the browser.
///
/// JS constructs it once from `.reftx` bytes and reuses it across many HGVS /
/// VCF calls. Internally it owns a [`TxProvider`] over a WASM-safe in-memory
/// transcript store, the same concrete `TranscriptProvider` used on native — so
/// `c.`/`n.` (and gene-symbol → MANE) resolution in the browser is byte-for-byte
/// identical to the native path.
#[wasm_bindgen]
pub struct TranscriptStore {
    provider: TxProvider,
}

#[wasm_bindgen]
impl TranscriptStore {
    /// Build a transcript store from owned `.reftx` bytes.
    ///
    /// JS owns the async fetch and hands resident bytes in synchronously. The
    /// bytes are copied into an `Arc<Vec<u8>>` (the WASM-safe `TxBytes::InMemory`
    /// backend) and the `.reftx` MAGIC/version/header is validated up front.
    ///
    /// # Arguments
    /// * `reftx_bytes` - The full contents of a `.reftx` transcript store, as a
    ///   Uint8Array (already decompressed; the caller inflates any `.gz` in JS).
    ///
    /// # Errors
    /// Returns a `JsError` with a descriptive message if the bytes are not a
    /// valid `.reftx` image (bad magic, unsupported version, truncated header).
    #[wasm_bindgen(constructor)]
    pub fn new(reftx_bytes: &[u8]) -> Result<TranscriptStore, JsError> {
        let store = ReadonlyTxStore::from_bytes(reftx_bytes.to_vec())
            .map_err(|e| JsError::new(&format!("failed to load transcript store: {e}")))?;
        Ok(TranscriptStore {
            provider: TxProvider::new(Arc::new(store)),
        })
    }

    /// Number of transcripts in the store (for smoke testing / introspection).
    #[wasm_bindgen(js_name = "transcriptCount")]
    pub fn transcript_count(&self) -> usize {
        self.provider.store().len() as usize
    }
}

impl TranscriptStore {
    /// The transcript provider, for the hgvs/vcf bindings. Parallels
    /// [`RefgetStore::inner_readonly`](crate::refget::RefgetStore) /
    /// `name_to_digest`. Returned as `&dyn TranscriptProvider` so the readonly
    /// bridge can take it directly.
    pub(crate) fn provider(&self) -> &dyn TranscriptProvider {
        &self.provider
    }
}

#[cfg(test)]
#[cfg(target_arch = "wasm32")]
mod tests {
    use super::*;
    use gtars_refget::digest::digest_sequence;
    use gtars_refget::transcripts::{Exon, ManeStatus, Strand, Transcript};
    use wasm_bindgen_test::*;

    use crate::hgvs::hgvs_to_vrs_id_with_transcripts;

    // ── Synthetic .reftx builder (in-memory, no fs) ──────────────────────────
    //
    // The `.reftx` format constants (`MAGIC`, `VERSION`, `HEADER_SIZE`,
    // `NONE_SENTINEL`, the FNV-1a accession hash) are `pub(crate)` to
    // gtars-refget, so the wasm test replicates the minimal serializer here —
    // mirroring gtars-refget's own `build_reftx_bytes` core-test helper. This
    // keeps the test free of the `filesystem`-gated builder while still
    // exercising the real `ReadonlyTxStore::from_bytes` parse path.

    const MAGIC: &[u8; 4] = b"RFTX";
    const VERSION: u32 = 2;
    const HEADER_SIZE: usize = 40;
    const NONE_SENTINEL: u32 = 0xFFFF_FFFF;

    fn fnv1a_64(data: &[u8]) -> u64 {
        const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
        const FNV_PRIME: u64 = 0x100_0000_01b3;
        let mut hash = FNV_OFFSET;
        for &byte in data {
            hash ^= byte as u64;
            hash = hash.wrapping_mul(FNV_PRIME);
        }
        hash
    }

    fn serialize_record(tx: &Transcript) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.push(tx.accession.len() as u8);
        buf.extend_from_slice(tx.accession.as_bytes());
        buf.push(tx.gene.len() as u8);
        buf.extend_from_slice(tx.gene.as_bytes());
        buf.extend_from_slice(&tx.chrom_digest);
        buf.push(tx.strand.to_byte());
        buf.push(tx.mane.to_flags_byte());
        buf.extend_from_slice(&tx.cds_start.unwrap_or(NONE_SENTINEL).to_le_bytes());
        buf.extend_from_slice(&tx.cds_end.unwrap_or(NONE_SENTINEL).to_le_bytes());
        buf.extend_from_slice(&(tx.exons.len() as u16).to_le_bytes());
        for exon in &tx.exons {
            buf.extend_from_slice(&exon.start.to_le_bytes());
            buf.extend_from_slice(&exon.end.to_le_bytes());
        }
        buf
    }

    fn build_reftx_bytes(transcripts: &[Transcript]) -> Vec<u8> {
        let mut sorted: Vec<&Transcript> = transcripts.iter().collect();
        sorted.sort_by_key(|t| fnv1a_64(t.accession.as_bytes()));

        let mut out = vec![0u8; HEADER_SIZE];
        let mut index_entries: Vec<(u64, u64)> = Vec::new();
        let mut mane_entries: Vec<(u64, u64)> = Vec::new();
        let mut current_offset = HEADER_SIZE as u64;

        for tx in &sorted {
            let hash = fnv1a_64(tx.accession.as_bytes());
            index_entries.push((hash, current_offset));
            if tx.mane.mane_select {
                let gene_hash = fnv1a_64(tx.gene.to_uppercase().as_bytes());
                mane_entries.push((gene_hash, current_offset));
            }
            let rec = serialize_record(tx);
            current_offset += rec.len() as u64;
            out.extend_from_slice(&rec);
        }

        let index_offset = current_offset;
        for (hash, offset) in &index_entries {
            out.extend_from_slice(&hash.to_le_bytes());
            out.extend_from_slice(&offset.to_le_bytes());
        }

        let mane_index_offset = if mane_entries.is_empty() {
            0u64
        } else {
            mane_entries.sort_by_key(|(h, _)| *h);
            let off = out.len() as u64;
            out.extend_from_slice(&(mane_entries.len() as u64).to_le_bytes());
            for (hash, offset) in &mane_entries {
                out.extend_from_slice(&hash.to_le_bytes());
                out.extend_from_slice(&offset.to_le_bytes());
            }
            off
        };

        out[0..4].copy_from_slice(MAGIC);
        out[4..8].copy_from_slice(&VERSION.to_le_bytes());
        out[8..16].copy_from_slice(&(sorted.len() as u64).to_le_bytes());
        out[16..24].copy_from_slice(&index_offset.to_le_bytes());
        out[24..32].copy_from_slice(&mane_index_offset.to_le_bytes());
        out
    }

    // Forward-strand synthetic contig; "chr"-prefixed so the bridge treats the
    // name as an accession. 80 bases of ACGT repeats.
    const CHR_F_NAME: &str = "chrF";
    const CHR_F_BASES: &str =
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    /// A single forward-strand coding transcript on `chrF`. CDS starts at
    /// genomic interbase 10, single exon spanning the whole contig, so c.N maps
    /// to genomic interbase `10 + (N - 1)`.
    fn chr_f_transcript() -> Transcript {
        // chrom_digest is the 24 raw bytes of chrF's sha512t24u digest.
        let meta = digest_sequence(CHR_F_NAME, CHR_F_BASES.as_bytes());
        let digest_b64 = meta.metadata().sha512t24u.clone();
        let raw = base64_url::decode(&digest_b64).expect("decode chrF digest");
        let mut chrom_digest = [0u8; 24];
        chrom_digest.copy_from_slice(&raw[..24]);

        Transcript {
            accession: "NM_TESTF.1".to_string(),
            gene: "TESTF".to_string(),
            chrom_digest,
            strand: Strand::Forward,
            cds_start: Some(10),
            cds_end: Some(70),
            exons: vec![Exon { start: 0, end: 80 }],
            mane: ManeStatus {
                mane_select: true,
                mane_clinical: false,
            },
        }
    }

    #[wasm_bindgen_test]
    fn transcript_store_loads_from_bytes() {
        let bytes = build_reftx_bytes(&[chr_f_transcript()]);
        let tx = TranscriptStore::new(&bytes).expect("load .reftx bytes");
        assert_eq!(tx.transcript_count(), 1);
    }

    #[wasm_bindgen_test]
    fn transcript_store_rejects_bad_magic() {
        let mut bytes = build_reftx_bytes(&[chr_f_transcript()]);
        bytes[0] = b'X';
        assert!(
            TranscriptStore::new(&bytes).is_err(),
            "corrupt magic must be rejected"
        );
    }

    // Cross-surface parity anchor for the c. path. c.6 on a forward transcript
    // with cds_start interbase 10 → genomic interbase 15 (1-based pos 16). chrF
    // is ACGT-repeating, so base at 0-based 15 = 'T'. Sub T>A.
    //
    // This MUST equal the native bridge's result for the same .reftx + bases +
    // HGVS (see gtars-vrs bridge tests). It pins the wasm provider path to the
    // native TxProvider path; divergence fails exactly one side.
    const GOLDEN_CHRF_C6TA_VRS_ID: &str = "ga4gh:VA.cnTVylp8roGCzVzSN-SPsgkWZoGEaxHF";

    #[wasm_bindgen_test]
    fn c_variant_resolves_via_transcript_store() {
        let bytes = build_reftx_bytes(&[chr_f_transcript()]);
        let tx = TranscriptStore::new(&bytes).expect("load .reftx bytes");

        // chrF base at genomic interbase 15 (c.6) is 'T'; substitute T>A.
        let id = hgvs_to_vrs_id_with_transcripts(
            "NM_TESTF.1:c.6T>A",
            CHR_F_NAME,
            CHR_F_BASES,
            &tx,
        )
        .expect("c. variant should resolve");
        assert!(
            id.starts_with("ga4gh:VA."),
            "must be a well-formed VRS allele id, got {id}"
        );
        assert_eq!(id, GOLDEN_CHRF_C6TA_VRS_ID);
    }
}
