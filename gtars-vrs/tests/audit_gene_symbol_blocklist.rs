//! AUDIT TEST — finding #2 (HIGH): `looks_like_gene_symbol` blocklist false positives.
//!
//! `gtars_vrs::hgvs::bridge::looks_like_gene_symbol` (bridge.rs:552) classifies an
//! accession token as a *gene symbol* (eligible for MANE resolution via
//! `gene_to_mane_accession`) iff it does NOT start with one of a fixed list of
//! blocked prefixes. That list contains the TWO-LETTER prefixes `"MT"`, `"KI"`,
//! and `"GL"`, matched with `str::starts_with`. As a result, perfectly real
//! HGNC gene symbols whose names merely *begin* with those two letters —
//! `KIT`, `MTOR`, `MTHFR`, `GLI1`, ... — are misclassified as NON-gene-symbols.
//!
//! Consequence under test: an HGVS expression like `KIT:c.1A>G` is treated as if
//! "KIT" were itself a transcript accession. The MANE lookup
//! (`gene_to_mane_accession("KIT")`) is SKIPPED, the bridge asks the provider to
//! resolve the chromosome for accession "KIT" directly, that fails, and the call
//! errors (TranscriptNotFound / UnknownChrom) instead of resolving to the MANE
//! transcript and producing a VRS id.
//!
//! This file builds the same in-memory `RefgetStore` + `MockProvider` fixture as
//! `hgvs_bridge.rs`, but wires the mock's `gene_to_mane_accession` so that the
//! gene symbols "KIT" (blocked-prefix collision) and "BRAF" (clean control) both
//! map to the forward fixture transcript `NM_FWD.1`. A correctly-behaving bridge
//! resolves BOTH via MANE to the SAME VRS id. The blocklist bug makes only
//! "BRAF" resolve; "KIT" errors.
//!
//! The asserts encode the EXPECTED-CORRECT behavior ("KIT" resolves), so the
//! test FAILS while the bug is present — a failing run CONFIRMS the bug.

use std::collections::HashMap;

use gtars_refget::store::{FastaImportOptions, RefgetStore};
use gtars_vrs::hgvs::bridge::hgvs_str_to_vrs_id;
use gtars_vrs::models::SequenceLocation;
use gtars_vrs::provider::{ProviderError, TranscriptProvider};

// ── Fixture chromosome (copied from hgvs_bridge.rs) ─────────────────────────

const CHR_F_NAME: &str = "chrF";
const CHR_F_SEQ: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // len 40

const TX_FWD: &str = "NM_FWD.1";

/// MockProvider — trimmed copy of the `hgvs_bridge.rs` fixture mock, with
/// `gene_to_mane_accession` extended so the gene symbols "KIT" and "BRAF" both
/// resolve to the forward fixture transcript `NM_FWD.1`.
struct MockProvider {
    name_to_digest: HashMap<String, String>,
}

impl MockProvider {
    fn new(name_to_digest: HashMap<String, String>) -> Self {
        Self { name_to_digest }
    }

    fn chrom_sq(&self, name: &str) -> Result<String, ProviderError> {
        let d = self
            .name_to_digest
            .get(name)
            .ok_or_else(|| ProviderError::TranscriptNotFound(name.to_string()))?;
        Ok(format!("SQ.{}", d))
    }
}

impl TranscriptProvider for MockProvider {
    fn c_to_genomic(&self, accession: &str, c_pos: i64) -> Result<SequenceLocation, ProviderError> {
        self.c_to_genomic_full(accession, c_pos, 0, false)
    }

    fn n_to_genomic(&self, accession: &str, n_pos: u64) -> Result<SequenceLocation, ProviderError> {
        self.n_to_genomic_full(accession, n_pos as i64, 0)
    }

    fn get_chrom_accession(&self, accession: &str) -> Result<String, ProviderError> {
        match accession {
            TX_FWD => self.chrom_sq(CHR_F_NAME),
            _ => Err(ProviderError::TranscriptNotFound(accession.to_string())),
        }
    }

    fn get_strand(&self, accession: &str) -> Result<i8, ProviderError> {
        match accession {
            TX_FWD => Ok(1),
            _ => Err(ProviderError::TranscriptNotFound(accession.to_string())),
        }
    }

    fn c_to_genomic_full(
        &self,
        accession: &str,
        c_pos: i64,
        offset: i64,
        _is_cds_end: bool,
    ) -> Result<SequenceLocation, ProviderError> {
        let chrom_name = match accession {
            TX_FWD => CHR_F_NAME,
            _ => return Err(ProviderError::TranscriptNotFound(accession.to_string())),
        };
        let sq = self.chrom_sq(chrom_name)?;
        // Forward: cds_start = genomic 5 (1-based) = ib 4; c.1 -> ib 4.
        let ib = (4 + (c_pos - 1) + offset) as u64;
        Ok(SequenceLocation {
            sequence_reference: gtars_vrs::models::SequenceReference {
                refget_accession: sq,
            },
            start: ib,
            end: ib + 1,
        })
    }

    fn n_to_genomic_full(
        &self,
        accession: &str,
        n_pos: i64,
        offset: i64,
    ) -> Result<SequenceLocation, ProviderError> {
        self.c_to_genomic_full(accession, n_pos, offset, false)
    }

    fn gene_to_mane_accession(&self, gene: &str) -> Option<String> {
        // Real HGNC symbols. "KIT" and "GLI1" begin with blocked two-letter
        // prefixes ("KI", "GL"); "BRAF" does not (clean control).
        match gene {
            "KIT" => Some(TX_FWD.to_string()),
            "MTOR" => Some(TX_FWD.to_string()),
            "GLI1" => Some(TX_FWD.to_string()),
            "BRAF" => Some(TX_FWD.to_string()),
            _ => None,
        }
    }
}

fn build_store_and_provider() -> (RefgetStore, MockProvider, String, tempfile::TempDir) {
    let temp = tempfile::tempdir().unwrap();
    let fasta_path = temp.path().join("fix.fa");
    let fasta_content = format!(
        ">{}\n{}\n",
        CHR_F_NAME,
        std::str::from_utf8(CHR_F_SEQ).unwrap(),
    );
    std::fs::write(&fasta_path, fasta_content).unwrap();

    let mut store = RefgetStore::in_memory();
    store.disable_encoding();
    store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
        .unwrap();
    store.load_all_sequences().unwrap();

    let collection_digest = store
        .iter_collections()
        .next()
        .map(|c| c.metadata.digest.clone())
        .expect("collection");

    let coll = store.get_collection(&collection_digest).unwrap();
    let mut name_to_digest = HashMap::new();
    for r in &coll.sequences {
        let m = r.metadata();
        name_to_digest.insert(m.name.clone(), m.sha512t24u.clone());
    }

    let provider = MockProvider::new(name_to_digest);
    (store, provider, collection_digest, temp)
}

// ── Positive control: a real symbol that does NOT collide with the blocklist ─

/// `BRAF` does not start with any blocked prefix, so it is (correctly) treated
/// as a gene symbol, resolved via MANE to NM_FWD.1, and produces a VRS id.
/// This shows the gene-symbol → MANE path itself works end-to-end in the
/// fixture, isolating the failure below to the blocklist collision.
#[ignore = "audit reproduction (PR #256) — run with --ignored"]
#[test]
fn control_braf_gene_symbol_resolves_via_mane() {
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    // c.2A>... : c.2 -> ib 5, chrF[5] = 'C'. Use C>T (correct ref).
    let res = hgvs_str_to_vrs_id("BRAF:c.2C>T", &provider, &mut store, &coll);
    assert!(
        res.is_ok(),
        "control: BRAF (non-colliding symbol) should resolve via MANE, got {res:?}"
    );
    let id = res.unwrap().value;
    assert!(id.starts_with("ga4gh:VA."), "unexpected id form: {id}");
    eprintln!("control BRAF:c.2C>T -> {id}");
}

// ── BUG A: blocked two-letter prefixes shadow real gene symbols ──────────────

/// EXPECTED-CORRECT behavior: `KIT:c.2C>T` resolves via MANE (KIT -> NM_FWD.1)
/// to the SAME VRS id as the non-colliding control `BRAF:c.2C>T`.
///
/// ACTUAL (buggy) behavior: "KIT" starts with the blocked prefix "KI", so
/// `looks_like_gene_symbol("KIT")` returns false, MANE resolution is skipped,
/// and the bridge tries to resolve "KIT" as a transcript accession — which the
/// provider does not know — yielding an error.
///
/// This assertion (KIT == BRAF id) FAILS while the bug is present.
#[ignore = "audit reproduction (PR #256) — run with --ignored"]
#[test]
fn kit_gene_symbol_must_resolve_via_mane() {
    let (mut store, provider, coll, _temp) = build_store_and_provider();

    let braf = hgvs_str_to_vrs_id("BRAF:c.2C>T", &provider, &mut store, &coll)
        .expect("control BRAF must resolve")
        .value;

    let kit = hgvs_str_to_vrs_id("KIT:c.2C>T", &provider, &mut store, &coll);
    eprintln!("KIT:c.2C>T result = {kit:?}");

    assert!(
        kit.is_ok(),
        "BUG A: 'KIT' is a real gene symbol but is blocked by the two-letter \
         'KI' prefix in looks_like_gene_symbol, so MANE resolution is skipped \
         and the call errors instead of resolving. got {kit:?}"
    );
    assert_eq!(
        kit.unwrap().value,
        braf,
        "KIT (via MANE) should produce the same VRS id as the BRAF control"
    );
}

/// Same defect via a different colliding symbol family, to show it is the
/// two-letter-prefix rule and not specific to "KIT". `MTOR` collides with "MT",
/// `GLI1` collides with "GL". Both should resolve via MANE; both are blocked.
#[ignore = "audit reproduction (PR #256) — run with --ignored"]
#[test]
fn mtor_and_gli1_gene_symbols_must_resolve_via_mane() {
    let (mut store, provider, coll, _temp) = build_store_and_provider();

    let mtor = hgvs_str_to_vrs_id("MTOR:c.2C>T", &provider, &mut store, &coll);
    let gli1 = hgvs_str_to_vrs_id("GLI1:c.2C>T", &provider, &mut store, &coll);
    eprintln!("MTOR result = {mtor:?}");
    eprintln!("GLI1 result = {gli1:?}");

    assert!(
        mtor.is_ok(),
        "BUG A: 'MTOR' blocked by two-letter 'MT' prefix. got {mtor:?}"
    );
    assert!(
        gli1.is_ok(),
        "BUG A: 'GLI1' blocked by two-letter 'GL' prefix. got {gli1:?}"
    );
}
