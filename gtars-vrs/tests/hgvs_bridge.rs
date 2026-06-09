//! Integration tests for the HGVS-to-VRS bridge.
//!
//! These tests use a small in-memory `RefgetStore` plus a hand-rolled
//! `MockTranscriptProvider`. The contract under test: for any HGVS
//! expression with a known equivalent (chrom, pos, ref, alt) tuple,
//! `hgvs_str_to_vrs_id` produces the same VRS ID as the equivalent VCF
//! row would through the normalize+digest path.

use std::collections::HashMap;

use gtars_refget::store::{FastaImportOptions, RefgetStore};
use gtars_vrs::digest::DigestWriter;
use gtars_vrs::hgvs::bridge::{hgvs_str_to_vrs_id, hgvs_to_allele, BridgeError};
use gtars_vrs::models::SequenceLocation;
use gtars_vrs::normalize::normalize;
use gtars_vrs::provider::{ProviderError, TranscriptProvider};

// ── Fixture chromosome ──────────────────────────────────────────────────
//
// Two synthetic "chroms": a forward-strand chrom F (40 bp) and a
// reverse-strand transcript T_REV mapped onto chrom R (40 bp). Both
// chroms live in one collection. Coordinates are designed so we can
// hand-compute round-trip VRS IDs via VCF-equivalent tuples.

const CHR_F_NAME: &str = "chrF";
const CHR_F_SEQ: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // len 40

const CHR_R_NAME: &str = "chrR";
const CHR_R_SEQ: &[u8] = b"AAAAACCCCCGGGGGTTTTTACGTACGTACGTACGTACGT"; // len 40

const TX_FWD: &str = "NM_FWD.1";
const TX_REV: &str = "NM_REV.1";

/// MockTranscriptProvider — minimal implementation tailored to our
/// fixture. The forward transcript NM_FWD.1 is a 1:1 c.→g. mapping on
/// chrF starting at genomic position 5 (1-based). The reverse transcript
/// NM_REV.1 maps on chrR with its CDS start at genomic position 30
/// (1-based) on the reverse strand.
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
    fn c_to_genomic(
        &self,
        accession: &str,
        c_pos: i64,
    ) -> Result<SequenceLocation, ProviderError> {
        self.c_to_genomic_full(accession, c_pos, 0, false)
    }

    fn n_to_genomic(
        &self,
        accession: &str,
        n_pos: u64,
    ) -> Result<SequenceLocation, ProviderError> {
        self.n_to_genomic_full(accession, n_pos as i64, 0)
    }

    fn get_chrom_accession(&self, accession: &str) -> Result<String, ProviderError> {
        match accession {
            TX_FWD => self.chrom_sq(CHR_F_NAME),
            TX_REV => self.chrom_sq(CHR_R_NAME),
            _ => Err(ProviderError::TranscriptNotFound(accession.to_string())),
        }
    }

    fn get_strand(&self, accession: &str) -> Result<i8, ProviderError> {
        match accession {
            TX_FWD => Ok(1),
            TX_REV => Ok(-1),
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
            TX_REV => CHR_R_NAME,
            _ => return Err(ProviderError::TranscriptNotFound(accession.to_string())),
        };
        let sq = self.chrom_sq(chrom_name)?;
        // Forward: cds_start = genomic 5 (1-based) = ib 4
        // c.1 -> ib 4, c.2 -> ib 5, ...
        // c.1+1 -> ib 5 (intronic)
        // Reverse: cds_start = genomic 30 (1-based) = ib 29
        // c.1 -> ib 29, c.2 -> ib 28, ...
        // c.1+1 (intronic on -strand) -> ib 28
        let ib = match accession {
            TX_FWD => (4 + (c_pos - 1) + offset) as u64,
            TX_REV => (29 - (c_pos - 1) - offset) as u64,
            _ => unreachable!(),
        };
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
        // For our fixture we just delegate to c_to_genomic_full semantics.
        self.c_to_genomic_full(accession, n_pos, offset, false)
    }

    fn gene_to_mane_accession(&self, gene: &str) -> Option<String> {
        match gene {
            "GENE_FWD" => Some(TX_FWD.to_string()),
            "GENE_REV" => Some(TX_REV.to_string()),
            _ => None,
        }
    }
}

// ── Helpers ─────────────────────────────────────────────────────────────

fn build_store_and_provider() -> (RefgetStore, MockProvider, String, tempfile::TempDir) {
    let temp = tempfile::tempdir().unwrap();
    let fasta_path = temp.path().join("fix.fa");
    let fasta_content = format!(
        ">{}\n{}\n>{}\n{}\n",
        CHR_F_NAME,
        std::str::from_utf8(CHR_F_SEQ).unwrap(),
        CHR_R_NAME,
        std::str::from_utf8(CHR_R_SEQ).unwrap(),
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

/// Compute the equivalent VCF-derived VRS ID for comparison.
fn vcf_equiv_vrs_id(seq: &[u8], refget_accession: &str, pos_ib: u64, refb: &[u8], alt: &[u8]) -> String {
    let norm = normalize(seq, pos_ib, refb, alt).unwrap();
    let mut writer = DigestWriter::new();
    writer.allele_identifier_literal(
        refget_accession,
        norm.start,
        norm.end,
        std::str::from_utf8(&norm.allele).unwrap(),
    )
}

fn raw_digest(store: &mut RefgetStore, name: &str, collection_digest: &str) -> String {
    let coll = store.get_collection(collection_digest).unwrap();
    for r in &coll.sequences {
        let m = r.metadata();
        if m.name == name {
            return m.sha512t24u.clone();
        }
    }
    panic!("missing chrom {name}");
}

// ── Round-trip tests ────────────────────────────────────────────────────

#[test]
fn case_g_baseline_sub() {
    // Case 1: g. SNV on chrF at position 6 (1-based) → ib 5.
    // chrF[5] = 'C'. Sub C>T.
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let id = hgvs_str_to_vrs_id(
        &format!("{CHR_F_NAME}:g.6C>T"),
        &provider,
        &mut store,
        &coll,
    )
    .unwrap()
    .value;
    let raw = raw_digest(&mut store, CHR_F_NAME, &coll);
    let sq = format!("SQ.{raw}");
    let expected = vcf_equiv_vrs_id(CHR_F_SEQ, &sq, 5, b"C", b"T");
    assert_eq!(id, expected);
}

#[test]
fn case_c_forward_sub() {
    // c.2C>T on TX_FWD → genomic ib 5 (chrF[5]='C') → C>T
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let id = hgvs_str_to_vrs_id(
        &format!("{TX_FWD}:c.2C>T"),
        &provider,
        &mut store,
        &coll,
    )
    .unwrap()
    .value;
    let raw = raw_digest(&mut store, CHR_F_NAME, &coll);
    let sq = format!("SQ.{raw}");
    let expected = vcf_equiv_vrs_id(CHR_F_SEQ, &sq, 5, b"C", b"T");
    assert_eq!(id, expected);
}

#[test]
fn case_c_reverse_sub() {
    // c.1G>A on TX_REV → genomic ib 29.
    // chrR = "AAAAACCCCCGGGGGTTTTTACGTACGTACGTACGTACGT" (len 40)
    //         0         1         2         3
    //         0123456789012345678901234567890123456789
    // chrR[29] = 'C'. Reverse-strand: HGVS G → genomic C (revcomp). HGVS A (alt) → genomic T.
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let id = hgvs_str_to_vrs_id(
        &format!("{TX_REV}:c.1G>A"),
        &provider,
        &mut store,
        &coll,
    )
    .unwrap()
    .value;
    let raw = raw_digest(&mut store, CHR_R_NAME, &coll);
    let sq = format!("SQ.{raw}");
    let expected = vcf_equiv_vrs_id(CHR_R_SEQ, &sq, 29, b"C", b"T");
    assert_eq!(id, expected);
}

#[test]
fn case_gene_symbol_to_mane() {
    // GENE_REV → TX_REV → same as case_c_reverse_sub.
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let id =
        hgvs_str_to_vrs_id("GENE_REV:c.1G>A", &provider, &mut store, &coll)
            .unwrap()
            .value;
    let raw = raw_digest(&mut store, CHR_R_NAME, &coll);
    let sq = format!("SQ.{raw}");
    let expected = vcf_equiv_vrs_id(CHR_R_SEQ, &sq, 29, b"C", b"T");
    assert_eq!(id, expected);
}

#[test]
fn case_g_deletion_left_trim() {
    // chrF: ACGT ACGT ACGT...
    // delete 1 bp at g.6 (chrF[5]='C') → VCF equivalent: pos 5(0b)=4 (1b 5), ref="AC", alt="A".
    // Use HGVS del: NC...:g.6del
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let id = hgvs_str_to_vrs_id(
        &format!("{CHR_F_NAME}:g.6del"),
        &provider,
        &mut store,
        &coll,
    )
    .unwrap()
    .value;
    let raw = raw_digest(&mut store, CHR_F_NAME, &coll);
    let sq = format!("SQ.{raw}");
    // Bridge produces start_ib=5, ref="C", alt="" — then normalize.
    let expected = vcf_equiv_vrs_id(CHR_F_SEQ, &sq, 5, b"C", b"");
    assert_eq!(id, expected);
}

#[test]
fn case_g_insertion() {
    // Insert "TT" between chrF positions 6 and 7 (1-based).
    // → ib 6 (interbase): start=end=6, alt="TT"
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let id = hgvs_str_to_vrs_id(
        &format!("{CHR_F_NAME}:g.6_7insTT"),
        &provider,
        &mut store,
        &coll,
    )
    .unwrap()
    .value;
    let raw = raw_digest(&mut store, CHR_F_NAME, &coll);
    let sq = format!("SQ.{raw}");
    let expected = vcf_equiv_vrs_id(CHR_F_SEQ, &sq, 6, b"", b"TT");
    assert_eq!(id, expected);
}

#[test]
fn case_g_dup() {
    // dup at g.6 (chrF[5]='C') → ref="C", alt="CC".
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let id = hgvs_str_to_vrs_id(
        &format!("{CHR_F_NAME}:g.6dup"),
        &provider,
        &mut store,
        &coll,
    )
    .unwrap()
    .value;
    let raw = raw_digest(&mut store, CHR_F_NAME, &coll);
    let sq = format!("SQ.{raw}");
    let expected = vcf_equiv_vrs_id(CHR_F_SEQ, &sq, 5, b"C", b"CC");
    assert_eq!(id, expected);
}

#[test]
fn case_c_reverse_delins() {
    // c.1_2delinsAT on TX_REV.
    // c.1 → ib 29, c.2 → ib 28. Range → [28, 30) genomic. chrR[28..30] = "AC".
    // HGVS alt "AT" on -strand → genomic revcomp = "AT".
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let id = hgvs_str_to_vrs_id(
        &format!("{TX_REV}:c.1_2delinsAT"),
        &provider,
        &mut store,
        &coll,
    )
    .unwrap()
    .value;
    let raw = raw_digest(&mut store, CHR_R_NAME, &coll);
    let sq = format!("SQ.{raw}");
    let expected = vcf_equiv_vrs_id(CHR_R_SEQ, &sq, 28, b"AC", b"AT");
    assert_eq!(id, expected);
}

#[test]
fn case_c_reverse_ins() {
    // c.1_2insGCG on TX_REV.
    // c.1 → ib 29, c.2 → ib 28. Range adjacent → insertion site at hi=29.
    // -strand: HGVS alt "GCG" → genomic revcomp = "CGC".
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let id = hgvs_str_to_vrs_id(
        &format!("{TX_REV}:c.1_2insGCG"),
        &provider,
        &mut store,
        &coll,
    )
    .unwrap()
    .value;
    let raw = raw_digest(&mut store, CHR_R_NAME, &coll);
    let sq = format!("SQ.{raw}");
    let expected = vcf_equiv_vrs_id(CHR_R_SEQ, &sq, 29, b"", b"CGC");
    assert_eq!(id, expected);
}

#[test]
fn case_identity_noop() {
    // c.1= identity edge case. chrR[29]='C'. Identity → ref==alt.
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let id =
        hgvs_str_to_vrs_id(&format!("{TX_REV}:c.1="), &provider, &mut store, &coll)
            .unwrap()
            .value;
    // identity collapses on normalize to start==end, alt empty.
    let raw = raw_digest(&mut store, CHR_R_NAME, &coll);
    let sq = format!("SQ.{raw}");
    let expected = vcf_equiv_vrs_id(CHR_R_SEQ, &sq, 29, b"C", b"C");
    assert_eq!(id, expected);
}

// ── Negative tests ──────────────────────────────────────────────────────

#[test]
fn err_ref_mismatch() {
    // chrF[5]='C', not 'A'. Asking for g.6A>T must err RefMismatch.
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let err = hgvs_str_to_vrs_id(
        &format!("{CHR_F_NAME}:g.6A>T"),
        &provider,
        &mut store,
        &coll,
    )
    .unwrap_err();
    assert!(matches!(err, BridgeError::RefMismatch { .. }), "got {:?}", err);
}

#[test]
fn err_unsupported_protein() {
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    // p. always fails via UnsupportedReferenceType — but our parser may
    // not even support p. yet; use a bridge-level test against the AST.
    let result = hgvs_str_to_vrs_id(
        &format!("{TX_REV}:p.Val600Glu"),
        &provider,
        &mut store,
        &coll,
    );
    assert!(result.is_err());
    // Expected: either Parse error (parser rejects p.) or UnsupportedReferenceType.
    if let Err(e) = result {
        match e {
            BridgeError::UnsupportedReferenceType(_) | BridgeError::Parse(_) => {}
            other => panic!("unexpected error: {other:?}"),
        }
    }
}

#[test]
fn err_unsupported_inv() {
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let err = hgvs_str_to_vrs_id(
        &format!("{CHR_F_NAME}:g.6inv"),
        &provider,
        &mut store,
        &coll,
    );
    // Either parser rejects inv at this grammar position, or bridge errors UnsupportedEdit.
    match err {
        Err(BridgeError::UnsupportedEdit(_)) | Err(BridgeError::Parse(_)) => {}
        other => panic!("unexpected: {other:?}"),
    }
}

#[test]
fn err_out_of_bounds() {
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let err = hgvs_str_to_vrs_id(
        &format!("{CHR_F_NAME}:g.999A>T"),
        &provider,
        &mut store,
        &coll,
    )
    .unwrap_err();
    assert!(matches!(err, BridgeError::OutOfBounds { .. }), "got {:?}", err);
}

#[test]
fn err_unknown_gene() {
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let err =
        hgvs_str_to_vrs_id("BOGUSGENE:c.1A>T", &provider, &mut store, &coll)
            .unwrap_err();
    assert!(
        matches!(
            err,
            BridgeError::Provider(ProviderError::NoManeTranscript(_))
        ),
        "got {:?}",
        err
    );
}

// ── Allele construction (non-normalized) ────────────────────────────────

#[test]
fn allele_construction_only() {
    let (mut store, provider, coll, _temp) = build_store_and_provider();
    let s = format!("{CHR_F_NAME}:g.6C>T");
    let v = gtars_vrs::hgvs::parse(&s).unwrap();
    let allele = hgvs_to_allele(&v, &provider, &mut store, &coll).unwrap().value;
    assert_eq!(allele.location.start, 5);
    assert_eq!(allele.location.end, 6);
    match allele.state {
        gtars_vrs::models::AlleleState::LiteralSequenceExpression { sequence } => {
            assert_eq!(sequence, "T");
        }
        _ => panic!("expected literal sequence"),
    }
}
