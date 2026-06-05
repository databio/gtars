//! Cross-validation for the readonly HGVS->VRS path used by the WASM binding.
//!
//! The WASM entry point (`gtars-wasm::hgvs_to_vrs_id`) builds an in-memory
//! `ReadonlyRefgetStore` by inserting JS-provided bases as a `Full` record via
//! `digest_sequence`, builds a name->digest map, and calls
//! `hgvs_str_to_vrs_id_readonly`. These tests exercise that exact construction
//! in pure Rust (no WASM runtime) and assert the result matches the canonical
//! normalize+digest VRS id — the same value the native VCF/mutable path
//! produces (see `hgvs_bridge::case_g_baseline_sub`, which compares the mutable
//! path against this same `vcf_equiv_vrs_id` for `chrF:g.6C>T`).

use std::collections::HashMap;

use gtars_refget::digest::digest_sequence;
use gtars_refget::store::RefgetStore;
use gtars_vrs::digest::DigestWriter;
use gtars_vrs::hgvs::bridge::hgvs_str_to_vrs_id_readonly;
use gtars_vrs::normalize::normalize;
use gtars_vrs::provider::NoTranscriptProvider;

// "chr"-prefixed, so the bridge treats it as an accession (not a gene symbol).
const CHR_F_NAME: &str = "chrF";
const CHR_F_SEQ: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // len 40

/// Replicate the WASM `hgvs_to_vrs_id` route-(a) store construction and call
/// the readonly bridge — exactly what the browser binding does internally.
fn readonly_route_a(hgvs: &str, name: &str, bases: &[u8]) -> String {
    let record = digest_sequence(name, bases);
    let raw_digest = record.metadata().sha512t24u.clone();

    let mut store = RefgetStore::in_memory();
    store.add_sequence_record(record, true).unwrap();

    let mut name_to_digest = HashMap::new();
    name_to_digest.insert(name.to_string(), raw_digest);

    hgvs_str_to_vrs_id_readonly(hgvs, &NoTranscriptProvider, &store, &name_to_digest).unwrap()
}

/// Canonical VRS id via the normalize+digest path (VCF-equivalent), mirroring
/// `hgvs_bridge::vcf_equiv_vrs_id`.
fn vcf_equiv_vrs_id(seq: &[u8], name: &str, pos_ib: u64, refb: &[u8], alt: &[u8]) -> String {
    // The bridge's SequenceReference accession is the raw digest of the bases.
    let raw = digest_sequence(name, seq).metadata().sha512t24u.clone();
    let sq = format!("SQ.{raw}");
    let norm = normalize(seq, pos_ib, refb, alt).unwrap();
    let mut writer = DigestWriter::new();
    writer.allele_identifier_literal(
        &sq,
        norm.start,
        norm.end,
        std::str::from_utf8(&norm.allele).unwrap(),
    )
}

#[test]
fn readonly_g_substitution_matches_canonical() {
    // chrF[5] (1-based pos 6) = 'C'. Sub C>T at interbase 5.
    let id = readonly_route_a(&format!("{CHR_F_NAME}:g.6C>T"), CHR_F_NAME, CHR_F_SEQ);
    let expected = vcf_equiv_vrs_id(CHR_F_SEQ, CHR_F_NAME, 5, b"C", b"T");
    assert_eq!(id, expected);
    assert!(id.starts_with("ga4gh:VA."), "unexpected id form: {id}");
}

#[test]
fn readonly_g_deletion_matches_canonical() {
    // g.6del: delete chrF[5] ('C'). Equivalent VCF tuple: ref "C", alt "".
    let id = readonly_route_a(&format!("{CHR_F_NAME}:g.6del"), CHR_F_NAME, CHR_F_SEQ);
    let expected = vcf_equiv_vrs_id(CHR_F_SEQ, CHR_F_NAME, 5, b"C", b"");
    assert_eq!(id, expected);
}

#[test]
fn readonly_g_insertion_matches_canonical() {
    // g.6_7insTT: insert TT between chrF[5] and chrF[6]. Interbase 6, ref "", alt "TT".
    let id = readonly_route_a(&format!("{CHR_F_NAME}:g.6_7insTT"), CHR_F_NAME, CHR_F_SEQ);
    let expected = vcf_equiv_vrs_id(CHR_F_SEQ, CHR_F_NAME, 6, b"", b"TT");
    assert_eq!(id, expected);
}

#[test]
fn readonly_rejects_transcript_variant() {
    // c./n. are out of scope for v1 (NoTranscriptProvider). Must error, not panic.
    let record = digest_sequence(CHR_F_NAME, CHR_F_SEQ);
    let raw_digest = record.metadata().sha512t24u.clone();
    let mut store = RefgetStore::in_memory();
    store.add_sequence_record(record, true).unwrap();
    let mut name_to_digest = HashMap::new();
    name_to_digest.insert(CHR_F_NAME.to_string(), raw_digest);

    let res = hgvs_str_to_vrs_id_readonly(
        "NM_000546.6:c.215C>G",
        &NoTranscriptProvider,
        &store,
        &name_to_digest,
    );
    assert!(res.is_err(), "c. variant should be rejected in v1");
}
