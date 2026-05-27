//! Integration tests for the HGVS parser and TranscriptProvider trait.

use gtars_vrs::hgvs::{parse, Datum, Edit, LocationRange, ReferenceType};
use gtars_vrs::{NoTranscriptProvider, ProviderError, TranscriptProvider};

#[test]
fn test_parse_braf_v600e() {
    let v = parse("NM_004333.6:c.1799T>A").unwrap();
    assert_eq!(v.accession, "NM_004333.6");
    assert!(matches!(v.reference_type, ReferenceType::C));
    match v.posedit.pos {
        LocationRange::Single(p) => {
            assert_eq!(p.base, 1799);
            assert_eq!(p.offset, 0);
            assert!(matches!(p.datum, Datum::CdsStart));
        }
        _ => panic!("expected single position"),
    }
    match v.posedit.edit {
        Edit::Sub { reference, alternate } => {
            assert_eq!(reference, "T");
            assert_eq!(alternate, "A");
        }
        _ => panic!("expected substitution"),
    }
}

#[test]
fn test_parse_intronic_variant() {
    // Common splice-site clinical variant pattern.
    let v = parse("NM_004333.6:c.1798+5G>A").unwrap();
    match v.posedit.pos {
        LocationRange::Single(p) => {
            assert_eq!(p.base, 1798);
            assert_eq!(p.offset, 5);
            assert!(matches!(p.datum, Datum::CdsStart));
        }
        _ => panic!(),
    }
}

#[test]
fn test_parse_gene_symbol_only() {
    // MANE workflow: clinician supplies gene + c. coordinate.
    let v = parse("BRAF:c.1799T>A").unwrap();
    assert_eq!(v.accession, "BRAF");
}

#[test]
fn test_no_transcript_provider_rejects_c_variants() {
    let p = NoTranscriptProvider;
    let err = p.c_to_genomic("NM_004333.6", 1799).unwrap_err();
    assert!(matches!(err, ProviderError::TranscriptNotFound(_)));
}

#[test]
fn test_no_transcript_provider_default_full_methods() {
    let p = NoTranscriptProvider;
    // c. with offset must err since default-NoProvider doesn't have a real
    // mapping.
    assert!(p.c_to_genomic_full("NM_004333.6", 1798, 5, false).is_err());
    assert!(p.gene_to_mane_accession("BRAF").is_none());
}
