//! VRS 2.0 compliance tests.
//!
//! Test vectors from the official GA4GH VRS specification:
//! - https://github.com/ga4gh/vrs/blob/2.0/validation/models.yaml
//! - https://github.com/ga4gh/vrs-python/blob/main/tests/test_vrs.py

use gtars_vrs::digest::{allele_identifier, sequence_location_digest};
use gtars_vrs::models::{Allele, AlleleState, SequenceLocation, SequenceReference};

// ============================================================================
// sha512t24u primitive tests (already tested in gtars-refget, but sanity check)
// ============================================================================

#[test]
fn test_sha512t24u_empty_string() {
    let digest = gtars_refget::sha512t24u("");
    assert_eq!(digest, "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc");
}

#[test]
fn test_sha512t24u_acgt() {
    let digest = gtars_refget::sha512t24u("ACGT");
    assert_eq!(digest, "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2");
}

// ============================================================================
// VRS 2.0 SequenceLocation digest tests
// ============================================================================

#[test]
fn test_sequence_location_digest_chr19_rs7412() {
    // rs7412 on chr19 (NC_000019.10)
    // From VRS 2.0 validation/models.yaml
    let loc = SequenceLocation {
        sequence_reference: SequenceReference {
            refget_accession: "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl".to_string(),
        },
        start: 44908821,
        end: 44908822,
    };
    let digest = sequence_location_digest(&loc);
    assert_eq!(
        digest, "wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz",
        "SequenceLocation digest for rs7412 does not match VRS 2.0 spec"
    );
}

#[test]
fn test_sequence_location_digest_chr7() {
    // chr7 location (NC_000007.14)
    let loc = SequenceLocation {
        sequence_reference: SequenceReference {
            refget_accession: "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul".to_string(),
        },
        start: 44908821,
        end: 44908822,
    };
    let digest = sequence_location_digest(&loc);
    assert_eq!(
        digest, "4t6JnYWqHwYw9WzBT_lmWBb3tLQNalkT",
        "SequenceLocation digest for chr7 does not match VRS 2.0 spec"
    );
}

#[test]
fn test_sequence_location_digest_egfr() {
    // EGFR region on chr7
    // From vrs-python test_vrs.py
    let loc = SequenceLocation {
        sequence_reference: SequenceReference {
            refget_accession: "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul".to_string(),
        },
        start: 55181319,
        end: 55181320,
    };
    let digest = sequence_location_digest(&loc);
    assert_eq!(
        digest, "_G2K0qSioM74l_u3OaKR0mgLYdeTL7Xd",
        "SequenceLocation digest for EGFR region does not match VRS 2.0 spec"
    );
}

// ============================================================================
// VRS 2.0 Allele identifier tests
// ============================================================================

#[test]
fn test_allele_identifier_rs7412_snv() {
    // rs7412 C>T on chr19 — THE canonical VRS test vector
    // Expected: ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt
    let allele = Allele {
        location: SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl".to_string(),
            },
            start: 44908821,
            end: 44908822,
        },
        state: AlleleState::LiteralSequenceExpression {
            sequence: "T".to_string(),
        },
    };
    let id = allele_identifier(&allele);
    assert_eq!(
        id, "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt",
        "Allele identifier for rs7412 does not match VRS 2.0 spec"
    );
}

#[test]
fn test_allele_identifier_egfr() {
    // EGFR SNV on chr7
    // From vrs-python test_vrs.py
    let allele = Allele {
        location: SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul".to_string(),
            },
            start: 55181319,
            end: 55181320,
        },
        state: AlleleState::LiteralSequenceExpression {
            sequence: "T".to_string(),
        },
    };
    let id = allele_identifier(&allele);
    assert_eq!(
        id, "ga4gh:VA.Hy2XU_-rp4IMh6I_1NXNecBo8Qx8n0oE",
        "Allele identifier for EGFR variant does not match VRS 2.0 spec"
    );
}

#[test]
fn test_allele_identifier_clinvar_383650() {
    // ClinVar 383650
    // From vrs-python test_vrs.py
    let allele = Allele {
        location: SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: "SQ.KEO-4XBcm1cxeo_DIQ8_ofqGUkp4iZhI".to_string(),
            },
            start: 128325834,
            end: 128325835,
        },
        state: AlleleState::LiteralSequenceExpression {
            sequence: "T".to_string(),
        },
    };
    let id = allele_identifier(&allele);
    assert_eq!(
        id, "ga4gh:VA.SZIS2ua7AL-0YgUTAqyBsFPYK3vE8h_d",
        "Allele identifier for ClinVar 383650 does not match VRS 2.0 spec"
    );
}

#[test]
fn test_allele_identifier_reference_length_expression() {
    // ReferenceLengthExpression allele on chr1
    // From VRS 2.0 validation/models.yaml
    let allele = Allele {
        location: SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO".to_string(),
            },
            start: 40819438,
            end: 40819446,
        },
        state: AlleleState::ReferenceLengthExpression {
            length: 11,
            repeat_subunit_length: 3,
            sequence: None,
        },
    };
    let id = allele_identifier(&allele);
    assert_eq!(
        id, "ga4gh:VA.Oop4kjdTtKcg1kiZjIJAAR3bp7qi4aNT",
        "Allele identifier for RLE allele does not match VRS 2.0 spec"
    );
}
