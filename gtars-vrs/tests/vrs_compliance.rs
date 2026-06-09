//! VRS compliance tests — GA4GH VRS 2.0 hard-coded vectors + VRS 2.1 YAML.
//!
//! Merged from the former `test_vrs_compliance` (2.0 hand-coded vectors) and
//! `test_vrs_2_1_validation` (YAML-driven against the vrs-spec repo) binaries.
//! Test-fn names are disjoint across the two suites; no renames needed.

use gtars_refget::sha512t24u;
use gtars_vrs::digest::{allele_digest, allele_identifier, sequence_location_digest};
use gtars_vrs::models::{Allele, AlleleState, SequenceLocation, SequenceReference};

// ============================================================================
// VRS 2.0 hard-coded vectors
// ============================================================================
//
// Test vectors from the official GA4GH VRS specification:
// - https://github.com/ga4gh/vrs/blob/2.0/validation/models.yaml
// - https://github.com/ga4gh/vrs-python/blob/main/tests/test_vrs.py

// ── sha512t24u primitive (also tested in gtars-refget; sanity check) ─────

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

// ── VRS 2.0 SequenceLocation digests ─────────────────────────────────────

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

// ── VRS 2.0 Allele identifiers ───────────────────────────────────────────

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

// ============================================================================
// VRS 2.1 YAML validation
// ============================================================================
//
// Loads YAML validation files from the cloned vrs-spec repo at
// `vrs-spec/validation/` and checks our implementation against them.
//
// Currently tests:
//   - sha512t24u hash function
//   - SequenceLocation digest/identify (exact coordinates only)
//   - Allele digest/identify (LiteralSequenceExpression + ReferenceLengthExpression)
//
// Types not yet implemented in gtars-vrs are reported as skipped. The tests
// degrade to a no-op `eprintln!("SKIP")` when the vrs-spec repo isn't present
// (via `require_validation_dir!`).
//
// Run with:
//   cargo test -p gtars-vrs --test vrs_compliance -- --nocapture

use std::fs;
use std::path::{Path, PathBuf};

use serde_json::Value;

/// Find the validation directory relative to the workspace root.
/// Returns None if the vrs-spec repo is not available.
fn validation_dir() -> Option<PathBuf> {
    // Check VRS_SPEC_DIR env var first, then fall back to sibling repo convention
    if let Ok(dir) = std::env::var("VRS_SPEC_DIR") {
        let p = PathBuf::from(dir).join("validation");
        if p.exists() {
            return Some(p);
        }
    }

    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo_root = manifest.parent().unwrap(); // gtars repo root
    // Try sibling to gtars (workspace repos/ layout)
    let sibling = repo_root.parent().unwrap().join("vrs-spec/validation");
    // Try inside gtars (standalone checkout)
    let inside = repo_root.join("vrs-spec/validation");

    if sibling.exists() {
        Some(sibling)
    } else if inside.exists() {
        Some(inside)
    } else {
        None
    }
}

/// Skip test if validation data is unavailable.
macro_rules! require_validation_dir {
    () => {
        match validation_dir() {
            Some(dir) => dir,
            None => {
                eprintln!("SKIP: vrs-spec validation data not available");
                return;
            }
        }
    };
}

fn load_yaml(dir: &Path, filename: &str) -> Value {
    let path = dir.join(filename);
    let content =
        fs::read_to_string(&path).unwrap_or_else(|e| panic!("Failed to read {path:?}: {e}"));
    serde_yaml::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {path:?}: {e}"))
}

/// Try to parse a SequenceLocation from the YAML input.
/// Returns None if the location uses Range coordinates (arrays) which we don't support yet.
fn parse_sequence_location(input: &Value) -> Option<SequenceLocation> {
    let refget_accession = input
        .get("sequenceReference")
        .and_then(|sr| sr.get("refgetAccession"))
        .and_then(|v| v.as_str())?;

    // Only handle simple integer coordinates, not ranges
    let start = input.get("start").and_then(|v| v.as_u64())?;
    let end = input.get("end").and_then(|v| v.as_u64())?;

    Some(SequenceLocation {
        sequence_reference: SequenceReference {
            refget_accession: refget_accession.to_string(),
        },
        start,
        end,
    })
}

/// Try to parse an Allele from the YAML input.
/// Returns None if the allele uses unsupported features or has incomplete data.
fn parse_allele(input: &Value) -> Option<Allele> {
    let location_input = input.get("location")?;
    let location = parse_sequence_location(location_input)?;

    let state_input = input.get("state")?;
    let state_type = state_input.get("type").and_then(|v| v.as_str())?;

    let state = match state_type {
        // Handle trailing comma in YAML (serde_yaml may parse "Type," as "Type,")
        t if t.starts_with("LiteralSequenceExpression") => {
            let sequence = state_input.get("sequence").and_then(|v| v.as_str())?;
            if sequence.is_empty() {
                return None;
            }
            AlleleState::LiteralSequenceExpression {
                sequence: sequence.to_string(),
            }
        }
        t if t.starts_with("ReferenceLengthExpression") => {
            let length = state_input.get("length").and_then(|v| v.as_u64())?;
            let repeat_subunit_length = state_input
                .get("repeatSubunitLength")
                .and_then(|v| v.as_u64())
                .unwrap_or(0);
            let sequence = state_input
                .get("sequence")
                .and_then(|v| v.as_str())
                .map(|s| s.to_string());
            AlleleState::ReferenceLengthExpression {
                length,
                repeat_subunit_length,
                sequence,
            }
        }
        _ => return None,
    };

    Some(Allele { location, state })
}

struct TestResult {
    passed: usize,
    skipped: usize,
    failed: usize,
}

impl TestResult {
    fn new() -> Self {
        Self {
            passed: 0,
            skipped: 0,
            failed: 0,
        }
    }
}

// ── sha512t24u function tests (functions.yaml) ───────────────────────────

#[test]
fn test_sha512t24u_from_validation() {
    let dir = require_validation_dir!();
    let data = load_yaml(&dir, "functions.yaml");
    let cases = data["sha512t24u"]
        .as_array()
        .expect("sha512t24u should be an array");

    for case in cases {
        let blob = case["in"]["blob"].as_str().unwrap();
        let expected = case["out"].as_str().unwrap();
        let actual = sha512t24u(blob);
        assert_eq!(
            actual, expected,
            "sha512t24u({:?}) = {:?}, expected {:?}",
            blob, actual, expected
        );
    }
    eprintln!("sha512t24u: {}/{} passed", cases.len(), cases.len());
}

// ── SequenceLocation tests (models.yaml) ─────────────────────────────────

#[test]
fn test_sequence_location_from_validation() {
    let dir = require_validation_dir!();
    let data = load_yaml(&dir, "models.yaml");
    let cases = data["SequenceLocation"]
        .as_array()
        .expect("SequenceLocation should be an array");

    let mut r = TestResult::new();

    for case in cases {
        let name = case
            .get("name")
            .and_then(|v| v.as_str())
            .unwrap_or("unnamed");
        let input = &case["in"];
        let output = &case["out"];

        let expected_digest = output.get("ga4gh_digest").and_then(|v| v.as_str());
        let expected_identify = output.get("ga4gh_identify").and_then(|v| v.as_str());

        match parse_sequence_location(input) {
            Some(loc) => {
                let actual_digest = sequence_location_digest(&loc);

                if let Some(exp) = expected_digest {
                    if actual_digest == exp {
                        r.passed += 1;
                        eprintln!("  PASS: {name} (digest)");
                    } else {
                        r.failed += 1;
                        eprintln!("  FAIL: {name} (digest): got {actual_digest}, expected {exp}");
                    }
                }

                if let Some(exp) = expected_identify {
                    let actual_identify = format!("ga4gh:SL.{actual_digest}");
                    if actual_identify == exp {
                        r.passed += 1;
                        eprintln!("  PASS: {name} (identify)");
                    } else {
                        r.failed += 1;
                        eprintln!(
                            "  FAIL: {name} (identify): got {actual_identify}, expected {exp}"
                        );
                    }
                }
            }
            None => {
                r.skipped += 1;
                eprintln!("  SKIP: {name} (Range coordinates not yet supported)");
            }
        }
    }

    eprintln!(
        "SequenceLocation: {} passed, {} skipped, {} failed",
        r.passed, r.skipped, r.failed
    );
    assert_eq!(r.failed, 0, "{} SequenceLocation test(s) failed", r.failed);
}

// ── Allele tests (models.yaml - canonical reference vectors) ─────────────

#[test]
fn test_alleles_from_models_yaml() {
    let dir = require_validation_dir!();
    let data = load_yaml(&dir, "models.yaml");
    let cases = data["Allele"]
        .as_array()
        .expect("Allele should be an array in models.yaml");

    let r = run_allele_cases(cases, "models.yaml");
    // models.yaml contains the canonical reference vectors - these MUST pass
    assert_eq!(
        r.failed, 0,
        "{} Allele test(s) failed from models.yaml (canonical vectors)",
        r.failed
    );
}

// ── Allele tests (alleles.yaml - extended test cases, some WIP) ──────────

// IGNORED: `alleles.yaml` is sourced from an in-progress VRS spec branch
// (2.1.0-connect_2026_#10) whose expected digests are not yet finalized — some
// vectors have known issues (start > end, optional-field serialization
// differences). It previously ran but only `eprintln!`'d failures and never
// asserted, so a real regression would have passed silently. Until the WIP
// spec stabilizes this test is `#[ignore]`d; when it does, drop the attribute
// and have `run_allele_cases` failures hard-fail via `assert_eq!(r.failed, 0)`.
#[test]
#[ignore = "alleles.yaml expected digests are from an unfinalized WIP VRS spec branch"]
fn test_alleles_from_alleles_yaml() {
    let dir = require_validation_dir!();
    let data = load_yaml(&dir, "alleles.yaml");
    let cases = data["Allele"]
        .as_array()
        .expect("Allele should be an array in alleles.yaml");

    let r = run_allele_cases(cases, "alleles.yaml");

    // When run explicitly (`--ignored`), surface the WIP failure count but do
    // not hard-fail: the expected vectors themselves are not yet final.
    if r.failed > 0 {
        eprintln!(
            "\n  NOTE: {} alleles.yaml case(s) produced different digests than expected.\n  \
             This may indicate spec changes on the WIP branch (2.1.0-connect_2026_#10)\n  \
             or differences in how optional fields (sequence, repeatSubunitLength) are\n  \
             included in serialization. Investigate individually.\n",
            r.failed
        );
    }
}

fn run_allele_cases(cases: &[Value], source: &str) -> TestResult {
    let mut r = TestResult::new();

    for case in cases {
        let name = case
            .get("name")
            .and_then(|v| v.as_str())
            .unwrap_or("unnamed");
        let input = &case["in"];
        let output = &case["out"];

        let expected_digest = output.get("ga4gh_digest").and_then(|v| v.as_str());
        let expected_identify = output.get("ga4gh_identify").and_then(|v| v.as_str());

        // Skip cases with "tbd" or missing expected values
        if expected_digest.is_none_or(|d| d == "tbd" || d.is_empty()) {
            r.skipped += 1;
            eprintln!("  SKIP: {name} (expected value is tbd or missing)");
            continue;
        }

        match parse_allele(input) {
            Some(allele) => {
                let actual_digest = allele_digest(&allele);
                let actual_identify = allele_identifier(&allele);

                if let Some(exp) = expected_digest {
                    if actual_digest == exp {
                        r.passed += 1;
                        eprintln!("  PASS: {name} (digest)");
                    } else {
                        r.failed += 1;
                        eprintln!("  FAIL: {name} (digest): got {actual_digest}, expected {exp}");
                    }
                }

                if let Some(exp) = expected_identify {
                    if actual_identify == exp {
                        r.passed += 1;
                        eprintln!("  PASS: {name} (identify)");
                    } else {
                        r.failed += 1;
                        eprintln!(
                            "  FAIL: {name} (identify): got {actual_identify}, expected {exp}"
                        );
                    }
                }
            }
            None => {
                r.skipped += 1;
                eprintln!("  SKIP: {name} (unsupported or incomplete test case)");
            }
        }
    }

    eprintln!(
        "Allele ({source}): {} passed, {} skipped, {} failed",
        r.passed, r.skipped, r.failed
    );
    r
}

// ── Serialization tests (ga4gh_serialize from models.yaml) ───────────────

#[test]
fn test_serialization_from_models_yaml() {
    let dir = require_validation_dir!();
    let data = load_yaml(&dir, "models.yaml");
    let mut r = TestResult::new();

    // Test SequenceLocation serialization
    if let Some(cases) = data["SequenceLocation"].as_array() {
        for case in cases {
            let name = case
                .get("name")
                .and_then(|v| v.as_str())
                .unwrap_or("unnamed");
            let expected_serialize = case["out"].get("ga4gh_serialize").and_then(|v| v.as_str());
            if expected_serialize.is_none() || expected_serialize == Some("") {
                r.skipped += 1;
                continue;
            }
            let expected = expected_serialize.unwrap();

            match parse_sequence_location(&case["in"]) {
                Some(loc) => {
                    let json_val = serde_json::json!({
                        "end": loc.end,
                        "sequenceReference": {
                            "refgetAccession": loc.sequence_reference.refget_accession,
                            "type": "SequenceReference"
                        },
                        "start": loc.start,
                        "type": "SequenceLocation"
                    });
                    let actual = gtars_refget::digest::algorithms::canonicalize_json(&json_val);
                    if actual == expected {
                        r.passed += 1;
                        eprintln!("  PASS: SL serialize {name}");
                    } else {
                        r.failed += 1;
                        eprintln!("  FAIL: SL serialize {name}:");
                        eprintln!("    got:    {actual}");
                        eprintln!("    expect: {expected}");
                    }
                }
                None => {
                    r.skipped += 1;
                    eprintln!("  SKIP: SL serialize {name} (Range coordinates)");
                }
            }
        }
    }

    eprintln!(
        "Serialization: {} passed, {} skipped, {} failed",
        r.passed, r.skipped, r.failed
    );
    assert_eq!(r.failed, 0, "{} serialization test(s) failed", r.failed);
}

// ── Coverage report: types not yet implemented in gtars-vrs ──────────────

#[test]
fn report_coverage() {
    let dir = require_validation_dir!();
    let data = load_yaml(&dir, "models.yaml");
    let not_implemented = [
        "Adjacency",
        "CisPhasedBlock",
        "DerivativeMolecule",
        "Terminus",
        "CopyNumberCount",
        "CopyNumberChange",
    ];

    eprintln!("\n=== VRS 2.1 Validation Coverage ===");
    eprintln!("Implemented:");
    eprintln!("  sha512t24u            - fully tested");
    eprintln!("  SequenceReference     - serialization tested");
    eprintln!("  SequenceLocation      - digest + identify (exact coords)");
    eprintln!("  Allele (LSE)          - digest + identify");
    eprintln!("  Allele (RLE)          - digest + identify");
    eprintln!();
    eprintln!("Not yet implemented:");
    for typ in &not_implemented {
        if let Some(cases) = data.get(*typ).and_then(|v| v.as_array()) {
            eprintln!("  {typ:<24} - {} test case(s) available", cases.len());
        }
    }
    eprintln!();
    eprintln!("Not yet supported:");
    eprintln!("  SequenceLocation w/ Range coordinates (3 test cases)");
    eprintln!("===================================\n");
}
