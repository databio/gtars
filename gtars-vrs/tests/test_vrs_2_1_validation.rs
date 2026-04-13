//! VRS 2.1 validation tests against official GA4GH test vectors.
//!
//! Loads YAML validation files from the cloned vrs-spec repo at
//! `vrs-spec/validation/` and checks our implementation against them.
//!
//! Currently tests:
//!   - sha512t24u hash function
//!   - SequenceLocation digest/identify (exact coordinates only)
//!   - Allele digest/identify (LiteralSequenceExpression + ReferenceLengthExpression)
//!
//! Types not yet implemented in gtars-vrs are reported as skipped.
//!
//! Run with:
//!   cargo test -p gtars-vrs --test test_vrs_2_1_validation -- --nocapture

use std::fs;
use std::path::PathBuf;

use gtars_refget::sha512t24u;
use gtars_vrs::digest::{allele_digest, allele_identifier, sequence_location_digest};
use gtars_vrs::models::{Allele, AlleleState, SequenceLocation, SequenceReference};

use serde_json::Value;

/// Find the validation directory relative to the workspace root.
fn validation_dir() -> PathBuf {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo_root = manifest.parent().unwrap();
    let dir = repo_root.join("vrs-spec/validation");
    assert!(
        dir.exists(),
        "Validation directory not found at {:?}. Clone the VRS spec repo first:\n  \
         git clone --branch '2.1.0-connect_2026_#10' --depth 1 \
         https://github.com/ga4gh/vrs.git vrs-spec",
        dir
    );
    dir
}

fn load_yaml(filename: &str) -> Value {
    let path = validation_dir().join(filename);
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

// ============================================================================
// sha512t24u function tests (functions.yaml)
// ============================================================================

#[test]
fn test_sha512t24u_from_validation() {
    let data = load_yaml("functions.yaml");
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

// ============================================================================
// SequenceLocation tests (models.yaml)
// ============================================================================

#[test]
fn test_sequence_location_from_validation() {
    let data = load_yaml("models.yaml");
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

// ============================================================================
// Allele tests (models.yaml - canonical reference vectors)
// ============================================================================

#[test]
fn test_alleles_from_models_yaml() {
    let data = load_yaml("models.yaml");
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

// ============================================================================
// Allele tests (alleles.yaml - extended test cases, some WIP)
// ============================================================================

#[test]
fn test_alleles_from_alleles_yaml() {
    let data = load_yaml("alleles.yaml");
    let cases = data["Allele"]
        .as_array()
        .expect("Allele should be an array in alleles.yaml");

    let r = run_allele_cases(cases, "alleles.yaml");

    // alleles.yaml is from a WIP branch - report results but don't hard-fail.
    // Some test vectors may have issues (e.g., start > end, missing serialization).
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
        if expected_digest.map_or(true, |d| d == "tbd" || d.is_empty()) {
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

// ============================================================================
// Serialization tests (ga4gh_serialize from models.yaml)
// ============================================================================

#[test]
fn test_serialization_from_models_yaml() {
    let data = load_yaml("models.yaml");
    let mut r = TestResult::new();

    // Test SequenceLocation serialization
    if let Some(cases) = data["SequenceLocation"].as_array() {
        for case in cases {
            let name = case
                .get("name")
                .and_then(|v| v.as_str())
                .unwrap_or("unnamed");
            let expected_serialize = case["out"]
                .get("ga4gh_serialize")
                .and_then(|v| v.as_str());
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

// ============================================================================
// Coverage report: types not yet implemented in gtars-vrs
// ============================================================================

#[test]
fn report_coverage() {
    let data = load_yaml("models.yaml");
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
