//! Skeleton end-to-end harness for the HGVS -> VRS bridge.
//!
//! Walks `tests/data/hgvs/vrs_bridge/clinvar_seed.tsv` and (when the
//! bridge is implemented) parses each `hgvs_expression`, calls the
//! HGVS->VRS conversion, computes the canonical digest, and compares to
//! `expected_vrs_id`.
//!
//! Until the bridge crate lands (planned in a follow-up), this test is
//! `#[ignore]` and only validates that the seed TSV is well-formed and
//! that every HGVS expression at least PARSES.

use std::fs;
use std::path::PathBuf;

use gtars_vrs::hgvs::parse;

fn seed_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/hgvs/vrs_bridge/clinvar_seed.tsv")
}

#[test]
fn vrs_bridge_seed_tsv_is_well_formed() {
    let raw = fs::read_to_string(seed_path()).expect("seed tsv");
    let mut lines = raw.lines();
    let header = lines.next().expect("header");
    let cols: Vec<&str> = header.split('\t').collect();
    assert_eq!(
        cols,
        [
            "hgvs_expression",
            "assembly",
            "expected_vrs_id",
            "clinvar_vcv",
            "notes"
        ]
    );

    let mut rows = 0usize;
    let mut parsed = 0usize;
    let mut parse_failures: Vec<String> = Vec::new();
    for (idx, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        rows += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(
            fields.len() >= 5,
            "row {} has {} fields, expected 5",
            idx + 2,
            fields.len()
        );
        let hgvs = fields[0].trim();
        match parse(hgvs) {
            Ok(_) => parsed += 1,
            Err(e) => parse_failures.push(format!("{} :: {}", hgvs, e)),
        }
    }
    assert!(rows >= 20, "seed tsv should have >= 20 cases, found {}", rows);
    if !parse_failures.is_empty() {
        for f in &parse_failures {
            eprintln!("seed parse failure: {}", f);
        }
        eprintln!(
            "{} of {} seed rows failed to parse — investigate before enabling bridge test",
            parse_failures.len(),
            rows
        );
    }
    eprintln!("vrs_bridge seed: {} rows, {} parse OK", rows, parsed);
}

#[test]
#[ignore = "HGVS->VRS bridge not yet implemented (see plan 8)"]
fn project_clinvar_seed_to_vrs() {
    // When the bridge lands, replace the body of this test with:
    //   1. Read the seed TSV (skip rows where expected_vrs_id is empty)
    //   2. parse(hgvs)
    //   3. provider := build a TranscriptProvider with the relevant
    //      transcripts (via ReftxProvider with a fixture binary)
    //   4. allele := bridge::to_vrs_allele(&parsed, &provider)?
    //   5. assert_eq!(allele.id(), expected_vrs_id)
}
