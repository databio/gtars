//! External corpus tests from biocommons/hgvs and ferro-hgvs.
//!
//! These tests validate our HGVS parser against test cases extracted from
//! established HGVS implementations.
//!
//! Test fixtures are sourced from:
//! - biocommons/hgvs (Apache-2.0): https://github.com/biocommons/hgvs
//! - ferro-hgvs (MIT): https://github.com/varfish-org/hgvs-rs

use gtars_vrs::hgvs::parse;
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
struct TestCase {
    input: String,
    valid: bool,
    source: String,
    #[serde(default)]
    description: String,
}

fn load_corpus() -> Vec<TestCase> {
    let fixtures_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures");

    let mut all_cases = Vec::new();

    for filename in &["biocommons.json", "ferro_hgvs.json"] {
        let path = format!("{}/{}", fixtures_dir, filename);
        let content = fs::read_to_string(&path)
            .unwrap_or_else(|_| panic!("Failed to read {}", filename));
        let cases: Vec<TestCase> = serde_json::from_str(&content)
            .unwrap_or_else(|_| panic!("Failed to parse {}", filename));
        all_cases.extend(cases);
    }

    all_cases
}

#[test]
fn test_external_corpus_valid_cases() {
    let corpus = load_corpus();
    let valid_cases: Vec<_> = corpus.iter().filter(|c| c.valid).collect();

    let mut failures = Vec::new();

    for case in &valid_cases {
        let result = parse(&case.input);
        if result.is_err() {
            failures.push(format!(
                "SHOULD PARSE: '{}' (source: {}) - error: {:?}",
                case.input,
                case.source,
                result.err()
            ));
        }
    }

    if !failures.is_empty() {
        eprintln!("\n=== PARSER FAILURES (should parse but didn't) ===");
        for f in &failures {
            eprintln!("  {}", f);
        }
        eprintln!(
            "\nTotal: {} failures out of {} valid cases",
            failures.len(),
            valid_cases.len()
        );
        panic!(
            "{} valid cases failed to parse (see stderr for details)",
            failures.len()
        );
    }

    eprintln!("\nAll {} valid cases parsed successfully", valid_cases.len());
}

#[test]
fn test_external_corpus_invalid_cases() {
    let corpus = load_corpus();
    let invalid_cases: Vec<_> = corpus.iter().filter(|c| !c.valid).collect();

    let mut failures = Vec::new();

    for case in &invalid_cases {
        let result = parse(&case.input);
        if result.is_ok() {
            failures.push(format!(
                "SHOULD REJECT: '{}' (source: {})",
                case.input, case.source
            ));
        }
    }

    if !failures.is_empty() {
        eprintln!("\n=== PARSER FAILURES (should reject but parsed) ===");
        for f in &failures {
            eprintln!("  {}", f);
        }
        eprintln!(
            "\nTotal: {} failures out of {} invalid cases",
            failures.len(),
            invalid_cases.len()
        );
        panic!(
            "{} invalid cases were incorrectly accepted (see stderr for details)",
            failures.len()
        );
    }

    eprintln!("\nAll {} invalid cases correctly rejected", invalid_cases.len());
}

/// Run corpus and report detailed statistics.
#[test]
fn test_external_corpus_stats() {
    let corpus = load_corpus();

    let mut by_source: std::collections::HashMap<String, (usize, usize)> = std::collections::HashMap::new();

    for case in &corpus {
        let entry = by_source.entry(case.source.clone()).or_insert((0, 0));
        let result = parse(&case.input);
        let passed = if case.valid {
            result.is_ok()
        } else {
            result.is_err()
        };
        entry.0 += 1; // total
        if passed {
            entry.1 += 1; // passed
        }
    }

    eprintln!("\n=== CORPUS STATISTICS ===");
    let mut total = 0;
    let mut total_passed = 0;
    for (source, (count, passed)) in by_source.iter() {
        let pct = (*passed as f64 / *count as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", source, passed, count, pct);
        total += count;
        total_passed += passed;
    }
    let pct = (total_passed as f64 / total as f64) * 100.0;
    eprintln!("  TOTAL: {}/{} ({:.1}%)", total_passed, total, pct);
}
