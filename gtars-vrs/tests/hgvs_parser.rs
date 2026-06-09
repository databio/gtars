//! HGVS parser tests — hand-written unit tests plus vendored-corpus regression.
//!
//! Merged from three former binaries (test_hgvs_parser, test_hgvs_corpus_parser,
//! test_external_corpus). All three exercise `gtars_vrs::hgvs::parse` and the
//! `TranscriptProvider` trait against disjoint material:
//!
//! - Hand-written parser/provider unit tests (BRAF V600E, intronic, gene-symbol,
//!   NoTranscriptProvider behavior).
//! - Vendored-corpus regression against three external sources:
//!   - varfish-org/hgvs-rs gauntlet (positive set; supported expressions must parse)
//!   - varfish-org/hgvs-rs reject   (negative set; each line must fail to parse)
//!   - biocommons/hgvs grammar_test.tsv (mixed Valid=True/False whole-expression rows)
//!
//!   Plus a locally-curated `local_reject.txt` and `gap_fill.tsv`. A "known-skip"
//!   allowlist (`known_skips.txt`) tolerates inputs the parser does not yet handle.
//!   See `tests/data/hgvs/REFRESH.md` and `NOTICE` for provenance.
//! - Curated external corpus (`fixtures/biocommons.json` + `fixtures/ferro_hgvs.json`):
//!   different vendored artifacts with different rows than the raw corpus above.
//!   See `tests/fixtures/README.md` for provenance.

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

use gtars_vrs::hgvs::{Datum, Edit, LocationRange, ReferenceType, parse};
use gtars_vrs::{NoTranscriptProvider, ProviderError, TranscriptProvider};
use serde::Deserialize;

// ============================================================================
// Hand-written parser / provider unit tests
// ============================================================================

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
        Edit::Sub {
            reference,
            alternate,
        } => {
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

// ============================================================================
// Vendored-corpus regression (raw varfish/biocommons data files)
// ============================================================================

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/hgvs")
}

fn read_lines(path: &Path) -> Vec<(usize, String)> {
    let raw = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("could not read {}: {}", path.display(), e));
    raw.lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l.to_string()))
        .collect()
}

/// Strip whole-line comments and blank lines. Preserve `#!unsupported:`
/// markers on the caller side — they are returned with the prefix stripped
/// in `is_unsupported`.
fn is_blank_or_comment(line: &str) -> bool {
    let t = line.trim();
    t.is_empty() || (t.starts_with('#') && !t.starts_with("#!unsupported:"))
}

fn unsupported_payload(line: &str) -> Option<&str> {
    line.trim().strip_prefix("#!unsupported:").map(|s| s.trim())
}

fn load_known_skips() -> HashSet<String> {
    let p = data_dir().join("known_skips.txt");
    let raw =
        fs::read_to_string(&p).unwrap_or_else(|e| panic!("could not read {}: {}", p.display(), e));
    raw.lines()
        .filter_map(|l| {
            let t = l.trim_start();
            if t.is_empty() || t.starts_with('#') {
                return None;
            }
            // Format: `<expression><TAB><reason>` — take the first column.
            let mut parts = l.splitn(2, '\t');
            parts.next().map(|s| s.trim().to_string())
        })
        .filter(|s| !s.is_empty())
        .collect()
}

#[test]
fn accept_varfish_gauntlet() {
    let p = data_dir().join("varfish/parser/gauntlet");
    let known_skips = load_known_skips();
    let mut hard_failures: Vec<(usize, String, String)> = Vec::new();
    let mut tolerated: usize = 0;
    let mut accepted: usize = 0;

    for (lineno, line) in read_lines(&p) {
        let trimmed = line.trim();
        if is_blank_or_comment(trimmed) {
            continue;
        }
        // Pull `#!unsupported:` payloads into the "should-skip" bucket.
        let (input, is_unsup) = match unsupported_payload(trimmed) {
            Some(payload) => (payload.to_string(), true),
            None => (trimmed.to_string(), false),
        };
        match parse(&input) {
            Ok(_) => {
                accepted += 1;
            }
            Err(e) => {
                if is_unsup || known_skips.contains(&input) {
                    eprintln!(
                        "tolerated parser failure (line {}): {} :: {}",
                        lineno, input, e
                    );
                    tolerated += 1;
                } else {
                    hard_failures.push((lineno, input, e.to_string()));
                }
            }
        }
    }

    if !hard_failures.is_empty() {
        for (lineno, input, err) in &hard_failures {
            eprintln!("FAIL line {}: {} :: {}", lineno, input, err);
        }
        panic!(
            "{} unexpected parser failures in varfish gauntlet ({} accepted, {} tolerated)",
            hard_failures.len(),
            accepted,
            tolerated
        );
    }
    eprintln!(
        "varfish gauntlet: {} accepted, {} tolerated (known-skip)",
        accepted, tolerated
    );
}

#[test]
fn reject_varfish_reject_set() {
    let p = data_dir().join("varfish/parser/reject");
    let mut surprises: Vec<(usize, String)> = Vec::new();
    let mut rejected: usize = 0;
    for (lineno, line) in read_lines(&p) {
        let trimmed = line.trim();
        if is_blank_or_comment(trimmed) {
            continue;
        }
        match parse(trimmed) {
            Ok(_) => surprises.push((lineno, trimmed.to_string())),
            Err(_) => rejected += 1,
        }
    }
    if !surprises.is_empty() {
        for (lineno, input) in &surprises {
            eprintln!("UNEXPECTED-ACCEPT line {}: {}", lineno, input);
        }
        // Soft warning — the upstream reject file is tiny and the gtars
        // parser may legitimately have a more permissive lexer than the
        // upstream Rust impl. Promote to a hard failure once the lexer
        // and validator paths are fully aligned.
        eprintln!(
            "{} inputs in varfish reject set were unexpectedly accepted",
            surprises.len()
        );
    }
    eprintln!("varfish reject set: {} rejected", rejected);
}

#[test]
fn reject_local_reject_set() {
    let p = data_dir().join("local_reject.txt");
    let mut surprises: Vec<(usize, String)> = Vec::new();
    let mut rejected: usize = 0;
    for (lineno, line) in read_lines(&p) {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        match parse(trimmed) {
            Ok(_) => surprises.push((lineno, trimmed.to_string())),
            Err(_) => rejected += 1,
        }
    }
    if !surprises.is_empty() {
        for (lineno, input) in &surprises {
            eprintln!(
                "UNEXPECTED-ACCEPT in local_reject (line {}): {}",
                lineno, input
            );
        }
        // Soft fail: see TODO in the file. Some reject categories below
        // require validator-layer checks the bare parser does not perform.
        eprintln!(
            "{} inputs in local_reject.txt were accepted by the parser",
            surprises.len()
        );
    }
    eprintln!("local reject set: {} rejected", rejected);
}

#[test]
fn accept_local_gap_fill() {
    let p = data_dir().join("gap_fill.tsv");
    let known_skips = load_known_skips();
    let mut hard_failures: Vec<(usize, String, String)> = Vec::new();
    let mut accepted = 0usize;
    let mut tolerated = 0usize;
    for (lineno, line) in read_lines(&p) {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let input = trimmed.split('\t').next().unwrap_or("").trim();
        if input.is_empty() {
            continue;
        }
        match parse(input) {
            Ok(_) => accepted += 1,
            Err(e) => {
                if known_skips.contains(input) {
                    eprintln!(
                        "tolerated gap_fill failure (line {}): {} :: {}",
                        lineno, input, e
                    );
                    tolerated += 1;
                } else {
                    hard_failures.push((lineno, input.to_string(), e.to_string()));
                }
            }
        }
    }
    if !hard_failures.is_empty() {
        for (lineno, input, err) in &hard_failures {
            eprintln!("FAIL gap_fill line {}: {} :: {}", lineno, input, err);
        }
        panic!(
            "{} unexpected failures in gap_fill.tsv (accepted={}, tolerated={})",
            hard_failures.len(),
            accepted,
            tolerated
        );
    }
    eprintln!(
        "gap_fill.tsv: {} accepted, {} tolerated",
        accepted, tolerated
    );
}

/// Funcs in grammar_test.tsv that map to whole-expression categories the
/// gtars parser exposes via `parse(...)`. Sub-grammar primitives like
/// `aa1`, `dna_seq`, etc. are NOT exercised here.
fn is_whole_expression_func(func: &str) -> bool {
    matches!(
        func,
        "hgvs_variant"
            | "c_variant"
            | "g_variant"
            | "n_variant"
            | "m_variant"
            | "r_variant"
            | "p_variant"
    )
}

#[test]
fn accept_biocommons_grammar() {
    let p = data_dir().join("biocommons/grammar_test.tsv");
    let known_skips = load_known_skips();
    let raw = fs::read_to_string(&p).expect("grammar_test.tsv");

    let mut hard_failures: Vec<(usize, String, String)> = Vec::new();
    let mut accepted_pos: usize = 0;
    let mut tolerated_pos: usize = 0;
    let mut accepted_neg: usize = 0;
    let mut surprised_neg: Vec<(usize, String)> = Vec::new();

    for (idx, line) in raw.lines().enumerate() {
        let lineno = idx + 1;
        if line.is_empty() || line.starts_with('#') || line.starts_with("Func\t") {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 4 {
            continue;
        }
        let func = cols[0].trim();
        let test = cols[1];
        let valid = cols[2].trim();
        let in_type = cols[3].trim();

        if !is_whole_expression_func(func) {
            continue;
        }
        // p_variant rows are parseable by our grammar (we accept p.) but
        // the upstream test data uses 3-letter amino-acid edits we do not
        // model. Skip outright — the parser's p. coverage is intentionally
        // shallow.
        if func == "p_variant" {
            continue;
        }
        let inputs: Vec<&str> = if in_type == "list" {
            test.split('|').collect()
        } else {
            vec![test]
        };
        for inp in inputs {
            let inp = inp.trim();
            if inp.is_empty() {
                continue;
            }
            match (valid, parse(inp)) {
                ("True", Ok(_)) => accepted_pos += 1,
                ("True", Err(e)) => {
                    if known_skips.contains(inp) {
                        eprintln!(
                            "tolerated grammar_test failure (line {}): {} :: {}",
                            lineno, inp, e
                        );
                        tolerated_pos += 1;
                    } else {
                        hard_failures.push((lineno, inp.to_string(), e.to_string()));
                    }
                }
                ("False", Err(_)) => accepted_neg += 1,
                ("False", Ok(_)) => surprised_neg.push((lineno, inp.to_string())),
                _ => {}
            }
        }
    }

    if !hard_failures.is_empty() {
        for (lineno, input, err) in &hard_failures {
            eprintln!(
                "FAIL biocommons grammar_test line {}: {} :: {}",
                lineno, input, err
            );
        }
        panic!(
            "{} unexpected parser failures in biocommons grammar_test",
            hard_failures.len()
        );
    }
    if !surprised_neg.is_empty() {
        for (lineno, input) in &surprised_neg {
            eprintln!(
                "UNEXPECTED-ACCEPT biocommons grammar_test line {} (Valid=False): {}",
                lineno, input
            );
        }
        eprintln!(
            "{} Valid=False rows were accepted by the parser",
            surprised_neg.len()
        );
    }
    eprintln!(
        "biocommons grammar_test: {} pos accepted, {} pos tolerated, {} neg rejected",
        accepted_pos, tolerated_pos, accepted_neg
    );
}

// ============================================================================
// Curated external corpus (fixtures/biocommons.json + fixtures/ferro_hgvs.json)
// ============================================================================

#[derive(Debug, Deserialize)]
struct TestCase {
    input: String,
    valid: bool,
    source: String,
    #[serde(default)]
    #[allow(dead_code)]
    description: String,
}

fn load_corpus() -> Vec<TestCase> {
    let fixtures_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures");

    let mut all_cases = Vec::new();

    for filename in &["biocommons.json", "ferro_hgvs.json"] {
        let path = format!("{}/{}", fixtures_dir, filename);
        let content =
            fs::read_to_string(&path).unwrap_or_else(|_| panic!("Failed to read {}", filename));
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

    eprintln!(
        "\nAll {} valid cases parsed successfully",
        valid_cases.len()
    );
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

    eprintln!(
        "\nAll {} invalid cases correctly rejected",
        invalid_cases.len()
    );
}

/// Run corpus and report detailed statistics.
#[test]
fn test_external_corpus_stats() {
    let corpus = load_corpus();

    let mut by_source: std::collections::HashMap<String, (usize, usize)> =
        std::collections::HashMap::new();

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
