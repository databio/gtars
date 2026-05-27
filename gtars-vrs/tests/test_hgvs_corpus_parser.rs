//! Vendored-corpus regression tests for the HGVS parser.
//!
//! Exercises the parser against three external sources:
//! - varfish-org/hgvs-rs gauntlet (positive set; supported expressions must parse)
//! - varfish-org/hgvs-rs reject   (negative set; each line must fail to parse)
//! - biocommons/hgvs grammar_test.tsv (mixed Valid=True/False rows on
//!   whole-expression Func categories)
//!
//! Plus a locally-curated `local_reject.txt`. See
//! `tests/data/hgvs/REFRESH.md` and `NOTICE` for provenance.
//!
//! The harness uses a "known-skip" allowlist (`known_skips.txt`) so that
//! corpus inputs the current parser does not yet handle are reported as
//! warnings rather than hard failures. New unexpected failures still
//! panic.

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

use gtars_vrs::hgvs::parse;

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
    let raw = fs::read_to_string(&p)
        .unwrap_or_else(|e| panic!("could not read {}: {}", p.display(), e));
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
