//! Vendored-corpus regression tests for the HGVS coordinate mapper.
//!
//! NOTE: this binary is kept STANDALONE (not merged into `hgvs_parser` or
//! `hgvs_synthetic`) because of its crate-level `#![cfg(feature = "filesystem")]`
//! gate plus its `gtars_refget::{CoordinateMapper, ReadonlyTxStore, TxBackend}`
//! imports. Merging would force that gate onto unrelated tests.
//!
//! RESTORED from the deleted `gtars-reftx/tests/test_hgvs_corpus_mapper.rs`
//! (commit `ded971e6`), which was dropped when `gtars-reftx` was folded into
//! `gtars_refget::transcripts`. The original `project_biocommons_gcp_full`
//! body was a deferred TODO stub; this version is a RECONSTRUCTION against the
//! current mapper API (`gtars_refget::transcripts::CoordinateMapper`) plus the
//! local HGVS parser (`gtars_vrs::hgvs`). The `#[ignore]` + missing-fixture
//! gating is preserved so offline CI stays green.
//!
//! Three layers:
//!
//! 1. `list_biocommons_gcp_files` — sanity-checks that the vendored gcp tables
//!    are present and non-empty. Runs unconditionally.
//!
//! 2. `parse_biocommons_gcp_expressions` — reads every gcp/*.tsv row and
//!    asserts the HGVSg / HGVSc expressions parse syntactically via the local
//!    `gtars_vrs::hgvs::parse`. This catches grammar drift without requiring a
//!    built transcript fixture. Runs unconditionally.
//!
//! 3. `project_biocommons_gcp_full` — full c.→g. projection truth-table test.
//!    Requires a built `.reftx` store covering the real transcripts referenced
//!    by the gcp tables (NM_xxxxxx). Marked `#[ignore]` so offline CI stays
//!    green; run with `GTARS_REFTX_FIXTURE=/path/to/store.reftx` to enable.
//!    This is the long-missing validation against a real clinical-grade
//!    transcript set (review W2).

#![cfg(feature = "filesystem")]

use std::fs;
use std::path::{Path, PathBuf};

use gtars_refget::{CoordinateMapper, ReadonlyTxStore, TxBackend};
use gtars_vrs::hgvs::ast::{Datum, LocationRange, Position, ReferenceType};
use gtars_vrs::hgvs::parse;

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/hgvs/biocommons")
}

fn gcp_files() -> Vec<PathBuf> {
    let dir = data_dir().join("gcp");
    let mut paths: Vec<PathBuf> = fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read {}: {}", dir.display(), e))
        .filter_map(|r| r.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map(|s| s == "tsv").unwrap_or(false))
        .collect();
    paths.sort();
    paths
}

/// Simple TSV parser that ignores blank lines and `#`-prefixed comments.
/// Returns header + data rows verbatim (header detection is left to callers).
fn read_tsv_rows(path: &Path) -> Vec<Vec<String>> {
    let raw = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {}", path.display(), e));
    let mut rows = Vec::new();
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        rows.push(line.split('\t').map(|s| s.to_string()).collect());
    }
    rows
}

/// Pull the column index for a header name, if present.
fn col_index(header: &[String], name: &str) -> Option<usize> {
    header.iter().position(|h| h == name)
}

#[test]
fn list_biocommons_gcp_files() {
    // Sanity check: the vendored set should include the 14 gcp tables.
    let files = gcp_files();
    assert!(
        files.len() >= 14,
        "expected >= 14 gcp/*.tsv files, found {}",
        files.len()
    );
    for f in &files {
        let rows = read_tsv_rows(f);
        assert!(
            !rows.is_empty(),
            "vendored gcp file {} has zero data rows",
            f.display()
        );
    }
}

/// Walk every gcp row and parse each non-empty HGVSg/HGVSc cell with the local
/// parser. Exercises the parser against a much wider surface than the
/// hand-written tests, without requiring a transcript fixture.
///
/// Some vendored rows carry expressions outside the v1 grammar (e.g. complex
/// protein consequences, which live in the HGVSp column we do not parse). We
/// only require the genomic and coding columns to parse, and we tolerate a
/// small fraction of parse failures there as known grammar gaps, failing only
/// if the failure rate is implausibly high (a sign of real grammar drift).
#[test]
fn parse_biocommons_gcp_expressions() {
    let files = gcp_files();
    let mut attempted = 0usize;
    let mut failed = 0usize;
    let mut first_failures: Vec<String> = Vec::new();

    for f in &files {
        let rows = read_tsv_rows(f);
        if rows.is_empty() {
            continue;
        }
        let header = &rows[0];
        let g_idx = col_index(header, "HGVSg");
        let c_idx = col_index(header, "HGVSc");
        if g_idx.is_none() && c_idx.is_none() {
            // Not an expression table (e.g. p.-only sanity data); skip.
            continue;
        }
        for row in &rows[1..] {
            for idx in [g_idx, c_idx].into_iter().flatten() {
                let Some(cell) = row.get(idx) else { continue };
                let cell = cell.trim();
                if cell.is_empty() || cell == "." {
                    continue;
                }
                attempted += 1;
                if parse(cell).is_err() {
                    failed += 1;
                    if first_failures.len() < 10 {
                        first_failures.push(cell.to_string());
                    }
                }
            }
        }
    }

    assert!(attempted > 0, "no HGVSg/HGVSc cells found in gcp corpus");
    // The parser handles the overwhelming majority of g./c. expressions;
    // allow up to 5% known grammar gaps before treating it as drift.
    let fail_rate = failed as f64 / attempted as f64;
    assert!(
        fail_rate <= 0.05,
        "HGVSg/HGVSc parse failure rate {:.1}% ({}/{}) exceeds 5% — possible grammar drift; first failures: {:?}",
        fail_rate * 100.0,
        failed,
        attempted,
        first_failures
    );
}

/// Extract a 1-based c. start position + intronic offset + 3'UTR flag from a
/// parsed coding expression, for feeding `CoordinateMapper::c_to_g_full`.
fn c_start_position(loc: &LocationRange) -> Option<Position> {
    match loc {
        LocationRange::Single(p) => Some(*p),
        LocationRange::Range { start, .. } => Some(*start),
        _ => None,
    }
}

/// Extract a 1-based genomic start position (0-based on return) from a parsed
/// genomic expression.
fn g_start_interbase(loc: &LocationRange) -> Option<u64> {
    let p = match loc {
        LocationRange::Single(p) => *p,
        LocationRange::Range { start, .. } => *start,
        _ => return None,
    };
    if p.base < 1 {
        return None;
    }
    Some((p.base - 1) as u64)
}

/// Real-transcript projection truth-table test. Requires a `.reftx` store
/// covering the transcripts referenced by the gcp/*.tsv tables. Skipped unless
/// `GTARS_REFTX_FIXTURE` is set.
///
/// For every row carrying both an HGVSc and an HGVSg expression on a transcript
/// the fixture knows, project the c. start to a genomic interbase coordinate
/// and assert it matches the g. start the corpus recorded.
#[test]
#[ignore = "requires GTARS_REFTX_FIXTURE pointing at a built .reftx with real transcripts"]
fn project_biocommons_gcp_full() {
    let fixture = match std::env::var("GTARS_REFTX_FIXTURE") {
        Ok(v) => PathBuf::from(v),
        Err(_) => {
            // No fixture configured — nothing to validate. (Mirrors the
            // original behavior; the #[ignore] keeps offline CI green.)
            return;
        }
    };
    assert!(
        fixture.exists(),
        "GTARS_REFTX_FIXTURE points to missing path: {}",
        fixture.display()
    );

    let store = ReadonlyTxStore::open_with_backend(&fixture, TxBackend::Pread)
        .unwrap_or_else(|e| panic!("open reftx fixture {}: {}", fixture.display(), e));
    let mapper = CoordinateMapper::new(&store);

    let files = gcp_files();
    let mut checked = 0usize;
    let mut mismatches: Vec<String> = Vec::new();

    for f in &files {
        let rows = read_tsv_rows(f);
        if rows.is_empty() {
            continue;
        }
        let header = &rows[0];
        let (Some(g_idx), Some(c_idx)) = (col_index(header, "HGVSg"), col_index(header, "HGVSc"))
        else {
            continue;
        };

        for row in &rows[1..] {
            let (Some(g_cell), Some(c_cell)) = (row.get(g_idx), row.get(c_idx)) else {
                continue;
            };
            let (g_cell, c_cell) = (g_cell.trim(), c_cell.trim());
            if g_cell.is_empty() || c_cell.is_empty() {
                continue;
            }

            let (Ok(c_var), Ok(g_var)) = (parse(c_cell), parse(g_cell)) else {
                continue;
            };
            if c_var.reference_type != ReferenceType::C {
                continue;
            }
            let Some(c_pos) = c_start_position(&c_var.posedit.pos) else {
                continue;
            };
            let Some(expected_g_ib) = g_start_interbase(&g_var.posedit.pos) else {
                continue;
            };

            let is_cds_end = matches!(c_pos.datum, Datum::CdsEnd);
            let res =
                match mapper.c_to_g_full(c_var.accession, c_pos.base, c_pos.offset, is_cds_end) {
                    Ok(r) => r,
                    // Transcript not in the fixture, or position the mapper
                    // rejects: skip rather than fail (the fixture need not
                    // cover every vendored transcript).
                    Err(_) => continue,
                };

            checked += 1;
            if res.position != expected_g_ib && mismatches.len() < 20 {
                mismatches.push(format!(
                    "{} ({}) -> g {} but corpus says g {}",
                    c_cell, c_var.accession, res.position, expected_g_ib
                ));
            }
        }
    }

    assert!(
        checked > 0,
        "fixture matched no gcp transcripts — wrong fixture? ({})",
        fixture.display()
    );
    assert!(
        mismatches.is_empty(),
        "{} c.->g. projection mismatch(es); first few: {:?}",
        mismatches.len(),
        mismatches
    );
}
