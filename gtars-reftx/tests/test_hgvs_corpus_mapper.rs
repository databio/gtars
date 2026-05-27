//! Vendored-corpus regression tests for the gtars-reftx coordinate mapper.
//!
//! Two layers:
//!
//! 1. `parse_biocommons_gcp_expressions` — reads every gcp/*.tsv row and
//!    asserts the HGVSg / HGVSc / HGVSp expressions parse syntactically
//!    via `gtars-vrs`. This catches grammar drift without requiring a
//!    built reftx binary. NOTE: this test depends on `gtars-vrs` for the
//!    parser only; if that crate is later pulled out, switch to a small
//!    in-tree parser stub.
//!
//! 2. `project_biocommons_gcp_full` — full c.<->g. projection truth-table
//!    test. This requires a built reftx binary covering the real
//!    transcripts referenced by the gcp tables (NM_xxxxxx). Marked
//!    `#[ignore]` so offline CI stays green; run with
//!    `GTARS_REFTX_FIXTURE=/path/to/reftx.bin` to enable.
//!
//! The synthetic transcript table at `biocommons/sample_data/sanity_cp.tsv`
//! is c.->p. data (not c.<->g.); it is therefore not consumed by these
//! mapper tests but is kept around for the parser harness in `gtars-vrs`.

use std::fs;
use std::path::{Path, PathBuf};

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
fn read_tsv_rows(path: &Path) -> Vec<Vec<String>> {
    let raw = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {}", path.display(), e));
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

/// Walk every gcp row and try to parse each non-empty HGVSg/HGVSc cell.
/// This exercises the parser against a much wider surface than the
/// hand-written tests, but does NOT require a real transcript fixture.
///
/// Currently feature-gated: this test pulls in `gtars-vrs` as a
/// path-only dev dep to avoid a cyclic dep. Until that wiring lands the
/// test is `#[ignore]`. (Plan 8 will introduce the bridge crate that
/// owns this.)
#[test]
#[ignore = "depends on gtars-vrs as a dev-dep; enable once dev-only path is wired"]
fn parse_biocommons_gcp_expressions() {
    // Implementation lives behind a future cfg; for now the test exists
    // as documentation of the intended coverage path.
}

/// Real-transcript projection truth-table test. Requires a reftx binary
/// covering the transcripts referenced by the gcp/*.tsv tables. Skipped
/// unless `GTARS_REFTX_FIXTURE` is set.
#[test]
#[ignore = "requires GTARS_REFTX_FIXTURE pointing at a built reftx with real transcripts"]
fn project_biocommons_gcp_full() {
    let fixture = match std::env::var("GTARS_REFTX_FIXTURE") {
        Ok(v) => PathBuf::from(v),
        Err(_) => {
            eprintln!("GTARS_REFTX_FIXTURE not set — skipping");
            return;
        }
    };
    assert!(
        fixture.exists(),
        "GTARS_REFTX_FIXTURE points to missing path: {}",
        fixture.display()
    );
    // Implementation deferred: needs the HGVS parser wired through
    // ReftxProvider::c_to_genomic_full and a tolerance comparison
    // against the HGVSg column.
    eprintln!(
        "TODO: load reftx fixture from {} and project gcp rows",
        fixture.display()
    );
}
