//! CLI smoke tests for `gtars refget export -c NAMESPACE:ALIAS`.
//!
//! Covers resolving a collection alias registered via `build
//! --collection-alias`, the not-found error listing known aliases, the
//! malformed-value rejection, and a bare-digest regression guard.

use std::path::{Path, PathBuf};
use std::process::Command;

/// A unique scratch dir under the system temp dir, cleaned up on drop.
struct Scratch {
    path: PathBuf,
}

impl Scratch {
    fn new(tag: &str) -> Self {
        let mut p = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        p.push(format!("gtars_exportalias_{}_{}_{}", tag, std::process::id(), nanos));
        std::fs::create_dir_all(&p).unwrap();
        Scratch { path: p }
    }
    fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

fn write_fasta(path: &Path, body: &str) {
    std::fs::write(path, body).unwrap();
}

fn run_build(args: &[&str]) -> std::process::Output {
    let bin = env!("CARGO_BIN_EXE_gtars");
    Command::new(bin)
        .arg("refget")
        .arg("build")
        .args(args)
        .output()
        .expect("failed to run gtars refget build")
}

fn run_export(args: &[&str]) -> std::process::Output {
    let bin = env!("CARGO_BIN_EXE_gtars");
    Command::new(bin)
        .arg("refget")
        .arg("export")
        .args(args)
        .output()
        .expect("failed to run gtars refget export")
}

/// Build a store with one FASTA registered under collection alias `ucsc:test`,
/// returning the store dir.
fn build_aliased_store(scratch: &Scratch) -> PathBuf {
    let fasta = scratch.path().join("a.fa");
    write_fasta(&fasta, ">chr1\nAAAACCCC\n>chr2\nGGGGTTTT\n");
    let out_dir = scratch.path().join("store");

    let output = run_build(&[
        fasta.to_str().unwrap(),
        "-o",
        out_dir.to_str().unwrap(),
        "--collection-alias",
        "ucsc:test",
    ]);
    assert!(
        output.status.success(),
        "build should succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    out_dir
}

#[test]
fn export_resolves_collection_alias() {
    let scratch = Scratch::new("ok");
    let store_dir = build_aliased_store(&scratch);
    let out_fa = scratch.path().join("out.fa");

    let output = run_export(&[
        "-s",
        store_dir.to_str().unwrap(),
        "-c",
        "ucsc:test",
        "-o",
        out_fa.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "export should succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let contents = std::fs::read_to_string(&out_fa).unwrap();
    assert!(contents.contains("AAAACCCC"), "output should contain chr1's sequence, got: {}", contents);
    assert!(contents.contains("GGGGTTTT"), "output should contain chr2's sequence, got: {}", contents);
}

#[test]
fn export_unknown_alias_lists_known_aliases() {
    let scratch = Scratch::new("unknown");
    let store_dir = build_aliased_store(&scratch);
    let out_fa = scratch.path().join("out.fa");

    let output = run_export(&[
        "-s",
        store_dir.to_str().unwrap(),
        "-c",
        "ucsc:nope",
        "-o",
        out_fa.to_str().unwrap(),
    ]);
    assert!(!output.status.success(), "unknown alias must be rejected");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Collection alias 'ucsc:nope' not found"),
        "error should name the missing alias, got: {}",
        stderr
    );
    assert!(
        stderr.contains("ucsc:test"),
        "error should list the known alias, got: {}",
        stderr
    );
}

#[test]
fn export_rejects_malformed_alias() {
    let scratch = Scratch::new("badval");
    let store_dir = build_aliased_store(&scratch);
    let out_fa = scratch.path().join("out.fa");

    let output = run_export(&[
        "-s",
        store_dir.to_str().unwrap(),
        "-c",
        ":bad",
        "-o",
        out_fa.to_str().unwrap(),
    ]);
    assert!(!output.status.success(), "malformed alias must be rejected");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("NAMESPACE:ALIAS"),
        "error should explain the expected form, got: {}",
        stderr
    );
}

#[test]
fn export_bare_digest_still_works() {
    let scratch = Scratch::new("digest");
    let store_dir = build_aliased_store(&scratch);

    // The store's one collection file is named by its digest.
    let digest = std::fs::read_dir(store_dir.join("collections"))
        .unwrap()
        .filter_map(|e| e.ok())
        .filter_map(|e| e.file_name().to_str()?.strip_suffix(".rgsi").map(str::to_string))
        .next()
        .expect("expected one collection file");

    let out_fa = scratch.path().join("out.fa");
    let output = run_export(&[
        "-s",
        store_dir.to_str().unwrap(),
        "-c",
        &digest,
        "-o",
        out_fa.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "export by bare digest should still work: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(out_fa.exists());
}

#[test]
fn export_streams_to_stdout_with_labels() {
    let scratch = Scratch::new("stdout");
    let store_dir = build_aliased_store(&scratch);

    let output = run_export(&["-s", store_dir.to_str().unwrap(), "--all", "-l", "{col:ucsc}/{name}"]);
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(stdout, ">test/chr1\nAAAACCCC\n>test/chr2\nGGGGTTTT\n");

    let output = run_export(&["-s", store_dir.to_str().unwrap(), "--all"]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(stdout.matches(">SQ.").count(), 2);

    let output = run_export(&["-s", store_dir.to_str().unwrap(), "--all", "-l", "{bad}"]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("unknown field"));
}
