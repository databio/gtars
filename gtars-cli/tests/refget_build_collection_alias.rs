//! CLI smoke tests for `gtars refget build --collection-alias`.
//!
//! Covers the happy path (a single FASTA named `ucsc:hg38` produces a
//! collection alias TSV on disk), the malformed-value rejection, and the
//! multi-file rejection whose error message teaches the per-file workaround.

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
        p.push(format!("gtars_collalias_{}_{}_{}", tag, std::process::id(), nanos));
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

#[test]
fn build_with_collection_alias_writes_alias_tsv() {
    let scratch = Scratch::new("ok");
    let fasta = scratch.path().join("a.fa");
    write_fasta(&fasta, ">chr1\nAAAACCCC\n>chr2\nGGGGTTTT\n");
    let out_dir = scratch.path().join("store");

    let output = run_build(&[
        fasta.to_str().unwrap(),
        "-o",
        out_dir.to_str().unwrap(),
        "--collection-alias",
        "ucsc:hg38",
    ]);
    assert!(
        output.status.success(),
        "build should succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let tsv = out_dir.join("aliases").join("collections").join("ucsc.tsv");
    assert!(tsv.exists(), "collection alias TSV should exist at {:?}", tsv);
    let contents = std::fs::read_to_string(&tsv).unwrap();
    assert!(
        contents.contains("hg38\t"),
        "TSV should contain the hg38 alias, got: {}",
        contents
    );
}

#[test]
fn build_rejects_collection_alias_without_colon() {
    let scratch = Scratch::new("badval");
    let fasta = scratch.path().join("a.fa");
    write_fasta(&fasta, ">chr1\nAAAACCCC\n");
    let out_dir = scratch.path().join("store");

    let output = run_build(&[
        fasta.to_str().unwrap(),
        "-o",
        out_dir.to_str().unwrap(),
        "--collection-alias",
        "hg38",
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
fn build_rejects_collection_alias_with_multiple_fastas() {
    let scratch = Scratch::new("multi");
    let fasta_a = scratch.path().join("a.fa");
    let fasta_b = scratch.path().join("b.fa");
    write_fasta(&fasta_a, ">chr1\nAAAACCCC\n");
    write_fasta(&fasta_b, ">chr1\nGGGGTTTT\n");
    let out_dir = scratch.path().join("store");

    let output = run_build(&[
        fasta_a.to_str().unwrap(),
        fasta_b.to_str().unwrap(),
        "-o",
        out_dir.to_str().unwrap(),
        "--collection-alias",
        "ucsc:hg38",
    ]);
    assert!(!output.status.success(), "multi-file alias must be rejected");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("add_collection_alias"),
        "error should teach the per-file workaround, got: {}",
        stderr
    );
}
