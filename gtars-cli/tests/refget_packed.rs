//! CLI tests for the packed layout: `gtars refget compact`, and the
//! standalone file limit on `build` that points users at `--packed`.

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
        p.push(format!("gtars_packed_{}_{}_{}", tag, std::process::id(), nanos));
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

fn run_refget(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_gtars"))
        .arg("refget")
        .args(args)
        .output()
        .expect("failed to run gtars refget")
}

/// A mistyped store path is an error and does not create a directory.
#[test]
fn compact_missing_store_fails_without_creating_it() {
    let scratch = Scratch::new("missing");
    let typo = scratch.path().join("no_such_store");
    let output = run_refget(&["compact", typo.to_str().unwrap()]);
    assert!(!output.status.success(), "compact of a missing store must fail");
    assert!(!typo.exists(), "compact must not create {}", typo.display());
}

/// Compact opens an existing packed store and leaves it readable.
#[test]
fn compact_existing_packed_store() {
    let scratch = Scratch::new("packed");
    let fasta = scratch.path().join("a.fa");
    std::fs::write(&fasta, ">chr1\nAAAACCCC\n>chr2\nGGGGTTTT\n").unwrap();
    let store = scratch.path().join("store");
    let build = run_refget(&[
        "build",
        fasta.to_str().unwrap(),
        "-o",
        store.to_str().unwrap(),
        "--packed",
    ]);
    assert!(build.status.success(), "{}", String::from_utf8_lossy(&build.stderr));

    let compact = run_refget(&["compact", store.to_str().unwrap()]);
    assert!(compact.status.success(), "{}", String::from_utf8_lossy(&compact.stderr));

    // Pack block cache off (the default), then on.
    for cache_mb in ["0", "4"] {
        let out = scratch.path().join(format!("out{cache_mb}.fa"));
        let export = run_refget(&[
            "export",
            "-s",
            store.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "-l",
            "{name}",
            "--pack-cache-mb",
            cache_mb,
        ]);
        assert!(export.status.success(), "{}", String::from_utf8_lossy(&export.stderr));
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            ">chr1\nAAAACCCC\n>chr2\nGGGGTTTT\n",
            "--pack-cache-mb {cache_mb}"
        );
    }
}

/// `compact --prune-unreferenced` deletes a pack file that no index entry
/// points to; plain `compact` keeps it. The store reads back afterwards.
#[test]
fn compact_prune_unreferenced_deletes_unlisted_packs() {
    let scratch = Scratch::new("prune");
    let fasta = scratch.path().join("a.fa");
    std::fs::write(&fasta, ">chr1\nAAAACCCC\n>chr2\nGGGGTTTT\n").unwrap();
    let store = scratch.path().join("store");
    let store_arg = store.to_str().unwrap();
    let build = run_refget(&["build", fasta.to_str().unwrap(), "-o", store_arg, "--packed"]);
    assert!(build.status.success(), "{}", String::from_utf8_lossy(&build.stderr));

    let packs: Vec<PathBuf> = std::fs::read_dir(store.join("packs"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "pack"))
        .collect();
    assert_eq!(packs.len(), 1, "the build should commit one pack: {packs:?}");
    let unlisted = store.join("packs").join("0123456789abcdef0123456789abcdef.pack");
    std::fs::copy(&packs[0], &unlisted).unwrap();

    let plain = run_refget(&["compact", store_arg]);
    assert!(plain.status.success(), "{}", String::from_utf8_lossy(&plain.stderr));
    assert!(unlisted.exists(), "compact without the flag deleted an unlisted pack");

    let pruned = run_refget(&["compact", store_arg, "--prune-unreferenced"]);
    let stderr = String::from_utf8_lossy(&pruned.stderr);
    assert!(pruned.status.success(), "{stderr}");
    assert!(!unlisted.exists(), "the unlisted pack survived --prune-unreferenced");
    assert!(stderr.contains("Pruned 1 unreferenced pack file(s)."), "{stderr}");

    let out = scratch.path().join("out.fa");
    let out_arg = out.to_str().unwrap();
    let export = run_refget(&["export", "-s", store_arg, "-o", out_arg, "-l", "{name}"]);
    assert!(export.status.success(), "{}", String::from_utf8_lossy(&export.stderr));
    assert_eq!(std::fs::read_to_string(&out).unwrap(), ">chr1\nAAAACCCC\n>chr2\nGGGGTTTT\n");
}

/// `--standalone-file-limit` stops a standalone build that would go over it,
/// `--force` does not lift it, and `--packed` is not subject to it.
#[test]
fn build_honors_standalone_file_limit() {
    let scratch = Scratch::new("limit");
    let fasta = scratch.path().join("a.fa");
    std::fs::write(&fasta, ">chr1\nAAAACCCC\n>chr2\nGGGGTTTT\n").unwrap();
    let fasta = fasta.to_str().unwrap();
    let store = scratch.path().join("store");
    let store = store.to_str().unwrap();

    let over = run_refget(&["build", fasta, "-o", store, "--standalone-file-limit", "1", "--force"]);
    assert!(!over.status.success(), "a build over the limit must fail");
    let stderr = String::from_utf8_lossy(&over.stderr);
    assert!(stderr.contains("standalone_file_limit"), "{stderr}");

    let ok = run_refget(&["build", fasta, "-o", store, "--standalone-file-limit", "2"]);
    assert!(ok.status.success(), "{}", String::from_utf8_lossy(&ok.stderr));

    let packed = scratch.path().join("packed");
    let packed = run_refget(&[
        "build", fasta, "-o", packed.to_str().unwrap(), "--packed", "--standalone-file-limit", "1",
    ]);
    assert!(packed.status.success(), "{}", String::from_utf8_lossy(&packed.stderr));
}
