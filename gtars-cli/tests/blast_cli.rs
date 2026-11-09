//! `gtars blast`: the same commands as the `rgmakeblastdb` and `rgblast` binaries, through the
//! `gtars` binary (the crate's own tests cover the commands themselves). Stores are built with
//! `gtars refget build --collection-alias`.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const FASTA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/data/fasta");

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
        p.push(format!(
            "gtars_blast_{}_{}_{}",
            tag,
            std::process::id(),
            nanos
        ));
        std::fs::create_dir_all(&p).unwrap();
        Scratch { path: p }
    }
    fn path(&self) -> &Path {
        &self.path
    }
    fn s(&self, name: &str) -> String {
        self.path.join(name).to_str().unwrap().to_string()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

fn gtars(sc: &Scratch, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_gtars"))
        .args(args)
        .env("GTARS_RGBLAST_CACHE", sc.path().join("nstats-cache"))
        .env("SOURCE_DATE_EPOCH", "0")
        .output()
        .expect("run gtars")
}

fn blast(sc: &Scratch, args: &[&str]) -> Output {
    let mut v = vec!["blast"];
    v.extend_from_slice(args);
    gtars(sc, &v)
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

fn ok(o: Output) -> String {
    assert!(o.status.success(), "failed: {}", stderr(&o));
    String::from_utf8_lossy(&o.stdout).to_string()
}

/// A store holding `base.fa` (ucsc:base), `different_names.fa` (ucsc:other, the same sequences
/// under other names) and a collection with sequences of its own (ucsc:uniq).
fn make_store(sc: &Scratch) -> String {
    let store = sc.s("store");
    for (f, alias) in [
        ("base.fa", "ucsc:base"),
        ("different_names.fa", "ucsc:other"),
    ] {
        // build from copies: `refget build` writes an .rgsi next to its input
        let copy = sc.path().join(f);
        std::fs::copy(format!("{FASTA}/{f}"), &copy).unwrap();
        let o = gtars(
            sc,
            &[
                "refget",
                "build",
                copy.to_str().unwrap(),
                "-o",
                &store,
                "--collection-alias",
                alias,
            ],
        );
        assert!(o.status.success(), "{}", stderr(&o));
    }
    let uniq = sc.path().join("uniq.fa");
    std::fs::write(&uniq, ">u1\nACGTACGTACGTTTGACCA\n>u2\nCCCGGGATATATTT\n").unwrap();
    let o = gtars(
        sc,
        &[
            "refget",
            "build",
            uniq.to_str().unwrap(),
            "-o",
            &store,
            "--collection-alias",
            "ucsc:uniq",
        ],
    );
    assert!(o.status.success(), "{}", stderr(&o));
    store
}

#[test]
fn help_lists_all_verbs() {
    let sc = Scratch::new("help");
    let out = ok(gtars(&sc, &["blast", "--help"]));
    for v in [
        "build", "verify", "view", "nstats", "search", "report", "origin", "reps",
    ] {
        assert!(out.contains(v), "{v} missing from help:\n{out}");
    }
}

#[test]
fn build_verify_and_a_usage_error() {
    let sc = Scratch::new("smoke");
    let store = make_store(&sc);
    let db = sc.s("db");
    ok(blast(
        &sc,
        &["build", "--store", &store, "--all", "--out", &db, "-q"],
    ));
    assert!(Path::new(&format!("{db}.nin")).is_file());
    let o = blast(&sc, &["verify", "--db", &db, "--store", &store]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(stderr(&o).contains("OK"));
    // a usage error exits 2
    let o = blast(&sc, &["build", "--store", &store, "--out", &db]);
    assert_eq!(o.status.code(), Some(2), "{}", stderr(&o));
}
