//! The `rgmakeblastdb` binary: usage errors, messages, behavior.

mod common;

use common::*;
use std::path::Path;
use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rgmakeblastdb"))
        .args(args)
        .env("REFGET_CACHE_DIR", nstats_cache())
        .output()
        .unwrap()
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).into_owned()
}

fn code(o: &Output) -> i32 {
    o.status.code().unwrap_or(-1)
}

fn p(x: &Path) -> String {
    x.display().to_string()
}

#[test]
fn all_and_collection_conflict() {
    let o = run(&[
        "build",
        "--store",
        "/s",
        "--all",
        "--collection",
        "x",
        "--out",
        "/o",
    ]);
    assert_eq!(code(&o), 2, "{}", stderr(&o));
    assert!(stderr(&o).contains("cannot be used with"), "{}", stderr(&o));
    let o = run(&["build", "--store", "/s", "--out", "/o"]);
    assert_eq!(code(&o), 2);
}

#[test]
fn fetch_flags_with_local_are_errors() {
    for flags in [&["--fetch-jobs", "4"][..], &["--resume"]] {
        let flag = flags[0];
        let mut a = vec!["build", "--store", "/s", "--all", "--out", "/o"];
        a.extend(flags);
        let o = run(&a);
        assert_eq!(code(&o), 2, "{flag}: {}", stderr(&o));
        assert!(
            stderr(&o).contains(&format!(
                "{flag} only applies to remote stores (--store http(s)://…)"
            )),
            "{}",
            stderr(&o)
        );
    }
}

#[test]
fn removed_cache_flags_fail() {
    for (flag, val) in [("--cache", "/c"), ("--fetch-log", "/l")] {
        let o = run(&["build", "--store", "/s", "--all", "--out", "/o", flag, val]);
        assert_eq!(code(&o), 2, "{flag}: {}", stderr(&o));
        assert!(stderr(&o).contains("unexpected argument"), "{}", stderr(&o));
    }
}

#[test]
fn quiet_and_verbose_conflict() {
    let o = run(&["build", "--store", "/s", "--all", "--out", "/o", "-q", "-v"]);
    assert_eq!(code(&o), 2, "{}", stderr(&o));
    assert!(stderr(&o).contains("cannot be used with"));
}

#[test]
fn removed_fasta_subcommand_fails() {
    let o = run(&["fasta", "--store", "/s", "--all", "--out", "/o.fa"]);
    assert_eq!(code(&o), 2, "{}", stderr(&o));
}

#[test]
fn removed_ids_flag_fails() {
    for extra in [
        ["--ids", "refget"],
        ["--blastdb-version", "5"],
        ["--max-len", "10"],
        ["--ids-out", "x.tsv"],
    ] {
        let mut a = vec!["build", "--store", "/s", "--all", "--out", "/o"];
        a.extend(extra);
        let o = run(&a);
        assert_eq!(code(&o), 2, "{extra:?}: {}", stderr(&o));
        assert!(stderr(&o).contains("unexpected argument"), "{}", stderr(&o));
    }
}

#[test]
fn missing_store_message() {
    let d = tempdir();
    let empty = d.path().join("empty");
    std::fs::create_dir(&empty).unwrap();
    let o = run(&[
        "build",
        "--store",
        &p(&empty),
        "--all",
        "--out",
        &p(&d.path().join("db")),
    ]);
    assert_eq!(code(&o), 1);
    assert_eq!(
        stderr(&o).trim_end(),
        format!(
            "rgmakeblastdb: error: no RefgetStore at {} (rgstore.json not found)",
            p(&empty)
        )
    );
    let o = run(&[
        "build",
        "--store",
        "/nonexistent/x",
        "--all",
        "--out",
        &p(&d.path().join("db")),
    ]);
    assert_eq!(code(&o), 1);
    assert!(
        stderr(&o).contains("no RefgetStore at /nonexistent/x (no such directory)"),
        "{}",
        stderr(&o)
    );
}

#[test]
fn unknown_collection_lists_namespaces() {
    let d = tempdir();
    let o = run(&[
        "build",
        "--store",
        &encoded().spec(),
        "--collection",
        "nope",
        "--out",
        &p(&d.path().join("db")),
    ]);
    assert_eq!(code(&o), 1);
    let e = stderr(&o);
    assert!(e.contains("unknown collection 'nope'"), "{e}");
    assert!(
        e.contains("collection alias namespaces: name, refseq"),
        "{e}"
    );
}

#[test]
fn collection_value_may_start_with_dash() {
    let f = encoded();
    let dd = f.digest("asmD");
    assert!(dd.starts_with('-'));
    for spelling in [
        vec!["--collection".to_string(), dd.clone()],
        vec![format!("--collection={dd}")],
    ] {
        let d = tempdir();
        let out = p(&d.path().join("db"));
        let mut a: Vec<&str> = vec![
            "build",
            "--store",
            &f.path.to_str().unwrap(),
            "--out",
            &out,
            "-q",
        ];
        a.extend(spelling.iter().map(|s| s.as_str()));
        let o = run(&a);
        assert_eq!(code(&o), 0, "{spelling:?}: {}", stderr(&o));
        assert!(d.path().join("db.nsq").exists());
    }
}

#[test]
fn dry_run_writes_nothing() {
    let d = tempdir();
    let o = run(&[
        "build",
        "--store",
        &encoded().spec(),
        "--all",
        "--out",
        &p(&d.path().join("db")),
        "--dry-run",
    ]);
    assert_eq!(code(&o), 0, "{}", stderr(&o));
    let e = stderr(&o);
    assert!(e.contains("database: 25 sequences"), "{e}");
    assert!(e.contains("exact +"), "{e}");
    assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 0);
}

#[test]
fn force_replaces() {
    let d = tempdir();
    let out = p(&d.path().join("db"));
    let s = encoded().spec();
    let o = run(&[
        "build",
        "--store",
        &s,
        "--collection",
        "asmB",
        "--out",
        &out,
        "-q",
    ]);
    assert_eq!(code(&o), 0, "{}", stderr(&o));
    let o = run(&["build", "--store", &s, "--all", "--out", &out, "-q"]);
    assert_eq!(code(&o), 1);
    assert!(stderr(&o).contains("already exists"), "{}", stderr(&o));
    assert!(
        stderr(&o).contains("hint: use --force to replace it"),
        "{}",
        stderr(&o)
    );
    let o = run(&["build", "--store", &s, "--all", "--out", &out, "--force"]);
    assert_eq!(code(&o), 0, "{}", stderr(&o));
    let e = stderr(&o);
    assert!(e.contains("wrote "), "{e}");
    assert!(e.contains("search it with: blastn -db"), "{e}");
}

#[test]
fn verify_subcommand_passes_and_fails() {
    let d = tempdir();
    let out = d.path().join("db");
    let s = encoded().spec();
    let o = run(&["build", "--store", &s, "--all", "--out", &p(&out), "-q"]);
    assert_eq!(code(&o), 0, "{}", stderr(&o));
    let o = run(&["verify", "--db", &p(&out), "--store", &s]);
    assert_eq!(code(&o), 0, "{}", stderr(&o));
    assert!(stderr(&o).contains("OK"), "{}", stderr(&o));
    let nsq = out.with_extension("nsq");
    let mut b = std::fs::read(&nsq).unwrap();
    b[20] ^= 0x0F;
    std::fs::write(&nsq, b).unwrap();
    let o = run(&["verify", "--db", &p(&out), "--store", &s]);
    assert_eq!(code(&o), 1, "{}", stderr(&o));
    assert!(stderr(&o).contains("1 md5 mismatches"), "{}", stderr(&o));
}

#[test]
fn selection_conflicts_exit_2_with_clap_wording() {
    let o = run(&[
        "build",
        "--store",
        "/s",
        "--namespace",
        "refseq",
        "--all",
        "--out",
        "/o",
    ]);
    assert_eq!(code(&o), 2, "{}", stderr(&o));
    assert!(
        stderr(&o).contains("the argument '--namespace <NS>' cannot be used with '--all'"),
        "{}",
        stderr(&o)
    );
    let o = run(&["build", "--store", "/s", "--out", "/o"]);
    assert_eq!(code(&o), 2);
    assert!(
        stderr(&o).contains("choose collections with"),
        "{}",
        stderr(&o)
    );
}

#[test]
fn view_nstats_subcommands() {
    let f = taxa_encoded();
    let d = tempdir();
    let out = p(&d.path().join("db"));
    let o = run(&[
        "build",
        "--store",
        &f.spec(),
        "--collection",
        "human_a",
        "--out",
        &out,
    ]);
    assert_eq!(code(&o), 0, "{}", stderr(&o));
    let e = stderr(&o);
    assert!(
        e.contains("selection: 1 of 5 collections; taxa: Homo sapiens 1"),
        "{e}"
    );
    assert!(
        e.contains("taxids written for 1 taxon. Filter exactly with: blastn -taxids 9606 -no_taxid_expansion"),
        "{e}"
    );
    // a named view
    let o = run(&[
        "view",
        "--db",
        &out,
        "--name",
        "humans",
        "--taxon",
        "Homo_sapiens",
    ]);
    assert_eq!(code(&o), 0, "{}", stderr(&o));
    assert!(
        stderr(&o).contains("view humans: 4 sequences, 1,040 bases"),
        "{}",
        stderr(&o)
    );
    // nstats prints one line per collection, highest N fraction first
    let o = run(&["nstats", "--store", &f.spec(), "--all"]);
    assert_eq!(code(&o), 0, "{}", stderr(&o));
    let lines: Vec<String> = String::from_utf8_lossy(&o.stdout)
        .lines()
        .map(String::from)
        .collect();
    assert_eq!(lines.len(), 6);
    assert!(
        lines[1].contains("\thuman_masked\tHomo sapiens\t2\t580\t198\t0.3414"),
        "{}",
        lines[1]
    );
}

#[test]
fn rgblast_engine_rule() {
    let rgblast = |args: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_rgblast"))
            .args(args)
            .output()
            .unwrap()
    };
    let selects = "select collections when building the database";
    for (args, msg) in [
        (
            &["search", "--engine", "blastn", "--", "-query", "q.fa"][..],
            "--engine blastn needs --db BASE",
        ),
        (
            &[
                "search", "--engine", "blastn", "--db", "d", "--taxon", "9606",
            ][..],
            selects,
        ),
        (
            &["search", "--engine", "blastn", "--db", "d", "--store", "/s"][..],
            selects,
        ),
        (
            &["search", "--store", "/s", "--all", "--db", "d"][..],
            "--db and --blastn go with --engine blastn",
        ),
        (
            &["search", "--all", "--", "-query", "q.fa"][..],
            "--store PATH|URL is required",
        ),
        (
            &[
                "origin", "--engine", "blastn", "--db", "d", "--staged", "--query", "q.fa",
            ][..],
            "needs the default engine, rgblastn",
        ),
        (
            &[
                "origin", "--engine", "blastn", "--db", "d", "--rep", "x", "--query", "q.fa",
            ][..],
            "needs the default engine, rgblastn",
        ),
        (
            &["search", "--store", "/s", "--all", "--", "-db", "x"][..],
            "-db is not allowed",
        ),
    ] {
        let o = rgblast(args);
        assert_eq!(code(&o), 2, "{args:?}: {}", stderr(&o));
        assert!(stderr(&o).contains(msg), "{args:?}: {}", stderr(&o));
    }
    // removed: the cache verb, text format, taxon report, report --store
    for args in [
        &["cache", "list"][..],
        &["report", "h.tsv", "--format", "text"][..],
        &["report", "h.tsv", "--taxon-report", "t.tsv"][..],
        &["report", "h.tsv", "--store", "/s"][..],
    ] {
        let o = rgblast(args);
        assert_eq!(code(&o), 2, "{args:?}: {}", stderr(&o));
    }
}
