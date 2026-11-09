//! The `rgblastn` binary without BLAST+: flag handling and output through the binary on the
//! fixture store.

mod common;

use common::*;
use std::path::Path;
use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rgblastn"))
        .args(args)
        .output()
        .unwrap()
}

fn out(o: &Output) -> String {
    String::from_utf8_lossy(&o.stdout).into_owned()
}

fn err(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).into_owned()
}

fn seq_of(fasta: &str, name: &str) -> String {
    let text = std::fs::read_to_string(Path::new(FIXTURES).join(fasta)).unwrap();
    let mut cur = None;
    let mut s = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if cur.is_some() {
                break;
            }
            if h.split_whitespace().next() == Some(name) {
                cur = Some(());
            }
        } else if cur.is_some() {
            s.push_str(line.trim());
        }
    }
    s
}

#[test]
fn flags_in_blastn_spelling_and_refusals() {
    let f = encoded();
    let d = tempdir();
    let q = d.path().join("q.fa");
    std::fs::write(
        &q,
        format!(">q1\n{}\n", &seq_of("asmA.fa", "chr1")[100..400]),
    )
    .unwrap();
    let qs = q.to_str().unwrap();
    let base = |extra: &[&str]| -> Output {
        let mut v = vec!["--store", f.spec().leak() as &str, "--all", "-query", qs];
        v.extend(extra);
        run(&v)
    };
    // underscore flags are blastn's spelling
    let o = base(&[
        "-outfmt",
        "6 qseqid sseqid pident",
        "-max_target_seqs",
        "2",
        "-perc_identity",
        "99",
        "-soft_masking",
        "false",
        "-word_size",
        "20",
    ]);
    assert!(o.status.success(), "{}", err(&o));
    let first = out(&o).lines().next().unwrap().to_string();
    assert!(
        first.starts_with("q1\tSQ.") && first.ends_with("\t100.000"),
        "{first}"
    );
    // unknown -outfmt tokens are ignored with a note, as blastn ignores them
    let o = base(&["-outfmt", "6 std bogus"]);
    assert!(o.status.success());
    assert!(
        err(&o).contains("'bogus' is not a format specifier; ignored"),
        "{}",
        err(&o)
    );
    // taxonomy names, -mt_mode and unported flags are refused with exit 2
    for (extra, msg) in [
        (vec!["-outfmt", "6 std ssciname"], "ssciname"),
        (vec!["-mt_mode", "1"], "-mt_mode"),
        (vec!["-num_threads", "0"], "-num_threads must be at least 1"),
        (vec!["-sum_stats", "true"], "-sum_stats"),
        (vec!["-outfmt", "7"], "-outfmt 7"),
        (vec!["-outfmt", "10"], "-outfmt 10"),
        (vec!["-task", "dc-megablast"], "-task dc-megablast"),
        (
            vec!["-dbsize", "1000"],
            "-dbsize is not supported by rgblastn",
        ),
        (vec!["-taxids", "9606"], "rgmakeblastdb"),
        (vec!["-culling_limit", "1"], "rgmakeblastdb"),
        (vec!["-reward", "1"], "rgmakeblastdb"),
        (vec!["-lcase_masking"], "-lcase_masking"),
        (vec!["--dump-stats"], "--dump-stats"),
        (vec!["--stats-nseq", "3"], "must be given together"),
    ] {
        let o = base(&extra);
        assert_eq!(o.status.code(), Some(2), "{extra:?}: {}", err(&o));
        assert!(err(&o).contains(msg), "{extra:?}: {}", err(&o));
        assert!(o.stdout.is_empty());
    }
}
