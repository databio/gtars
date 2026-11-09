//! `rgblast origin`, `search` and `report` on the synthetic origin store (feature `golden`):
//! rgblastn's staged stages, and real blastn (BLAST+ 2.17.0 from `GTARS_RGBLAST_BIN`, PATH, or
//! a test-only shim that runs docker `GTARS_RGBLAST_DOCKER`, default `ncbi/blast:2.17.0`) on a
//! database built from the store, against the hits recorded in `origin_store_hits/`, which the
//! default tests replay; plus the BLAST behaviors `rgblast` depends on. `GOLDEN_UPDATE=1`
//! re-records the hits.

mod common;

use common::blast::Blast;
use common::origin::{HITS, ORIGIN_FIXTURES, origin_store};
use gtars_rgblast::api::{self, OriginRequest, ReportRequest, SearchRequest};
use gtars_rgblast::origin::{ModeChoice, OriginOptions};
use gtars_rgblast::report::output;
use gtars_rgblast::report::parse::DEFAULT_FIELDS;
use gtars_rgblast::run::{
    BlastRunner, BlastnRunner, EngineChoice, EngineOutput, RunSpec, SearchTarget,
};
use gtars_rgblast::select::Selection;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// The database built from the origin store, and a blastn for the runner.
struct Db {
    dir: tempfile::TempDir,
    db: PathBuf,
    blastn: PathBuf,
}

fn db() -> Db {
    let dir = common::tempdir();
    let db = dir.path().join("db");
    let mut b = gtars_rgblast::BuildRequest {
        store: origin_store().store.display().to_string(),
        cache_dir: None,
        selection: Selection::all(),
        options: Default::default(),
    };
    b.options.out = db.clone();
    gtars_rgblast::build(&b, &mut |_| true).unwrap();
    let blastn = Blast::find().blastn_path(dir.path());
    let r = BlastnRunner::resolve(Some(&blastn)).unwrap();
    assert_eq!(
        r.version, "2.17.0+",
        "the golden tests are pinned to BLAST+ 2.17.0"
    );
    Db { dir, db, blastn }
}

fn blastn(d: &Db) -> EngineChoice {
    EngineChoice::Blastn {
        db: d.db.clone(),
        blastn: Some(d.blastn.clone()),
    }
}

fn origin_req(engine: EngineChoice, mode: ModeChoice, keep: &Path) -> OriginRequest {
    let f = origin_store();
    OriginRequest {
        store: match engine {
            EngineChoice::Rgblastn => f.store.display().to_string(),
            _ => String::new(),
        },
        selection: Selection::all(),
        engine,
        cache_dir: None,
        query: f.queries.clone(),
        args: vec!["-num_threads".into(), "2".into()],
        reps: Vec::new(),
        options: OriginOptions {
            mode,
            ..Default::default()
        },
        taxdump: Some(Path::new(ORIGIN_FIXTURES).join("taxdump_mini")),
        keep_hits: Some(keep.to_path_buf()),
        fetch: Default::default(),
    }
}

/// The first `n` columns of each line.
fn cols(text: &str, n: usize) -> Vec<Vec<&str>> {
    text.lines()
        .map(|l| l.split('\t').take(n).collect())
        .collect()
}

fn expected(name: &str) -> String {
    std::fs::read_to_string(Path::new(ORIGIN_FIXTURES).join(name)).unwrap_or_else(|_| {
        panic!("missing {name}; run the default origin tests with UPDATE_EXPECTED=1")
    })
}

#[test]
fn golden_origin_search_report() {
    let d = db();
    let update = std::env::var_os("GOLDEN_UPDATE").is_some();
    if update {
        std::fs::create_dir_all(HITS).unwrap();
    }
    // staged on the store (rgblastn, under the full selection's statistics), exhaustive on the
    // database (blastn): the stage hits equal the recorded ones
    for (engine, mode, name, stages) in [
        (
            EngineChoice::Rgblastn,
            ModeChoice::Staged,
            "staged",
            &["stage1", "stage2"][..],
        ),
        (blastn(&d), ModeChoice::Auto, "exhaustive", &["full"][..]),
    ] {
        let keep = d.dir.path().join(name);
        let o = api::origin(&origin_req(engine, mode, &keep), &mut |_| true).unwrap();
        for s in stages {
            let got = std::fs::read_to_string(keep.join(format!("{s}.tsv"))).unwrap();
            let rec = Path::new(HITS).join(format!("{s}.tsv"));
            if update {
                std::fs::write(&rec, &got).unwrap();
            }
            let want = std::fs::read_to_string(&rec).unwrap();
            // stage2.tsv was recorded with blastn through an SQ. seqid list, which keeps only
            // the SQ. defline; rgblastn keeps every defline, so compare the alignments
            let n = if *s == "stage2" { 14 } else { usize::MAX };
            assert_eq!(cols(&got, n), cols(&want, n), "{name}: {s}.tsv");
        }
        if update {
            continue;
        }
        assert_eq!(o.meta.mode, name);
        assert_eq!(
            output::tsv(&o.report),
            expected(&format!("expected_origin_{name}.tsv")),
            "{name}"
        );
    }
    if update {
        return;
    }
    // search + report give the exhaustive report
    let f = origin_store();
    let out = d.dir.path().join("search.tsv");
    let s = api::search(
        &SearchRequest {
            engine: blastn(&d),
            extra_fields: true,
            args: vec![
                "-query".into(),
                f.queries.display().to_string(),
                "-outfmt".into(),
                "6".into(),
                "-out".into(),
                out.display().to_string(),
            ],
            meta: Some(PathBuf::from(format!("{}.meta.json", out.display()))),
            ..Default::default()
        },
        &mut |_| true,
    )
    .unwrap();
    assert_eq!(s.exit_code, 0);
    assert_eq!(s.db.as_deref(), Some(d.db.as_path()));
    let r = api::report(&ReportRequest {
        hits: vec![out],
        query: Some(f.queries.clone()),
        fields: DEFAULT_FIELDS.into(),
        taxdump: Some(Path::new(ORIGIN_FIXTURES).join("taxdump_mini")),
        ..Default::default()
    })
    .unwrap();
    assert_eq!(output::tsv(&r), expected("expected_origin_exhaustive.tsv"));
}

/// Run blastn on the database through the runner: (exit code, output, stderr).
fn run(d: &Db, extra: &[&str]) -> (i32, String, String) {
    let runner = BlastnRunner::resolve(Some(&d.blastn)).unwrap();
    let mut args: Vec<String> = vec![
        "-query".into(),
        origin_store().queries.display().to_string(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    let spec = RunSpec {
        stage: "full".into(),
        target: SearchTarget::Db(d.db.clone()),
        args,
        stats: None,
        query: None,
        output: EngineOutput::Capture,
        threads_added: false,
    };
    let mut err = Vec::new();
    let o = runner
        .run(&spec, None, &mut |e| {
            if let gtars_rgblast::ProgressEvent::EngineStderr { line } = e {
                err.push(line.clone());
            }
            true
        })
        .unwrap();
    (o.exit_code, o.output.unwrap_or_default(), err.join("\n"))
}

const F6: &str = "6 qaccver sseqid evalue bitscore sallseqid";

#[test]
fn golden_blast_behaviors() {
    let d = db();
    let s = common::origin::seqs();
    let digests: HashMap<&str, String> = [("s1", &s.s1)]
        .into_iter()
        .map(|(k, v)| (k, gtars_refget::digest::sha512t24u(v.as_bytes())))
        .collect();
    let s1 = format!("SQ.{}", digests["s1"]);
    // one sallseqid id per defline: SQ. plus one gnl|<db_tag>| per collection holding chr1
    let (code, out, _) = run(&d, &["-outfmt", F6]);
    assert_eq!(code, 0);
    let line = out
        .lines()
        .find(|l| l.starts_with("q_s1\t") && l.contains(&s1))
        .unwrap();
    let all: Vec<&str> = line.split('\t').nth(4).unwrap().split(';').collect();
    assert_eq!(
        all,
        vec![s1.as_str(), "gnl|A1|chr1", "gnl|A2|chr1", "gnl|A3|chr1"]
    );
    // a gnl| seqid list keeps only its deflines; an SQ. list keeps only SQ.
    let ids = d.dir.path().join("ids.txt");
    let ids_s = ids.display().to_string();
    std::fs::write(&ids, "gnl|A2|chr1\n").unwrap();
    let (_, out, _) = run(&d, &["-outfmt", F6, "-seqidlist", &ids_s]);
    let l: Vec<&str> = out.lines().next().unwrap().split('\t').collect();
    assert_eq!((l[1], l[4]), ("gnl|A2|chr1", "gnl|A2|chr1"));
    std::fs::write(&ids, format!("{s1}\n")).unwrap();
    let (_, out, _) = run(&d, &["-outfmt", F6, "-seqidlist", &ids_s]);
    let l: Vec<&str> = out.lines().next().unwrap().split('\t').collect();
    assert_eq!((l[1], l[4]), (s1.as_str(), s1.as_str()));
    // -taxids without -no_taxid_expansion and no taxonomy4blast.sqlite3: unfiltered
    let (code, out, err) = run(&d, &["-outfmt", F6, "-taxids", "10090"]);
    assert_eq!(code, 0);
    assert!(err.contains("requires additional data files"), "{err}");
    assert!(
        out.lines().any(|l| l.contains("gnl|A1|")),
        "-taxids 10090 still returns human hits"
    );
    let (_, out, _) = run(
        &d,
        &["-outfmt", F6, "-taxids", "10090", "-no_taxid_expansion"],
    );
    assert!(!out.lines().any(|l| l.contains("gnl|A1|")));
    assert!(out.lines().any(|l| l.contains("gnl|B1|")));
}
