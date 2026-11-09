//! `rgblast report` on committed jungle hits (`tests/data/fixtures/origin`, made by
//! `scripts/make_origin_fixtures.sh`) and on hand-written hits.
//! `UPDATE_EXPECTED=1` rewrites the `expected_*` files.

use gtars_rgblast::api::{ReportRequest, report};
use gtars_rgblast::report::output;
use gtars_rgblast::report::parse::DEFAULT_FIELDS;
use gtars_rgblast::report::{RankOptions, Report};
use std::path::{Path, PathBuf};

const DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/fixtures/origin");

fn fx(name: &str) -> PathBuf {
    Path::new(DIR).join(name)
}

fn req(hits: &str) -> ReportRequest {
    ReportRequest {
        hits: vec![fx(hits)],
        db: Some(fx("small_jungle")),
        query: Some(fx("small_jungle_queries.fa")),
        fields: DEFAULT_FIELDS.into(),
        ..Default::default()
    }
}

fn run(r: &ReportRequest) -> Report {
    let mut rep = report(r).unwrap();
    rep.meta.db = None; // an absolute path
    rep
}

fn check(name: &str, got: &str) {
    let p = fx(name);
    if std::env::var_os("UPDATE_EXPECTED").is_some() {
        std::fs::write(&p, got).unwrap();
    }
    let want = std::fs::read_to_string(&p).unwrap_or_else(|_| panic!("missing {name}"));
    assert_eq!(got, want, "{name} differs (UPDATE_EXPECTED=1 rewrites it)");
}

#[test]
fn small_jungle_report() {
    let r = run(&req("small_jungle_hits.tsv"));
    let tsv = output::tsv(&r);
    // the measured tie counts: 37 human genomes, 27 mouse, 13 for EBV, and a none line
    let lines: Vec<&str> = tsv.lines().collect();
    assert_eq!(lines.len(), 5);
    assert!(lines[1].starts_with("q_human_chrM_2001\ttaxon\t9606\tHomo sapiens\t-\tGRCh37-igenomes-ensembl\tMT\t100.000\t100.0\t1847\t0.0\t37\t1\t79/79\t-"));
    assert!(lines[2].contains("\t27\t1\t79/79\t-"));
    assert!(lines[3].contains("\tchrEBV\t") && lines[3].contains("\t13\t1\t79/79\t-"));
    assert_eq!(
        lines[4],
        "q_random\tnone\t-\t-\t-\t-\t-\t-\t-\t-\t-\t0\t0\t79/79\t-"
    );
    check("expected_summary.tsv", &tsv);
    check("expected_detail.tsv", &output::detail(&r, None));
    check("expected.json", &output::json(&r));
    let human = &r.queries[0];
    assert_eq!((human.taxa[0].hit, human.taxa[1].hit), (49, 28));
    assert_eq!(human.taxa[1].best_pident, "75.766");
    // the JSON document reads back
    let back: Report = serde_json::from_str(&output::json(&r)).unwrap();
    assert_eq!(back.queries, r.queries);
}

fn rows_without_evalue(r: &Report) -> Vec<String> {
    output::detail(r, None)
        .lines()
        .map(|l| {
            let mut f: Vec<&str> = l.split('\t').collect();
            f.remove(13);
            f.join("\t")
        })
        .collect()
}

#[test]
fn views_give_the_same_genome_rows() {
    let full = run(&req("small_jungle_hits.tsv"));
    // an SQ. seqid list keeps only SQ. ids: genomes come from the digest; E-values differ
    let sq = run(&req("small_jungle_hits.sqview.tsv"));
    assert_eq!(rows_without_evalue(&full), rows_without_evalue(&sq));
    assert!(output::tsv(&sq).contains("\t37\t1\t79/79"));
    // a gnl| view keeps only the view's deflines: the same rows for its genome
    let mm10 = run(&req("small_jungle_hits.mm10view.tsv"));
    let pick = |r: &Report| -> Vec<String> {
        rows_without_evalue(r)
            .into_iter()
            .filter(|l| l.contains("\tmm10-ucsc-initial-soft-masked\t"))
            .map(|l| l.replace("\tyes\t", "\t").replace("\tno\t", "\t"))
            .collect()
    };
    assert_eq!(pick(&full), pick(&mm10));
    assert_eq!(pick(&mm10).len(), 2);
}

#[test]
fn missing_column_is_named() {
    let mut r = req("small_jungle_hits.tsv");
    r.hits = vec![fx("small_jungle_queries.fa")];
    r.fields = "std sallseqid".into();
    let e = report(&r).unwrap_err();
    assert!(
        e.to_string()
            .contains("report needs columns: qlen (add it to -outfmt"),
        "{e}"
    );
}

const HUMAN: &str = "SQ.k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct";
const MOUSE: &str = "SQ.1BOrY0gjxARFX2lPLmYnmsXK8yfyWLAX";

/// A hit line in the default fields.
fn line(q: &str, subject: &str, pident: &str, bits: &str, score: u32) -> String {
    format!(
        "{q}\t{subject}\t{pident}\t1000\t0\t0\t1\t1000\t1\t1000\t0.0\t{bits}\t1000\t{score}\t{subject}\tt\t0\n"
    )
}

fn hand(lines: &[String], rank: RankOptions, taxdump: bool) -> Report {
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("hits.tsv");
    std::fs::write(&p, lines.concat()).unwrap();
    let mut r = req("small_jungle_hits.tsv");
    r.hits = vec![p];
    r.query = None;
    r.rank = rank;
    r.taxdump = taxdump.then(|| fx("taxdump_mini"));
    run(&r)
}

#[test]
fn ties_lca_and_multi() {
    let tie = vec![
        line("q", HUMAN, "100.000", "1847", 1000),
        line("q", MOUSE, "100.000", "1847", 1000),
    ];
    let r = hand(&tie, RankOptions::default(), true);
    let row = &r.summary_rows()[0];
    assert_eq!(
        (
            row.status.as_str(),
            row.taxid.as_str(),
            row.taxon.as_str(),
            row.rank.as_str()
        ),
        ("lca", "314146", "Euarchontoglires", "superorder")
    );
    // every human genome with that chrM and every mouse genome with that chrM
    assert_eq!((row.tied_genomes, row.tied_taxa), (37 + 24, 2));
    let r = hand(&tie, RankOptions::default(), false);
    let row = &r.summary_rows()[0];
    // more tied human genomes first
    assert_eq!(
        (
            row.status.as_str(),
            row.taxid.as_str(),
            row.taxon.as_str(),
            row.rank.as_str()
        ),
        ("multi", "9606,10090", "Homo sapiens,Mus musculus", "-")
    );
    // raw score decides over equal printed bits
    let r = hand(
        &[
            line("q", HUMAN, "100.000", "1847", 1000),
            line("q", MOUSE, "100.000", "1847", 999),
        ],
        RankOptions::default(),
        false,
    );
    assert_eq!(r.summary_rows()[0].status, "taxon");
    assert_eq!(r.summary_rows()[0].taxid, "9606");
}

#[test]
fn tie_percent_pulls_in_the_next_tier() {
    let mut r = req("small_jungle_hits.tsv");
    r.rank.tie_percent = 1.0;
    let rep = run(&r);
    // 37 identical + the hg18/hg19 UCSC genomes at 99.700% (score 991 >= 990): 8 genomes hold
    // that sequence, one of which also holds the identical one
    assert_eq!(rep.queries[0].tied.genomes, 44);
    assert_eq!(rep.queries[0].status, "taxon");
}

#[test]
fn group_rank_folds_subspecies() {
    let d = tempfile::tempdir().unwrap();
    let base = d.path().join("db");
    // one mouse genome recorded as Mus musculus castaneus (10091)
    let colls = std::fs::read_to_string(fx("small_jungle.collections.tsv")).unwrap();
    let mut changed = false;
    let colls: String = colls
        .lines()
        .map(|l| {
            if !changed && l.contains("\t10090\t") && l.contains("mm10-ucsc-initial-soft-masked") {
                changed = true;
                l.replace("\t10090\t", "\t10091\t") + "\n"
            } else {
                format!("{l}\n")
            }
        })
        .collect();
    assert!(changed);
    std::fs::write(base.with_extension("collections.tsv"), colls).unwrap();
    std::fs::copy(
        fx("small_jungle.members.tsv"),
        base.with_extension("members.tsv"),
    )
    .unwrap();
    let mut r = req("small_jungle_hits.tsv");
    r.db = Some(base);
    r.taxdump = Some(fx("taxdump_mini"));
    let q = &run(&r).queries[1];
    assert_eq!((q.status.as_str(), q.tied.taxa), ("lca", 2));
    assert_eq!(q.call.as_ref().unwrap().name, "Mus musculus");
    r.rank.group_rank = Some("species".into());
    let q = &run(&r).queries[1];
    assert_eq!((q.status.as_str(), q.tied.taxa), ("taxon", 1));
    assert_eq!(q.call.as_ref().unwrap().rank.as_deref(), Some("species"));
}
