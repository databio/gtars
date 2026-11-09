//! `rgblast origin` orchestration on the synthetic origin store, with a fake runner that
//! replays `origin_store_hits/` (recorded by the golden test) and records every `RunSpec`.
//! `UPDATE_EXPECTED=1` rewrites `expected_origin_*.tsv`.

mod common;

use common::origin::{HITS, ORIGIN_FIXTURES, origin_store, seqs};
use gtars_rgblast::api::{OriginRequest, origin_with};
use gtars_rgblast::origin::{ModeChoice, OriginOptions};
use gtars_rgblast::report::output;
use gtars_rgblast::run::{
    BlastRunner, EngineChoice, EngineInfo, EngineOutput, RunOutcome, RunSpec, SearchTarget,
};
use gtars_rgblast::select::Selection;
use gtars_rgblast::{Error, OpenedStore, ProgressEvent};
use std::path::Path;
use std::sync::Mutex;

#[derive(Default)]
struct FakeRunner {
    specs: Mutex<Vec<RunSpec>>,
}

impl BlastRunner for FakeRunner {
    fn describe(&self) -> String {
        "fake blastn (replay)".into()
    }
    fn info(&self) -> EngineInfo {
        EngineInfo {
            name: "fake".into(),
            version: "0".into(),
            via: "replay".into(),
        }
    }
    fn run(
        &self,
        spec: &RunSpec,
        _store: Option<&mut OpenedStore>,
        _progress: &mut dyn FnMut(&ProgressEvent) -> bool,
    ) -> Result<RunOutcome, Error> {
        let EngineOutput::File(out) = &spec.output else {
            panic!("origin writes stage output to a file")
        };
        std::fs::copy(Path::new(HITS).join(format!("{}.tsv", spec.stage)), out).unwrap();
        self.specs.lock().unwrap().push(spec.clone());
        Ok(RunOutcome::default())
    }
}

fn req(mode: ModeChoice) -> OriginRequest {
    let f = origin_store();
    OriginRequest {
        store: f.store.display().to_string(),
        selection: Selection::all(),
        engine: Default::default(),
        cache_dir: None,
        query: f.queries.clone(),
        args: vec!["-num_threads".into(), "2".into()],
        reps: Vec::new(),
        options: OriginOptions {
            mode,
            ..Default::default()
        },
        taxdump: Some(Path::new(ORIGIN_FIXTURES).join("taxdump_mini")),
        keep_hits: None,
        fetch: Default::default(),
    }
}

/// The subjects of a store search.
fn subjects(sp: &RunSpec) -> Option<Vec<String>> {
    match &sp.target {
        SearchTarget::Store { subjects, .. } => subjects.as_ref().map(|d| d.to_vec()),
        x => panic!("{x:?}"),
    }
}

fn check(name: &str, got: &str) {
    let p = Path::new(ORIGIN_FIXTURES).join(name);
    if std::env::var_os("UPDATE_EXPECTED").is_some() {
        std::fs::write(&p, got).unwrap();
    }
    assert_eq!(got, std::fs::read_to_string(&p).unwrap(), "{name}");
}

fn digest(s: &str) -> String {
    gtars_refget::digest::sha512t24u(s.as_bytes())
}

fn sorted(v: &[&String]) -> Vec<String> {
    let mut v: Vec<String> = v.iter().map(|s| digest(s)).collect();
    v.sort();
    v
}

fn row<'a>(tsv: &'a str, q: &str) -> Vec<&'a str> {
    tsv.lines()
        .find(|l| l.starts_with(&format!("{q}\t")))
        .unwrap()
        .split('\t')
        .collect()
}

#[test]
fn staged_stages_and_scopes() {
    let s = seqs();
    let fake = FakeRunner::default();
    let o = origin_with(&req(ModeChoice::Staged), &fake, &mut |_| true).unwrap();
    let specs = fake.specs.lock().unwrap();
    assert_eq!(
        specs.iter().map(|s| s.stage.as_str()).collect::<Vec<_>>(),
        ["stage1", "stage2"]
    );
    // every stage searches under the full selection's statistics
    for sp in specs.iter() {
        let st = sp.stats.as_ref().unwrap();
        assert_eq!((st.num_seqs, st.total_length), (10, 23000), "{}", sp.stage);
        assert!(
            sp.args
                .windows(2)
                .any(|w| w[0] == "-outfmt" && w[1].starts_with("6 qaccver"))
        );
    }
    // stage 1: the representatives A1, B1, C1
    assert_eq!(
        subjects(&specs[0]).unwrap(),
        sorted(&[&s.s1, &s.s2, &s.sm, &s.b1, &s.bm, &s.b2, &s.u])
    );
    // stage 2: the candidates' sequences minus the representatives
    assert_eq!(subjects(&specs[1]).unwrap(), sorted(&[&s.s2p, &s.d, &s.v]));
    let tsv = output::tsv(&o.report);
    check("expected_origin_staged.tsv", &tsv);
    // identical copies tie; the SNP copy is found only by stage 2
    assert_eq!(row(&tsv, "q_s1")[11], "3");
    let snp = row(&tsv, "q_snp");
    assert_eq!((snp[5], snp[7], snp[11]), ("A2", "100.000", "1"));
    // the decoy is in no representative
    let decoy = row(&tsv, "q_decoy");
    assert_eq!(
        (decoy[1], decoy[13], decoy[14]),
        ("none", "3/6", "reps-only")
    );
    // the chimera's stage 2 hit on vec1 fans out to its own candidate (A3), not to B2
    let chim = &o.report.queries[5];
    assert_eq!(chim.query, "q_chim");
    assert_eq!(chim.status, "taxon");
    assert!(
        chim.taxa
            .iter()
            .all(|t| t.genomes.iter().all(|g| g.genome != "B2"))
    );
    assert!(
        chim.taxa[0]
            .genomes
            .iter()
            .any(|g| g.genome == "A3" && g.stage == "fanout")
    );
    let c = o
        .meta
        .candidates
        .iter()
        .find(|c| c.query == "q_chim")
        .unwrap();
    assert_eq!(c.groups, ["Homo sapiens"]);
    // the collection with no FHR is its own group
    let u = row(&tsv, "q_u");
    assert_eq!((u[1], u[2], u[3]), ("taxon", "unknown:C1", "C1"));
    assert_eq!(row(&tsv, "q_s1")[13], "5/6");
    assert!(o.meta.notes[0].starts_with("origin: staged; stage 1: 3 representatives"));
}

#[test]
fn exhaustive_on_the_store_and_on_a_database() {
    let fake = FakeRunner::default();
    let o = origin_with(&req(ModeChoice::Exhaustive), &fake, &mut |_| true).unwrap();
    {
        let specs = fake.specs.lock().unwrap();
        assert_eq!(specs.len(), 1);
        assert_eq!(subjects(&specs[0]), None);
        assert!(specs[0].stats.is_none());
    }
    let tsv = output::tsv(&o.report);
    check("expected_origin_exhaustive.tsv", &tsv);
    assert_eq!(row(&tsv, "q_decoy")[1], "taxon");
    // exhaustive sees B2's vec1 for the chimera: human and mouse tie
    assert_eq!(row(&tsv, "q_chim")[1], "lca");
    let vec = row(&tsv, "q_vec");
    assert_eq!(
        (vec[1], vec[2], vec[3], vec[4]),
        ("lca", "314146", "Euarchontoglires", "superorder")
    );

    // blastn on a database built from the same store: one search, the same report
    let t = common::tempdir();
    let db = t.path().join("db");
    let mut b = gtars_rgblast::BuildRequest {
        store: origin_store().store.display().to_string(),
        cache_dir: None,
        selection: Selection::all(),
        options: Default::default(),
    };
    b.options.out = db.clone();
    gtars_rgblast::build(&b, &mut |_| true).unwrap();
    let mut r = req(ModeChoice::Auto);
    r.store = String::new();
    r.engine = EngineChoice::Blastn {
        db: db.clone(),
        blastn: None,
    };
    let fake = FakeRunner::default();
    let o = origin_with(&r, &fake, &mut |_| true).unwrap();
    {
        let specs = fake.specs.lock().unwrap();
        assert_eq!(specs.len(), 1);
        assert!(matches!(&specs[0].target, SearchTarget::Db(p) if *p == db));
    }
    assert_eq!(o.meta.mode, "exhaustive");
    assert!(o.meta.reps.is_none());
    assert_eq!(output::tsv(&o.report), tsv);
    // staged origin needs rgblastn
    r.options.mode = ModeChoice::Staged;
    let e = origin_with(&r, &FakeRunner::default(), &mut |_| true).unwrap_err();
    assert!(
        e.to_string().contains("needs the default engine, rgblastn"),
        "{e}"
    );
}

#[test]
fn auto_mode() {
    // representatives plus the largest remainder are the whole of this tiny store
    let fake = FakeRunner::default();
    let o = origin_with(&req(ModeChoice::Auto), &fake, &mut |_| true).unwrap();
    assert_eq!(o.meta.mode, "exhaustive");
    assert!(
        o.meta
            .notes
            .iter()
            .any(|n| n.contains("running one full search"))
    );
}

#[test]
fn flag_spellings_resolve_alike() {
    use gtars_rgblast::cli::SelectionArgs;
    let f = origin_store();
    let mut st = OpenedStore::open(&f.store.display().to_string(), None).unwrap();
    let key = |a: SelectionArgs, st: &mut OpenedStore| {
        let mut c = a.resolve(st).unwrap().collections;
        c.sort();
        c
    };
    let by_taxon = key(
        SelectionArgs {
            taxon: vec!["9606".into()],
            ..Default::default()
        },
        &mut st,
    );
    let by_list = key(
        SelectionArgs {
            collection: vec!["A3".into(), "name:A1".into(), "A2".into()],
            ..Default::default()
        },
        &mut st,
    );
    assert_eq!(by_taxon, by_list);
    let by_name = key(
        SelectionArgs {
            taxon: vec!["homo_sapiens".into()],
            ..Default::default()
        },
        &mut st,
    );
    assert_eq!(by_taxon, by_name);
}

#[test]
fn rgblastn_search_then_report_reads_the_store_from_the_meta() {
    use gtars_rgblast::api::{ReportRequest, SearchRequest, report, search};
    let f = origin_store();
    let t = common::tempdir();
    let out = t.path().join("hits.tsv");
    let s = search(
        &SearchRequest {
            store: f.store.display().to_string(),
            selection: Selection::all(),
            extra_fields: true,
            args: vec![
                "-query".into(),
                f.queries.display().to_string(),
                "-outfmt".into(),
                "6".into(),
                "-out".into(),
                out.display().to_string(),
            ],
            meta: Some(t.path().join("hits.tsv.meta.json")),
            ..Default::default()
        },
        &mut |_| true,
    )
    .unwrap();
    assert_eq!(s.engine.name, "rgblastn");
    assert_eq!(s.store.as_ref().unwrap().collections.len(), 6);
    let r = report(&ReportRequest {
        hits: vec![out],
        query: Some(f.queries.clone()),
        fields: gtars_rgblast::report::parse::DEFAULT_FIELDS.into(),
        taxdump: Some(Path::new(ORIGIN_FIXTURES).join("taxdump_mini")),
        ..Default::default()
    })
    .unwrap();
    assert_eq!(
        output::tsv(&r),
        std::fs::read_to_string(Path::new(ORIGIN_FIXTURES).join("expected_origin_exhaustive.tsv"))
            .unwrap()
    );
}
