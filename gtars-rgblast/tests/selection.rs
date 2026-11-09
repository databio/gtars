//! Selection flags, their rules, and resolution against the `taxa` fixture store.

mod common;

use clap::Parser;
use common::*;
use gtars_refget::store::StorageMode;
use gtars_rgblast::api::ViewRequest;
use gtars_rgblast::cli::{NstatsArgs, SelectionArgs, ViewArgs};
use gtars_rgblast::select::{CollectionSet, Selection};
use gtars_rgblast::{Error, ErrorKind, NstatsRequest, OpenedStore, ResolvedSelection};

#[derive(Parser, Debug)]
struct S {
    #[command(flatten)]
    a: SelectionArgs,
}
#[derive(Parser, Debug)]
struct V {
    #[command(flatten)]
    a: ViewArgs,
}
#[derive(Parser, Debug)]
struct N {
    #[command(flatten)]
    a: NstatsArgs,
}

fn sel(args: &[&str]) -> Result<SelectionArgs, clap::Error> {
    let mut v = vec!["x"];
    v.extend(args);
    S::try_parse_from(v).map(|s| s.a)
}

/// One conflicting pair: the flags (in the order the rule table lists them) and the clap
/// wording both paths must report.
fn conflicts() -> Vec<(Vec<&'static str>, SelectionArgs, &'static str)> {
    let d = SelectionArgs::default;
    vec![
        (
            vec!["--collection", "a", "--all"],
            SelectionArgs {
                collection: vec!["a".into()],
                all: true,
                ..d()
            },
            "the argument '--collection <SEL>' cannot be used with '--all'",
        ),
        (
            vec!["--collection", "a", "--namespace", "refseq"],
            SelectionArgs {
                collection: vec!["a".into()],
                namespace: vec!["refseq".into()],
                ..d()
            },
            "the argument '--collection <SEL>' cannot be used with '--namespace <NS>'",
        ),
        (
            vec!["--collection", "a", "--taxon", "9606"],
            SelectionArgs {
                collection: vec!["a".into()],
                taxon: vec!["9606".into()],
                ..d()
            },
            "the argument '--collection <SEL>' cannot be used with '--taxon <TAXID|NAME>'",
        ),
        (
            vec!["--collections-file", "f", "--all"],
            SelectionArgs {
                collections_file: vec!["f".into()],
                all: true,
                ..d()
            },
            "the argument '--collections-file <FILE>' cannot be used with '--all'",
        ),
        (
            vec!["--collections-file", "f", "--namespace", "refseq"],
            SelectionArgs {
                collections_file: vec!["f".into()],
                namespace: vec!["refseq".into()],
                ..d()
            },
            "the argument '--collections-file <FILE>' cannot be used with '--namespace <NS>'",
        ),
        (
            vec!["--collections-file", "f", "--taxon", "9606"],
            SelectionArgs {
                collections_file: vec!["f".into()],
                taxon: vec!["9606".into()],
                ..d()
            },
            "the argument '--collections-file <FILE>' cannot be used with '--taxon <TAXID|NAME>'",
        ),
        (
            vec!["--namespace", "refseq", "--all"],
            SelectionArgs {
                namespace: vec!["refseq".into()],
                all: true,
                ..d()
            },
            "the argument '--namespace <NS>' cannot be used with '--all'",
        ),
        (
            vec!["--taxon", "9606", "--all"],
            SelectionArgs {
                taxon: vec!["9606".into()],
                all: true,
                ..d()
            },
            "the argument '--taxon <TAXID|NAME>' cannot be used with '--all'",
        ),
        (
            vec!["--all", "--max-n-fraction", "0.3", "--exclude-masked"],
            SelectionArgs {
                all: true,
                max_n_fraction: Some(0.3),
                exclude_masked: true,
                ..d()
            },
            "the argument '--max-n-fraction <F>' cannot be used with '--exclude-masked'",
        ),
    ]
}

#[test]
fn every_conflict_pair_same_wording_through_clap_and_tryfrom() {
    for (argv, direct, want) in conflicts() {
        let e = sel(&argv).unwrap_err();
        assert_eq!(
            e.kind(),
            clap::error::ErrorKind::ArgumentConflict,
            "{argv:?}"
        );
        let text = e.to_string();
        assert!(
            text.starts_with(&format!("error: {want}")),
            "{argv:?}: {text}"
        );
        let e = Selection::try_from(&direct).unwrap_err();
        assert_eq!(e.kind(), ErrorKind::InvalidInput);
        assert_eq!(e.to_string(), want, "{argv:?}");
    }
    // allowed combinations
    for ok in [
        vec!["--collection", "a", "--collection", "b", "--exclude-masked"],
        vec![
            "--namespace",
            "refseq",
            "--taxon",
            "9606",
            "--max-n-fraction",
            "0.5",
        ],
        vec!["--all", "--sequences-file", "f", "--exclude-masked"],
        vec!["--taxon", "9606", "--sequences-file", "f"],
    ] {
        sel(&ok).unwrap_or_else(|e| panic!("{ok:?}: {e}"));
    }
}

#[test]
fn every_violation_reported_at_once() {
    let a = SelectionArgs {
        collection: vec!["x".into()],
        namespace: vec!["refseq".into()],
        all: true,
        max_n_fraction: Some(0.3),
        exclude_masked: true,
        ..Default::default()
    };
    let e = Selection::try_from(&a).unwrap_err().to_string();
    assert_eq!(
        e,
        "the argument '--collection <SEL>' cannot be used with '--namespace <NS>'; \
         the argument '--collection <SEL>' cannot be used with '--all'; \
         the argument '--namespace <NS>' cannot be used with '--all'; \
         the argument '--max-n-fraction <F>' cannot be used with '--exclude-masked'"
    );
}

#[test]
fn direct_struct_all_plus_collection_is_invalid() {
    let a = SelectionArgs {
        all: true,
        collection: vec!["x".into()],
        ..Default::default()
    };
    assert!(matches!(
        Selection::try_from(&a),
        Err(Error::InvalidOptions(_))
    ));
}

#[test]
fn no_mode_and_bad_fraction() {
    let e = Selection::try_from(&SelectionArgs::default()).unwrap_err();
    assert_eq!(
        e.to_string(),
        "choose collections with --collection, --collections-file, --namespace, --taxon or --all"
    );
    let e = sel(&["--all", "--max-n-fraction", "0"])
        .unwrap_err()
        .to_string();
    assert!(
        e.contains(
            "invalid value '0' for '--max-n-fraction <F>': must be greater than 0 and at most 1"
        ),
        "{e}"
    );
    let a = SelectionArgs {
        all: true,
        max_n_fraction: Some(1.5),
        ..Default::default()
    };
    assert_eq!(
        Selection::try_from(&a).unwrap_err().to_string(),
        "invalid value '1.5' for '--max-n-fraction <F>': must be greater than 0 and at most 1"
    );
    // --exclude-masked is exactly 0.25
    let s = Selection::try_from(&sel(&["--all", "--exclude-masked"]).unwrap()).unwrap();
    assert_eq!(s.max_n_fraction, Some(0.25));
}

/// Parse with clap, then convert; either error as text.
fn view(args: &[&str]) -> Result<ViewRequest, String> {
    let a = V::try_parse_from([&["x"], args].concat()).map_err(|e| e.to_string())?;
    ViewRequest::try_from(&a.a).map_err(|e| e.to_string())
}

fn nstats(args: &[&str]) -> Result<NstatsRequest, String> {
    let a =
        N::try_parse_from([&["x", "--store", "/s"], args].concat()).map_err(|e| e.to_string())?;
    NstatsRequest::try_from(&a.a).map_err(|e| e.to_string())
}

#[test]
fn view_nstats_rules() {
    // exactly one kind (clap's group)
    for args in [
        &["--db", "d"][..],
        &["--db", "d", "--per-collection", "--name", "n"],
    ] {
        assert!(view(args).is_err(), "{args:?}");
    }
    // an all-default selection: allowed for per-collection views
    assert_eq!(
        view(&["--db", "d", "--per-collection"]).unwrap().selection,
        None
    );
    let e = view(&["--db", "d", "--name", "n"]).unwrap_err();
    assert!(e.contains("--name requires a selection"), "{e}");
    let e = view(&[
        "--db",
        "d",
        "--per-collection",
        "--all",
        "--collection",
        "x",
    ])
    .unwrap_err();
    assert!(e.contains("cannot be used with"), "{e}");
    let e = view(&["--db", "d", "--name", "../x", "--all"]).unwrap_err();
    assert!(e.contains("may hold only letters"), "{e}");

    // nstats: a mode is required; --resume only with a remote store
    let e = nstats(&[]).unwrap_err();
    assert!(e.contains("choose collections with"), "{e}");
    let e = nstats(&["--resume", "--all"]).unwrap_err();
    assert!(e.contains("--resume only applies to remote stores"), "{e}");
}

// ---------------------------------------------------------------------------
// resolution
// ---------------------------------------------------------------------------

#[allow(clippy::result_large_err)]
fn resolve(s: &Selection) -> Result<ResolvedSelection, Error> {
    let mut st = OpenedStore::open(&taxa_encoded().spec(), None).unwrap();
    gtars_rgblast::api::resolve_in_store(&mut st, s, &resolve_opts(), &mut |_| true).map(|r| r.0)
}

fn names(r: &ResolvedSelection) -> Vec<String> {
    let f = taxa_encoded();
    let mut v: Vec<String> = r
        .collections
        .iter()
        .map(|d| {
            f.digests
                .iter()
                .find(|(_, x)| *x == d)
                .map(|(k, _)| k.clone())
                .unwrap()
        })
        .collect();
    v.sort();
    v
}

#[test]
fn explicit_selectors() {
    let f = taxa_encoded();
    let ha = f.digest("human_a");
    // digest, NS:ALIAS, bare alias, and a duplicate are kept once; output sorted bytewise
    let r = resolve(&Selection::explicit(&[
        ha.as_str(),
        "name:human_b",
        "mouse_a",
        "human_a",
        "ucsc:hg_test",
        "insdc:GCA_TEST.1",
    ]))
    .unwrap();
    assert_eq!(names(&r), vec!["human_a", "human_b", "mouse_a"]);
    let mut sorted = r.collections.clone();
    sorted.sort();
    assert_eq!(r.collections, sorted);
    assert!(!r.include_orphans);
    // NS:ALIAS disambiguates an alias two namespaces use
    assert_eq!(
        names(&resolve(&Selection::explicit(&["ucsc:dup"])).unwrap()),
        vec!["human_a"]
    );
    assert_eq!(
        names(&resolve(&Selection::explicit(&["refseq:dup"])).unwrap()),
        vec!["mouse_a"]
    );
    let e = resolve(&Selection::explicit(&["dup"])).unwrap_err();
    assert!(matches!(e, Error::InvalidOptions(_)), "{e:?}");
    assert!(
        e.to_string()
            .starts_with("collection alias 'dup' names more than one collection"),
        "{e}"
    );
    // unknown alias in a known namespace, unknown namespace, unknown bare alias
    let e = resolve(&Selection::explicit(&["name:nope"])).unwrap_err();
    assert!(matches!(e, Error::UnknownCollection { .. }), "{e:?}");
    assert!(e.to_string().contains("'name:nope'"), "{e}");
    let e = resolve(&Selection::explicit(&["ensembl:x"])).unwrap_err();
    assert!(matches!(&e, Error::NotFound(m) if m.contains("'ensembl'") && m.contains("ucsc")));
    assert_eq!(e.kind(), ErrorKind::NotFound);
    let e = resolve(&Selection::explicit(&["nope"])).unwrap_err();
    assert!(
        e.to_string()
            .contains("collection alias namespaces: insdc, name, refseq, ucsc"),
        "{e}"
    );
}

#[test]
fn attribute_and_or() {
    let r = |ns: &[&str], tx: &[&str]| names(&resolve(&Selection::attributes(ns, tx)).unwrap());
    assert_eq!(r(&["refseq"], &[]), vec!["mouse_a"]);
    assert_eq!(r(&["insdc"], &[]), vec!["human_b"]);
    // values OR'd
    assert_eq!(r(&["insdc", "refseq"], &[]), vec!["human_b", "mouse_a"]);
    assert_eq!(
        r(&[], &["9606"]),
        vec!["human_a", "human_b", "human_masked"]
    );
    assert_eq!(
        r(&[], &["9606", "10090"]),
        vec!["human_a", "human_b", "human_masked", "mouse_a"]
    );
    // names: case and '_' ignored
    assert_eq!(r(&[], &["mus_MUSCULUS"]), vec!["mouse_a"]);
    assert_eq!(
        r(&[], &["homo sapiens"]),
        vec!["human_a", "human_b", "human_masked"]
    );
    // namespace AND taxon
    assert_eq!(r(&["insdc"], &["9606"]), vec!["human_b"]);
    assert_eq!(r(&["ucsc", "refseq"], &["Homo_sapiens"]), vec!["human_a"]);
    // nothing matches both
    let e = resolve(&Selection::attributes(&["insdc"], &["10090"])).unwrap_err();
    assert!(
        e.to_string() == "the selection is empty (no collection matches)",
        "{e:?}"
    );
}

#[test]
fn unknown_taxon_lists_taxa_and_unknown_namespace() {
    let e = resolve(&Selection::attributes(&[], &["Danio rerio"])).unwrap_err();
    assert_eq!(
        e.to_string(),
        "unknown taxon 'Danio rerio'; taxa in this store: 9606 Homo sapiens (3 collections), 10090 Mus musculus (1 collection)"
    );
    assert!(matches!(
        resolve(&Selection::attributes(&[], &["7955"])),
        Err(Error::NotFound(_))
    ));
    let e = resolve(&Selection::attributes(&["ensembl"], &[])).unwrap_err();
    assert!(matches!(e, Error::NotFound(_)), "{e:?}");
}

#[test]
fn n_filter_excludes_or_refuses() {
    // all mode: human_masked (198 N of 580 bp) leaves, the rest stay
    let r = resolve(&Selection::all().with_max_n_fraction(0.25)).unwrap();
    assert_eq!(names(&r), vec!["human_a", "human_b", "mouse_a", "nofhr"]);
    assert_eq!(r.excluded.len(), 1);
    // excluded sequences must not return as orphans
    assert!(!r.include_orphans);
    assert_eq!(r.excluded[0].display, "human_masked");
    assert!((r.excluded[0].n_fraction - 198.0 / 580.0).abs() < 1e-12);
    assert!(
        r.summary(5)
            .contains("1 excluded: N \u{2265} 0.25, e.g. human_masked 34.1%"),
        "{}",
        r.summary(5)
    );
    // attribute mode too
    let r = resolve(&Selection::attributes(&[], &["9606"]).with_max_n_fraction(0.25)).unwrap();
    assert_eq!(names(&r), vec!["human_a", "human_b"]);
    // a higher threshold keeps it
    let r = resolve(&Selection::all().with_max_n_fraction(0.5)).unwrap();
    assert_eq!(r.collections.len(), 5);
    // explicit mode: an error naming the collection
    let e = resolve(&Selection::explicit(&["human_a", "human_masked"]).with_max_n_fraction(0.25))
        .unwrap_err();
    assert_eq!(
        e.to_string(),
        "name:human_masked is 34.1% N, at or above --max-n-fraction 0.25; remove it or drop the N filter"
    );
    assert_eq!(e.kind(), ErrorKind::InvalidInput);
    // everything excluded
    let e = resolve(&Selection::explicit(&["human_masked"]).with_max_n_fraction(0.3)).unwrap_err();
    assert!(
        e.to_string().contains("at or above --max-n-fraction"),
        "{e}"
    );
    // every candidate excluded: human_masked alone under its own taxid (a --taxid-map entry)
    let f = taxa_encoded();
    let mut st = OpenedStore::open(&f.spec(), None).unwrap();
    let d = tempdir();
    let map = d.path().join("map.tsv");
    std::fs::write(&map, "human_masked\t63221\tHomo sapiens neanderthalensis\n").unwrap();
    let mut ro = resolve_opts();
    ro.taxid_map_file = Some(map);
    let e = gtars_rgblast::api::resolve_in_store(
        &mut st,
        &Selection::attributes(&[], &["63221"]).with_max_n_fraction(0.25),
        &ro,
        &mut |_| true,
    )
    .unwrap_err();
    assert_eq!(
        e.to_string(),
        "the selection is empty: 1 collection excluded (N fraction \u{2265} 0.25)"
    );
    assert_eq!(e.kind(), ErrorKind::InvalidInput);
}

#[test]
fn n_fraction_identical_across_storage_modes() {
    let mut got = Vec::new();
    for mode in [StorageMode::Encoded, StorageMode::Raw, StorageMode::Zstd] {
        let f = taxa(mode);
        let mut st = OpenedStore::open(&f.spec(), None).unwrap();
        let (r, _) = gtars_rgblast::api::resolve_in_store(
            &mut st,
            &Selection::all().with_max_n_fraction(1.0),
            &resolve_opts(),
            &mut |_| true,
        )
        .unwrap();
        got.push(r.info.values().map(|i| i.n_fraction).collect::<Vec<_>>());
    }
    assert_eq!(got[0], got[1]);
    assert_eq!(got[0], got[2]);
}

#[test]
fn sequences_file_rules() {
    let f = taxa_encoded();
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    let chr1 = plan
        .subjects
        .iter()
        .find(|s| s.occurrences.iter().any(|o| o.name == "chr1"))
        .unwrap()
        .digest
        .clone();
    let u1 = plan
        .subjects
        .iter()
        .find(|s| s.occurrences.iter().any(|o| o.name == "u1"))
        .unwrap()
        .digest
        .clone();
    // SQ. optional, deduplicated, sorted
    let r = resolve(
        &Selection::explicit(&["human_a"]).with_sequences(&[format!("SQ.{chr1}"), chr1.clone()]),
    )
    .unwrap();
    assert_eq!(r.sequences, Some(vec![chr1.clone()]));
    // a digest no selected collection holds
    let e = resolve(&Selection::explicit(&["human_a"]).with_sequences(std::slice::from_ref(&u1)))
        .unwrap_err();
    assert_eq!(
        e.to_string(),
        format!("--sequences-file: SQ.{u1} is in no selected collection")
    );
    // orphans count as selected under --all
    let r = resolve(&Selection::all().with_sequences(std::slice::from_ref(&f.orphan))).unwrap();
    assert_eq!(r.sequences, Some(vec![f.orphan.clone()]));
    let e =
        resolve(&Selection::explicit(&["nofhr"]).with_sequences(std::slice::from_ref(&f.orphan)))
            .unwrap_err();
    assert!(
        e.to_string().contains("is in no selected collection"),
        "{e}"
    );
    // through a file and the plan: subjects are only the listed ones, with every defline
    let d = tempdir();
    let list = d.path().join("seqs.txt");
    std::fs::write(&list, format!("# chr1\nSQ.{chr1}\textra\n\n")).unwrap();
    let a = SelectionArgs {
        all: true,
        sequences_file: Some(list),
        ..Default::default()
    };
    let s = Selection::try_from(&a).unwrap();
    assert_eq!(s.sequences, Some(vec![chr1.clone()]));
    let (_, p) = plan_for(&f.spec(), &s);
    assert_eq!(p.subjects.len(), 1);
    assert_eq!(
        p.subjects[0].occurrences.len(),
        2,
        "human_a and human_b chr1"
    );
}

#[test]
fn collections_file_and_canonical_json() {
    let f = taxa_encoded();
    let d = tempdir();
    let list = d.path().join("c.txt");
    std::fs::write(&list, "# my list\nhuman_a  first\nrefseq:dup\n\n").unwrap();
    let a = SelectionArgs {
        collections_file: vec![list],
        collection: vec![f.digest("human_b")],
        ..Default::default()
    };
    let s = Selection::try_from(&a).unwrap();
    match &s.set {
        CollectionSet::Explicit(v) => assert_eq!(v.len(), 3),
        o => panic!("{o:?}"),
    }
    let r = resolve(&s).unwrap();
    assert_eq!(names(&r), vec!["human_a", "human_b", "mouse_a"]);
    // the canonical form stores digests, so the same set spelled by digest has the same digest
    let ds: Vec<String> = ["human_a", "human_b", "mouse_a"]
        .iter()
        .map(|n| f.digest(n))
        .collect();
    let r2 = resolve(&Selection::explicit(&ds)).unwrap();
    assert_eq!(r.canonical, r2.canonical);
    assert_eq!(r.digest, r2.digest);
    assert_eq!(r.canonical["mode"]["explicit"].as_array().unwrap().len(), 3);
    // a taxon given by name is stored as its taxid
    let t1 = resolve(&Selection::attributes(&[], &["Homo sapiens"])).unwrap();
    let t2 = resolve(&Selection::attributes(&[], &["9606"])).unwrap();
    assert_eq!(t1.digest, t2.digest);
    // warnings name the collection with no FHR taxon
    let all = resolve(&Selection::all()).unwrap();
    assert!(all.include_orphans);
    assert_eq!(
        all.warnings,
        vec!["1 collection has no FHR taxon (taxid 0): nofhr".to_string()]
    );
}
