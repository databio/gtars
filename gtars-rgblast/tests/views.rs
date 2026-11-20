//! `rgmakeblastdb view`: files, statistics, names, the index and the manifest record.

mod common;

use common::*;
use gtars_rgblast::api::{ViewRequest, view};
use gtars_rgblast::blastdb::alias;
use gtars_rgblast::blastdb::output::with_ext;
use gtars_rgblast::blastdb::sidecars::{self, ViewKind};
use gtars_rgblast::select::Selection;
use gtars_rgblast::views::{DbCatalog, read_index};
use gtars_rgblast::{Error, ErrorKind};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

fn taxa_db(dir: &Path) -> PathBuf {
    let out = dir.join("db");
    build(&taxa_encoded().spec(), Selection::all(), &out);
    out
}

fn req(db: &Path, kind: ViewKind) -> ViewRequest {
    ViewRequest {
        db: db.to_path_buf(),
        kind,
        name: None,
        selection: None,
        out_dir: None,
        title: None,
        full_db_stats: false,
    }
}

/// (NSEQ, LENGTH) a view of these collections must have: distinct OIDs and their lengths.
fn expected(db: &Path, colls: &[String]) -> (u64, u64) {
    let cat = DbCatalog::open(db).unwrap();
    let mut oids = BTreeSet::new();
    for c in colls {
        for m in cat.members_of(c) {
            oids.insert(m.oid);
        }
    }
    let len = oids.iter().map(|&o| cat.oids[o as usize].length).sum();
    (oids.len() as u64, len)
}

fn nal_value(nal: &Path, key: &str) -> Option<String> {
    std::fs::read_to_string(nal)
        .unwrap()
        .lines()
        .find_map(|l| l.strip_prefix(&format!("{key} ")).map(String::from))
}

#[test]
fn per_collection_views() {
    let f = taxa_encoded();
    let d = tempdir();
    let db = taxa_db(d.path());
    let s = view(&req(&db, ViewKind::PerCollection), &mut |_| true).unwrap();
    assert_eq!(s.out_dir, d.path().join("db.views"));
    let names: Vec<&str> = s.views.iter().map(|v| v.name.as_str()).collect();
    assert_eq!(
        names,
        vec!["human_a", "human_b", "human_masked", "mouse_a", "nofhr"]
    );
    for v in &s.views {
        let (nseq, len) = expected(&db, &v.collections);
        assert_eq!((v.nseq, v.length), (nseq, len), "{}", v.name);
        assert_eq!(nal_value(&v.nal, "NSEQ").unwrap(), nseq.to_string());
        assert_eq!(nal_value(&v.nal, "LENGTH").unwrap(), len.to_string());
        assert_eq!(nal_value(&v.nal, "DBLIST").unwrap(), "../db");
        assert_eq!(
            nal_value(&v.nal, "SEQIDLIST").unwrap(),
            format!("{}.bsl", v.name)
        );
        assert!(nal_value(&v.nal, "STATS_NSEQ").is_none());
        // the list holds <collection digest>:<id_name> ids only, never SQ.
        let b = alias::read_bsl(&v.bsl).unwrap();
        assert_eq!(b.title, v.name);
        assert!(
            b.ids
                .iter()
                .all(|i| i.starts_with(&format!("{}:", v.collections[0])))
        );
        assert!(!b.ids.iter().any(|i| i.starts_with("SQ.")));
    }
    let ha = &s.views[0];
    assert_eq!(ha.collections, vec![f.digest("human_a")]);
    assert_eq!(ha.taxid, Some(9606));
    assert_eq!((ha.nseq, ha.length), (4, 300 + 280 + 260 + 200));
    assert_eq!(
        nal_value(&ha.nal, "TITLE").unwrap(),
        format!("{TITLE}: Homo sapiens human_a")
    );
    assert_eq!(
        nal_value(&s.views[4].nal, "TITLE").unwrap(),
        format!("{TITLE}: nofhr")
    );
    // views.tsv: one row per view
    let idx = read_index(&s.index).unwrap();
    assert_eq!(idx.len(), 5);
    assert_eq!(
        idx[0][..2],
        ["human_a".to_string(), "per-collection".to_string()]
    );
}

#[test]
fn named_views_share_the_index() {
    let f = taxa_encoded();
    let d = tempdir();
    let db = taxa_db(d.path());
    let s = view(&req(&db, ViewKind::PerCollection), &mut |_| true).unwrap();
    // a named view: an unmasked subset with no rebuild, N fractions from collections.tsv
    let mut r = req(&db, ViewKind::Named);
    r.name = Some("unmasked".into());
    r.selection = Some(Selection::all().with_max_n_fraction(0.25));
    let n = view(&r, &mut |_| true).unwrap();
    assert_eq!(n.views.len(), 1);
    let keep: Vec<String> = ["human_a", "human_b", "mouse_a", "nofhr"]
        .iter()
        .map(|x| f.digest(x))
        .collect();
    assert_eq!((n.views[0].nseq, n.views[0].length), expected(&db, &keep));
    assert_eq!(n.views[0].taxid, None, "two taxa");
    // both batches are in one index
    let idx = read_index(&s.index).unwrap();
    assert_eq!(idx.len(), 5 + 1);
    // rerunning a kind replaces its rows only
    view(&req(&db, ViewKind::PerCollection), &mut |_| true).unwrap();
    assert_eq!(read_index(&s.index).unwrap().len(), 6);
    // a selection limits per-collection views; --sequences-file limits a named view's rows
    let mut r = req(&db, ViewKind::PerCollection);
    r.selection = Some(Selection::attributes(&["refseq"], &[]));
    r.out_dir = Some(d.path().join("sub"));
    let s = view(&r, &mut |_| true).unwrap();
    assert_eq!(s.views.len(), 1);
    assert_eq!(nal_value(&s.views[0].nal, "DBLIST").unwrap(), "../db");
    let cat = DbCatalog::open(&db).unwrap();
    let chr1 = cat
        .members_of(&f.digest("human_a"))
        .into_iter()
        .find(|m| m.name == "chr1")
        .unwrap()
        .digest
        .clone();
    let mut r = req(&db, ViewKind::Named);
    r.name = Some("chr1".into());
    r.selection = Some(Selection::attributes(&[], &["9606"]).with_sequences(&[chr1]));
    let s = view(&r, &mut |_| true).unwrap();
    assert_eq!((s.views[0].nseq, s.views[0].length), (1, 300));
    let b = alias::read_bsl(&s.views[0].bsl).unwrap();
    assert_eq!(b.ids.len(), 2, "human_a:chr1 and human_b:chr1: {:?}", b.ids);
}

#[test]
fn full_db_stats_lines() {
    let d = tempdir();
    let db = taxa_db(d.path());
    let m = sidecars::read_manifest(&db).unwrap();
    let mut r = req(&db, ViewKind::PerCollection);
    r.full_db_stats = true;
    let s = view(&r, &mut |_| true).unwrap();
    let nal = &s.views[0].nal;
    assert_eq!(
        nal_value(nal, "STATS_NSEQ").unwrap(),
        m.totals.num_seqs.to_string()
    );
    assert_eq!(
        nal_value(nal, "STATS_TOTLEN").unwrap(),
        m.totals.total_length.to_string()
    );
    assert_eq!(nal_value(nal, "NSEQ").unwrap(), s.views[0].nseq.to_string());
    assert!(m.totals.num_seqs as u64 > s.views[0].nseq);
}

#[test]
fn view_errors() {
    let d = tempdir();
    let db = taxa_db(d.path());
    // a collection that is not in the database
    let small = d.path().join("small");
    build(
        &taxa_encoded().spec(),
        Selection::explicit(&["human_a"]),
        &small,
    );
    let mut r = req(&small, ViewKind::PerCollection);
    r.selection = Some(Selection::explicit(&["mouse_a"]));
    let e = view(&r, &mut |_| true).unwrap_err();
    assert_eq!(
        e.to_string(),
        "collection mouse_a is not in database small; to add it, rebuild the database with rgmakeblastdb build --force"
    );
    // a database without a manifest
    let bare = d.path().join("bare");
    std::fs::copy(with_ext(&db, "nin"), with_ext(&bare, "nin")).unwrap();
    let e = view(&req(&bare, ViewKind::PerCollection), &mut |_| true).unwrap_err();
    assert!(
        matches!(&e, Error::NotFound(m) if m.contains("not a database written by rgmakeblastdb")),
        "{e:?}"
    );
    assert_eq!(e.kind(), ErrorKind::NotFound);
    assert!(e.to_string().contains("rgmakeblastdb build --force"));
    // a named view needs a name and a selection
    let r = req(&db, ViewKind::Named);
    assert!(matches!(
        view(&r, &mut |_| true),
        Err(Error::InvalidOptions(_))
    ));
}

#[test]
fn view_file_name_collisions() {
    use gtars_refget::store::{FastaImportOptions, RefgetStore};
    let d = tempdir();
    let path = d.path().join("st");
    let mut s = RefgetStore::on_disk(&path).unwrap();
    s.set_quiet(true);
    let src = d.path().join("src");
    fixture_copies(&src, &["asmB.fa", "asmD.fa"]);
    let (b, _) = s
        .add_sequence_collection_from_fasta(src.join("asmB.fa"), FastaImportOptions::new())
        .unwrap();
    let (dd, _) = s
        .add_sequence_collection_from_fasta(src.join("asmD.fa"), FastaImportOptions::new())
        .unwrap();
    // "a:b" and "a_b" are distinct db_tags but the same file name
    s.add_collection_alias("name", "a:b", &b.digest).unwrap();
    s.add_collection_alias("name", "a_b", &dd.digest).unwrap();
    s.write().unwrap();
    drop(s);
    let db = d.path().join("db");
    build(&path.display().to_string(), Selection::all(), &db);
    let v = view(&req(&db, ViewKind::PerCollection), &mut |_| true).unwrap();
    let names: Vec<&str> = v.views.iter().map(|x| x.name.as_str()).collect();
    // the second gets its digest's first 8 characters, made safe (a leading '-' becomes '_')
    assert!(dd.digest.starts_with('-'));
    assert_eq!(
        names,
        vec!["a_b".to_string(), format!("a_b._{}", &dd.digest[1..8])]
    );
}
