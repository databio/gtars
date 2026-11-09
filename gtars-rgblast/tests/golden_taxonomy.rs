//! Byte-identity canary for taxonomy and seqid lists, pinned to BLAST+ 2.17.0 (feature
//! `golden`): `-taxid_map` parity, makeblastdb as an oracle for the taxonomy files of plans
//! `-taxid_map` can't express, and `.bsl` parity with `blastdb_aliastool`.

mod common;

use common::blast::Blast;
use common::*;
use gtars_rgblast::blastdb::output::with_ext;
use gtars_rgblast::blastdb::reader::{DbReader, dump_lmdb, parse_nin};
use gtars_rgblast::blastdb::{alias, asn, njs, oidfiles};
use gtars_rgblast::plan::{Defline, SubjectPlan};
use gtars_rgblast::reference_fasta;
use gtars_rgblast::select::Selection;
use gtars_rgblast::store::StoreFiles;
use std::path::Path;

fn read(p: &Path) -> Vec<u8> {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

#[test]
fn golden_taxonomy_requires_blast_2_17() {
    let d = tempdir();
    assert_eq!(Blast::find().version(d.path()), "2.17.0");
}

#[test]
fn taxid_map_parity() {
    let f = taxa_encoded();
    let blast = Blast::find();
    let d = tempdir();
    let work = d.path();
    // one taxon only, so every SQ. defline agrees and every bare tag has one taxid
    let sel = Selection::explicit(&["human_a", "human_b"]);
    build(&f.spec(), sel.clone(), &work.join("ours/db"));
    let (st, plan) = plan_for(&f.spec(), &sel);
    let fetch = StoreFiles::new(&st);
    reference_fasta::write(&plan, &fetch, &work.join("ref.fa")).unwrap();
    reference_fasta::write_taxid_map(&plan, &work.join("map.txt")).unwrap();
    let map = std::fs::read_to_string(work.join("map.txt")).unwrap();
    assert!(
        map.lines()
            .any(|l| l.starts_with("SQ.") && l.ends_with(" 9606")),
        "{map}"
    );
    assert!(map.lines().any(|l| l == "chr1 9606"), "{map}");
    std::fs::create_dir_all(work.join("mk")).unwrap();
    blast.ok(
        work,
        "makeblastdb",
        &[
            "-in",
            "ref.fa",
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-blastdb_version",
            "5",
            "-title",
            TITLE,
            "-taxid_map",
            "map.txt",
            "-out",
            "mk/db",
        ],
    );
    let (o, m) = (work.join("ours/db"), work.join("mk/db"));
    for ext in ["nhr", "nsq", "not", "nto", "nos"] {
        assert!(
            read(&with_ext(&o, ext)) == read(&with_ext(&m, ext)),
            ".{ext} differs"
        );
    }
    for ext in ["ndb", "ntf"] {
        assert_eq!(
            dump_lmdb(&with_ext(&o, ext)).unwrap(),
            dump_lmdb(&with_ext(&m, ext)).unwrap(),
            ".{ext}"
        );
    }
    let (mut ho, ..) = parse_nin(&with_ext(&o, "nin")).unwrap();
    let (mut hm, ..) = parse_nin(&with_ext(&m, "nin")).unwrap();
    ho.date.clear();
    hm.date.clear();
    assert_eq!(ho, hm);
    let (jo, jm) = (
        njs::read(&with_ext(&o, "njs")).unwrap(),
        njs::read(&with_ext(&m, "njs")).unwrap(),
    );
    assert_eq!(jo.number_of_taxids, Some(1));
    assert_eq!(jo.number_of_taxids, jm.number_of_taxids);
    // plans -taxid_map can't express: the decoy's mixed SQ. defline, chr1 in human and mouse
    let (_, all) = plan_for(&f.spec(), &Selection::all());
    let e = reference_fasta::write_taxid_map(&all, &work.join("x.txt")).unwrap_err();
    assert!(e.to_string().contains("can't express"), "{e}");
    let (_, mixed_tag) = plan_for(&f.spec(), &Selection::explicit(&["human_a", "mouse_a"]));
    let e = reference_fasta::write_taxid_map(&mixed_tag, &work.join("y.txt")).unwrap_err();
    assert!(e.to_string().contains("'chr1'"), "{e}");
}

/// Our deflines with every omitted taxid turned into 0, as makeblastdb writes them.
fn with_zero(plan: &SubjectPlan, oid: usize) -> Vec<u8> {
    let d: Vec<Defline> = plan.subjects[oid]
        .deflines
        .iter()
        .map(|x| Defline {
            taxid: Some(x.taxid.unwrap_or(0)),
            ..x.clone()
        })
        .collect();
    let mut b = Vec::new();
    asn::encode(&d, &mut b);
    b
}

#[test]
fn taxonomy_oracle_input_blastdb() {
    let f = taxa_encoded();
    let blast = Blast::find();
    let d = tempdir();
    let work = d.path();
    build(&f.spec(), Selection::all(), &work.join("ours/db"));
    std::fs::create_dir_all(work.join("mk")).unwrap();
    blast.ok(
        work,
        "makeblastdb",
        &[
            "-in",
            "ours/db",
            "-input_type",
            "blastdb",
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-blastdb_version",
            "5",
            "-title",
            TITLE,
            "-out",
            "mk/db",
        ],
    );
    let (o, m) = (work.join("ours/db"), work.join("mk/db"));
    assert!(
        read(&with_ext(&o, "nsq")) == read(&with_ext(&m, "nsq")),
        ".nsq differs"
    );
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    let rm = DbReader::open(&m).unwrap();
    let mut mixed = 0;
    for oid in 0..plan.subjects.len() {
        assert_eq!(
            rm.header_bytes(oid as u32).unwrap(),
            with_zero(&plan, oid),
            "OID {oid}"
        );
        if plan.subjects[oid].deflines[0].taxid.is_none() {
            mixed += 1;
        }
    }
    assert_eq!(mixed, 1, "the decoy");
    // NCBI's taxonomy files equal ours once each mixed OID also carries 0
    let sets: Vec<Vec<i32>> = plan
        .subjects
        .iter()
        .map(|s| {
            let mut v = gtars_rgblast::plan::taxonomy::taxid_set(&s.deflines);
            if s.deflines[0].taxid.is_none() {
                v.push(0);
                v.sort();
            }
            v
        })
        .collect();
    assert!(
        oidfiles::not_bytes(&sets) == read(&with_ext(&m, "not")),
        ".not"
    );
    let (nto, offsets) = oidfiles::nto_bytes(&sets);
    assert!(nto == read(&with_ext(&m, "nto")), ".nto");
    assert_eq!(
        gtars_rgblast::blastdb::lmdb::read_ntf(&with_ext(&m, "ntf")).unwrap(),
        offsets,
        ".ntf"
    );
    // and ours are exactly the same files without that 0
    let ours: Vec<Vec<i32>> = oidfiles::taxid_sets(&plan.subjects);
    assert!(oidfiles::not_bytes(&ours) == read(&with_ext(&o, "not")));
    assert!(oidfiles::nto_bytes(&ours).0 == read(&with_ext(&o, "nto")));
}

#[test]
fn bsl_parity() {
    let blast = Blast::find();
    let d = tempdir();
    let work = d.path();
    let f = taxa_encoded();
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    let mut ids: Vec<String> = plan
        .subjects
        .iter()
        .map(|s| format!("SQ.{}", s.digest))
        .collect();
    for (oid, s) in plan.subjects.iter().enumerate() {
        for d in &s.deflines[1..] {
            ids.push(d.ids[0].fasta());
        }
        if oid == 0 {
            ids.push(format!("{}:{}", f.digest("mouse_a"), "7"));
        }
    }
    ids.push(format!("x:{}", "y".repeat(298)));
    std::fs::write(work.join("ids.txt"), ids.join("\n") + "\n").unwrap();
    blast.ok(
        work,
        "blastdb_aliastool",
        &[
            "-seqid_file_in",
            "ids.txt",
            "-seqid_title",
            "parity",
            "-seqid_file_out",
            "ref.bsl",
        ],
    );
    let reference = read(&work.join("ref.bsl"));
    let date = alias::parse_bsl(&reference).unwrap().date;
    alias::write_bsl_dated(&ids, &work.join("parity.bsl"), &date).unwrap();
    assert!(
        read(&work.join("parity.bsl")) == reference,
        "write_bsl differs from blastdb_aliastool"
    );
    // the tool reads ours
    alias::write_bsl(&ids, &work.join("ours.bsl")).unwrap();
    let info = blast.ok(work, "blastdb_aliastool", &["-seqid_file_info", "ours.bsl"]);
    assert!(
        info.contains(&format!(
            "Num of Ids: {}",
            alias::read_bsl(&work.join("ours.bsl")).unwrap().ids.len()
        )),
        "{info}"
    );
    assert!(info.contains("Title: ours"), "{info}");
}
