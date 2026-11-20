//! The functional contract on small committed fixtures: no network, docker or BLAST+.

mod common;

use common::*;
use gtars_refget::store::StorageMode;
use gtars_rgblast::blastdb::output::{OutputSet, with_ext};
use gtars_rgblast::blastdb::reader::DbReader;
use gtars_rgblast::blastdb::verify::{lookup_key, verify_db};
use gtars_rgblast::blastdb::{self, BlastDate, WriteOptions, njs, sidecars};
use gtars_rgblast::plan::{OnOversize, PlanOptions, SeqId, SkipReason, SubjectPlan};
use gtars_rgblast::select::Selection;
use gtars_rgblast::store::fetch::{SeqFetch, StoredBytes};
use gtars_rgblast::store::{OpenedStore, StoreFiles};
use gtars_rgblast::{Error, ProgressEvent, SUBJECT_PLAN_VERSION};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::atomic::Ordering;

fn all() -> Selection {
    Selection::all()
}

fn colls(keys: &[&str]) -> Selection {
    Selection::explicit(keys)
}

/// A manifest's JSON without its `store` field.
fn manifest_sans_store(b: &[u8]) -> serde_json::Value {
    let mut v: serde_json::Value = serde_json::from_slice(b).unwrap();
    v.as_object_mut().unwrap().remove("store");
    v
}

/// Build the whole encoded fixture store into `<tmp>/db` once per test that asks.
fn build_all(dir: &Path) -> std::path::PathBuf {
    let out = dir.join("db");
    build(&encoded().spec(), all(), &out);
    out
}

// ---------------------------------------------------------------------------
// storage modes and remote
// ---------------------------------------------------------------------------

#[test]
fn db_identical_across_storage_modes() {
    let mut outputs = Vec::new();
    for mode in [StorageMode::Encoded, StorageMode::Raw, StorageMode::Zstd] {
        let f = fixture(mode);
        let d = tempdir();
        let s = build(&f.spec(), all(), &d.path().join("db"));
        assert_eq!(s.preflight.mode, format!("{mode:?}"));
        outputs.push((mode, files_in(d.path()), d));
    }
    let (_, first, _) = &outputs[0];
    assert!(
        first.len() >= 16,
        "{:?}",
        first.iter().map(|f| &f.0).collect::<Vec<_>>()
    );
    for (mode, files, _) in &outputs[1..] {
        assert_eq!(files.len(), first.len(), "{mode:?}");
        for ((n1, b1), (n2, b2)) in first.iter().zip(files) {
            assert_eq!(n1, n2);
            if n1.ends_with(".manifest.json") {
                // records where the store is (and its digests)
                assert_eq!(manifest_sans_store(b1), manifest_sans_store(b2), "{mode:?}");
                continue;
            }
            assert!(b1 == b2, "{mode:?}: {n1} differs from the Encoded build");
        }
    }
}

#[test]
fn remote_build_equals_local_build() {
    let f = encoded();
    let server = http::serve(f.path.parent().unwrap());
    let url = format!("{}/fixture", server.url);
    let local = tempdir();
    build(&f.spec(), all(), &local.path().join("db"));
    let remote = tempdir();
    let cache = tempdir();
    let mut req = request(&url, all(), &remote.path().join("db"));
    req.cache_dir = Some(cache.path().join("c"));
    req.options.remote.fetch_jobs = 4;
    let s = build_req(&req);
    assert!(s.preflight.remote);
    let net = s.network.expect("network report");
    assert_eq!(net.sequences_fetched, s.preflight.subjects);
    same_database(local.path(), remote.path());
    // the selection was pulled into the cache and read from there
    assert!(!no_files_under(&cache.path().join("c/sequences")));
}

/// Every database file identical; manifests equal apart from where the store is.
fn same_database(a: &Path, b: &Path) {
    let (a, b) = (files_in(a), files_in(b));
    assert_eq!(a.len(), b.len());
    for ((n1, b1), (n2, b2)) in a.iter().zip(&b) {
        assert_eq!(n1, n2);
        if n1.ends_with(".manifest.json") {
            let (ma, mb) = (manifest_sans_store(b1), manifest_sans_store(b2));
            assert_eq!(ma, mb, "manifest beyond the store location");
            continue;
        }
        assert!(b1 == b2, "{n1} differs between the two builds");
    }
}

/// True when `dir` holds no files (it may be missing or hold empty folders).
fn no_files_under(dir: &Path) -> bool {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return true;
    };
    rd.flatten()
        .all(|e| e.path().is_dir() && no_files_under(&e.path()))
}

#[test]
fn remote_resume_downloads_nothing_twice() {
    let f = encoded();
    let server = http::serve(f.path.parent().unwrap());
    let url = format!("{}/fixture", server.url);
    let local = tempdir();
    build(&f.spec(), all(), &local.path().join("db"));
    let cache = tempdir();
    let run = || {
        let out = tempdir();
        let mut req = request(&url, all(), &out.path().join("db"));
        req.cache_dir = Some(cache.path().to_path_buf());
        req.options.remote.fetch_jobs = 4;
        let s = build_req(&req);
        same_database(local.path(), out.path());
        s
    };
    let first = run();
    let n1 = first.network.expect("network report");
    assert!(n1.bytes > 0 && n1.requests > 0, "{n1:?}");
    assert_eq!(n1.sequences_fetched, first.preflight.subjects);
    assert_eq!(n1.cache_hits, 0);
    assert!(!no_files_under(&cache.path().join("sequences")));
    let before = server.requests.load(Ordering::SeqCst);
    let second = run();
    let n2 = second.network.expect("network report");
    assert_eq!(
        (n2.bytes, n2.requests, n2.sequences_fetched),
        (0, 0, 0),
        "the second build downloads no sequence: {n2:?}"
    );
    assert_eq!(n2.cache_hits, second.preflight.subjects);
    assert!(n2.line().contains(&format!("cache hits {}", n2.cache_hits)));
    // only metadata went over the wire (the manifest, FHR sidecars)
    let meta_requests = server.requests.load(Ordering::SeqCst) - before;
    assert!(
        meta_requests < second.preflight.subjects,
        "{meta_requests} requests on the cached rerun"
    );
}

#[test]
fn remote_missing_store_is_not_found_or_denied() {
    let f = encoded();
    let server = http::serve(f.path.parent().unwrap());
    let cache = tempdir();
    let e = OpenedStore::open(&format!("{}/no-such-store", server.url), Some(cache.path()))
        .err()
        .unwrap();
    assert!(matches!(&e, Error::NotFound(m) if m.contains("/no-such-store/rgstore.json")));
    assert!(
        e.to_string()
            .contains("not found or access denied (HTTP 403)"),
        "{e}"
    );
}

// ---------------------------------------------------------------------------
// the functional contract
// ---------------------------------------------------------------------------

#[test]
fn md5_of_every_entry() {
    let d = tempdir();
    let out = build_all(d.path());
    let (st, plan) = plan_for(&encoded().spec(), &all());
    let r = DbReader::open(&out).unwrap();
    assert_eq!(r.num_oids() as usize, plan.subjects.len());
    let fetch = StoreFiles::new(&st);
    for (oid, s) in plan.subjects.iter().enumerate() {
        let seq = r.sequence(oid as u32).unwrap();
        assert_eq!(seq.len(), s.length as usize);
        if md5_hex(&seq) != s.md5 {
            // BLAST has no U: the store sequence with U -> T must match
            let stored = fetch.read(&s.digest).unwrap();
            let ascii = gtars_rgblast::store::decode::to_ascii(
                st.mode(),
                s.alphabet,
                (*stored).as_ref(),
                s.length,
            )
            .unwrap();
            assert!(ascii.contains(&b'U'), "OID {oid}: md5 mismatch without U");
            let t: Vec<u8> = ascii
                .iter()
                .map(|&b| if b == b'U' { b'T' } else { b })
                .collect();
            assert_eq!(md5_hex(&seq), md5_hex(&t), "OID {oid}");
        }
    }
    // and through verify
    let rep = verify_db(&out, &st.store, &fetch, 2, &mut |_| true).unwrap();
    assert!(rep.ok(), "{rep:?}");
    assert_eq!(rep.sequences_checked as usize, plan.subjects.len());
}

#[test]
fn every_listed_id_resolves() {
    let d = tempdir();
    let out = build_all(d.path());
    let ids = sidecars::read_ids(&with_ext(&out, "ids.tsv")).unwrap();
    let ndb = DbReader::open(&out).unwrap().ndb().unwrap();
    assert!(ids.len() > 100);
    for (id, oid, _) in &ids {
        assert_eq!(ndb.lookup(&lookup_key(id)).unwrap(), vec![*oid], "{id}");
    }
}

#[test]
fn every_lmdb_key_resolves() {
    let d = tempdir();
    let out = build_all(d.path());
    let (_, plan) = plan_for(&encoded().spec(), &all());
    let ndb = DbReader::open(&out).unwrap().ndb().unwrap();
    let mut want: HashSet<(String, u32)> = HashSet::new();
    for (oid, s) in plan.subjects.iter().enumerate() {
        for dl in &s.deflines {
            for id in &dl.ids {
                let mut keys = Vec::new();
                id.lmdb_keys(&mut keys);
                for (k, _) in keys {
                    assert_eq!(ndb.lookup(&k).unwrap(), vec![oid as u32], "{k}");
                    want.insert((k, oid as u32));
                }
            }
        }
        // every lookup form, including the dash-leading db:tag ones left out of ids.tsv
        for f in plan.lookup_forms(oid as u32) {
            assert_eq!(
                ndb.lookup(&lookup_key(&f.id)).unwrap(),
                vec![oid as u32],
                "{}",
                f.id
            );
        }
    }
    let got: HashSet<(String, u32)> = ndb.all_keys().unwrap().into_iter().collect();
    assert_eq!(got, want, "acc2oid holds exactly the planned keys");
}

#[test]
fn nhr_deflines_match_plan() {
    let d = tempdir();
    let out = build_all(d.path());
    let (_, plan) = plan_for(&encoded().spec(), &all());
    let r = DbReader::open(&out).unwrap();
    for (oid, s) in plan.subjects.iter().enumerate() {
        assert_eq!(r.deflines(oid as u32).unwrap(), s.deflines, "OID {oid}");
    }
}

#[test]
fn dash_digest_ids() {
    let f = encoded();
    let dd = f.digest("asmD");
    assert!(
        dd.starts_with('-'),
        "fixture asmD digest {dd} must start with '-'"
    );
    let d = tempdir();
    let out = build_all(d.path());
    let (_, plan) = plan_for(&f.spec(), &all());
    let ci = plan
        .collections
        .iter()
        .position(|c| c.digest == dd)
        .unwrap();
    assert_eq!(plan.collections[ci].db_tag, "asmD");
    let (oid, s) = plan
        .subjects
        .iter()
        .enumerate()
        .find(|(_, s)| s.occurrences.iter().any(|o| o.collection as usize == ci))
        .unwrap();
    let name = &s.occurrences[0].name;
    let gnl = format!("gnl|{dd}|{name}");
    let colon = format!("{dd}:{name}");
    let ids = std::fs::read_to_string(with_ext(&out, "ids.tsv")).unwrap();
    assert!(
        ids.lines()
            .any(|l| l.starts_with(&format!("{gnl}\t{oid}\t"))),
        "ids.tsv lists {gnl}"
    );
    assert!(
        !ids.lines().any(|l| l.starts_with(&colon)),
        "ids.tsv must not list {colon}"
    );
    let ndb = DbReader::open(&out).unwrap().ndb().unwrap();
    assert_eq!(ndb.lookup(&lookup_key(&gnl)).unwrap(), vec![oid as u32]);
    assert_eq!(ndb.lookup(&colon).unwrap(), vec![oid as u32]);
    // a database of only asmD (by its dash-leading digest) works too
    let d2 = tempdir();
    let s2 = build(&f.spec(), colls(&[&dd]), &d2.path().join("db"));
    assert_eq!(s2.preflight.dash_leading_collections, 1);
    assert!(
        s2.warnings
            .iter()
            .any(|w| w.starts_with("1 collection digests start with '-'"))
    );
}

// ---------------------------------------------------------------------------
// the plan
// ---------------------------------------------------------------------------

fn plan_snapshot(plan: &SubjectPlan) -> String {
    let mut s = format!(
        "# SUBJECT_PLAN_VERSION={SUBJECT_PLAN_VERSION}\n#oid\tdigest\tdefline\tids\ttitle\n"
    );
    for (oid, sub) in plan.subjects.iter().enumerate() {
        for (i, d) in sub.deflines.iter().enumerate() {
            let ids: Vec<String> = d.ids.iter().map(|x| x.fasta()).collect();
            s += &format!(
                "{oid}\t{}\t{i}\t{}\t{}\n",
                sub.digest,
                ids.join(" "),
                d.title
            );
        }
    }
    s
}

#[test]
fn oid_order_is_canonical() {
    let (_, plan) = plan_for(&encoded().spec(), &all());
    let got = plan_snapshot(&plan);
    let path = Path::new(EXPECTED).join("plan_fixture.tsv");
    if std::env::var_os("GTARS_RGBLAST_UPDATE_EXPECTED").is_some() {
        std::fs::write(&path, &got).unwrap();
    }
    let want = std::fs::read_to_string(&path).expect("tests/data/expected/plan_fixture.tsv");
    let first = want.lines().next().unwrap();
    assert_eq!(
        first,
        format!("# SUBJECT_PLAN_VERSION={SUBJECT_PLAN_VERSION}"),
        "a deliberate plan change updates the snapshot and bumps SUBJECT_PLAN_VERSION together"
    );
    assert_eq!(
        got, want,
        "OID order, ids or titles changed; if deliberate, bump SUBJECT_PLAN_VERSION and update the snapshot (GTARS_RGBLAST_UPDATE_EXPECTED=1)"
    );
}

#[test]
fn first_id_invariants() {
    let d = tempdir();
    let out = build_all(d.path());
    let (_, plan) = plan_for(&encoded().spec(), &all());
    let r = DbReader::open(&out).unwrap();
    for (oid, s) in plan.subjects.iter().enumerate() {
        for deflines in [s.deflines.clone(), r.deflines(oid as u32).unwrap()] {
            assert_eq!(
                deflines[0].ids[0],
                SeqId::Local(format!("SQ.{}", s.digest)),
                "OID {oid}"
            );
            assert_eq!(deflines.len(), s.occurrences.len() + 1);
            for (dl, o) in deflines[1..].iter().zip(&s.occurrences) {
                let c = &plan.collections[o.collection as usize];
                assert_eq!(
                    dl.ids[0],
                    SeqId::General {
                        db: c.db_tag.clone(),
                        tag: o.id_name.clone()
                    },
                    "OID {oid}"
                );
            }
            // only lcl and gnl ids
            for dl in &deflines {
                for id in &dl.ids {
                    assert!(id.fasta().starts_with("lcl|") || id.fasta().starts_with("gnl|"));
                }
            }
        }
    }
}

#[test]
fn db_tag_independent_of_selection() {
    let f = encoded();
    let mut tags: HashMap<String, String> = HashMap::new();
    for sel in [colls(&["asmA"]), colls(&["asmA", "asmB"]), all()] {
        let (_, plan) = plan_for(&f.spec(), &sel);
        for c in &plan.collections {
            let prev = tags
                .entry(c.digest.clone())
                .or_insert_with(|| c.db_tag.clone());
            assert_eq!(
                *prev, c.db_tag,
                "db_tag of {} changed with the selection",
                c.digest
            );
        }
    }
    let upper: HashSet<String> = tags.values().map(|t| t.to_uppercase()).collect();
    assert_eq!(
        upper.len(),
        tags.len(),
        "db_tags are unique case-insensitively"
    );
    assert_eq!(tags[&f.digest("asmA")], "asmA");
    assert_eq!(tags[&f.digest("asmC")], "asmC");
}

#[test]
fn retain_keeps_full_deflines_and_recomputes_totals() {
    let f = encoded();
    let (st, plan) = plan_for(&f.spec(), &all());
    let keep: HashSet<String> = plan
        .subjects
        .iter()
        .step_by(3)
        .map(|s| s.digest.clone())
        .collect();
    let r = plan.retain(|d| keep.contains(d));
    assert_eq!(r.subjects.len(), keep.len());
    let mut last = None;
    for s in &r.subjects {
        let orig = plan
            .subjects
            .iter()
            .position(|o| o.digest == s.digest)
            .unwrap();
        assert_eq!(plan.subjects[orig].deflines, s.deflines);
        assert!(last.is_none_or(|l| orig > l), "relative order kept");
        last = Some(orig);
    }
    assert_eq!(r.totals.num_seqs as usize, keep.len());
    assert_eq!(
        r.totals.total_length,
        r.subjects.iter().map(|s| s.length as u64).sum::<u64>()
    );
    for (ci, c) in r.collections.iter().enumerate() {
        let rows: Vec<u64> = r
            .subjects
            .iter()
            .flat_map(|s| {
                s.occurrences
                    .iter()
                    .filter(|o| o.collection as usize == ci)
                    .map(move |_| s.length as u64)
            })
            .collect();
        assert_eq!(c.n_sequences as usize, rows.len());
        assert_eq!(c.total_length, rows.iter().sum::<u64>());
        assert_eq!(c.n_names, plan.collections[ci].n_names);
    }
    // a database written from the retained plan passes md5 and id checks
    let d = tempdir();
    let out = d.path().join("db");
    let fetch = StoreFiles::new(&st);
    blastdb::write_blastdb(&r, &fetch, &write_opts(&out), &mut |_| true).unwrap();
    let rep = verify_db(&out, &st.store, &fetch, 2, &mut |_| true).unwrap();
    assert!(rep.ok(), "{rep:?}");
    assert_eq!(rep.sequences_checked as usize, keep.len());
}

fn write_opts(out: &Path) -> WriteOptions {
    WriteOptions {
        out: out.to_path_buf(),
        title: TITLE.into(),
        max_file_size: 3_000_000_000,
        budget_bytes: 1 << 30,
        batch_bytes: 8 << 20,
        force: false,
        check_free_space: true,
        date: BlastDate::parse(DATE).unwrap(),
        progress_interval_ms: 0,
        threads: 2,
        store: "fixture".into(),
        store_label: "fixture".into(),
        remote: false,
        selection: String::new(),
        nstats: nstats_opts(),
        manifest: None,
    }
}

#[test]
fn plan_totals_need_no_sequence_bytes() {
    let f = encoded();
    let (_, plan) = plan_for(&f.spec(), &all());
    let d = tempdir();
    let copy = d.path().join("noseq");
    copy_dir(&f.path, &copy, &|p| {
        p.file_name().is_some_and(|n| n == "sequences")
    });
    assert!(!copy.join("sequences").exists());
    let (_, p2) = plan_for(&copy.display().to_string(), &all());
    assert_eq!(p2.totals, plan.totals);
    assert_eq!(p2.subjects.len(), plan.subjects.len());
}

#[test]
fn store_digests() {
    let f = encoded();
    let st = OpenedStore::open(&f.spec(), None).unwrap();
    let m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(f.path.join("rgstore.json")).unwrap())
            .unwrap();
    let got = st.digests();
    for k in [
        "sequences_digest",
        "collections_digest",
        "aliases_digest",
        "fhr_digest",
    ] {
        assert_eq!(got[k], m[k].as_str().unwrap(), "{k}");
    }
}

#[test]
fn collection_rows_match_rgsi() {
    let f = encoded();
    let mut st = OpenedStore::open(&f.spec(), None).unwrap();
    for alias in COLLECTIONS {
        let d = f.digest(alias);
        let rows = st.collection_rows(&d).unwrap();
        let text = std::fs::read_to_string(f.path.join(format!("collections/{d}.rgsi"))).unwrap();
        let want: Vec<(String, u64, String)> = text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
            .map(|l| {
                let x: Vec<&str> = l.split('\t').collect();
                (x[0].to_string(), x[1].parse().unwrap(), x[3].to_string())
            })
            .collect();
        let got: Vec<(String, u64, String)> = rows
            .into_iter()
            .map(|r| (r.name, r.length, r.digest))
            .collect();
        assert_eq!(got, want, "{alias}");
    }
}

#[test]
fn sidecar_readers_roundtrip() {
    let d = tempdir();
    let out = build_all(d.path());
    let (_, plan) = plan_for(&encoded().spec(), &all());
    let cm = sidecars::read_collections(&with_ext(&out, "collections.tsv")).unwrap();
    assert_eq!(cm.len(), plan.collections.len());
    for (m, c) in cm.iter().zip(&plan.collections) {
        assert_eq!(m.collection, c.digest);
        assert_eq!(m.db_tag, c.db_tag);
        assert_eq!(m.display, c.display);
        assert_eq!(m.aliases, c.aliases);
        assert_eq!(
            (m.n_sequences, m.total_length),
            (c.n_sequences, c.total_length)
        );
        // no fixture collection has an FHR taxon; asmC's assembly is its refseq alias
        assert_eq!((m.taxid, &m.species), (None, &None));
        let want_asm = (c.digest == encoded().digest("asmC")).then(|| "GCF_TEST.1".to_string());
        assert_eq!(m.assembly, want_asm, "{}", c.display);
        let f = m
            .n_fraction
            .expect("n_fraction filled from the build's counts");
        assert!((0.0..1.0).contains(&f), "{f}");
        // asmA chr2 (N runs) is also asmB chrUn_x
        let with_n = [encoded().digest("asmA"), encoded().digest("asmB")];
        assert_eq!(f > 0.0, with_n.contains(&c.digest), "{}: {f}", c.display);
    }
    let members = sidecars::read_members(&with_ext(&out, "members.tsv")).unwrap();
    let mut want = Vec::new();
    for (oid, s) in plan.subjects.iter().enumerate() {
        for o in &s.occurrences {
            want.push((
                oid as u32,
                s.digest.clone(),
                plan.collections[o.collection as usize].digest.clone(),
                o.name.clone(),
                o.id_name.clone(),
            ));
        }
    }
    let got: Vec<_> = members
        .into_iter()
        .map(|m| (m.oid, m.digest, m.collection, m.name, m.id_name))
        .collect();
    assert_eq!(got, want);
    assert_eq!(
        got.len() as u64,
        plan.collections
            .iter()
            .map(|c| c.n_sequences as u64)
            .sum::<u64>()
    );
}

#[test]
fn totals_match_plan() {
    let d = tempdir();
    let out = build_all(d.path());
    let (_, plan) = plan_for(&encoded().spec(), &all());
    let r = DbReader::open(&out).unwrap();
    assert_eq!(r.num_oids(), plan.totals.num_seqs);
    assert_eq!(r.total_letters(), plan.totals.total_length);
    assert_eq!(r.max_len(), plan.totals.max_len);
    let j = njs::read(&with_ext(&out, "njs")).unwrap();
    assert_eq!(j.number_of_sequences, plan.totals.num_seqs as u64);
    assert_eq!(j.number_of_letters, plan.totals.total_length);
    assert_eq!(j.number_of_volumes, 1);
    assert_eq!(j.dbname, "db");
    assert_eq!(j.description, TITLE);
    assert_eq!(j.last_updated, "2026-11-08T12:51:00");
    let listed: u64 = j
        .files
        .iter()
        .map(|f| std::fs::metadata(d.path().join(f)).unwrap().len())
        .sum();
    assert_eq!(j.bytes_total, listed);
    assert!(
        !j.files
            .iter()
            .any(|f| f.ends_with(".njs") || f.ends_with(".tsv"))
    );
}

#[test]
fn orphans_only_with_all() {
    let f = encoded();
    let (_, p_all) = plan_for(&f.spec(), &all());
    let last = p_all.subjects.last().unwrap();
    assert_eq!(last.digest, f.orphan);
    assert!(last.occurrences.is_empty());
    assert_eq!(last.deflines.len(), 1);
    assert!(
        last.deflines[0].title.starts_with("no collection; length:"),
        "{}",
        last.deflines[0].title
    );
    assert_eq!(last.deflines[0].taxid, Some(0));
    assert_eq!(p_all.n_orphans(), 1);
    let named: Vec<&str> = COLLECTIONS.to_vec();
    let (_, p_named) = plan_for(&f.spec(), &colls(&named));
    assert_eq!(p_named.n_orphans(), 0);
    assert_eq!(p_named.subjects.len() + 1, p_all.subjects.len());
}

// ---------------------------------------------------------------------------
// volumes
// ---------------------------------------------------------------------------

#[test]
fn multi_volume_build() {
    let f = encoded();
    let d = tempdir();
    let out = d.path().join("db");
    let mut req = request(&f.spec(), all(), &out);
    req.options.max_file_size = 4000;
    let s = build_req(&req);
    assert!(s.volumes >= 3, "{} volumes", s.volumes);
    let nal = std::fs::read_to_string(with_ext(&out, "nal")).unwrap();
    let list: Vec<String> = (0..s.volumes).map(|i| format!("db.{i:02}")).collect();
    assert!(
        nal.contains(&format!("DBLIST {}\n", list.join(" "))),
        "{nal}"
    );
    assert!(nal.contains(&format!("# Alias file created: {DATE}\n")));
    assert!(nal.contains(&format!("TITLE {TITLE}\n")));
    let r = DbReader::open(&out).unwrap();
    let vols = r.ndb().unwrap().volumes().unwrap();
    assert_eq!(vols.len() as u64, s.volumes);
    for (i, (name, n)) in vols.iter().enumerate() {
        assert_eq!(name, &list[i]);
        assert_eq!(*n, r.volumes[i].header.n_oids);
        assert_eq!(r.volumes[i].header.volume, i as u32);
    }
    let st = OpenedStore::open(&f.spec(), None).unwrap();
    let fetch = StoreFiles::new(&st);
    let rep = verify_db(&out, &st.store, &fetch, 2, &mut |_| true).unwrap();
    assert!(rep.ok(), "{rep:?}");
    let j = njs::read(&with_ext(&out, "njs")).unwrap();
    assert_eq!(j.number_of_volumes, s.volumes);
    assert!(j.files.contains(&"db.nal".to_string()));
}

// ---------------------------------------------------------------------------
// sidecars, skips, collisions
// ---------------------------------------------------------------------------

#[test]
fn sidecars_contents() {
    let d = tempdir();
    let out = build_all(d.path());
    let read = |s: &str| std::fs::read_to_string(with_ext(&out, s)).unwrap();
    let ids = read("ids.tsv");
    let mut l = ids.lines();
    assert_eq!(l.next(), Some("#id\toid\tsequence_digest"));
    assert_eq!(l.next(), Some(sidecars::IDS_HELP));
    assert!(ids.contains("gnl|asmA|ENA_KI270740_KI270740.1\t"));
    assert!(ids.contains("gnl|asmA|chrX_alt\t"));
    assert!(ids.contains("gnl|asmA|HLA-A*01:01:01:01\t"));
    assert!(ids.contains("\nCM000001.1\t"));
    let san = read("sanitized.tsv");
    assert!(san.starts_with("#kind\tcollection\toriginal\tsanitized\n"));
    let a = encoded().digest("asmA");
    assert!(san.contains(&format!(
        "name\t{a}\tENA|KI270740|KI270740.1\tENA_KI270740_KI270740.1\n"
    )));
    assert!(san.contains(&format!("name\t{a}\tchrX;alt\tchrX_alt\n")));
    let sk = read("skipped.tsv");
    assert!(sk.starts_with("#sequence_digest\tlength\talphabet\treason\tcollections\n"));
    assert!(sk.contains("\t59\tprotein\tnon_nucleotide\t1\n"), "{sk}");
    let c = read("collections.tsv");
    assert!(
        c.contains(&format!(
            "{}\tasmC\tasmC\tname:asmC,name:asmC-alt,refseq:GCF_TEST.1\t\t\tGCF_TEST.1\t2\t564\t0\n",
            encoded().digest("asmC")
        )),
        "{c}"
    );
    let m = read("members.tsv");
    assert!(m.starts_with("#oid\tdigest\tcollection\tname\tid_name\n"));
    assert!(m.contains(&format!("\t{a}\tchrX;alt\tchrX_alt\n")));
}

#[test]
fn non_nucleotide_skipped() {
    let f = encoded();
    let (_, plan) = plan_for(&f.spec(), &colls(&["asmC"]));
    assert_eq!(plan.skipped.len(), 1);
    assert_eq!(plan.skipped[0].reason, SkipReason::NonNucleotide);
    assert_eq!(plan.skipped[0].length, 59);
    assert_eq!(plan.subjects.len(), 2);
    let c = &plan.collections[0];
    assert_eq!((c.n_names, c.n_sequences, c.total_length), (3, 2, 564));
    assert!(
        !plan
            .subjects
            .iter()
            .any(|s| s.occurrences.iter().any(|o| o.name == "prot1"))
    );
}

#[test]
fn oversize_skip_and_fail() {
    let f = encoded();
    let o = PlanOptions {
        on_oversize: OnOversize::Skip,
        max_subject_len: 1001,
    };
    let (_, plan) = plan_with(&f.spec(), &colls(&["asmA", "asmB"]), &o);
    let over: Vec<u64> = plan
        .skipped
        .iter()
        .filter(|s| s.reason == SkipReason::Oversize)
        .map(|s| s.length)
        .collect();
    let mut over_sorted = over.clone();
    over_sorted.sort();
    assert_eq!(over_sorted, vec![1002, 1003, 12670]);
    assert!(plan.subjects.iter().all(|s| s.length <= 1001));
    // through the API: skip lists them, fail writes nothing
    let d = tempdir();
    let out = d.path().join("db");
    let mut req = request(&f.spec(), colls(&["asmA", "asmB"]), &out);
    req.options.max_subject_len = 1001;
    let s = build_req(&req);
    assert_eq!(s.preflight.skipped_oversize, 3);
    let sk = std::fs::read_to_string(with_ext(&out, "skipped.tsv")).unwrap();
    assert_eq!(sk.matches("\toversize\t").count(), 3);
    let d2 = tempdir();
    let out2 = d2.path().join("db");
    let mut req = request(&f.spec(), colls(&["asmA", "asmB"]), &out2);
    req.options.max_subject_len = 1001;
    req.options.on_oversize = OnOversize::Fail;
    let e = gtars_rgblast::build(&req, &mut |_| true).unwrap_err();
    assert!(matches!(e, Error::Oversize { count: 3, .. }), "{e:?}");
    assert!(OutputSet::new(&out2).existing().unwrap().is_empty());
}

#[test]
fn sanitizing_collision_is_error() {
    let d = tempdir();
    let store = d.path().join("collide");
    build_collide_store(&store);
    let out = d.path().join("db");
    let e = gtars_rgblast::build(
        &request(&store.display().to_string(), all(), &out),
        &mut |_| true,
    )
    .unwrap_err();
    match &e {
        Error::IdCollision { id, sanitized, .. } => {
            assert_eq!(id, "gnl|collide|a_b");
            assert!(*sanitized);
        }
        other => panic!("expected IdCollision, got {other:?}"),
    }
    assert_eq!(
        e.to_string(),
        "names 'a|b' and 'a_b' in collection collide both become id 'gnl|collide|a_b' after replacing characters BLAST can't store; rename one in the store"
    );
    assert!(OutputSet::new(&out).existing().unwrap().is_empty());
}

#[test]
fn alias_clash_falls_back_to_digest() {
    use gtars_refget::store::{FastaImportOptions, RefgetStore};
    let d = tempdir();
    let path = d.path().join("clash");
    let mut s = RefgetStore::on_disk(&path).unwrap();
    s.set_quiet(true);
    let src = d.path().join("src");
    fixture_copies(&src, &["asmA.fa", "asmB.fa"]);
    let fa = |n: &str| src.join(n);
    let (a, _) = s
        .add_sequence_collection_from_fasta(fa("asmA.fa"), FastaImportOptions::new())
        .unwrap();
    let (b, _) = s
        .add_sequence_collection_from_fasta(fa("asmB.fa"), FastaImportOptions::new())
        .unwrap();
    s.add_collection_alias("name", "X", &a.digest).unwrap();
    s.add_collection_alias("ucsc", "X", &b.digest).unwrap();
    s.write().unwrap();
    drop(s);
    let spec = path.display().to_string();
    let (st, plan) = plan_for(&spec, &all());
    for c in &plan.collections {
        assert_eq!(c.db_tag, c.digest);
    }
    assert_eq!(plan.ambiguous_aliases, vec!["X".to_string()]);
    for s in &plan.subjects {
        for dl in &s.deflines {
            assert!(
                !dl.ids
                    .iter()
                    .any(|id| matches!(id, SeqId::General { db, .. } if db == "X"))
            );
        }
    }
    assert!(matches!(
        st.resolve_collection("X"),
        Err(Error::InvalidOptions(m)) if m.contains("names more than one collection")
    ));
    let out = d.path().join("db");
    let sum = build(&spec, all(), &out);
    assert!(
        sum.warnings
            .iter()
            .any(|w| w.starts_with("1 collection aliases name more than one"))
    );
}

// ---------------------------------------------------------------------------
// preflight
// ---------------------------------------------------------------------------

#[test]
fn preflight_exact_parts_match_files() {
    let f = encoded();
    for max in [3_000_000_000u64, 4000] {
        let d = tempdir();
        let out = d.path().join("db");
        let mut req = request(&f.spec(), all(), &out);
        req.options.max_file_size = max;
        let pre = gtars_rgblast::preflight(&req).unwrap();
        assert!(
            OutputSet::new(&out).existing().unwrap().is_empty(),
            "preflight writes nothing"
        );
        let s = build_req(&req);
        let size = |p: &Path| std::fs::metadata(p).unwrap().len();
        let r = DbReader::open(&out).unwrap();
        let vols: Vec<std::path::PathBuf> = r.volumes.iter().map(|v| v.stem.clone()).collect();
        let sum = |ext: &str| vols.iter().map(|v| size(&with_ext(v, ext))).sum::<u64>();
        let amb: u64 = (0..r.num_oids())
            .map(|o| r.ambiguity_bytes(o).unwrap())
            .sum();
        assert_eq!(sum("nhr"), pre.exact.nhr);
        assert_eq!(size(&with_ext(&out, "nos")), pre.exact.nos);
        assert_eq!(size(&with_ext(&out, "not")), pre.exact.not);
        assert_eq!(size(&with_ext(&out, "nto")), pre.exact.nto);
        assert_eq!(
            sum("nsq") - amb,
            pre.exact.nsq_without_ambiguity,
            "max {max}"
        );
        if s.volumes == pre.predicted_volumes {
            assert_eq!(sum("nin"), pre.exact.nin);
        }
        let db_files: u64 = s
            .files
            .iter()
            .filter(|p| !p.to_string_lossy().ends_with(".njs"))
            .map(|p| size(p))
            .sum();
        assert!(
            db_files <= pre.exact_bytes + pre.allowance_bytes,
            "{db_files} > {} + {}",
            pre.exact_bytes,
            pre.allowance_bytes
        );
    }
}

// ---------------------------------------------------------------------------
// output safety
// ---------------------------------------------------------------------------

#[test]
fn existing_output_requires_force() {
    let f = encoded();
    let d = tempdir();
    let out = d.path().join("db");
    build(&f.spec(), colls(&["asmB"]), &out);
    let before = files_in(d.path());
    let e = gtars_rgblast::build(&request(&f.spec(), colls(&["asmA"]), &out), &mut |_| true)
        .unwrap_err();
    assert!(matches!(e, Error::OutputExists { .. }), "{e:?}");
    assert_eq!(files_in(d.path()), before, "nothing touched");
    // a stale volume file alone also counts
    let d2 = tempdir();
    std::fs::write(d2.path().join("db.07.nsq"), b"x").unwrap();
    let e = gtars_rgblast::build(
        &request(&f.spec(), all(), &d2.path().join("db")),
        &mut |_| true,
    )
    .unwrap_err();
    assert!(matches!(e, Error::OutputExists { .. }));
    let mut req = request(&f.spec(), colls(&["asmA"]), &out);
    req.options.force = true;
    build_req(&req);
    let after = files_in(d.path());
    assert_ne!(after, before);
    let mut req = request(&f.spec(), all(), &d2.path().join("db"));
    req.options.force = true;
    build_req(&req);
    assert!(
        !d2.path().join("db.07.nsq").exists(),
        "--force removes stale volume files"
    );
}

/// The store's files, failing at one subject.
struct FailAt<'a> {
    inner: StoreFiles<'a>,
    digest: String,
}

impl SeqFetch for FailAt<'_> {
    fn read(&self, digest: &str) -> Result<StoredBytes, Error> {
        if digest == self.digest {
            return Err(Error::Fetch {
                what: digest.to_string(),
                msg: "injected failure".into(),
            });
        }
        self.inner.read(digest)
    }
}

#[test]
fn partial_output_removed_on_error() {
    let f = encoded();
    let (st, plan) = plan_for(&f.spec(), &all());
    let d = tempdir();
    let out = d.path().join("db");
    let fetch = FailAt {
        inner: StoreFiles::new(&st),
        digest: plan.subjects[4].digest.clone(),
    };
    let mut wo = write_opts(&out);
    wo.max_file_size = 4000; // several volumes exist when it fails
    let e = blastdb::write_blastdb(&plan, &fetch, &wo, &mut |_| true).unwrap_err();
    assert!(matches!(e, Error::Fetch { .. }), "{e:?}");
    assert!(OutputSet::new(&out).existing().unwrap().is_empty());
    assert_eq!(
        std::fs::read_dir(d.path()).unwrap().count(),
        0,
        "no files left"
    );
}

#[test]
fn progress_events_and_cancel() {
    let f = encoded();
    let d = tempdir();
    let out = d.path().join("db");
    let mut events: Vec<ProgressEvent> = Vec::new();
    gtars_rgblast::build(&request(&f.spec(), all(), &out), &mut |e| {
        events.push(e.clone());
        true
    })
    .unwrap();
    let kinds: Vec<&str> = events
        .iter()
        .map(|e| match e {
            ProgressEvent::Opened { .. } => "opened",
            ProgressEvent::Planned(_) => "planned",
            ProgressEvent::Sequences { .. } => "sequences",
            ProgressEvent::Finishing { .. } => "finishing",
            ProgressEvent::NScan { .. } => "nscan",
            ProgressEvent::Warning(_) => "warning",
            ProgressEvent::Finished(_) => "finished",
            _ => "rgblast",
        })
        .collect();
    assert_eq!(kinds[0], "opened");
    assert_eq!(*kinds.last().unwrap(), "finished");
    let pos = |k: &str| kinds.iter().position(|x| *x == k).unwrap();
    // selection notes (no FHR taxon) come right after opening, before the plan
    for (i, e) in events.iter().enumerate().take(pos("planned")).skip(1) {
        assert!(matches!(e, ProgressEvent::Warning(_)), "event {i}: {e:?}");
    }
    assert!(pos("sequences") < pos("finishing"));
    let last_warning = kinds.iter().rposition(|x| *x == "warning").unwrap();
    assert!(pos("finishing") < last_warning);
    let mut last = 0;
    for e in &events {
        if let ProgressEvent::Sequences { done, total, .. } = e {
            assert!(*done >= last && done <= total);
            last = *done;
        }
    }
    let stages: Vec<String> = events
        .iter()
        .filter_map(|e| {
            if let ProgressEvent::Finishing { stage } = e {
                Some(stage.clone())
            } else {
                None
            }
        })
        .collect();
    // `.nos` reads the acc2oid keys in OID order, before the LMDB sorts them
    assert_eq!(stages, ["oid files", "lmdb", "sidecars", "njs"]);
    // cancel during the sequence phase
    let d2 = tempdir();
    let out2 = d2.path().join("db");
    let mut n = 0;
    let e = gtars_rgblast::build(&request(&f.spec(), all(), &out2), &mut |e| {
        if matches!(e, ProgressEvent::Sequences { .. }) {
            n += 1;
            return n < 3;
        }
        true
    })
    .unwrap_err();
    assert!(matches!(e, Error::Cancelled));
    assert_eq!(e.kind(), gtars_rgblast::ErrorKind::Cancelled);
    assert!(OutputSet::new(&out2).existing().unwrap().is_empty());
}

// ---------------------------------------------------------------------------
// verify
// ---------------------------------------------------------------------------

#[test]
fn verify_detects_corruption() {
    let f = encoded();
    let d = tempdir();
    let out = build_all(d.path());
    let st = OpenedStore::open(&f.spec(), None).unwrap();
    let fetch = StoreFiles::new(&st);
    assert!(
        verify_db(&out, &st.store, &fetch, 2, &mut |_| true)
            .unwrap()
            .ok()
    );
    let nsq = with_ext(&out, "nsq");
    let mut b = std::fs::read(&nsq).unwrap();
    b[10] ^= 0x55;
    std::fs::write(&nsq, &b).unwrap();
    let rep = verify_db(&out, &st.store, &fetch, 2, &mut |_| true).unwrap();
    assert!(!rep.ok());
    assert_eq!(rep.md5_mismatches, 1, "{rep:?}");
    assert!(rep.examples[0].contains("md5"), "{:?}", rep.examples);
}
