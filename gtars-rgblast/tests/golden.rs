//! Byte-identity canary pinned to BLAST+ 2.17.0 (feature `golden`). A newer BLAST+ that
//! breaks byte identity but passes `blastplus` is a canary failure to investigate, not a
//! broken product: the functional checks are the contract.

mod common;

use common::blast::Blast;
use common::golden::{Lcg, push_query};
use common::txome::txome;
use common::*;
use gtars_refget::store::StorageMode;
use gtars_rgblast::blastdb::njs;
use gtars_rgblast::blastdb::output::with_ext;
use gtars_rgblast::blastdb::reader::{DbReader, dump_lmdb, parse_nin};
use gtars_rgblast::reference_fasta;
use gtars_rgblast::select::Selection;
use gtars_rgblast::store::StoreFiles;
use std::path::Path;

/// Build ours at `<work>/ours/db` and makeblastdb's at `<work>/mk/db` from the same plan.
fn build_pair(store: &str, work: &Path, max_file_size: Option<u64>, threads: usize) -> Blast {
    build_pair_with(store, work, max_file_size, threads, false)
}

/// [`build_pair`], with `taxid_map` handing makeblastdb the plan's taxa as `-taxid_map`.
fn build_pair_with(
    store: &str,
    work: &Path,
    max_file_size: Option<u64>,
    threads: usize,
    taxid_map: bool,
) -> Blast {
    let blast = Blast::find();
    let mut req = request(store, Selection::all(), &work.join("ours/db"));
    req.options.threads = threads;
    if let Some(m) = max_file_size {
        req.options.max_file_size = m;
    }
    build_req(&req);
    let (st, plan) = plan_for(store, &Selection::all());
    let fetch = StoreFiles::new(&st);
    reference_fasta::write(&plan, &fetch, &work.join("ref.fa")).unwrap();
    if taxid_map {
        reference_fasta::write_taxid_map(&plan, &work.join("map.txt")).unwrap();
    }
    std::fs::create_dir_all(work.join("mk")).unwrap();
    let mut args = vec![
        "-in",
        "ref.fa",
        "-dbtype",
        "nucl",
        "-parse_seqids",
        "-blastdb_version",
        "5",
        "-title",
        TITLE,
        "-out",
        "mk/db",
    ];
    let sz = max_file_size.map(|m| format!("{m}B"));
    if let Some(s) = &sz {
        args.extend(["-max_file_sz", s.as_str()]);
    }
    if taxid_map {
        args.extend(["-taxid_map", "map.txt"]);
    }
    blast.ok(work, "makeblastdb", &args);
    blast
}

/// Every difference between our database and makeblastdb's, except the documented ones
/// (date strings, LMDB page layout, `.nog`).
fn compare(ours: &Path, mk: &Path) -> Vec<String> {
    let mut diffs = Vec::new();
    let read = |p: &Path| std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    let (ro, rm) = (DbReader::open(ours).unwrap(), DbReader::open(mk).unwrap());
    if ro.volumes.len() != rm.volumes.len() {
        diffs.push(format!(
            "volumes: ours {} makeblastdb {}",
            ro.volumes.len(),
            rm.volumes.len()
        ));
        return diffs;
    }
    for (vo, vm) in ro.volumes.iter().zip(&rm.volumes) {
        for ext in ["nhr", "nsq"] {
            let (a, b) = (
                read(&with_ext(&vo.stem, ext)),
                read(&with_ext(&vm.stem, ext)),
            );
            if a != b {
                let first = a
                    .iter()
                    .zip(&b)
                    .position(|(x, y)| x != y)
                    .unwrap_or(a.len().min(b.len()));
                diffs.push(format!(
                    "{}.{ext}: sizes {} {}, first difference at byte {first}",
                    vo.stem.display(),
                    a.len(),
                    b.len()
                ));
            }
        }
        let (mut ho, oo, so, ao) = parse_nin(&with_ext(&vo.stem, "nin")).unwrap();
        let (mut hm, om, sm, am) = parse_nin(&with_ext(&vm.stem, "nin")).unwrap();
        ho.date.clear();
        hm.date.clear();
        if (ho, oo, so, ao) != (hm, om, sm, am) {
            diffs.push(format!(
                "{}.nin: fields differ (beyond the date)",
                vo.stem.display()
            ));
        }
    }
    for ext in ["nos", "not", "nto"] {
        if read(&with_ext(ours, ext)) != read(&with_ext(mk, ext)) {
            diffs.push(format!(".{ext} differs"));
        }
    }
    for ext in ["ndb", "ntf"] {
        let (a, b) = (
            dump_lmdb(&with_ext(ours, ext)).unwrap(),
            dump_lmdb(&with_ext(mk, ext)).unwrap(),
        );
        if a != b {
            let only_ours: Vec<&String> = a.iter().filter(|x| !b.contains(x)).take(5).collect();
            let only_mk: Vec<&String> = b.iter().filter(|x| !a.contains(x)).take(5).collect();
            diffs.push(format!(
                ".{ext} LMDB content differs: only ours {only_ours:?}; only makeblastdb {only_mk:?}"
            ));
        }
    }
    let nal = |p: &Path| -> Option<Vec<String>> {
        std::fs::read_to_string(with_ext(p, "nal")).ok().map(|t| {
            t.lines()
                .filter(|l| !l.starts_with("# Alias file created"))
                .map(String::from)
                .collect()
        })
    };
    if nal(ours) != nal(mk) {
        diffs.push(format!(".nal differs: {:?} vs {:?}", nal(ours), nal(mk)));
    }
    let mut jo = njs::read(&with_ext(ours, "njs")).unwrap();
    let mut jm = njs::read(&with_ext(mk, "njs")).unwrap();
    jo.last_updated.clear();
    jm.last_updated.clear();
    // makeblastdb also writes an OID-to-GI file per volume (.nog); we do not. Sizes that
    // depend on the date text (.nin padding, .nal) or the LMDB page layout (.ndb, .ntf) are
    // left out of the byte totals.
    let size = |p: std::path::PathBuf| std::fs::metadata(p).unwrap().len();
    let (dir_o, dir_m) = (ours.parent().unwrap(), mk.parent().unwrap());
    for f in jm.files.iter().filter(|f| f.ends_with(".nog")) {
        jm.bytes_total -= size(dir_m.join(f));
    }
    jm.files.retain(|f| !f.ends_with(".nog"));
    for (j, dir) in [(&mut jo, dir_o), (&mut jm, dir_m)] {
        for f in j.files.clone() {
            let n = size(dir.join(&f));
            if [".nin", ".nal", ".ndb", ".ntf"]
                .iter()
                .any(|e| f.ends_with(e))
            {
                j.bytes_total -= n;
            }
            if f.ends_with(".nin") {
                j.bytes_to_cache -= n;
            }
        }
    }
    if jo != jm {
        diffs.push(format!(".njs differs:\n ours {jo:?}\n mk   {jm:?}"));
    }
    diffs
}

#[test]
fn golden_requires_blast_2_17() {
    let d = tempdir();
    let v = Blast::find().version(d.path());
    assert_eq!(
        v, "2.17.0",
        "golden canary is pinned to BLAST+ 2.17.0, found {v}"
    );
}

#[test]
fn byte_identity_single_volume() {
    let d = tempdir();
    build_pair(&encoded().spec(), d.path(), None, 2);
    let diffs = compare(&d.path().join("ours/db"), &d.path().join("mk/db"));
    assert!(diffs.is_empty(), "{diffs:#?}");
}

#[test]
fn byte_identity_multi_volume() {
    let d = tempdir();
    build_pair(&encoded().spec(), d.path(), Some(4000), 2);
    let ours = d.path().join("ours/db");
    assert!(with_ext(&ours, "nal").exists());
    assert!(DbReader::open(&ours).unwrap().volumes.len() >= 3);
    let diffs = compare(&ours, &d.path().join("mk/db"));
    assert!(diffs.is_empty(), "{diffs:#?}");
}

#[test]
fn byte_identity_over_100_volumes() {
    // more than 100 volumes: makeblastdb renames db.00..db.99 to db.000..db.099
    let d = tempdir();
    let mut rng = Lcg(7);
    let mut fa = Vec::new();
    for i in 0..105 {
        let s: Vec<u8> = (0..400)
            .map(|_| b"ACGT"[(rng.next() % 4) as usize])
            .collect();
        push_query(&mut fa, &format!("s{i}"), &s);
    }
    let fasta = d.path().join("many.fa");
    std::fs::write(&fasta, fa).unwrap();
    let store = d.path().join("store");
    build_store_from_fasta(&store, &fasta, "many", StorageMode::Encoded);
    build_pair(&store.display().to_string(), d.path(), Some(150), 2);
    let ours = d.path().join("ours/db");
    assert!(DbReader::open(&ours).unwrap().volumes.len() > 100);
    assert!(with_ext(&d.path().join("ours/db.000"), "nsq").exists());
    let diffs = compare(&ours, &d.path().join("mk/db"));
    assert!(diffs.is_empty(), "{diffs:#?}");
}

#[test]
fn numeric_and_dash_ids_match_makeblastdb() {
    let d = tempdir();
    let f = encoded();
    build_pair(&f.spec(), d.path(), None, 2);
    let (ro, rm) = (
        DbReader::open(&d.path().join("ours/db")).unwrap(),
        DbReader::open(&d.path().join("mk/db")).unwrap(),
    );
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    let dd = f.digest("asmD");
    let interesting = |oid: usize| {
        let s = &plan.subjects[oid];
        s.deflines.iter().flat_map(|x| x.ids.iter()).any(|id| {
            let t = id.fasta();
            t.contains("|numeric|") || t.contains(&dd) || t == "lcl|123456" || t == "lcl|0099"
        })
    };
    let mut n = 0;
    for oid in (0..plan.subjects.len()).filter(|&o| interesting(o)) {
        n += 1;
        assert_eq!(
            ro.header_bytes(oid as u32).unwrap(),
            rm.header_bytes(oid as u32).unwrap(),
            "OID {oid}: {:?} vs makeblastdb {:?}",
            ro.deflines(oid as u32).unwrap(),
            rm.deflines(oid as u32).unwrap()
        );
    }
    assert!(n >= 13, "{n} OIDs checked");
    let (ko, km) = (
        ro.ndb().unwrap().all_keys().unwrap(),
        rm.ndb().unwrap().all_keys().unwrap(),
    );
    assert_eq!(ko, km, "acc2oid keys");
}

/// The same comparison on the jungle mirror (about 20 min):
/// `GTARS_RGBLAST_STORE_JUNGLE=/path/to/jungle cargo test -p gtars-rgblast --features golden -- --ignored golden_jungle`
#[test]
#[ignore]
fn golden_jungle() {
    let store = std::env::var("GTARS_RGBLAST_STORE_JUNGLE")
        .expect("set GTARS_RGBLAST_STORE_JUNGLE to a local jungle mirror");
    let work = std::env::var("GTARS_RGBLAST_WORK")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| tempdir().keep());
    build_pair(&store, &work, None, 8);
    let diffs = compare(&work.join("ours/db"), &work.join("mk/db"));
    assert!(diffs.is_empty(), "{diffs:#?}");
}

/// The `txome` fixture (2,826 short subjects in 42 collections, 40 of them bins with their own
/// taxids): byte identity, and an `.ndb` of hundreds of pages, so the LMDB comparison exercises
/// a real B-tree. Also split into volumes by the header bound.
#[test]
fn txome_golden_equals_makeblastdb() {
    let f = txome(StorageMode::Encoded);
    let d = tempdir();
    build_pair_with(&f.spec(), d.path(), None, 4, true);
    let ours = d.path().join("ours/db");
    let pages = std::fs::metadata(with_ext(&ours, "ndb")).unwrap().len() / 4096;
    assert!(pages >= 100, ".ndb has {pages} pages");
    let diffs = compare(&ours, &d.path().join("mk/db"));
    assert!(diffs.is_empty(), "{diffs:#?}");

    let d = tempdir();
    build_pair_with(&f.spec(), d.path(), Some(400_000), 3, true);
    let ours = d.path().join("ours/db");
    assert!(DbReader::open(&ours).unwrap().volumes.len() >= 3);
    let diffs = compare(&ours, &d.path().join("mk/db"));
    assert!(diffs.is_empty(), "{diffs:#?}");
}

/// The GENCODE v47 transcript store of `scripts/bench/txome_store.sh` (384,354 subjects). Leaves
/// `ref.fa` (the makeblastdb input) and both databases in the work folder for the benchmarks:
/// `GTARS_RGBLAST_STORE_TXOME=~/work/gtars-rgblast/bench/txome/standalone GTARS_RGBLAST_WORK=~/work/gtars-rgblast/bench/txome/golden cargo test -p gtars-rgblast --release --features golden --test golden -- --ignored golden_txome_gencode`
#[test]
#[ignore]
fn golden_txome_gencode() {
    let store = std::env::var("GTARS_RGBLAST_STORE_TXOME")
        .expect("set GTARS_RGBLAST_STORE_TXOME to the store of scripts/bench/txome_store.sh");
    let work = std::env::var("GTARS_RGBLAST_WORK")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| tempdir().keep());
    std::fs::create_dir_all(&work).unwrap();
    for d in ["ours", "mk"] {
        let _ = std::fs::remove_dir_all(work.join(d));
    }
    build_pair(&store, &work, None, 16);
    let diffs = compare(&work.join("ours/db"), &work.join("mk/db"));
    assert!(diffs.is_empty(), "{diffs:#?}");
}
