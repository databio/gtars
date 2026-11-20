//! Taxids, titles, the taxonomy index files, `--taxid-map`, alias ranking, `.bsl`, and N
//! counts, on the `taxa` fixture (and plan 1's fixture for the no-FHR case).

mod common;

use common::*;
use gtars_rgblast::api::NstatsRequest;
use gtars_rgblast::blastdb::alias;
use gtars_rgblast::blastdb::output::with_ext;
use gtars_rgblast::blastdb::reader::DbReader;
use gtars_rgblast::blastdb::{lmdb, njs, sidecars};
use gtars_rgblast::plan::{SeqId, SubjectPlan};
use gtars_rgblast::select::Selection;
use gtars_rgblast::select::nstats::{NCount, NstatsCache, NstatsOptions};
use gtars_rgblast::{Error, ProgressEvent};
use std::collections::BTreeMap;
use std::path::Path;

fn subject_named<'a>(
    plan: &'a SubjectPlan,
    coll: &str,
    name: &str,
) -> (usize, &'a gtars_rgblast::plan::Subject) {
    let f = taxa_encoded();
    let ci = plan
        .collections
        .iter()
        .position(|c| c.digest == f.digest(coll))
        .unwrap() as u32;
    plan.subjects
        .iter()
        .enumerate()
        .find(|(_, s)| {
            s.occurrences
                .iter()
                .any(|o| o.collection == ci && o.name == name)
        })
        .unwrap_or_else(|| panic!("{coll} {name}"))
}

#[test]
fn per_defline_taxid_rule() {
    let f = taxa_encoded();
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    let taxids =
        |s: &gtars_rgblast::plan::Subject| s.deflines.iter().map(|d| d.taxid).collect::<Vec<_>>();
    // collections agree -> Some(t)
    let (_, chr1) = subject_named(&plan, "human_a", "chr1");
    assert_eq!(taxids(chr1), vec![Some(9606), Some(9606), Some(9606)]);
    // mixed -> the SQ. defline omits it
    let (_, decoy) = subject_named(&plan, "mouse_a", "decoy");
    assert_eq!(taxids(decoy), vec![None, Some(9606), Some(10090)]);
    // unknown -> 0
    let (_, u1) = subject_named(&plan, "nofhr", "u1");
    assert_eq!(taxids(u1), vec![Some(0), Some(0)]);
    // orphan -> 0
    let o = plan.subjects.iter().find(|s| s.digest == f.orphan).unwrap();
    assert_eq!(taxids(o), vec![Some(0)]);
    // the collection info reached the plan
    let c = |n: &str| {
        plan.collections
            .iter()
            .find(|c| c.digest == f.digest(n))
            .unwrap()
    };
    assert_eq!(
        (c("human_a").taxid, c("human_a").species.as_deref()),
        (Some(9606), Some("Homo sapiens"))
    );
    assert_eq!(
        c("human_masked").species.as_deref(),
        Some("Homo sapiens"),
        "trimmed"
    );
    assert_eq!((c("nofhr").taxid, c("nofhr").species.clone()), (None, None));
}

#[test]
fn exact_titles() {
    let f = taxa_encoded();
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    let (_, chr1) = subject_named(&plan, "human_a", "chr1");
    assert_eq!(
        chr1.deflines[0].title,
        format!(
            "Homo sapiens chr1; 2 assemblies; length:300; md5:{}",
            chr1.md5
        )
    );
    assert_eq!(chr1.deflines[1].title, "Homo sapiens human_a (GCF_TEST.1)");
    assert_eq!(chr1.deflines[2].title, "Homo sapiens human_b (GCA_TEST.1)");
    let (_, decoy) = subject_named(&plan, "mouse_a", "decoy");
    assert_eq!(
        decoy.deflines[0].title,
        format!(
            "Homo sapiens, Mus musculus decoy; 2 assemblies; length:250; md5:{}",
            decoy.md5
        )
    );
    assert_eq!(
        decoy.deflines[2].title,
        "Mus musculus mouse_a (GCF_000002.1)"
    );
    let (_, m1) = subject_named(&plan, "human_masked", "chr1");
    assert_eq!(m1.deflines[1].title, "Homo sapiens human_masked");
    let (_, u1) = subject_named(&plan, "nofhr", "u1");
    assert_eq!(
        u1.deflines[0].title,
        format!("unknown species u1; 1 assembly; length:180; md5:{}", u1.md5)
    );
    assert_eq!(u1.deflines[1].title, "nofhr");
    let o = plan.subjects.iter().find(|s| s.digest == f.orphan).unwrap();
    assert_eq!(
        o.deflines[0].title,
        format!("no collection; length:120; md5:{}", o.md5)
    );
    // no title repeats its defline's name tag
    for s in &plan.subjects {
        for (d, occ) in s.deflines[1..].iter().zip(&s.occurrences) {
            assert!(
                !d.title.split(' ').any(|w| w == occ.name),
                "{} {}",
                d.title,
                occ.name
            );
        }
    }
}

#[test]
fn alias_ranking_display_and_db_tag() {
    let f = taxa_encoded();
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    let c = |n: &str| {
        plan.collections
            .iter()
            .find(|c| c.digest == f.digest(n))
            .unwrap()
    };
    // `name` ranks first although `insdc` comes first in rgstore.json
    let hb = c("human_b");
    assert_eq!(hb.aliases[0], ("name".to_string(), "human_b".to_string()));
    assert_eq!(
        hb.aliases[1],
        ("insdc".to_string(), "GCA_TEST.1".to_string())
    );
    assert_eq!(
        (hb.display.as_str(), hb.db_tag.as_str()),
        ("human_b", "human_b")
    );
    // other namespaces before the accession ones
    let ha = c("human_a");
    let ns: Vec<&str> = ha.aliases.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(ns, vec!["name", "ucsc", "ucsc"]);
    let ma = c("mouse_a");
    let ns: Vec<&str> = ma.aliases.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(ns, vec!["name", "refseq", "refseq"]);
    // ids: db_tag first, unique aliases, the digest; `dup` (two collections) left out
    let (_, chr1) = subject_named(&plan, "human_a", "chr1");
    let g = |db: &str| SeqId::General {
        db: db.into(),
        tag: "chr1".into(),
    };
    assert_eq!(
        chr1.deflines[1].ids,
        vec![g("human_a"), g("hg_test"), g(&f.digest("human_a"))]
    );
    assert_eq!(
        chr1.deflines[2].ids,
        vec![g("human_b"), g("GCA_TEST.1"), g(&f.digest("human_b"))]
    );
    assert_eq!(plan.ambiguous_aliases, vec!["dup".to_string()]);
    // the sequence alias joins defline 0
    assert_eq!(chr1.deflines[0].ids[1], SeqId::Local("NC_TEST1.1".into()));
    // collection order follows display
    let order: Vec<&str> = plan
        .collections
        .iter()
        .map(|c| c.display.as_str())
        .collect();
    assert_eq!(
        order,
        vec!["human_a", "human_b", "human_masked", "mouse_a", "nofhr"]
    );
}

fn u64s(b: &[u8]) -> Vec<u64> {
    b.chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect()
}
fn u32s(b: &[u8]) -> Vec<u32> {
    b.chunks_exact(4)
        .map(|c| u32::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

/// (per-OID taxid sets from `.not`, taxid -> OIDs from `.nto`, `.nto` offsets).
fn read_taxonomy(base: &Path) -> (Vec<Vec<i32>>, BTreeMap<i32, Vec<u32>>) {
    let not = std::fs::read(with_ext(base, "not")).unwrap();
    let n = u64::from_le_bytes(not[..8].try_into().unwrap()) as usize;
    let ends = u64s(&not[8..8 + 8 * n]);
    let vals: Vec<i32> = not[8 + 8 * n..]
        .chunks_exact(4)
        .map(|c| i32::from_le_bytes(c.try_into().unwrap()))
        .collect();
    let mut sets = Vec::new();
    let mut start = 0usize;
    for e in ends {
        sets.push(vals[start..e as usize].to_vec());
        start = e as usize;
    }
    assert_eq!(start, vals.len());
    let nto = u32s(&std::fs::read(with_ext(base, "nto")).unwrap());
    let ntf = lmdb::read_ntf(&with_ext(base, "ntf")).unwrap();
    let mut by = BTreeMap::new();
    for (t, off) in ntf {
        let i = off as usize / 4;
        let k = nto[i] as usize;
        by.insert(t, nto[i + 1..i + 1 + k].to_vec());
    }
    (sets, by)
}

#[test]
fn taxonomy_files_read_back() {
    let f = taxa_encoded();
    let d = tempdir();
    let out = d.path().join("db");
    let s = build(&f.spec(), Selection::all(), &out);
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    let (sets, by) = read_taxonomy(&out);
    assert_eq!(sets.len(), plan.subjects.len());
    let mut want_by: BTreeMap<i32, Vec<u32>> = BTreeMap::new();
    for (oid, s) in plan.subjects.iter().enumerate() {
        let want = gtars_rgblast::plan::taxonomy::taxid_set(&s.deflines);
        assert_eq!(sets[oid], want, "OID {oid}");
        for t in want {
            want_by.entry(t).or_default().push(oid as u32);
        }
    }
    assert_eq!(by, want_by);
    assert_eq!(by.keys().copied().collect::<Vec<_>>(), vec![0, 9606, 10090]);
    // the decoy OID is in both taxa and not in 0
    let (doid, _) = subject_named(&plan, "mouse_a", "decoy");
    assert_eq!(sets[doid], vec![9606, 10090]);
    let j = njs::read(&with_ext(&out, "njs")).unwrap();
    assert_eq!(j.number_of_taxids, Some(3));
    // deflines read back from .nhr, the omitted taxid included
    let r = DbReader::open(&out).unwrap();
    for (oid, s) in plan.subjects.iter().enumerate() {
        assert_eq!(r.deflines(oid as u32).unwrap(), s.deflines, "OID {oid}");
    }
    // the build reports its taxa (for the CLI's -taxids hint)
    let taxa: Vec<(i32, Option<String>, u64)> = s
        .taxa
        .iter()
        .map(|t| (t.taxid, t.species.clone(), t.collections))
        .collect();
    assert_eq!(
        taxa,
        vec![
            (9606, Some("Homo sapiens".into()), 3),
            (10090, Some("Mus musculus".into()), 1)
        ]
    );
    let hint = gtars_rgblast::progress::taxid_hint(&s.taxa).unwrap();
    assert!(hint.starts_with("taxids written for 2 taxa. Filter exactly with: blastn -taxids 9606 -no_taxid_expansion\n"), "{hint}");
    assert!(hint.contains("taxdb.tar.gz"));
}

#[test]
fn store_without_fhr_taxa_keeps_taxid_zero() {
    let f = encoded();
    let d = tempdir();
    let out = d.path().join("db");
    let s = build(&f.spec(), Selection::all(), &out);
    assert!(s.taxa.is_empty());
    assert!(gtars_rgblast::progress::taxid_hint(&s.taxa).is_none());
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    for sub in &plan.subjects {
        assert!(sub.deflines.iter().all(|d| d.taxid == Some(0)));
    }
    let n = plan.subjects.len();
    // the same bytes as makeblastdb writes for an all-zero database
    let not = std::fs::read(with_ext(&out, "not")).unwrap();
    let mut want = (n as u64).to_le_bytes().to_vec();
    for i in 1..=n as u64 {
        want.extend_from_slice(&i.to_le_bytes());
    }
    for _ in 0..n {
        want.extend_from_slice(&0i32.to_le_bytes());
    }
    assert_eq!(not, want);
    let nto = std::fs::read(with_ext(&out, "nto")).unwrap();
    let mut want = (n as u32).to_le_bytes().to_vec();
    for i in 0..n as u32 {
        want.extend_from_slice(&i.to_le_bytes());
    }
    assert_eq!(nto, want);
    assert_eq!(
        lmdb::read_ntf(&with_ext(&out, "ntf")).unwrap(),
        vec![(0, 0)]
    );
    assert_eq!(
        njs::read(&with_ext(&out, "njs")).unwrap().number_of_taxids,
        None
    );
}

#[test]
fn collections_tsv_filled() {
    let f = taxa_encoded();
    let d = tempdir();
    let out = d.path().join("db");
    build(&f.spec(), Selection::all(), &out);
    let cm = sidecars::read_collections(&with_ext(&out, "collections.tsv")).unwrap();
    let get = |n: &str| cm.iter().find(|c| c.collection == f.digest(n)).unwrap();
    let ha = get("human_a");
    assert_eq!(
        (
            ha.taxid,
            ha.species.as_deref(),
            ha.assembly.as_deref(),
            ha.n_fraction
        ),
        (
            Some(9606),
            Some("Homo sapiens"),
            Some("GCF_TEST.1"),
            Some(0.0)
        )
    );
    assert_eq!(get("human_b").assembly.as_deref(), Some("GCA_TEST.1"));
    assert_eq!(get("mouse_a").assembly.as_deref(), Some("GCF_000002.1"));
    let hm = get("human_masked");
    assert!((hm.n_fraction.unwrap() - 198.0 / 580.0).abs() < 1e-12);
    assert_eq!(hm.assembly, None);
    let nf = get("nofhr");
    assert_eq!(
        (nf.taxid, nf.species.clone(), nf.assembly.clone()),
        (None, None, None)
    );
    assert_eq!(nf.n_fraction, Some(0.0));
    // the manifest records the selection and the oids sidecar every OID
    let m = sidecars::read_manifest(&out).unwrap();
    assert_eq!(m.selection.given, Selection::all());
    let oids = sidecars::read_oids(&with_ext(&out, "oids.tsv")).unwrap();
    assert_eq!(oids.len() as u32, m.totals.num_seqs);
    let masked = oids.iter().find(|r| r.n == 198).expect("the masked chr1");
    assert_eq!((masked.length, masked.ambiguous), (300, 198));
}

#[test]
fn taxid_map_overrides_fhr() {
    let f = taxa_encoded();
    let d = tempdir();
    let map = d.path().join("map.tsv");
    std::fs::write(
        &map,
        "# collection\ttaxid\tspecies\nnofhr\t32630\tsynthetic construct\nname:human_b\t9598\nrefseq:dup\t10090\n",
    )
    .unwrap();
    let out = d.path().join("db");
    let mut req = request(&f.spec(), Selection::explicit(&["human_b", "nofhr"]), &out);
    req.options.taxid_map = Some(map.clone());
    build_req(&req);
    let cm = sidecars::read_collections(&with_ext(&out, "collections.tsv")).unwrap();
    let get = |n: &str| cm.iter().find(|c| c.collection == f.digest(n)).unwrap();
    assert_eq!(
        (get("nofhr").taxid, get("nofhr").species.as_deref()),
        (Some(32630), Some("synthetic construct"))
    );
    // a different taxid without a species name drops FHR's name
    assert_eq!(
        (get("human_b").taxid, get("human_b").species.clone()),
        (Some(9598), None)
    );
    // the line for mouse_a (not selected) is ignored, but kept in the manifest
    assert_eq!(cm.len(), 2);
    let m = sidecars::read_manifest(&out).unwrap();
    assert_eq!(m.taxid_map.len(), 3);
    // an unknown selector is an error naming the line
    std::fs::write(&map, "nofhr\t1\nnope\t2\n").unwrap();
    let mut req = request(&f.spec(), Selection::all(), &d.path().join("db2"));
    req.options.taxid_map = Some(map);
    let e = gtars_rgblast::build(&req, &mut |_| true).unwrap_err();
    assert!(matches!(e, Error::InvalidOptions(_)), "{e:?}");
    assert!(e.to_string().contains("line 2"), "{e}");
}

#[test]
fn bsl_matches_blastdb_aliastool() {
    let dir = Path::new(EXPECTED).join("bsl");
    let reference = std::fs::read(dir.join("reference.bsl")).unwrap();
    let ids: Vec<String> = std::fs::read_to_string(dir.join("ids.txt"))
        .unwrap()
        .lines()
        .map(String::from)
        .collect();
    let parsed = alias::parse_bsl(&reference).unwrap();
    assert_eq!(parsed.ids.len(), 9);
    let ours = alias::bsl_bytes(&ids, "reference", &parsed.date);
    assert_eq!(ours, reference);
    // write_bsl takes the title from the file stem
    let d = tempdir();
    let p = d.path().join("reference.bsl");
    alias::write_bsl_dated(&ids, &p, &parsed.date).unwrap();
    assert_eq!(std::fs::read(&p).unwrap(), reference);
    let back = alias::read_bsl(&p).unwrap();
    assert_eq!(back, parsed);
}

#[test]
fn nstats_by_products_cache_and_report() {
    let f = taxa_encoded();
    let d = tempdir();
    let cache_dir = d.path().join("cache");
    let opts = NstatsOptions {
        cache_dir: Some(cache_dir.clone()),
        ..Default::default()
    };
    // a build appends its conversion counts for 3- and 4-bit sequences
    let out = d.path().join("db");
    let mut req = request(&f.spec(), Selection::all(), &out);
    req.options.nstats = opts.clone();
    build_req(&req);
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    let (_, m1) = subject_named(&plan, "human_masked", "chr1");
    let (_, iu) = subject_named(&plan, "mouse_a", "chrUn_iupac");
    let cache = NstatsCache {
        dir: cache_dir.join("blast").join("nstats"),
    };
    let got = cache.read(&[m1.digest.as_str(), iu.digest.as_str()]);
    assert_eq!(
        got[&m1.digest][0],
        NCount {
            length: 300,
            n: 198,
            ambiguous: 198
        }
    );
    assert_eq!(got[&iu.digest][0].n, 0);
    assert!(got[&iu.digest][0].ambiguous > 0);
    // the report reads them without a scan
    let nreq = NstatsRequest {
        store: f.spec(),
        cache_dir: None,
        fetch: Default::default(),
        selection: Selection::all(),
        nstats: opts.clone(),
        taxid_map: None,
        threads: 2,
    };
    let mut scans = 0;
    let s = gtars_rgblast::nstats(&nreq, &mut |e| {
        if matches!(e, ProgressEvent::NScan { .. }) {
            scans += 1;
        }
        true
    })
    .unwrap();
    assert_eq!(scans, 0, "everything was cached");
    assert_eq!(s.collections[0].display, "human_masked");
    assert!((s.collections[0].n_fraction.unwrap() - 198.0 / 580.0).abs() < 1e-12);
    assert_eq!(
        (
            s.collections[0].rows,
            s.collections[0].bases,
            s.collections[0].n
        ),
        (2, 580, 198)
    );
    assert_eq!(s.collections.len(), 5);
    // with an empty cache: announced scan, same answer
    let empty = NstatsOptions {
        cache_dir: Some(d.path().join("cache2")),
        ..Default::default()
    };
    let mut events = Vec::new();
    let s2 = gtars_rgblast::nstats(
        &NstatsRequest {
            nstats: empty.clone(),
            ..nreq.clone()
        },
        &mut |e| {
            if let ProgressEvent::NScan { done, total, .. } = e {
                events.push((*done, *total));
            }
            true
        },
    )
    .unwrap();
    assert_eq!(
        events.first(),
        Some(&(0, 2)),
        "announced: two sequences need reading"
    );
    assert_eq!(events.last(), Some(&(2, 2)));
    assert_eq!(
        gtars_rgblast::progress::describe(&ProgressEvent::NScan {
            done: 0,
            total: 2334,
            bytes_read: 0,
            bytes_total: 15_800_000_000
        }),
        "computing N counts for 2,334 sequences, 15.80 GB to read"
    );
    let fr = |s: &gtars_rgblast::api::NstatsSummary| {
        s.collections
            .iter()
            .map(|c| (c.digest.clone(), c.n_fraction))
            .collect::<Vec<_>>()
    };
    assert_eq!(fr(&s), fr(&s2));
}

#[test]
fn n_fraction_counts_built_rows_whatever_the_cache_holds() {
    // a subset build reports its own rows' N fraction, whether or not the N-count cache
    // already holds the whole collection's counts
    let f = taxa_encoded();
    let d = tempdir();
    let (_, plan) = plan_for(&f.spec(), &Selection::all());
    let (_, masked) = subject_named(&plan, "human_masked", "chr1");
    let ci = plan
        .collections
        .iter()
        .position(|c| c.digest == f.digest("human_masked"))
        .unwrap() as u32;
    let clean = plan
        .subjects
        .iter()
        .find(|s| s.digest != masked.digest && s.occurrences.iter().any(|o| o.collection == ci))
        .unwrap();
    let subset = Selection::all().with_sequences(&[clean.digest.as_str()]);
    let n_fraction = |cache: &Path, name: &str| {
        let out = d.path().join(name);
        let mut req = request(&f.spec(), subset.clone(), &out);
        req.options.nstats = NstatsOptions {
            cache_dir: Some(cache.to_path_buf()),
            ..Default::default()
        };
        build_req(&req);
        sidecars::read_collections(&with_ext(&out, "collections.tsv"))
            .unwrap()
            .into_iter()
            .find(|c| c.collection == f.digest("human_masked"))
            .unwrap()
            .n_fraction
    };
    let cold = n_fraction(&d.path().join("cold"), "cold");
    // warm the cache with every sequence of every collection
    let warm_dir = d.path().join("warm");
    let mut full = request(&f.spec(), Selection::all(), &d.path().join("full"));
    full.options.nstats = NstatsOptions {
        cache_dir: Some(warm_dir.clone()),
        ..Default::default()
    };
    build_req(&full);
    let warm = n_fraction(&warm_dir, "warm");
    assert_eq!(cold, Some(0.0), "the built row has no N");
    assert_eq!(warm, cold);
}
