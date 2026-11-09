//! Many small subjects (plan 8): the `txome` fixture, batch sizes, exact preflight sizes,
//! the sidecar policy, and thousands of collections.

mod common;

use common::txome::{Layout, TxomeSpec, build_txome_store, txome, write_txome_sources};
use common::*;
use gtars_refget::store::StorageMode;
use gtars_rgblast::api::{ViewRequest, view};
use gtars_rgblast::blastdb::output::with_ext;
use gtars_rgblast::blastdb::sidecars::{self, ViewKind};
use gtars_rgblast::blastdb::{asn, lmdb};
use gtars_rgblast::select::Selection;
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

const EXPECTED_JSON: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/fixtures/txome.expected.json"
);

/// Records and md5 of every generated FASTA, and the totals of a plan of the whole store.
fn fixture_summary() -> serde_json::Value {
    let d = tempdir();
    let mut fastas = BTreeMap::new();
    for s in write_txome_sources(d.path(), &TxomeSpec::STANDARD) {
        let bytes = std::fs::read(&s.fasta).unwrap();
        let records = bytes.iter().filter(|&&b| b == b'>').count();
        let rel = s
            .fasta
            .strip_prefix(d.path())
            .unwrap()
            .display()
            .to_string();
        fastas.insert(
            rel,
            serde_json::json!({"records": records, "md5": md5_hex(&bytes)}),
        );
    }
    let (_, plan) = plan_for(&txome(StorageMode::Encoded).spec(), &Selection::all());
    serde_json::json!({
        "collections": plan.collections.len(),
        "names": plan.n_names(),
        "totals": plan.totals,
        "fastas": fastas,
    })
}

/// The generator is deterministic and the store it gives is the one recorded in
/// `txome.expected.json` (regenerate with `GTARS_RGBLAST_UPDATE_EXPECTED=1`).
#[test]
fn txome_fixture_matches_expected() {
    let got = fixture_summary();
    if std::env::var("GTARS_RGBLAST_UPDATE_EXPECTED").is_ok() {
        let mut t = serde_json::to_string_pretty(&got).unwrap();
        t.push('\n');
        std::fs::write(EXPECTED_JSON, t).unwrap();
    }
    let want: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(EXPECTED_JSON).unwrap()).unwrap();
    assert_eq!(got, want);
    let n = got["totals"]["num_seqs"].as_u64().unwrap();
    assert!(n > 2_500, "{n} subjects");
}

/// `asn::encoded_len` equals the length `asn::encode` writes, and `SeqId::lmdb_key_lens` the
/// keys `lmdb_keys` gives, for every subject of every fixture store.
#[test]
fn exact_sizes_match_encoding_on_every_fixture() {
    let stores = [
        encoded().spec(),
        taxa_encoded().spec(),
        txome(StorageMode::Encoded).spec(),
    ];
    let mut n = 0;
    let mut h = Vec::new();
    for st in &stores {
        let (_, plan) = plan_for(st, &Selection::all());
        for s in &plan.subjects {
            asn::encode(&s.deflines, &mut h);
            assert_eq!(
                asn::encoded_len(&s.deflines) as usize,
                h.len(),
                "{}",
                s.digest
            );
            for id in s.deflines.iter().flat_map(|d| &d.ids) {
                let mut keys = Vec::new();
                id.lmdb_keys(&mut keys);
                let mut lens = Vec::new();
                id.lmdb_key_lens(|l, save| lens.push((l, save)));
                let want: Vec<(usize, bool)> = keys.iter().map(|(k, s)| (k.len(), *s)).collect();
                assert_eq!(lens, want, "{id:?}");
            }
            n += 1;
        }
    }
    assert!(n > 2_500, "{n} subjects checked");
}

/// The preflight's exact sizes are the sizes written.
#[test]
fn txome_preflight_sizes_are_exact() {
    let f = txome(StorageMode::Encoded);
    let d = tempdir();
    let out = d.path().join("db");
    let s = build(&f.spec(), Selection::all(), &out);
    let len = |e: &str| std::fs::metadata(with_ext(&out, e)).unwrap().len();
    let x = &s.preflight.exact;
    assert_eq!(x.nhr, len("nhr"));
    assert_eq!(x.nin, len("nin"));
    assert_eq!(x.nos, len("nos"));
    assert_eq!(x.not, len("not"));
    assert_eq!(x.nto, len("nto"));
    assert!(x.nsq_without_ambiguity <= len("nsq"));
}

/// Every file of a build is byte-identical whatever the batch size, thread count or storage mode.
#[test]
fn batch_size_does_not_change_output() {
    let mut first: Option<Vec<(String, Vec<u8>)>> = None;
    // (mode, threads, batch bytes): the default, one subject per batch, and batches that end
    // inside isoform groups and bins
    let runs = [
        (StorageMode::Encoded, 1, 8u64 << 20),
        (StorageMode::Encoded, 4, 1),
        (StorageMode::Encoded, 3, 64 << 10),
        (StorageMode::Raw, 4, 1),
        (StorageMode::Zstd, 2, 4 << 10),
    ];
    for (mode, threads, batch) in runs {
        let d = tempdir();
        let mut req = request(&txome(mode).spec(), Selection::all(), &d.path().join("db"));
        req.options.threads = threads;
        let mut p = gtars_rgblast::api::plan_request(&req, &mut |_| true).unwrap();
        p.write.batch_bytes = batch;
        let fetch = gtars_rgblast::store::StoreFiles::new(&p.store);
        gtars_rgblast::blastdb::write_blastdb(&p.plan, &fetch, &p.write, &mut |_| true).unwrap();
        let files: Vec<(String, Vec<u8>)> = files_in(d.path())
            .into_iter()
            .map(|(n, b)| {
                if n.ends_with(".manifest.json") {
                    (n, manifest_sans_store(&b))
                } else {
                    (n, b)
                }
            })
            .collect();
        match &first {
            None => first = Some(files),
            Some(want) => {
                assert_eq!(files.len(), want.len());
                for ((n1, b1), (n2, b2)) in files.iter().zip(want) {
                    assert_eq!(n1, n2);
                    assert!(
                        b1 == b2,
                        "{n1} differs ({mode:?}, {threads} threads, batch {batch})"
                    );
                }
            }
        }
    }
}

/// The manifest without the store location (stores of different modes live in different folders).
fn manifest_sans_store(b: &[u8]) -> Vec<u8> {
    let mut v: serde_json::Value = serde_json::from_slice(b).unwrap();
    v["store"]["location"] = serde_json::Value::Null;
    serde_json::to_vec(&v).unwrap()
}

#[test]
fn every_sidecar_is_written() {
    let f = encoded();
    let d = tempdir();
    let full = d.path().join("full");
    let s = build(&f.spec(), Selection::all(), &full);
    let names: BTreeSet<String> = s
        .sidecars
        .iter()
        .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
        .collect();
    for n in [
        "ids.tsv",
        "collections.tsv",
        "members.tsv",
        "sanitized.tsv",
        "skipped.tsv",
        "oids.tsv",
        "manifest.json",
    ] {
        assert!(names.contains(&format!("full.{n}")), "{n}");
    }
}

/// Thousands of collections: a store of `n` tiny bins, each with its own alias and taxid.
/// Returns the plan build time.
fn many_collections(n: usize) -> f64 {
    let d = tempdir();
    let store = d.path().join("bins");
    let spec = TxomeSpec::many_bins(n);
    let digests = build_txome_store(&store, StorageMode::Encoded, Layout::Standalone, &spec);
    assert_eq!(digests.len(), n);
    let st = store.display().to_string();
    let t = Instant::now();
    let (_, plan) = plan_for(&st, &Selection::all());
    let plan_s = t.elapsed().as_secs_f64();
    assert_eq!(plan.collections.len(), n);
    let out = d.path().join("db");
    build(&st, Selection::all(), &out);

    // collections.tsv: one row per bin, its alias and taxid
    let taxid_of = |alias: &str| 100_000 + 7 * (alias[3..].parse::<i32>().unwrap() - 1);
    let rows = sidecars::read_collections(&with_ext(&out, "collections.tsv")).unwrap();
    assert_eq!(rows.len(), n);
    let mut n_by_taxid = BTreeMap::new();
    for r in &rows {
        assert_eq!(r.aliases.len(), 1);
        let (ns, alias) = &r.aliases[0];
        assert_eq!(ns, "name");
        assert_eq!(digests[alias], r.collection);
        assert_eq!(r.taxid, Some(taxid_of(alias)));
        n_by_taxid.insert(taxid_of(alias), r.n_sequences);
    }
    // .nto and .ntf: every taxid with the OIDs of its bin (random contigs: no sharing)
    let ntf = lmdb::read_ntf(&with_ext(&out, "ntf")).unwrap();
    assert_eq!(
        ntf.iter().map(|(t, _)| *t).collect::<Vec<_>>(),
        n_by_taxid.keys().copied().collect::<Vec<_>>()
    );
    let nto = std::fs::read(with_ext(&out, "nto")).unwrap();
    let u32_at =
        |off: u64| u32::from_le_bytes(nto[off as usize..off as usize + 4].try_into().unwrap());
    for (t, off) in &ntf {
        assert_eq!(u32_at(*off), n_by_taxid[t], "taxid {t}");
    }
    // per-collection views: one per bin, each with its bin's sequences
    let s = view(
        &ViewRequest {
            db: out.clone(),
            kind: ViewKind::PerCollection,
            name: None,
            selection: None,
            out_dir: None,
            title: None,
            full_db_stats: false,
        },
        &mut |_| true,
    )
    .unwrap();
    assert_eq!(s.views.len(), n);
    let by_digest: BTreeMap<&str, &sidecars::CollectionMeta> =
        rows.iter().map(|r| (r.collection.as_str(), r)).collect();
    for v in &s.views {
        assert_eq!(v.collections.len(), 1);
        let r = by_digest[v.collections[0].as_str()];
        assert_eq!(
            (v.nseq, v.length),
            (r.n_sequences as u64, r.total_length),
            "{}",
            v.name
        );
    }
    plan_s
}

#[test]
fn many_collections_reduced() {
    many_collections(200);
}

/// The full count, and plan build time linear in it:
/// `cargo test -p gtars-rgblast --test txome -- --ignored many_collections_full`
#[test]
#[ignore]
fn many_collections_full() {
    let small = many_collections(200);
    let big = many_collections(2_000);
    assert!(
        big < 25.0 * small.max(0.01),
        "plan build {big:.3}s for 2,000 collections vs {small:.3}s for 200"
    );
}
