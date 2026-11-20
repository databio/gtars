//! Search, origin and the byte readers on packed stores (`RefgetStore::on_disk_packed`) give
//! the same output as the same store standalone. No BLAST+ needed: every packed run is compared
//! with the standalone run, which the golden tests compare with blastn.

mod common;

use common::golden::*;
use common::http;
use common::origin::{ORIGIN_FIXTURES, OriginFixture, origin_store, origin_store_in};
use common::txome::{txome, txome_in};
use common::*;
use gtars_refget::store::StorageMode;
use gtars_rgblast::api::{self, OriginRequest};
use gtars_rgblast::engine::options::Task;
use gtars_rgblast::origin::{ModeChoice, OriginOptions};
use gtars_rgblast::report::output;
use gtars_rgblast::run::EngineChoice;
use gtars_rgblast::select::nstats::NstatsOptions;
use gtars_rgblast::select::{ResolveOptions, Selection};
use gtars_rgblast::store::access::MappedFiles;
use gtars_rgblast::{NstatsRequest, OpenedStore, VerifyRequest};
use std::path::{Path, PathBuf};
use std::sync::atomic::Ordering;

const MODES: [StorageMode; 3] = [StorageMode::Encoded, StorageMode::Raw, StorageMode::Zstd];

/// rgblastn in process on `store`: rendered output and hit count.
fn search(store: &str, q: &Path, outfmt: &str, args: &[&str]) -> (String, u64) {
    let o = run(&rgblastn_req(store, q, outfmt, args));
    let hits = o.results.summary.hits as u64;
    (render(&o), hits)
}

/// The fixture queries (as `threads.rs` builds them).
fn fixture_queries(dir: &Path) -> PathBuf {
    let fixtures = Path::new(FIXTURES);
    let q = dir.join("q_fixture.fa");
    let fastas = ["asmA.fa", "asmB.fa", "asmD.fa"].map(|f| fixtures.join(f));
    std::fs::write(&q, queries_from(&fastas, 6, 42)).unwrap();
    q
}

/// A fresh N-count cache: a shared one would already hold the counts and hide the reads.
fn fresh_nstats(dir: &Path) -> NstatsOptions {
    NstatsOptions {
        cache_dir: Some(dir.to_path_buf()),
        ..Default::default()
    }
}

#[test]
fn packed_rgblastn_equals_standalone() {
    let dir = tempdir();
    let q = fixture_queries(dir.path());
    for mode in MODES {
        let packed = fixture_in(mode, Layout::Packed);
        assert_packed(&packed.path);
        let standalone = fixture(mode);
        for (task, mts) in [
            (Task::Megablast, 500),
            (Task::Megablast, 1),
            (Task::Blastn, 500),
            (Task::BlastnShort, 500),
        ] {
            // btop and qseq/sseq read subject letters back (`blastna_window`)
            for outfmt in ["6 std stitle btop", "6 std qseq sseq"] {
                for threads in [1, 4] {
                    let (mts, nt) = (mts.to_string(), threads.to_string());
                    let args = [
                        "-task",
                        task.name(),
                        "-max_target_seqs",
                        &mts,
                        "-num_threads",
                        &nt,
                    ];
                    let (want, hits) = search(&standalone.spec(), &q, outfmt, &args);
                    assert!(hits > 0, "{mode:?} {task:?}: vacuous");
                    let (got, _) = search(&packed.spec(), &q, outfmt, &args);
                    assert!(
                        got == want,
                        "{mode:?} {task:?} mts {mts} {outfmt:?} {threads} threads: packed differs"
                    );
                }
            }
        }
    }
}

/// Long subjects: the first search stage reads them piece by piece at their pack offset, the
/// traceback maps the span at an unaligned offset, and decoding once (`decode_whole`) reads
/// them whole.
#[test]
fn packed_long_subjects_equal_standalone() {
    let (pd, sd) = (tempdir(), tempdir());
    let (packed, q) = long_subject_store(pd.path(), Layout::Packed);
    let (standalone, _) = long_subject_store(sd.path(), Layout::Standalone);
    assert_packed(Path::new(&packed));
    // some long record sits at an offset that is not page aligned
    let (st, plan) = plan_for(&packed, &Selection::all());
    let offsets: Vec<(u64, Option<u64>)> = plan
        .subjects
        .iter()
        .map(|s| {
            let span = st.store.stored_span(&s.digest).unwrap();
            (span.offset, span.len)
        })
        .collect();
    assert!(
        offsets
            .iter()
            .any(|&(o, l)| o % 4096 != 0 && l.unwrap() >= MappedFiles::MMAP_MIN),
        "{offsets:?}"
    );
    drop(st);
    let stores = |root: &str| {
        let root = Path::new(root).parent().unwrap().to_path_buf();
        ["store", "raw", "zstd"].map(|n| root.join(n).display().to_string())
    };
    for (p, s) in stores(&packed).iter().zip(stores(&standalone).iter()) {
        assert_packed(Path::new(p));
        for task in [Task::Megablast, Task::Blastn] {
            for threads in [1, 3] {
                let nt = threads.to_string();
                let args = ["-task", task.name(), "-num_threads", &nt];
                let (want, hits) = search(s, &q, "6 std stitle btop", &args);
                assert!(hits > 0, "{s} {task:?}");
                let (got, _) = search(p, &q, "6 std stitle btop", &args);
                assert!(got == want, "{p} {task:?} {threads} threads");
            }
            // decoded once for all query batches
            let one = |store: &str| {
                let args = ["-task", task.name(), "-num_threads", "4"];
                let mut req = rgblastn_req(store, &q, "6 std stitle btop", &args);
                req.params.query_batch_size = Some(400);
                req.params.decode_cache_bytes = Some(u64::MAX);
                render(&run(&req))
            };
            assert!(one(p) == one(s), "{p} {task:?}: decode_whole");
        }
    }
}

/// Thousands of small subjects, read with one positioned read each from shared pack handles.
#[test]
fn packed_txome_equals_standalone() {
    let packed = txome_in(StorageMode::Encoded, Layout::Packed);
    assert_packed(&packed.path);
    let standalone = txome(StorageMode::Encoded);
    let src = packed.path.with_extension("src");
    let mut fastas: Vec<PathBuf> = std::fs::read_dir(&src)
        .unwrap()
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|e| e == "fa"))
        .collect();
    fastas.sort();
    fastas.push(src.join("bins/bin01.fa"));
    let dir = tempdir();
    let q = dir.path().join("q_txome.fa");
    std::fs::write(&q, queries_from(&fastas, 12, 7)).unwrap();
    for (task, mts) in [
        (Task::Megablast, 500),
        (Task::Megablast, 5),
        (Task::Blastn, 500),
    ] {
        let mts = mts.to_string();
        let args = [
            "-task",
            task.name(),
            "-max_target_seqs",
            &mts,
            "-num_threads",
            "4",
        ];
        let (want, hits) = search(&standalone.spec(), &q, "6 std stitle btop", &args);
        assert!(hits > 0, "{task:?}");
        let (got, _) = search(&packed.spec(), &q, "6 std stitle btop", &args);
        assert!(got == want, "{task:?} mts {mts}: packed differs");
    }
}

fn origin_req(f: &OriginFixture, mode: ModeChoice, cache: &Path) -> OriginRequest {
    OriginRequest {
        store: f.store.display().to_string(),
        selection: Selection::all(),
        engine: EngineChoice::Rgblastn,
        cache_dir: Some(cache.to_path_buf()),
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

#[test]
fn packed_origin_equals_standalone() {
    let packed = origin_store_in(Layout::Packed);
    assert_packed(&packed.store);
    let standalone = origin_store();
    let expected_exhaustive =
        std::fs::read_to_string(Path::new(ORIGIN_FIXTURES).join("expected_origin_exhaustive.tsv"))
            .unwrap();
    for mode in [ModeChoice::Staged, ModeChoice::Exhaustive] {
        let (pc, sc) = (tempdir(), tempdir());
        let got = api::origin(&origin_req(packed, mode, pc.path()), &mut |_| true).unwrap();
        let want = api::origin(&origin_req(standalone, mode, sc.path()), &mut |_| true).unwrap();
        assert!(!got.report.queries.is_empty(), "{mode:?}");
        assert_eq!(
            output::tsv(&got.report),
            output::tsv(&want.report),
            "{mode:?}"
        );
        if mode == ModeChoice::Exhaustive {
            assert_eq!(output::tsv(&got.report), expected_exhaustive);
        }
    }
    // search + report read the store from the search's meta
    use gtars_rgblast::api::{ReportRequest, SearchRequest, report};
    let t = tempdir();
    let out = t.path().join("hits.tsv");
    api::search(
        &SearchRequest {
            store: packed.store.display().to_string(),
            selection: Selection::all(),
            extra_fields: true,
            args: vec![
                "-query".into(),
                packed.queries.display().to_string(),
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
    let r = report(&ReportRequest {
        hits: vec![out],
        query: Some(packed.queries.clone()),
        fields: gtars_rgblast::report::parse::DEFAULT_FIELDS.into(),
        taxdump: Some(Path::new(ORIGIN_FIXTURES).join("taxdump_mini")),
        ..Default::default()
    })
    .unwrap();
    assert_eq!(output::tsv(&r), expected_exhaustive);
}

/// The N-count scan (`--max-n-fraction`, `nstats`), the database writer and `verify` read
/// whole records (`StoreFiles::read`) from pack spans.
#[test]
fn packed_selection_nstats_verify() {
    let packed = taxa_in(StorageMode::Encoded, Layout::Packed);
    assert_packed(&packed.path);
    let standalone = taxa(StorageMode::Encoded);
    let resolve = |f: &Fixture| {
        let cache = tempdir();
        let mut st = OpenedStore::open(&f.spec(), None).unwrap();
        let ro = ResolveOptions {
            nstats: fresh_nstats(cache.path()),
            threads: 2,
            ..Default::default()
        };
        let sel = Selection::all().with_max_n_fraction(0.25);
        api::resolve_in_store(&mut st, &sel, &ro, &mut |_| true)
            .unwrap()
            .0
    };
    let (got, want) = (resolve(packed), resolve(standalone));
    assert_eq!(got.collections, want.collections);
    assert_eq!(got.collections.len(), 4);
    assert_eq!(got.excluded.len(), 1);
    assert_eq!(got.excluded[0].display, "human_masked");
    assert!((got.excluded[0].n_fraction - 198.0 / 580.0).abs() < 1e-12);
    assert_eq!(got.excluded, want.excluded);

    let nstats = |f: &Fixture| {
        let cache = tempdir();
        api::nstats(
            &NstatsRequest {
                store: f.spec(),
                cache_dir: None,
                fetch: Default::default(),
                selection: Selection::all(),
                nstats: fresh_nstats(cache.path()),
                taxid_map: None,
                threads: 2,
            },
            &mut |_| true,
        )
        .unwrap()
    };
    let (got, want) = (nstats(packed), nstats(standalone));
    assert_eq!(got.collections, want.collections);
    assert_eq!(got.sequences, want.sequences);
    assert!(got.collections.iter().any(|c| c.n > 0));

    let (pd, sd) = (tempdir(), tempdir());
    for (f, d) in [(packed, &pd), (standalone, &sd)] {
        let ncache = tempdir();
        let mut req = request(&f.spec(), Selection::all(), &d.path().join("db"));
        req.options.nstats = fresh_nstats(ncache.path());
        build_req(&req);
    }
    same_database(pd.path(), sd.path());
    let rep = api::verify(
        &VerifyRequest {
            db: pd.path().join("db"),
            store: packed.spec(),
            cache_dir: None,
            threads: 2,
            fetch: Default::default(),
        },
        &mut |_| true,
    )
    .unwrap();
    assert!(rep.ok(), "{rep:?}");
    assert!(rep.sequences_checked > 0);
    // the stored letters, through StoreFiles::read
    let (st, plan) = plan_for(&packed.spec(), &Selection::all());
    drop(st);
    for s in plan.subjects.iter().take(5) {
        let n = s.length as usize;
        assert_eq!(
            slice(&packed.spec(), &s.digest, 1, n),
            slice(&standalone.spec(), &s.digest, 1, n)
        );
    }
}

/// A remote packed store: the subjects' pack spans come down with ranged GETs into the shared
/// cache as `.seq` files, and the search maps those.
#[test]
fn packed_remote_rgblastn_equals_local() {
    let packed = fixture_in(StorageMode::Encoded, Layout::Packed);
    let dir = tempdir();
    let q = dir.path().join("q.fa");
    let fixtures = Path::new(FIXTURES);
    std::fs::write(
        &q,
        queries_from(&[fixtures.join("asmA.fa"), fixtures.join("asmB.fa")], 6, 99),
    )
    .unwrap();
    let local = ours(
        &fixture(StorageMode::Encoded).spec(),
        &q,
        "6 std stitle btop",
        &["-num_threads", "2"],
    );
    assert!(local.lines().count() > 4);
    let server = http::serve(packed.path.parent().unwrap());
    let url = format!("{}/fixture", server.url);
    let cache = dir.path().join("cache");
    let remote_run = || {
        let mut req = rgblastn_req(&url, &q, "6 std stitle btop", &["-num_threads", "2"]);
        req.cache_dir = Some(cache.clone());
        req.remote.fetch_jobs = 4;
        render(&run(&req))
    };
    let first = remote_run();
    assert!(
        first == local,
        "remote packed output differs from the local run"
    );
    let mut seqs = 0;
    let mut stack = vec![cache.join("sequences")];
    while let Some(d) = stack.pop() {
        for e in std::fs::read_dir(&d).into_iter().flatten().flatten() {
            let p = e.path();
            if p.is_dir() {
                stack.push(p);
            } else if p.extension().is_some_and(|x| x == "seq") {
                seqs += 1;
            }
        }
    }
    assert!(seqs > 0, "no .seq files in the cache");
    assert!(server.ranged.load(Ordering::SeqCst) > 0, "no ranged GETs");
    let requests_after_first = server.requests.load(Ordering::SeqCst);
    assert!(remote_run() == local);
    let delta = server.requests.load(Ordering::SeqCst) - requests_after_first;
    assert!(
        delta < requests_after_first / 2,
        "the rerun fetched sequences again: {delta} requests after {requests_after_first}"
    );
}

#[test]
fn packed_rgblastn_binary() {
    let packed = fixture_in(StorageMode::Encoded, Layout::Packed);
    let dir = tempdir();
    let q = fixture_queries(dir.path());
    let want = ours(
        &fixture(StorageMode::Encoded).spec(),
        &q,
        "6 std stitle",
        &[],
    );
    assert!(!want.is_empty());
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_rgblastn"))
        .args([
            "--store",
            &packed.spec(),
            "--all",
            "-query",
            q.to_str().unwrap(),
            "-outfmt",
            "6 std stitle",
            "--cache-dir",
            dir.path().join("cache").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(String::from_utf8_lossy(&out.stdout) == want);
}
