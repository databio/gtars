//! Threads, work items and remote access (plan 4 steps 25 and 26) cannot change the output.
//! No BLAST+ needed: every run is compared with the one-thread, default-batch run of the same
//! search; the golden tests compare that run with blastn.

mod common;

use common::golden::*;
use common::http;
use common::txome::txome;
use common::*;
use gtars_refget::store::StorageMode;
use gtars_rgblast::engine::options::Task;
use gtars_rgblast::engine::search::{build_work_items, debug_shuffle, default_batch_letters};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::Ordering;

/// Queries over sequences of a FASTA set: exact slices, substitution mutants and reverse
/// complements of up to `take` sequences, plus one random no-hit query.
fn queries_from(fastas: &[PathBuf], take: usize, seed: u64) -> Vec<u8> {
    let mut rng = Lcg(seed);
    let mut fa = Vec::new();
    let mut n = 0;
    let mut picked = 0;
    for p in fastas {
        for (name, s) in read_fasta(p) {
            if picked >= take {
                break;
            }
            if s.len() < 300 {
                continue;
            }
            picked += 1;
            let a = (rng.next() as usize) % (s.len() / 3);
            let b = (a + 150 + (rng.next() as usize) % 250).min(s.len());
            n += 1;
            push_query(&mut fa, &format!("q{n:03}_{name}_exact"), &s[a..b]);
            let mut m = s[a..b].to_vec();
            for _ in 0..(m.len() / 25) {
                let p = (rng.next() as usize) % m.len();
                m[p] = b"ACGT"[(rng.next() % 4) as usize];
            }
            n += 1;
            push_query(&mut fa, &format!("q{n:03}_{name}_mut"), &m);
            if picked % 2 == 0 {
                n += 1;
                push_query(
                    &mut fa,
                    &format!("q{n:03}_{name}_rc"),
                    &rc_letters(&s[a..b]),
                );
            }
        }
    }
    let random: Vec<u8> = (0..400)
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect();
    push_query(&mut fa, "q_random no hit expected", &random);
    fa
}

/// The shuffle switch is process-wide; runs that set it are serialised here.
static SHUFFLE: Mutex<()> = Mutex::new(());

/// One search: `threads`, the letters per work item (`SearchParams::batch_letters`, None =
/// default), and whether the work lists are shuffled (`RGBLASTN_DEBUG_SHUFFLE`).
fn run(
    store: &str,
    query: &Path,
    task: Task,
    max_target_seqs: u32,
    threads: usize,
    batch_letters: Option<u64>,
    shuffle: Option<u64>,
) -> (String, gtars_rgblast::engine::search::SearchSummary) {
    let (mts, nt) = (max_target_seqs.to_string(), threads.to_string());
    let args = [
        "-task",
        task.name(),
        "-max_target_seqs",
        &mts,
        "-num_threads",
        &nt,
    ];
    let mut req = rgblastn_req(store, query, "6 std stitle btop", &args);
    req.params.batch_letters = batch_letters;
    let _g = SHUFFLE.lock().unwrap_or_else(|e| e.into_inner());
    // SAFETY: the test binary's only environment writes happen under SHUFFLE; the engine
    // reads the variable once per query batch on this thread.
    match shuffle {
        Some(seed) => unsafe { std::env::set_var("RGBLASTN_DEBUG_SHUFFLE", seed.to_string()) },
        None => unsafe { std::env::remove_var("RGBLASTN_DEBUG_SHUFFLE") },
    }
    let out = common::golden::run(&req);
    unsafe { std::env::remove_var("RGBLASTN_DEBUG_SHUFFLE") };
    let summary = out.results.summary.clone();
    (render(&out), summary)
}

/// Thread counts, batch sizes and work order to vary against the one-thread default run.
const VARIANTS: [(usize, Option<u64>, Option<u64>); 7] = [
    (4, None, None),
    (16, None, None),
    (1, Some(250_000), None),
    (4, Some(250_000), Some(7)),
    (16, Some(1_024_000_000), Some(11)),
    (3, Some(1), Some(3)),
    (1, Some(1), Some(5)),
];

/// `-max_target_seqs` below the number of hitting subjects engages the hit-list heap
/// (`Blast_HitListUpdate`'s replacement path) in the preliminary and the traceback stage.
fn check_store(label: &str, store: &str, query: &Path, tasks: &[(Task, u32)]) {
    for (task, mts) in tasks {
        let (base, summary) = run(store, query, *task, *mts, 1, None, None);
        assert!(
            summary.hits > 0,
            "{label} {task:?}: the base run found nothing; the test is vacuous"
        );
        for (threads, batch, shuffle) in VARIANTS {
            let (mine, s) = run(store, query, *task, *mts, threads, batch, shuffle);
            assert!(
                mine == base,
                "{label} {task:?} mts {mts}: {threads} threads, batch {batch:?}, shuffle {shuffle:?} changed the output ({} vs {} hits)",
                s.hits,
                summary.hits
            );
            assert_eq!(
                s.good_init_extends, summary.good_init_extends,
                "{label} {task:?}"
            );
            assert_eq!(s.batches, summary.batches, "{label} {task:?}");
        }
    }
}

/// `batch_size_does_not_change_output` (plan 4 step 25): threads {1, 3, 4, 16}, letters per
/// work item {1, 250,000, default, 1,024,000,000} and shuffled work-item and traceback order
/// give the same bytes on the fixture and the `txome` stores.
#[test]
fn threads_and_batches_do_not_change_output() {
    let f = encoded();
    let dir = tempdir();
    let fixtures = Path::new(FIXTURES);
    let q = dir.path().join("q_fixture.fa");
    std::fs::write(
        &q,
        queries_from(
            &[
                fixtures.join("asmA.fa"),
                fixtures.join("asmB.fa"),
                fixtures.join("asmD.fa"),
            ],
            6,
            42,
        ),
    )
    .unwrap();
    check_store(
        "fixture",
        &f.spec(),
        &q,
        &[
            (Task::Megablast, 500),
            (Task::Megablast, 1),
            (Task::Blastn, 500),
            (Task::BlastnShort, 500),
        ],
    );

    let t = txome(StorageMode::Encoded);
    let src = t.path.with_extension("src");
    let mut fastas: Vec<PathBuf> = std::fs::read_dir(&src)
        .unwrap()
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|e| e == "fa"))
        .collect();
    fastas.sort();
    let mut bins: Vec<PathBuf> = std::fs::read_dir(src.join("bins"))
        .unwrap()
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|e| e == "fa"))
        .collect();
    bins.sort();
    fastas.extend(bins.into_iter().take(2));
    let q = dir.path().join("q_txome.fa");
    std::fs::write(&q, queries_from(&fastas, 12, 7)).unwrap();
    check_store(
        "txome",
        &t.spec(),
        &q,
        &[
            (Task::Megablast, 500),
            (Task::Megablast, 1),
            (Task::Megablast, 5),
            (Task::Blastn, 500),
            (Task::Blastn, 2),
        ],
    );
}

/// A store with subjects longer than `MAX_DBSEQ_LEN` (5 Mbp), so the search decodes them chunk
/// by chunk: a 12.4 Mbp sequence (3 chunks) with planted copies of the query material across
/// the chunk boundaries and N runs (so it is 3-bit, its second chunk starts mid-group and the
/// random fill carries over from chunk to chunk), a 5.2 Mbp one (2 chunks) and a short one; the
/// same sequences also as raw and zstd stores beside it (`raw/`, `zstd/`).
fn long_subject_store(dir: &Path) -> (String, PathBuf) {
    let mut rng = Lcg(2024);
    let mut seq =
        |n: usize| -> Vec<u8> { (0..n).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect() };
    let unit: Vec<u8> = seq(900);
    let mut long = seq(12_400_000);
    let step = 5_000_000 - 100;
    // copies inside each chunk, one straddling every chunk boundary and one in the overlap
    for at in [
        1_000,
        2_500_000,
        step - 450,
        step + 20,
        7_000_000,
        2 * step - 450,
        9_900_000,
        12_398_000,
    ] {
        let end = (at + unit.len()).min(long.len());
        long[at..end].copy_from_slice(&unit[..end - at]);
    }
    for (at, n) in [
        (500, 300),
        (3_000_003, 5_000),
        (6_000_001, 100),
        (11_000_000, 7),
    ] {
        long[at..at + n].fill(b'N');
    }
    let mut mid = seq(5_200_000);
    for at in [3_000, step - 300, 5_100_000] {
        mid[at..at + unit.len()].copy_from_slice(&unit);
    }
    let mut short = seq(20_000);
    short[5_000..5_900].copy_from_slice(&unit);
    let mut fa = Vec::new();
    push_query(&mut fa, "longA planted copies", &long);
    push_query(&mut fa, "midB two chunks", &mid);
    push_query(&mut fa, "shortC", &short);
    let fasta = dir.join("long.fa");
    std::fs::write(&fasta, &fa).unwrap();
    let store = dir.join("store");
    build_store_from_fasta(&store, &fasta, "long", StorageMode::Encoded);
    for (name, mode) in [("raw", StorageMode::Raw), ("zstd", StorageMode::Zstd)] {
        build_store_from_fasta(&dir.join(name), &fasta, "long", mode);
    }
    let mut qfa = Vec::new();
    push_query(&mut qfa, "u_exact", &unit);
    let mut m = unit.clone();
    for p in (0..m.len()).step_by(37) {
        m[p] = if m[p] == b'A' { b'C' } else { b'A' };
    }
    push_query(&mut qfa, "u_mut", &m);
    push_query(&mut qfa, "u_rc", &rc_letters(&unit));
    push_query(&mut qfa, "u_half", &unit[200..650]);
    // the boundary region itself, so a hit spans the overlap
    push_query(&mut qfa, "boundary", &long[step - 700..step + 300]);
    let q = dir.join("q_long.fa");
    std::fs::write(&q, &qfa).unwrap();
    (store.display().to_string(), q)
}

#[test]
fn chunked_subjects_do_not_depend_on_threads() {
    let dir = tempdir();
    let (store, q) = long_subject_store(dir.path());
    let (_, plan) = plan_for(&store, &gtars_rgblast::select::Selection::all());
    let lengths: Vec<u32> = plan.subjects.iter().map(|s| s.length).collect();
    let items = build_work_items(&lengths, default_batch_letters(17_620_000, 1));
    assert_eq!(items, vec![0..1, 1..2, 2..3]);
    check_store(
        "long",
        &store,
        &q,
        &[
            (Task::Megablast, 500),
            (Task::Megablast, 1),
            (Task::Blastn, 500),
        ],
    );
    let (base, s) = run(&store, &q, Task::Megablast, 500, 1, None, None);
    // every query but the random one hits all three subjects, and the boundary query's hit
    // on longA spans the chunk overlap
    assert!(s.hits >= 3 * 4, "{} hits\n{base}", s.hits);
    assert!(
        base.lines().any(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[0] == "boundary" && f[3] == "1000" && f[8] == "4999201" && f[12].contains("longA")
        }),
        "{base}"
    );
    // raw records are read piece by piece as packed ones are, zstd ones decompressed whole
    for name in ["raw", "zstd"] {
        let other = dir.path().join(name).display().to_string();
        let (mine, _) = run(&other, &q, Task::Megablast, 500, 3, None, None);
        assert!(mine == base, "the {name} store changed the output");
    }
}

/// `build_work_items`: grouping by letters and count, long subjects alone, longest first;
/// `debug_shuffle` only permutes.
#[test]
fn work_items_group_by_letters_and_split_long_subjects() {
    let lens = [100, 200, 300, 5_000_000, 5_000_001, 50, 10_000_000, 7, 8];
    let items = build_work_items(&lens, 600);
    assert_eq!(items, vec![6..7, 4..5, 3..4, 0..3, 5..6, 7..9]);
    // exactly 600 letters fit; 601 do not
    let items = build_work_items(&[300, 300, 1], 600);
    assert_eq!(items, vec![0..2, 2..3]);
    // the subject cap
    let items = build_work_items(&[1; 10_000], u64::MAX);
    assert_eq!(items, vec![0..4096, 4096..8192, 8192..10_000]);
    // the default: total / (16 * threads) within [1 Mbp, 64 Mbp]
    assert_eq!(default_batch_letters(3_088_286_401, 1), 64_000_000);
    assert_eq!(default_batch_letters(3_088_286_401, 4), 48_254_475);
    assert_eq!(default_batch_letters(3_088_286_401, 16), 12_063_618);
    assert_eq!(default_batch_letters(10_000_000, 16), 1_000_000);
    assert_eq!(default_batch_letters(0, 16), 1_000_000);
    // the shuffle is a permutation, and only with the variable set
    let mut v: Vec<u32> = (0..50).collect();
    let _g = SHUFFLE.lock().unwrap_or_else(|e| e.into_inner());
    unsafe { std::env::remove_var("RGBLASTN_DEBUG_SHUFFLE") };
    debug_shuffle(&mut v);
    assert_eq!(v, (0..50).collect::<Vec<_>>());
    unsafe { std::env::set_var("RGBLASTN_DEBUG_SHUFFLE", "9") };
    debug_shuffle(&mut v);
    unsafe { std::env::remove_var("RGBLASTN_DEBUG_SHUFFLE") };
    assert_ne!(v, (0..50).collect::<Vec<_>>());
    v.sort();
    assert_eq!(v, (0..50).collect::<Vec<_>>());
}

/// Remote `rgblastn` (plan 4 step 26) over the local HTTP fixture: the subjects are prefetched
/// into the shared cache with ranged GETs and the search maps the cached files, so the bytes
/// equal the local run's; a rerun on the same cache downloads no sequence again.
#[test]
fn remote_rgblastn_equals_local() {
    let f = encoded();
    let dir = tempdir();
    let q = dir.path().join("q.fa");
    let fixtures = Path::new(FIXTURES);
    std::fs::write(
        &q,
        queries_from(&[fixtures.join("asmA.fa"), fixtures.join("asmB.fa")], 6, 99),
    )
    .unwrap();
    let local = ours(&f.spec(), &q, "6 std stitle", &["-num_threads", "2"]);
    assert!(local.lines().count() > 4);
    let server = http::serve(f.path.parent().unwrap());
    let url = format!("{}/fixture", server.url);
    let cache = dir.path().join("cache");
    let remote_run = || {
        let mut req = rgblastn_req(&url, &q, "6 std stitle", &["-num_threads", "2"]);
        req.cache_dir = Some(cache.clone());
        req.remote.fetch_jobs = 4;
        render(&common::golden::run(&req))
    };
    let first = remote_run();
    assert!(first == local, "remote output differs from the local run");
    let requests_after_first = server.requests.load(Ordering::SeqCst);
    let second = remote_run();
    assert!(second == local);
    let delta = server.requests.load(Ordering::SeqCst) - requests_after_first;
    assert!(
        delta < requests_after_first / 2,
        "the rerun fetched sequences again: {delta} requests after {requests_after_first}"
    );
    // the binary: the same bytes through -num_threads
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_rgblastn"))
        .args([
            "--store",
            &url,
            "--all",
            "-query",
            q.to_str().unwrap(),
            "-outfmt",
            "6 std stitle",
            "-num_threads",
            "3",
            "--cache-dir",
            cache.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(String::from_utf8_lossy(&out.stdout) == local);
}
