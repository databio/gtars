//! rgblastn vs BLAST+ 2.17.0 (feature `golden`) for `-task blastn` and `-task blastn-short`:
//! byte-identical tabular output on the fixture store over the flag variants of the megablast
//! suite, the NA lookup table (long masked query, many short queries), the transcriptome
//! fixture, blastn's query batching, and the statistics footers of every task; `#[ignore]`
//! hg38 run when `GTARS_RGBLAST_STORE_JUNGLE` is set.

mod common;

use common::blast::Blast;
use common::golden::*;
use common::txome::txome;
use common::*;
use gtars_refget::store::StorageMode;
use gtars_rgblast::engine::search::SearchSummary;
use std::path::Path;

/// The megablast suite's query classes (same seed, same slices), plus heavier mutants that
/// only the gapped DP finds: 15 % substitutions and indel clusters.
fn fixture_queries() -> Vec<u8> {
    let seqs = read_fasta(&Path::new(FIXTURES).join("asmA.fa"));
    let mut rng = Lcg(42);
    let mut fa = Vec::new();
    let mut n = 0;
    let mut push = |fa: &mut Vec<u8>, name: &str, s: &[u8]| {
        n += 1;
        push_query(fa, &format!("q{n:02}_{name} test query"), s);
    };
    for (i, (name, s)) in seqs.iter().enumerate().take(5) {
        let len = s.len();
        if len < 120 {
            continue;
        }
        let a = (rng.next() as usize) % (len / 3);
        let b = (a + 150 + (rng.next() as usize) % 200).min(len);
        push(&mut fa, &format!("{name}_exact"), &s[a..b]);
        let mut m = s[a..b].to_vec();
        for _ in 0..(m.len() / 20) {
            let p = (rng.next() as usize) % m.len();
            m[p] = b"ACGT"[(rng.next() % 4) as usize];
        }
        push(&mut fa, &format!("{name}_mut"), &m);
        if i % 2 == 1 {
            push(&mut fa, &format!("{name}_rc"), &rc_letters(&s[a..b]));
        }
        if i == 2 {
            let mut withn = s[a..b].to_vec();
            for p in (10..withn.len().min(60)).step_by(7) {
                withn[p] = b'N';
            }
            push(&mut fa, &format!("{name}_n"), &withn);
            push(&mut fa, &format!("{name}_short"), &s[a..a + 40]);
        }
        // 15 % substitutions: below megablast's reach, within blastn's
        let mut m15 = s[a..b].to_vec();
        for _ in 0..(m15.len() * 15 / 100) {
            let p = (rng.next() as usize) % m15.len();
            m15[p] = b"ACGT"[(rng.next() % 4) as usize];
        }
        push(&mut fa, &format!("{name}_mut15"), &m15);
    }
    let random: Vec<u8> = (0..500)
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect();
    push(&mut fa, "random", &random);
    let c = read_fasta(&Path::new(FIXTURES).join("asmC.fa"));
    let c1 = &c.iter().find(|(n, _)| n == "c1").unwrap().1;
    push(&mut fa, "c1_exact", &c1[100..400]);
    let mut g = c1[50..450].to_vec();
    g.drain(100..103);
    g.insert(250, b'A');
    g.insert(250, b'C');
    g.drain(300..301);
    push(&mut fa, "c1_indel", &g);
    let chr1 = &seqs[0].1;
    let mut m = chr1[200..600].to_vec();
    for p in (0..m.len()).step_by(37) {
        m[p] = if m[p] == b'A' { b'C' } else { b'A' };
    }
    push(&mut fa, "chr1_mut_rc", &rc_letters(&m));
    // indel clusters every 60 bases on chr1 (chr2 of the fixture is almost all N)
    let mut ic = chr1[100..900].to_vec();
    let mut p = 60;
    while p + 70 < ic.len() {
        ic.drain(p..p + 2);
        ic.insert(p + 30, b'G');
        ic.insert(p + 30, b'T');
        ic.insert(p + 30, b'T');
        p += 60;
    }
    push(&mut fa, "chr2_indels", &ic);
    let lower: Vec<u8> = chr1[300..500].to_ascii_lowercase();
    push_query(&mut fa, "q_lower lower case", &lower);
    let u: Vec<u8> = chr1[400..620]
        .iter()
        .map(|&c| if c == b'T' { b'U' } else { c })
        .collect();
    push_query(&mut fa, "q_u RNA spelling", &u);
    push_query(
        &mut fa,
        "  q_spaced  defline | with | pipes ",
        &chr1[500..700],
    );
    push_query(&mut fa, "", &chr1[600..800]);
    fa
}

/// 40 kb of low complexity (dust masks it) followed by a real slice: few lookup words but
/// query offsets past 32767, so `BlastChooseNaLookupTable` picks the NA table for both tasks.
fn na_table_query() -> Vec<u8> {
    let seqs = read_fasta(&Path::new(FIXTURES).join("asmA.fa"));
    let chr1 = &seqs[0].1;
    let mut s = Vec::new();
    for _ in 0..20000 {
        s.extend_from_slice(b"AT");
    }
    s.extend_from_slice(&chr1[300..700]);
    let mut fa = Vec::new();
    push_query(&mut fa, "q_masked_long mostly dust", &s);
    push_query(&mut fa, "q_tail chr1 slice", &chr1[700..1000]);
    fa
}

/// 500 primer-sized slices (40 bases, both strands, a few with one mismatch): enough lookup
/// entries for `blastn-short` to need the NA table at width 7, stride 1.
fn short_queries() -> Vec<u8> {
    let seqs = read_fasta(&Path::new(FIXTURES).join("asmA.fa"));
    let mut rng = Lcg(7);
    let mut fa = Vec::new();
    let total: usize = seqs.iter().map(|(_, s)| s.len()).sum();
    for i in 0..500 {
        let mut pos = (rng.next() as usize) % (total - 40);
        let mut src = &seqs[0].1;
        for (_, s) in &seqs {
            if pos + 40 <= s.len() {
                src = s;
                break;
            }
            pos = pos.saturating_sub(s.len());
        }
        if pos + 40 > src.len() {
            pos = src.len() - 40;
        }
        let mut q = src[pos..pos + 40].to_vec();
        if i % 5 == 0 {
            let p = (rng.next() as usize) % 40;
            q[p] = if q[p] == b'A' { b'G' } else { b'A' };
        }
        if i % 3 == 0 {
            q = rc_letters(&q);
        }
        push_query(&mut fa, &format!("p{i:03}"), &q);
    }
    fa
}

/// rgblastn's output and search summary, with blastn's `BATCH_SIZE` as `query_batch_size`.
fn ours_with_summary(
    store: &str,
    query: &Path,
    args: &[&str],
    query_batch_size: Option<u32>,
) -> (String, SearchSummary) {
    let mut req = rgblastn_req(store, query, "6", args);
    req.params.query_batch_size = query_batch_size;
    let out = run(&req);
    (render(&out), out.results.summary.clone())
}

#[test]
fn blastn_fixture_byte_identical() {
    let f = encoded();
    let dir = tempdir();
    let work = dir.path();
    let blast = Blast::find();
    reference_db_all(&blast, &f.spec(), work);
    std::fs::write(work.join("q.fa"), fixture_queries()).unwrap();
    for (name, args) in [
        ("default", vec![]),
        ("evalue", vec!["-evalue", "1e-5"]),
        ("evalue1000", vec!["-evalue", "1000"]),
        ("plus", vec!["-strand", "plus"]),
        ("minus", vec!["-strand", "minus"]),
        ("nodust", vec!["-dust", "no"]),
        ("hard", vec!["-soft_masking", "false"]),
        ("mts1", vec!["-max_target_seqs", "1"]),
        ("w7", vec!["-word_size", "7"]),
        ("w16", vec!["-word_size", "16"]),
        ("pid", vec!["-perc_identity", "98"]),
        ("maxhsps", vec!["-max_hsps", "2"]),
    ] {
        let mut a = vec!["-task", "blastn"];
        a.extend(args.iter());
        let reference = check(&blast, &f.spec(), work, "q.fa", ALL_FIELDS, name, &a);
        assert!(
            !reference.is_empty() || name == "pid",
            "{name}: the reference produced no hits; the test is vacuous"
        );
    }
}

#[test]
fn blastn_short_fixture_byte_identical() {
    let f = encoded();
    let dir = tempdir();
    let work = dir.path();
    let blast = Blast::find();
    reference_db_all(&blast, &f.spec(), work);
    std::fs::write(work.join("q.fa"), fixture_queries()).unwrap();
    std::fs::write(work.join("short.fa"), short_queries()).unwrap();
    for (name, args) in [
        ("default", vec![]),
        ("evalue", vec!["-evalue", "1e-5"]),
        ("minus", vec!["-strand", "minus"]),
        ("mts1", vec!["-max_target_seqs", "1"]),
        ("dust", vec!["-dust", "yes"]),
        ("w11", vec!["-word_size", "11"]),
    ] {
        let mut a = vec!["-task", "blastn-short"];
        a.extend(args.iter());
        let reference = check(&blast, &f.spec(), work, "q.fa", ALL_FIELDS, name, &a);
        assert!(!reference.is_empty(), "{name}: vacuous");
    }
    // 500 short queries: the NA lookup table (width 7, stride 1)
    for (name, args) in [
        ("short_default", vec![]),
        ("short_mts1", vec!["-max_target_seqs", "1"]),
        ("short_w4", vec!["-word_size", "4"]),
    ] {
        let mut a = vec!["-task", "blastn-short"];
        a.extend(args.iter());
        let reference = check(
            &blast,
            &f.spec(),
            work,
            "short.fa",
            "6 std stitle",
            name,
            &a,
        );
        assert!(
            reference.lines().count() > if name == "short_default" { 500 } else { 100 },
            "{name}: too few hits to be meaningful"
        );
    }
}

/// The NA lookup table through the long masked query, for both tasks (width 8, stride 4
/// under blastn: the `s_BlastNaScanSubject_8_4` kernel; stride 21 under megablast).
#[test]
fn na_lookup_table_byte_identical() {
    let f = encoded();
    let dir = tempdir();
    let work = dir.path();
    let blast = Blast::find();
    reference_db_all(&blast, &f.spec(), work);
    std::fs::write(work.join("na.fa"), na_table_query()).unwrap();
    for (name, args) in [
        ("blastn", vec!["-task", "blastn"]),
        (
            "blastn_hard",
            vec!["-task", "blastn", "-soft_masking", "false"],
        ),
        ("megablast", vec!["-task", "megablast"]),
        ("blastn_w12", vec!["-task", "blastn", "-word_size", "12"]),
    ] {
        let reference = check(&blast, &f.spec(), work, "na.fa", ALL_FIELDS, name, &args);
        assert!(reference.lines().count() >= 2, "{name}: vacuous");
    }
}

/// Lambda, K, H and the search space of every query against blastn's `-outfmt 0` footer, for
/// every task.
#[test]
fn stats_footers() {
    let f = encoded();
    let dir = tempdir();
    let work = dir.path();
    let blast = Blast::find();
    reference_db_all(&blast, &f.spec(), work);
    std::fs::write(work.join("q.fa"), fixture_queries()).unwrap();
    for task in ["megablast", "blastn", "blastn-short"] {
        check_footer(&blast, &f.spec(), work, "q.fa", &["-task", task]);
    }
}

/// Plan 8's transcriptome fixture (2,826 short subjects in isoform groups): slices shared by
/// several transcripts under `-task blastn` with `-max_target_seqs` 1 and 10000, and
/// `blastn-short` on 30-base slices, byte for byte.
#[test]
fn blastn_txome_byte_identical() {
    let f = txome(StorageMode::Encoded);
    let dir = tempdir();
    let work = dir.path();
    let blast = Blast::find();
    reference_db_all(&blast, &f.spec(), work);
    let seqs = read_fasta(&work.join("mk/db.fa"));
    let mut fa = Vec::new();
    let mut short = Vec::new();
    let mut n = 0;
    for (i, (_, s)) in seqs.iter().enumerate() {
        if i % 97 != 3 || s.len() < 400 {
            continue;
        }
        n += 1;
        push_query(&mut fa, &format!("t{n:02}_slice"), &s[50..350]);
        let mut m = s[100..400].to_vec();
        for p in (0..m.len()).step_by(23) {
            m[p] = if m[p] == b'A' { b'G' } else { b'A' };
        }
        push_query(&mut fa, &format!("t{n:02}_mut"), &m);
        push_query(&mut short, &format!("s{n:02}"), &s[120..150]);
        if n % 2 == 0 {
            push_query(
                &mut short,
                &format!("s{n:02}_rc"),
                &rc_letters(&s[200..232]),
            );
        }
    }
    assert!(n >= 10, "{n} transcripts picked");
    std::fs::write(work.join("tx.fa"), &fa).unwrap();
    std::fs::write(work.join("txs.fa"), &short).unwrap();
    for (name, args) in [
        ("tx_default", vec!["-task", "blastn"]),
        ("tx_mts1", vec!["-task", "blastn", "-max_target_seqs", "1"]),
        (
            "tx_mts10000",
            vec!["-task", "blastn", "-max_target_seqs", "10000"],
        ),
        ("tx_evalue", vec!["-task", "blastn", "-evalue", "1e-20"]),
    ] {
        let reference = check(
            &blast,
            &f.spec(),
            work,
            "tx.fa",
            "6 std stitle",
            name,
            &args,
        );
        assert!(reference.lines().count() > n, "{name}: too few hits");
    }
    for (name, args) in [
        ("txs_default", vec!["-task", "blastn-short"]),
        (
            "txs_mts1",
            vec!["-task", "blastn-short", "-max_target_seqs", "1"],
        ),
    ] {
        let reference = check(
            &blast,
            &f.spec(),
            work,
            "txs.fa",
            "6 std stitle",
            name,
            &args,
        );
        assert!(!reference.is_empty(), "{name}: vacuous");
    }
}

/// blastn's adaptive query batching (`blastn_app.cpp`, `CBatchSizeMixer`). The store holds
/// four 10 kb slices of hg38 cut around weak q06 hits, at the same offset modulo 4 as in the
/// chromosome, and hg38's length gives hg38's cutoffs: `-dbsize 3088286401` for blastn,
/// `--stats-nseq 4 --stats-totlen 3088286401` for rgblastn (the same statistics: the store's
/// four sequences and hg38's length; the batch mixer reads the real total length under both). The ten hg38 queries read
/// 970 + 756 + 1089 + 330 + 1861 = 5006 letters for q01..q05, so the first batch (target
/// 1,000,000 hits / 200 = 5,000 letters) ends there and q06 starts the second. With everything
/// in one block (`BATCH_SIZE`) `-task blastn` reports four more q06 HSPs at the query's first
/// bases (`s_NuclUngappedExtend` scoring 4-base blocks across the sentinel before q06); the
/// default output does not, and rgblastn follows each mode byte for byte. The queries per
/// batch and the saved initial hits per batch (`CLocalBlast::GetNumExtensions`) are what an
/// instrumented blastn 2.17.0 printed in its batch loop. megablast runs the same batches; its
/// output does not depend on them here.
#[test]
fn query_batching_matches_blastn() {
    let blast = Blast::find();
    let dir = tempdir();
    let work = dir.path();
    let store = work.join("store");
    build_store_from_fasta(
        &store,
        &Path::new(FIXTURES).join("hg38_batch_slices.fa"),
        "hg38_slices",
        StorageMode::Encoded,
    );
    let spec = store.display().to_string();
    reference_db_all(&blast, &spec, work);
    std::fs::copy(
        Path::new(FIXTURES).join("hg38_batch_queries.fa"),
        work.join("q.fa"),
    )
    .unwrap();
    let q = work.join("q.fa");
    const DBSIZE: &str = "3088286401";
    for (task, hits) in [("blastn", [510, 385]), ("megablast", [14, 7])] {
        let args = ["-task", task, "--stats-nseq", "4", "--stats-totlen", DBSIZE];
        let bargs = [
            "-db", "mk/db", "-query", "q.fa", "-outfmt", "6", "-task", task, "-dbsize", DBSIZE,
        ];
        let reference = blast.ok(work, "blastn", &bargs);
        assert!(reference.lines().any(|l| l.starts_with("q06_")));
        let (mine, summary) = ours_with_summary(&spec, &q, &args, None);
        assert_same(
            &format!("{task} with blastn's query batches"),
            &reference,
            &mine,
        );
        assert_eq!(summary.batches, vec![5, 5], "{task}: queries per batch");
        assert_eq!(
            summary.good_init_extends, hits,
            "{task}: saved initial hits per batch"
        );
        let one_ref = blast.ok_env(work, "blastn", &bargs, &[("BATCH_SIZE", "10000000")]);
        let (one, s1) = ours_with_summary(&spec, &q, &args, Some(10_000_000));
        assert_same(&format!("{task} with one query batch"), &one_ref, &one);
        assert_eq!(s1.batches, vec![10]);
        if task == "blastn" {
            assert_ne!(reference, one, "the batch boundary before q06 must matter");
            let extra: Vec<&str> = one.lines().filter(|l| !reference.contains(l)).collect();
            assert_eq!(extra.len(), 4, "{extra:#?}");
            assert!(
                extra
                    .iter()
                    .all(|l| l.starts_with("q06_") && l.split('\t').nth(6) == Some("1")),
                "{extra:#?}"
            );
            assert!(reference.lines().all(|l| one.contains(l)));
        } else {
            assert_same(
                &format!("{task}: one batch equals the default"),
                &reference,
                &one,
            );
        }
    }
}

#[test]
#[ignore]
fn golden_search_blastn_hg38() {
    let Ok(jungle) = std::env::var("GTARS_RGBLAST_STORE_JUNGLE") else {
        return;
    };
    let blast = Blast::find();
    let dir = tempdir();
    let work = dir.path();
    let sel = gtars_rgblast::select::Selection::explicit(&["hg38-primary-refgenie"]);
    reference_db(&blast, &jungle, &sel, work, "mk/db", false);
    let q = Path::new(env!("CARGO_MANIFEST_DIR")).join("benchmarks/results/q_hg38.fa");
    std::fs::copy(&q, work.join("q.fa")).unwrap();
    let fmt = "6 std stitle";
    let reference = blast.ok(
        work,
        "blastn",
        &[
            "-db", "mk/db", "-query", "q.fa", "-task", "blastn", "-outfmt", fmt,
        ],
    );
    let mine = ours(
        &jungle,
        &work.join("q.fa"),
        fmt,
        &["--collection", "hg38-primary-refgenie", "-task", "blastn"],
    );
    assert_same("hg38 -task blastn", &reference, &mine);
}
