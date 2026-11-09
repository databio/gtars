//! rgblastn vs BLAST+ 2.17.0 (feature `golden`): byte-identical `-outfmt 6` (every field) for
//! megablast on the fixture stores (the footers are in `golden_search_blastn`); subject targets
//! (`--sequences-file`) and the `STATS_NSEQ`/`STATS_TOTLEN` override; `#[ignore]`
//! whole-collection runs on the jungle mirror when `GTARS_RGBLAST_STORE_JUNGLE` is set.

mod common;

use common::blast::Blast;
use common::golden::*;
use common::*;
use gtars_rgblast::blastdb::alias::{NalSpec, nal_text, write_bsl};
use gtars_rgblast::plan::SubjectPlan;
use gtars_rgblast::select::Selection;
use std::path::Path;
use std::process::{Command, Output};

const EVERY_FIELD: &str = "6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids staxid stitle salltitles sstrand qcovs qcovhsp qcovus";

/// Deterministic query classes over the fixture sequences: exact slices, substitution
/// mutants, reverse complements, a slice with Ns, a short one, an asmC slice (a defline with
/// several `gnl|` ids), lower-case and `U` spellings, odd deflines, and a random (no-hit) one.
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
    }
    let random: Vec<u8> = (0..500)
        .map(|_| b"ACGT"[(rng.next() % 4) as usize])
        .collect();
    push(&mut fa, "random", &random);
    // asmC c1: its occurrence defline carries three gnl| ids (asmC, asmC-alt, GCF_TEST.1)
    let c = read_fasta(&Path::new(FIXTURES).join("asmC.fa"));
    let c1 = &c.iter().find(|(n, _)| n == "c1").unwrap().1;
    push(&mut fa, "c1_exact", &c1[100..400]);
    // a gapped mutant of c1: deletions and an insertion
    let mut g = c1[50..450].to_vec();
    g.drain(100..103);
    g.insert(250, b'A');
    g.insert(250, b'C');
    g.drain(300..301);
    push(&mut fa, "c1_indel", &g);
    // the mutant of chr1 reverse complemented, so a minus-strand hit has mismatches
    let chr1 = &seqs[0].1;
    let mut m = chr1[200..600].to_vec();
    for p in (0..m.len()).step_by(37) {
        m[p] = if m[p] == b'A' { b'C' } else { b'A' };
    }
    push(&mut fa, "chr1_mut_rc", &rc_letters(&m));
    // spellings and deflines blastn's reader accepts
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

/// Query classes over the `taxa` fixture: human, mouse, the shared decoy, the taxid-0
/// collection, a 4-bit sequence and the N-heavy one.
fn taxa_queries() -> Vec<u8> {
    let taxa = Path::new(TAXA);
    let get = |file: &str, name: &str| -> Vec<u8> {
        read_fasta(&taxa.join(file))
            .into_iter()
            .find(|(n, _)| n == name)
            .unwrap()
            .1
    };
    let mut fa = Vec::new();
    push_query(
        &mut fa,
        "h_chr1 human_a chr1",
        &get("human_a.fa", "chr1")[20..260],
    );
    push_query(
        &mut fa,
        "h_chr2 human_a chr2",
        &get("human_a.fa", "chr2")[10..230],
    );
    push_query(
        &mut fa,
        "m_chr1 mouse_a chr1",
        &get("mouse_a.fa", "chr1")[30..270],
    );
    push_query(
        &mut fa,
        "m_7 mouse_a 7",
        &rc_letters(&get("mouse_a.fa", "7")[10..200]),
    );
    push_query(
        &mut fa,
        "decoy human_b and mouse_a",
        &get("human_b.fa", "decoy")[15..235],
    );
    push_query(&mut fa, "u1 nofhr", &get("nofhr.fa", "u1")[5..170]);
    push_query(
        &mut fa,
        "iupac mouse_a chrUn_iupac",
        &get("mouse_a.fa", "chrUn_iupac")[0..150],
    );
    push_query(
        &mut fa,
        "masked human_masked chr1",
        &get("human_masked.fa", "chr1")[0..300],
    );
    push_query(&mut fa, "orphan o1", &get("orphan.fa", "o1")[0..120]);
    fa
}

/// The `rgblastn` binary.
fn rgblastn_bin(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rgblastn"))
        .args(args)
        .output()
        .unwrap()
}

fn text(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

/// Digest of the subject `name` of the collection displayed as `coll` in `plan`.
fn digest_of(plan: &SubjectPlan, coll: &str, name: &str) -> String {
    plan.subjects
        .iter()
        .find(|s| {
            s.occurrences
                .iter()
                .any(|o| o.name == name && plan.collections[o.collection as usize].display == coll)
        })
        .unwrap_or_else(|| panic!("{coll}:{name}"))
        .digest
        .clone()
}

/// Megablast flag cells on the fixture store, `-outfmt "6 std stitle"`.
#[test]
fn megablast_fixture_byte_identical() {
    let f = encoded();
    let dir = tempdir();
    let work = dir.path();
    let blast = Blast::find();
    reference_db_all(&blast, &f.spec(), work);
    std::fs::write(work.join("q.fa"), fixture_queries()).unwrap();
    for (name, args) in [
        ("default", vec![]),
        ("evalue", vec!["-evalue", "1e-5"]),
        ("plus", vec!["-strand", "plus"]),
        ("minus", vec!["-strand", "minus"]),
        ("nodust", vec!["-dust", "no"]),
        ("hard", vec!["-soft_masking", "false"]),
        ("mts1", vec!["-max_target_seqs", "1"]),
        ("w16", vec!["-word_size", "16"]),
        ("pid", vec!["-perc_identity", "98"]),
    ] {
        let mut a = vec!["-task", "megablast"];
        a.extend(args.iter());
        let reference = check(&blast, &f.spec(), work, "q.fa", "6 std stitle", name, &a);
        assert!(
            !reference.is_empty() || name == "pid",
            "{name}: the reference produced no hits; the test is vacuous"
        );
    }
}

/// Every field, on the fixture and `taxa` stores (taxids through `-taxid_map`).
#[test]
fn megablast_fixture_all_fields() {
    let blast = Blast::find();
    for (store, queries, label) in [
        (encoded().spec(), fixture_queries(), "fixture"),
        (taxa_encoded().spec(), taxa_queries(), "taxa"),
    ] {
        let dir = tempdir();
        let work = dir.path();
        reference_db(&blast, &store, &Selection::all(), work, "mk/db", true);
        std::fs::write(work.join("q.fa"), &queries).unwrap();
        for (fmt, extra) in [
            (ALL_FIELDS, vec![]),
            (ALL_FIELDS, vec!["-strand", "minus"]),
            (EVERY_FIELD, vec![]),
            (EVERY_FIELD, vec!["-evalue", "1e-3"]),
            ("6 std bogus", vec![]),
        ] {
            let name = format!("{label} -outfmt '{fmt}' {extra:?}");
            let reference = check(&blast, &store, work, "q.fa", fmt, &name, &extra);
            assert!(!reference.is_empty(), "{name}: no hit lines; vacuous");
        }
    }
}

/// `sallseqid` prints one id per defline (the best, which is the first), not every id:
/// defline 0 of asmA chr1 holds `SQ.` plus two sequence aliases; asmC's occurrence defline
/// holds three `gnl|` ids.
#[test]
fn sallseqid_one_id_per_defline() {
    let f = encoded();
    let blast = Blast::find();
    let dir = tempdir();
    let work = dir.path();
    let plan = reference_db(&blast, &f.spec(), &Selection::all(), work, "mk/db", true);
    std::fs::write(work.join("q.fa"), fixture_queries()).unwrap();
    let fmt = "6 qseqid sseqid sallseqid sallacc salltitles staxids";
    let reference = blast.ok(
        work,
        "blastn",
        &["-db", "mk/db", "-query", "q.fa", "-outfmt", fmt],
    );
    let mine = ours(&f.spec(), &work.join("q.fa"), fmt, &[]);
    assert_same("sallseqid", &reference, &mine);
    let chr1 = digest_of(&plan, "asmA", "chr1");
    let c1 = digest_of(&plan, "asmC", "c1");
    let s = plan.subjects.iter().find(|s| s.digest == chr1).unwrap();
    assert_eq!(s.deflines[0].ids.len(), 3, "SQ. + CM000001.1 + NC_TEST.1");
    let sc = plan.subjects.iter().find(|s| s.digest == c1).unwrap();
    assert!(sc.deflines[1].ids.len() >= 3, "{:?}", sc.deflines[1].ids);
    let line = mine
        .lines()
        .find(|l| l.starts_with("q01_chr1_exact\t"))
        .unwrap();
    let f: Vec<&str> = line.split('\t').collect();
    assert_eq!(f[2], format!("SQ.{chr1};gnl|asmA|chr1;gnl|asmB|chr1"));
    assert_eq!(f[3], format!("SQ.{chr1};asmA:chr1;asmB:chr1"));
    let mut want: Vec<i32> = s.deflines.iter().filter_map(|d| d.taxid).collect();
    want.sort_unstable();
    want.dedup();
    assert_eq!(
        f[5],
        want.iter()
            .map(|t| t.to_string())
            .collect::<Vec<_>>()
            .join(";")
    );
    let line = mine
        .lines()
        .find(|l| l.contains(&format!("\tSQ.{c1}\t")))
        .unwrap();
    let f: Vec<&str> = line.split('\t').collect();
    assert_eq!(f[2], format!("SQ.{c1};gnl|asmC|c1"));
}

/// jungle chrM: one unique sequence held by 37 collections, searched as a one-sequence
/// `--sequences-file` target against the same database built by makeblastdb.
#[test]
#[ignore]
fn sallseqid_jungle_chrm() {
    let Ok(jungle) = std::env::var("GTARS_RGBLAST_STORE_JUNGLE") else {
        return;
    };
    let blast = Blast::find();
    let dir = tempdir();
    let work = dir.path();
    let chrm = "k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct";
    let sel = Selection::all().with_sequences(&[chrm]);
    let plan = reference_db(&blast, &jungle, &sel, work, "mk/db", true);
    assert_eq!(plan.subjects.len(), 1);
    assert!(
        plan.subjects[0].deflines.len() > 30,
        "{} deflines",
        plan.subjects[0].deflines.len()
    );
    let mut fa = Vec::new();
    push_query(&mut fa, "chrM_1001_2000", &slice(&jungle, chrm, 1001, 2000));
    push_query(
        &mut fa,
        "chrM_rc",
        &rc_letters(&slice(&jungle, chrm, 5001, 5800)),
    );
    std::fs::write(work.join("q.fa"), fa).unwrap();
    let ids = work.join("chrm.txt");
    std::fs::write(&ids, format!("SQ.{chrm}\n")).unwrap();
    let fmt = "6 sseqid sallseqid sallacc salltitles staxids";
    let reference = blast.ok(
        work,
        "blastn",
        &["-db", "mk/db", "-query", "q.fa", "-outfmt", fmt],
    );
    let mine = ours(
        &jungle,
        &work.join("q.fa"),
        fmt,
        &["--sequences-file", ids.to_str().unwrap()],
    );
    assert_same("jungle chrM", &reference, &mine);
    assert!(
        mine.lines()
            .next()
            .unwrap()
            .split('\t')
            .nth(1)
            .unwrap()
            .split(';')
            .count()
            > 30
    );
}

/// A wrapper alias `.nal` at `<work>/views/<name>.nal` over `db` (relative to `views/`).
fn wrapper(
    work: &Path,
    name: &str,
    db: &str,
    bsl: Option<&str>,
    stats: Option<(u64, u64)>,
) -> String {
    std::fs::create_dir_all(work.join("views")).unwrap();
    let text = nal_text(&NalSpec {
        title: name,
        dblist: &[db.to_string()],
        seqidlist: bsl,
        nseq: None,
        length: None,
        full_stats: stats,
        date: DATE,
    });
    std::fs::write(work.join("views").join(format!("{name}.nal")), text).unwrap();
    format!("views/{name}")
}

/// `STATS_NSEQ`/`STATS_TOTLEN` (the staged search of sub-plan 5): rgblastn on a
/// `--sequences-file` subset with the full database's statistics equals blastn under a
/// wrapper alias, for both reference shapes (the subset as its own database; the full
/// database with a `SEQIDLIST` view), in process and through the binary.
#[test]
fn stats_override_equals_nal_stats() {
    let f = encoded();
    let blast = Blast::find();
    let dir = tempdir();
    let work = dir.path();
    let full = reference_db(&blast, &f.spec(), &Selection::all(), work, "mk/db", false);
    let subset: Vec<String> = vec![
        digest_of(&full, "asmA", "chr1"),
        digest_of(&full, "asmB", "chrB2"),
        digest_of(&full, "asmC", "c1"),
    ];
    let sel = Selection::all().with_sequences(&subset);
    reference_db(&blast, &f.spec(), &sel, work, "mk/sub", false);
    std::fs::write(work.join("q.fa"), fixture_queries()).unwrap();
    let seqs = work.join("subset.txt");
    std::fs::write(
        &seqs,
        subset
            .iter()
            .map(|d| format!("SQ.{d}\n"))
            .collect::<String>(),
    )
    .unwrap();
    let (nseq, totlen) = (full.totals.num_seqs as u64, full.totals.total_length);
    let (ns, tl) = (nseq.to_string(), totlen.to_string());
    let seqs_arg = seqs.to_str().unwrap();
    // shape A: the subset database under STATS_*
    let a = wrapper(work, "a", "../mk/sub", None, Some((nseq, totlen)));
    // shape B: the full database with a SEQIDLIST view under STATS_* (only SQ. deflines show)
    let ids: Vec<String> = subset.iter().map(|d| format!("SQ.{d}")).collect();
    write_bsl(&ids, &work.join("views/b.bsl")).unwrap();
    let b = wrapper(work, "b", "../mk/db", Some("b.bsl"), Some((nseq, totlen)));
    let fmt_a = "6 std stitle sallseqid qlen slen";
    let fmt_b = "6 std qlen slen";
    let qpath = work.join("q.fa");
    for task_args in [vec![], vec!["-evalue", "1e-5"]] {
        for (shape, db, fmt) in [("A", &a, fmt_a), ("B", &b, fmt_b)] {
            let mut args = vec!["-db", db, "-query", "q.fa", "-task", "megablast"];
            args.extend(["-outfmt", fmt]);
            args.extend(task_args.iter());
            let reference = blast.ok(work, "blastn", &args);
            assert!(
                reference.lines().count() >= 3,
                "shape {shape}: too few hits ({})",
                reference.lines().count()
            );
            let mut ours_args = vec![
                "--sequences-file",
                seqs_arg,
                "--stats-nseq",
                &ns,
                "--stats-totlen",
                &tl,
            ];
            ours_args.extend(task_args.iter());
            let mine = ours(&f.spec(), &qpath, fmt, &ours_args);
            assert_same(
                &format!("stats shape {shape} {task_args:?}"),
                &reference,
                &mine,
            );
            let spec = f.spec();
            let mut bin = vec!["--store", spec.as_str(), "--all", "-query"];
            bin.extend([qpath.to_str().unwrap(), "-outfmt", fmt]);
            bin.extend(ours_args.iter());
            let o = rgblastn_bin(&bin);
            assert!(o.status.success(), "{}", text(&o.stderr));
            assert_same(
                &format!("stats shape {shape} binary"),
                &reference,
                &text(&o.stdout),
            );
        }
    }
    // the subset's own statistics differ from the full database's (the override matters)
    let plain = ours(&f.spec(), &qpath, fmt_a, &["--sequences-file", seqs_arg]);
    let with = ours(
        &f.spec(),
        &qpath,
        fmt_a,
        &[
            "--sequences-file",
            seqs_arg,
            "--stats-nseq",
            &ns,
            "--stats-totlen",
            &tl,
        ],
    );
    assert_ne!(plain, with, "the statistics override changed nothing");
    // --sequences-file alone equals a database of just those subjects
    let reference = blast.ok(
        work,
        "blastn",
        &[
            "-db",
            "mk/sub",
            "-query",
            "q.fa",
            "-task",
            "megablast",
            "-outfmt",
            fmt_a,
        ],
    );
    assert_same("--sequences-file", &reference, &plain);
}

/// `rgblast search` (the default engine, rgblastn) and its runner (`-out`, refused flags).
#[test]
fn search_engine_rgblastn_matches_blastn() {
    let f = encoded();
    let blast = Blast::find();
    let dir = tempdir();
    let work = dir.path();
    reference_db(&blast, &f.spec(), &Selection::all(), work, "mk/db", true);
    std::fs::write(work.join("q.fa"), fixture_queries()).unwrap();
    let q = work.join("q.fa");
    let rgblast = |args: &[&str]| -> Output {
        Command::new(env!("CARGO_BIN_EXE_rgblast"))
            .args(["search", "--store", &f.spec(), "--all"])
            .args(args)
            .output()
            .unwrap()
    };
    let fmt = "6 std qlen score sallseqid salltitles staxids";
    let reference = blast.ok(
        work,
        "blastn",
        &[
            "-db", "mk/db", "-query", "q.fa", "-outfmt", fmt, "-evalue", "1e-3",
        ],
    );
    let o = rgblast(&[
        "--quiet",
        "-query",
        q.to_str().unwrap(),
        "-outfmt",
        "6",
        "-evalue",
        "1e-3",
    ]);
    assert!(o.status.success(), "{}", text(&o.stderr));
    assert_same("rgblast search", &reference, &text(&o.stdout));
    // -out through the runner
    let out = work.join("hits.tsv");
    let o = rgblast(&[
        "--no-extra-fields",
        "-query",
        q.to_str().unwrap(),
        "-outfmt",
        "6 std stitle",
        "-out",
        out.to_str().unwrap(),
    ]);
    assert!(o.status.success(), "{}", text(&o.stderr));
    let reference = blast.ok(
        work,
        "blastn",
        &["-db", "mk/db", "-query", "q.fa", "-outfmt", "6 std stitle"],
    );
    assert_same(
        "runner -out",
        &reference,
        &std::fs::read_to_string(&out).unwrap(),
    );
    // a flag rgblastn does not implement is refused (exit 2) with its wording
    let o = rgblast(&["-query", "q.fa", "-sum_stats", "true"]);
    assert_eq!(o.status.code(), Some(2), "{}", text(&o.stderr));
    assert!(
        text(&o.stderr).contains("-sum_stats is not supported"),
        "{}",
        text(&o.stderr)
    );
}

#[test]
#[ignore]
fn golden_search_hg38() {
    let Ok(jungle) = std::env::var("GTARS_RGBLAST_STORE_JUNGLE") else {
        return;
    };
    let blast = Blast::find();
    let dir = tempdir();
    let work = dir.path();
    let sel = Selection::explicit(&["hg38-primary-refgenie"]);
    reference_db(&blast, &jungle, &sel, work, "mk/db", false);
    let q = Path::new(env!("CARGO_MANIFEST_DIR")).join("benchmarks/results/q_hg38.fa");
    std::fs::copy(&q, work.join("q.fa")).unwrap();
    let fmt = "6 std stitle";
    let reference = blast.ok(
        work,
        "blastn",
        &[
            "-db",
            "mk/db",
            "-query",
            "q.fa",
            "-num_threads",
            "4",
            "-outfmt",
            fmt,
        ],
    );
    let mine = ours(
        &jungle,
        &work.join("q.fa"),
        fmt,
        &["--collection", "hg38-primary-refgenie", "-num_threads", "4"],
    );
    assert_same("hg38", &reference, &mine);
}
