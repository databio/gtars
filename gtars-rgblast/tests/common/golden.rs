//! Helpers for the rgblastn tests: query building, rgblastn through its real argument parser,
//! and (feature `reference-fasta`) makeblastdb reference databases and blastn's statistics
//! footer.

use super::store::{Layout, build_store_from_fasta};
use clap::Parser;
use gtars_refget::store::StorageMode;
use gtars_rgblast::api::rgblastn::{QuerySource, RgblastnOutcome, RgblastnRequest};
use gtars_rgblast::cli::RgblastnArgs;
use gtars_rgblast::cli::rgblastn::normalize_argv;
use gtars_rgblast::engine::encoding::{ASCII_TO_BLASTNA, reverse_complement};
use std::path::{Path, PathBuf};

pub const ALL_FIELDS: &str =
    "6 std stitle sallseqid salltitles staxid staxids sallacc qlen slen sstrand btop qseq sseq";

/// (name, sequence) records; the name is the defline's first word.
pub fn read_fasta(p: &Path) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(p).unwrap();
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            out.push((
                h.split_whitespace().next().unwrap_or("").to_string(),
                Vec::new(),
            ));
        } else if let Some(last) = out.last_mut() {
            last.1.extend(line.trim().bytes());
        }
    }
    out
}

/// A seeded linear congruential generator, so query sets are the same on every run.
pub struct Lcg(pub u64);

impl Lcg {
    pub fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

/// Append one FASTA record, 70 letters per line.
pub fn push_query(fa: &mut Vec<u8>, defline: &str, s: &[u8]) {
    fa.extend_from_slice(format!(">{defline}\n").as_bytes());
    for chunk in s.chunks(70) {
        fa.extend_from_slice(chunk);
        fa.push(b'\n');
    }
}

/// The reverse complement of IUPAC letters.
pub fn rc_letters(s: &[u8]) -> Vec<u8> {
    let letters: Vec<u8> = s.iter().map(|&c| ASCII_TO_BLASTNA[c as usize]).collect();
    reverse_complement(&letters)
        .iter()
        .map(|&c| b"ACGTRYMKWSBDHVN-"[c as usize])
        .collect()
}

/// An rgblastn request from blastn-style flags (`-task blastn`, `-dust no`, ...), through the
/// command line's own parser. `--all` is added unless `args` select collections.
pub fn rgblastn_req(store: &str, query: &Path, outfmt: &str, args: &[&str]) -> RgblastnRequest {
    let mut argv: Vec<String> = ["rgblastn", "--store", store, "-query"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    argv.push(query.display().to_string());
    argv.extend(["-outfmt".to_string(), outfmt.to_string()]);
    let selects = args.iter().any(|a| {
        matches!(
            *a,
            "--all" | "--collection" | "--namespace" | "--taxon" | "--collections-file"
        )
    });
    if !selects {
        argv.push("--all".into());
    }
    argv.extend(args.iter().map(|s| s.to_string()));
    let a = RgblastnArgs::try_parse_from(normalize_argv(argv)).unwrap();
    RgblastnRequest::try_from(&a).unwrap()
}

/// `rgblastn_req` with a query from memory.
pub fn rgblastn_req_bytes(
    store: &str,
    fasta: &[u8],
    outfmt: &str,
    args: &[&str],
) -> RgblastnRequest {
    let mut req = rgblastn_req(store, Path::new("unused.fa"), outfmt, args);
    req.query = QuerySource::Bytes(fasta.to_vec());
    req
}

pub fn run(req: &RgblastnRequest) -> RgblastnOutcome {
    gtars_rgblast::api::rgblastn::rgblastn(req, &mut |_| true).unwrap()
}

/// The tabular output, as the binary writes it.
pub fn render(o: &RgblastnOutcome) -> String {
    let mut buf = Vec::new();
    gtars_rgblast::api::rgblastn::write_tabular(o, &mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

/// rgblastn in process: `rgblastn_req`, run, render.
pub fn ours(store: &str, query: &Path, outfmt: &str, args: &[&str]) -> String {
    render(&run(&rgblastn_req(store, query, outfmt, args)))
}

/// The upper-case letters `start..=end` (1-based) of the stored sequence `digest`.
pub fn slice(store: &str, digest: &str, start: usize, end: usize) -> Vec<u8> {
    use gtars_rgblast::store::{SeqFetch, StoreFiles, decode};
    let mut st = gtars_rgblast::OpenedStore::open(store, None).unwrap();
    st.store.load_sequence_index().unwrap();
    let m = st.store.get_sequence_metadata(digest).unwrap().clone();
    let stored = StoreFiles::new(&st).read(digest).unwrap();
    let ascii = decode::to_ascii(st.mode(), m.alphabet, &stored, m.length as u32).unwrap();
    ascii[start - 1..end].to_vec()
}

/// Queries over sequences of a FASTA set: exact slices, substitution mutants and reverse
/// complements of up to `take` sequences, plus one random no-hit query.
pub fn queries_from(fastas: &[PathBuf], take: usize, seed: u64) -> Vec<u8> {
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

/// A store with subjects longer than `MAX_DBSEQ_LEN` (5 Mbp), so the search decodes them chunk
/// by chunk: a 12.4 Mbp sequence (3 chunks) with planted copies of the query material across
/// the chunk boundaries and N runs (so it is 3-bit, its second chunk starts mid-group and the
/// random fill carries over from chunk to chunk), a 5.2 Mbp one (2 chunks) and a short one; the
/// same sequences also as raw and zstd stores beside it (`raw/`, `zstd/`), all in `layout`
/// (default pack cap: the 2-bit `midB` record then sits at the unaligned offset 4,650,000,
/// behind `longA`'s 3-bit one).
pub fn long_subject_store(dir: &Path, layout: Layout) -> (String, PathBuf) {
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
    build_store_from_fasta(&store, &fasta, "long", StorageMode::Encoded, layout);
    for (name, mode) in [("raw", StorageMode::Raw), ("zstd", StorageMode::Zstd)] {
        build_store_from_fasta(&dir.join(name), &fasta, "long", mode, layout);
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

pub fn assert_same(name: &str, reference: &str, mine: &str) {
    assert!(
        reference == mine,
        "{name}: output differs\n--- blastn\n{reference}\n--- rgblastn\n{mine}"
    );
}

#[cfg(feature = "reference-fasta")]
#[allow(unused_imports)]
pub use reference::*;

#[cfg(feature = "reference-fasta")]
mod reference {
    use super::*;
    use crate::common::blast::Blast;
    use crate::common::store::{TITLE, build, plan_for};
    use gtars_rgblast::SubjectPlan;
    use gtars_rgblast::reference_fasta;
    use gtars_rgblast::select::Selection;
    use gtars_rgblast::store::StoreFiles;

    /// The makeblastdb database of `sel` at `<work>/<name>`, from `reference_fasta::write`,
    /// and the plan it was built from. With `taxonomy`, the plan's taxids go in through
    /// `-taxid_map` when it can express them; otherwise the `rgmakeblastdb` database itself
    /// is the reference (makeblastdb would give every defline a taxid).
    pub fn reference_db(
        blast: &Blast,
        store: &str,
        sel: &Selection,
        work: &Path,
        name: &str,
        taxonomy: bool,
    ) -> SubjectPlan {
        let (st, plan) = plan_for(store, sel);
        let fa = format!("{name}.fa");
        std::fs::create_dir_all(work.join(name).parent().unwrap()).unwrap();
        let map = format!("{name}.taxid_map");
        let mut args = vec![
            "-in",
            &fa,
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-blastdb_version",
            "5",
            "-title",
            TITLE,
            "-out",
            name,
        ];
        if taxonomy {
            if reference_fasta::write_taxid_map(&plan, &work.join(&map)).is_ok() {
                args.extend(["-taxid_map", &map]);
            } else {
                build(store, sel.clone(), &work.join(name));
                return plan;
            }
        }
        let fetch = StoreFiles::new(&st);
        reference_fasta::write(&plan, &fetch, &work.join(&fa)).unwrap();
        blast.ok(work, "makeblastdb", &args);
        plan
    }

    /// The whole store's reference database at `<work>/mk/db`.
    pub fn reference_db_all(blast: &Blast, store: &str, work: &Path) {
        reference_db(blast, store, &Selection::all(), work, "mk/db", false);
    }

    /// blastn on `<work>/mk/db` and rgblastn on `store`, `query` (a file in `work`) with
    /// `args`, byte for byte; returns blastn's output.
    pub fn check(
        blast: &Blast,
        store: &str,
        work: &Path,
        query: &str,
        outfmt: &str,
        name: &str,
        args: &[&str],
    ) -> String {
        let mut a = vec!["-db", "mk/db", "-query", query, "-outfmt", outfmt];
        a.extend(args.iter());
        let reference = blast.ok(work, "blastn", &a);
        let mine = ours(store, &work.join(query), outfmt, args);
        assert_same(name, &reference, &mine);
        reference
    }

    /// One query's statistics in blastn's pairwise footer: ungapped and gapped Lambda, K, H,
    /// and the effective search space.
    type Footer = (f64, f64, f64, f64, f64, f64, i64);

    fn parse_footer(text: &str) -> Vec<Footer> {
        let mut out = Vec::new();
        let lines: Vec<&str> = text.lines().collect();
        let nums =
            |l: &str| -> Vec<f64> { l.split_whitespace().map(|x| x.parse().unwrap()).collect() };
        let mut i = 0;
        while i < lines.len() {
            if lines[i].starts_with("Lambda") {
                let a = nums(lines[i + 1]);
                let mut j = i + 2;
                while !lines[j].starts_with("Lambda") {
                    j += 1;
                }
                let b = nums(lines[j + 1]);
                let mut k = j + 2;
                while !lines[k].starts_with("Effective search space used:") {
                    k += 1;
                }
                let ss: i64 = lines[k].rsplit(' ').next().unwrap().parse().unwrap();
                out.push((a[0], a[1], a[2], b[0], b[1], b[2], ss));
                i = k;
            }
            i += 1;
        }
        out
    }

    /// Three significant digits, as the footer prints them.
    fn sig3(x: f64) -> f64 {
        format!("{:.2e}", x).parse().unwrap()
    }

    /// blastn's pairwise footer (`-outfmt 0`) against rgblastn's per-context statistics for
    /// the same flags: Lambda, K and H to the footer's three digits, the search space exactly.
    /// blastn prints no block for a query whose statistics could not be computed (all `N`),
    /// so only valid contexts are compared, the plus strand of each query.
    pub fn check_footer(blast: &Blast, store: &str, work: &Path, query: &str, args: &[&str]) {
        let mut a = vec!["-db", "mk/db", "-query", query, "-outfmt", "0"];
        a.extend(args.iter());
        let footer = parse_footer(&blast.ok(work, "blastn", &a));
        let out = run(&rgblastn_req(store, &work.join(query), "6", args));
        let stats = &out.results.stats;
        let valid: Vec<_> = stats.iter().filter(|c| c.is_valid).collect();
        assert_eq!(valid.len(), footer.len() * 2, "{args:?}: valid contexts");
        for (qi, f) in footer.iter().enumerate() {
            let c = valid[qi * 2];
            let (u, g) = (c.kbp_ungapped.unwrap(), c.kbp_gapped.unwrap());
            let same = |want: f64, got: f64| (sig3(want) - sig3(got)).abs() < 1e-9;
            assert!(
                same(f.0, u.lambda) && same(f.1, u.k) && same(f.2, u.h),
                "{args:?} query {qi} ungapped: footer {f:?} vs {u:?}"
            );
            assert!(
                same(f.3, g.lambda) && same(f.4, g.k) && same(f.5, g.h),
                "{args:?} query {qi} gapped: footer {f:?} vs {g:?}"
            );
            assert_eq!(f.6, c.eff_searchsp, "{args:?} query {qi} search space");
        }
    }
}
