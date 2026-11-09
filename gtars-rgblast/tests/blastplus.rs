//! The functional contract through real BLAST+ tools (feature `blastplus`): blastdbcmd sees
//! every sequence and id, and blastn gives the same tabular output on our database as on a
//! makeblastdb database built from the reference FASTA.

mod common;

use common::blast::Blast;
use common::*;
use gtars_rgblast::blastdb::output::with_ext;
use gtars_rgblast::blastdb::{njs, sidecars};
use gtars_rgblast::plan::SubjectPlan;
use gtars_rgblast::reference_fasta;
use gtars_rgblast::select::Selection;
use gtars_rgblast::store::decode::to_ascii;
use gtars_rgblast::store::{SeqFetch, StoreFiles};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::OnceLock;

pub struct Setup {
    _dir: tempfile::TempDir,
    pub work: PathBuf,
    pub plan: SubjectPlan,
    pub blast: Blast,
    /// OID -> sequence with U written as T (what BLAST stores).
    pub seqs: Vec<Vec<u8>>,
}

/// Our database at `ours/db`, the reference FASTA at `ref.fa` and makeblastdb's database
/// from it at `mk/db`, built once per test binary.
pub fn setup() -> &'static Setup {
    static S: OnceLock<Setup> = OnceLock::new();
    S.get_or_init(|| {
        let f = encoded();
        let blast = Blast::find();
        let dir = tempdir();
        let work = dir.path().to_path_buf();
        build(&f.spec(), Selection::all(), &work.join("ours/db"));
        let (st, plan) = plan_for(&f.spec(), &Selection::all());
        let fetch = StoreFiles::new(&st);
        reference_fasta::write(&plan, &fetch, &work.join("ref.fa")).unwrap();
        std::fs::create_dir_all(work.join("mk")).unwrap();
        blast.ok(
            &work,
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
                "-out",
                "mk/db",
            ],
        );
        let seqs = plan
            .subjects
            .iter()
            .map(|s| {
                let stored = fetch.read(&s.digest).unwrap();
                let stored = &stored[..];
                let a = to_ascii(plan.store_mode, s.alphabet, stored, s.length).unwrap();
                a.into_iter()
                    .map(|b| if b == b'U' { b'T' } else { b })
                    .collect()
            })
            .collect();
        Setup {
            _dir: dir,
            work,
            plan,
            blast,
            seqs,
        }
    })
}

#[test]
fn blastdbcmd_md5_every_entry() {
    let s = setup();
    let out = s.blast.ok(
        &s.work,
        "blastdbcmd",
        &["-db", "ours/db", "-entry", "all", "-outfmt", "%o %s"],
    );
    let mut seen: HashMap<u32, String> = HashMap::new();
    for line in out.lines() {
        let (oid, seq) = line.split_once(' ').unwrap();
        seen.entry(oid.parse().unwrap())
            .or_insert_with(|| md5_hex(seq.to_ascii_uppercase().as_bytes()));
    }
    assert_eq!(seen.len(), s.plan.subjects.len());
    for (oid, sub) in s.plan.subjects.iter().enumerate() {
        let got = &seen[&(oid as u32)];
        // the store md5, or (for sequences with U) the md5 with U -> T
        assert!(
            *got == sub.md5 || *got == md5_hex(&s.seqs[oid]),
            "OID {oid} SQ.{}",
            sub.digest
        );
    }
}

#[test]
fn blastdbcmd_every_id_resolves() {
    let s = setup();
    let ids = sidecars::read_ids(&with_ext(&s.work.join("ours/db"), "ids.tsv")).unwrap();
    let batch: String = ids.iter().map(|(id, _, _)| format!("{id}\n")).collect();
    std::fs::write(s.work.join("ids.batch"), batch).unwrap();
    let out = s.blast.ok(
        &s.work,
        "blastdbcmd",
        &[
            "-db",
            "ours/db",
            "-entry_batch",
            "ids.batch",
            "-outfmt",
            "%o",
        ],
    );
    let got: Vec<u32> = out.lines().map(|l| l.trim().parse().unwrap()).collect();
    let want: Vec<u32> = ids.iter().map(|(_, oid, _)| *oid).collect();
    assert_eq!(got, want);
}

#[test]
fn blastdbcmd_dash_digest_entry() {
    let s = setup();
    let dd = encoded().digest("asmD");
    let (oid, sub) = s
        .plan
        .subjects
        .iter()
        .enumerate()
        .find(|(_, x)| {
            x.occurrences
                .iter()
                .any(|o| s.plan.collections[o.collection as usize].digest == dd)
        })
        .unwrap();
    let entry = format!("gnl|{dd}|{}", sub.occurrences[0].id_name);
    let out = s.blast.ok(
        &s.work,
        "blastdbcmd",
        &["-db", "ours/db", "-entry", &entry, "-outfmt", "%o"],
    );
    assert_eq!(out.trim(), oid.to_string());
}

#[test]
fn blastdbcmd_metadata_equals_njs() {
    let s = setup();
    let out = s
        .blast
        .ok(&s.work, "blastdbcmd", &["-db", "ours/db", "-metadata"]);
    let mut got: njs::Njs = serde_json::from_str(&out).unwrap();
    let ours = njs::read(&s.work.join("ours/db.njs")).unwrap();
    // blastdbcmd recomputes the metadata with .njs itself present, so it lists it too
    let njs_size = std::fs::metadata(s.work.join("ours/db.njs")).unwrap().len();
    assert!(got.files.contains(&"db.njs".to_string()));
    got.files.retain(|f| f != "db.njs");
    got.bytes_total -= njs_size;
    assert_eq!(got, ours);
}

#[test]
fn blastdbcmd_info_totals() {
    let s = setup();
    let out = s
        .blast
        .ok(&s.work, "blastdbcmd", &["-db", "ours/db", "-info"]);
    let t = &s.plan.totals;
    let want = format!(
        "{} sequences; {} total bases",
        gtars_rgblast::progress::commas(t.num_seqs as u64),
        gtars_rgblast::progress::commas(t.total_length)
    );
    assert!(out.contains(&want), "{out}\nexpected: {want}");
    assert!(out.contains(&format!("Database: {TITLE}")), "{out}");
}

/// Deterministic queries (ported from the prototype's gen_queries.py): exact slices, slices
/// with 2-8 substitutions, reverse complements, slices overlapping N runs, and a random query.
fn queries(s: &Setup) -> String {
    let mut x: u64 = 0x5eed_1234_abcd_0001;
    let mut next = move |n: u64| {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x % n
    };
    let rc = |q: &[u8]| -> Vec<u8> {
        q.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                o => o,
            })
            .collect()
    };
    let mut out = String::new();
    let big: Vec<usize> = (0..s.seqs.len())
        .filter(|&i| s.seqs[i].len() >= 120)
        .collect();
    for k in 0..24 {
        let i = big[next(big.len() as u64) as usize];
        let sq = &s.seqs[i];
        let len = (60 + next(240) as usize).min(sq.len());
        let a = next((sq.len() - len + 1) as u64) as usize;
        let mut q = sq[a..a + len].to_vec();
        let mut tag = String::new();
        if k % 3 == 1 {
            let n = 2 + next(7) as usize;
            for _ in 0..n {
                let p = next(q.len() as u64) as usize;
                q[p] = match q[p] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }
            tag += &format!("_mut{n}");
        }
        if k % 4 == 2 {
            q = rc(&q);
            tag += "_rc";
        }
        out += &format!(
            ">q{k:02}_oid{i}_{a}_{len}{tag}\n{}\n",
            String::from_utf8_lossy(&q)
        );
    }
    // slices overlapping the N runs of chr2 (the sequence with the most N)
    let (i, sq) = s
        .seqs
        .iter()
        .enumerate()
        .max_by_key(|(_, q)| q.iter().filter(|&&b| b == b'N').count())
        .unwrap();
    for (k, start) in [40usize, 105, 160, 4300].iter().enumerate() {
        let q = &sq[*start..(*start + 150).min(sq.len())];
        out += &format!(">n{k}_oid{i}_{start}\n{}\n", String::from_utf8_lossy(q));
    }
    let random: String = (0..200)
        .map(|_| b"ACGT"[next(4) as usize] as char)
        .collect();
    out += &format!(">nohit\n{random}\n");
    out
}

const FMT: &str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sallseqid salltitles staxids";

#[test]
fn blastn_tabular_equals_makeblastdb() {
    let s = setup();
    std::fs::write(s.work.join("q.fa"), queries(s)).unwrap();
    for task in ["megablast", "blastn"] {
        let mut outs = Vec::new();
        for db in ["ours/db", "mk/db"] {
            outs.push(s.blast.ok(
                &s.work,
                "blastn",
                &["-task", task, "-db", db, "-query", "q.fa", "-outfmt", FMT],
            ));
        }
        assert!(
            outs[0].lines().count() > 20,
            "{task}: too few hits:\n{}",
            outs[0]
        );
        assert_eq!(
            outs[0], outs[1],
            "{task}: blastn output differs between our database and makeblastdb's"
        );
    }
}

#[test]
fn sallseqid_shows_first_ids() {
    let s = setup();
    // one query per subject of at least 40 bp, from its first (up to) 100 bases
    let mut q = String::new();
    for (oid, sq) in s.seqs.iter().enumerate() {
        let head = &sq[..sq.len().min(100)];
        if head.len() >= 40 && !head.contains(&b'N') {
            q += &format!(">s{oid}\n{}\n", String::from_utf8_lossy(head));
        }
    }
    std::fs::write(s.work.join("first.fa"), q).unwrap();
    let out = s.blast.ok(
        &s.work,
        "blastn",
        &[
            "-task",
            "megablast",
            "-db",
            "ours/db",
            "-query",
            "first.fa",
            "-outfmt",
            "6 qseqid sseqid sallseqid pident length",
        ],
    );
    let mut checked = 0;
    for line in out.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let oid: usize = f[0][1..].parse().unwrap();
        let sub = &s.plan.subjects[oid];
        if f[1] != format!("SQ.{}", sub.digest) || f[3] != "100.000" {
            continue; // only self hits
        }
        let mut want = vec![format!("SQ.{}", sub.digest)];
        for o in &sub.occurrences {
            want.push(format!(
                "gnl|{}|{}",
                s.plan.collections[o.collection as usize].db_tag, o.id_name
            ));
        }
        assert_eq!(f[2], want.join(";"), "OID {oid}");
        checked += 1;
    }
    assert!(checked >= 20, "only {checked} self hits checked:\n{out}");
    assert!(
        out.contains(";gnl|numeric|22\t"),
        "integer tags print like string ones:\n{out}"
    );
    assert!(
        out.contains(";gnl|asmA|chrX_alt"),
        "the sanitized ';' name keeps the list intact:\n{out}"
    );
}
