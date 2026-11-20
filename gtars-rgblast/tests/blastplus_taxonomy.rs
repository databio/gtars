//! Taxonomy, titles and views through real BLAST+ (feature `blastplus`), on the `taxa` fixture:
//! `-taxids` with and without `-no_taxid_expansion` (the trap), `-negative_taxids`, `staxids`,
//! `-outfmt 0` and view statistics.

mod common;

use common::blast::Blast;
use common::*;
use gtars_rgblast::api::ViewRequest;
use gtars_rgblast::blastdb::output::with_ext;
use gtars_rgblast::blastdb::sidecars::{self, ViewKind};
use gtars_rgblast::select::Selection;
use gtars_rgblast::store::decode::to_ascii;
use gtars_rgblast::store::{SeqFetch, StoreFiles};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

struct Setup {
    _dir: tempfile::TempDir,
    work: PathBuf,
    blast: Blast,
    /// Every lookup id of the full database -> sequence digest.
    digest_of: HashMap<String, String>,
}

const TAB: &str = "6 qseqid sseqid evalue bitscore";

/// Queries from the fixture sequences: two exact slices and one with substitutions per
/// sequence of 100 bp or more (none from the masked chr1).
fn queries(store: &str) -> String {
    let (st, plan) = plan_for(store, &Selection::all());
    let fetch = StoreFiles::new(&st);
    let mut out = String::new();
    for (oid, s) in plan.subjects.iter().enumerate() {
        let stored = fetch.read(&s.digest).unwrap();
        let stored = &stored[..];
        let a = to_ascii(plan.store_mode, s.alphabet, stored, s.length).unwrap();
        if a.len() < 100 || a.iter().filter(|&&b| b == b'N').count() > 10 {
            continue;
        }
        let name = s
            .occurrences
            .first()
            .map(|o| o.name.clone())
            .unwrap_or_else(|| "orphan".into());
        let slice =
            |from: usize, len: usize| String::from_utf8_lossy(&a[from..from + len]).into_owned();
        out += &format!(">q{oid}_{name}_a\n{}\n", slice(5, 90));
        out += &format!(">q{oid}_{name}_b\n{}\n", slice(a.len() - 80, 70));
        let mut m: Vec<u8> = a[20..100].to_vec();
        for p in [7, 30, 55] {
            m[p] = if m[p] == b'A' { b'C' } else { b'A' };
        }
        out += &format!(">q{oid}_{name}_m\n{}\n", String::from_utf8_lossy(&m));
    }
    out
}

fn setup() -> &'static Setup {
    static S: OnceLock<Setup> = OnceLock::new();
    S.get_or_init(|| {
        let f = taxa_encoded();
        let blast = Blast::find();
        let dir = tempdir();
        let work = dir.path().to_path_buf();
        let db = work.join("ours/db");
        let r = request(&f.spec(), Selection::all(), &db);
        build_req(&r);
        // per-collection views, plus one named view per taxon
        let vr = |kind, name: Option<&str>, selection| gtars_rgblast::api::ViewRequest {
            db: db.clone(),
            kind,
            name: name.map(String::from),
            selection,
            out_dir: None,
            title: None,
            full_db_stats: false,
        };
        gtars_rgblast::api::view(&vr(ViewKind::PerCollection, None, None), &mut |_| true).unwrap();
        for t in ["9606", "10090"] {
            let sel = Selection::attributes::<&str>(&[], &[t]);
            let r = vr(ViewKind::Named, Some(&format!("taxon-{t}")), Some(sel));
            gtars_rgblast::api::view(&r, &mut |_| true).unwrap();
        }
        std::fs::write(work.join("q.fa"), queries(&f.spec())).unwrap();
        let digest_of = sidecars::read_ids(&with_ext(&db, "ids.tsv"))
            .unwrap()
            .into_iter()
            .map(|(id, _, d)| (id, d))
            .collect();
        Setup {
            _dir: dir,
            work,
            blast,
            digest_of,
        }
    })
}

fn blastn(s: &Setup, db: &str, extra: &[&str], fmt: &str) -> String {
    let mut a = vec![
        "-task",
        "megablast",
        "-db",
        db,
        "-query",
        "q.fa",
        "-outfmt",
        fmt,
    ];
    a.extend(extra);
    s.blast.ok(&s.work, "blastn", &a)
}

/// Tabular lines with `sseqid` (field 2) replaced by its digest, sorted.
fn by_digest(s: &Setup, out: &str) -> Vec<String> {
    let mut v: Vec<String> = out
        .lines()
        .map(|l| {
            let mut f: Vec<String> = l.split('\t').map(String::from).collect();
            f[1] = s
                .digest_of
                .get(&f[1])
                .cloned()
                .unwrap_or_else(|| panic!("unknown id {}", f[1]));
            f.join("\t")
        })
        .collect();
    v.sort();
    v
}

#[test]
fn taxids_filter_with_no_expansion() {
    let s = setup();
    let fmt = "6 qseqid sseqid staxids sallseqid";
    let all = blastn(s, "ours/db", &[], fmt);
    let mouse = blastn(
        s,
        "ours/db",
        &["-taxids", "10090", "-no_taxid_expansion"],
        fmt,
    );
    assert!(!mouse.is_empty());
    assert!(
        mouse.lines().count() < all.lines().count(),
        "filtered:\n{mouse}"
    );
    // Only mouse deflines are shown: the SQ. defline of a mouse-only sequence (its taxid is
    // 10090 too) and the mouse_a defline; never a human one.
    for l in mouse.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        assert_eq!(f[2], "10090", "{l}");
        for id in std::iter::once(f[1]).chain(f[3].split(';')) {
            assert!(
                id.starts_with("SQ.") || id.starts_with("gnl|mouse_a|"),
                "{l}"
            );
        }
    }
    assert!(
        mouse
            .lines()
            .any(|l| l.split('\t').nth(1).unwrap().starts_with("gnl|mouse_a|"))
    );
    // the shared decoy shows only its mouse defline
    let decoy: Vec<&str> = mouse.lines().filter(|l| l.contains("_decoy_")).collect();
    assert!(!decoy.is_empty(), "{mouse}");
    for l in decoy {
        assert!(l.ends_with("\tgnl|mouse_a|decoy"), "{l}");
    }
}

#[test]
fn taxids_trap_without_expansion() {
    // No taxonomy4blast.sqlite3 in the working directory or on BLASTDB, and no
    // -no_taxid_expansion: blastn 2.17 warns, exits 0 and does not filter at all (M5).
    let s = setup();
    assert!(!s.work.join("taxonomy4blast.sqlite3").exists());
    let all = blastn(s, "ours/db", &[], TAB);
    let o = s.blast.run(
        &s.work,
        "blastn",
        &[
            "-task",
            "megablast",
            "-db",
            "ours/db",
            "-query",
            "q.fa",
            "-outfmt",
            TAB,
            "-taxids",
            "10090",
        ],
    );
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(
        String::from_utf8_lossy(&o.stdout),
        all,
        "a BLAST+ that now filters here changed the documented trap"
    );
}

#[test]
fn negative_taxids() {
    let s = setup();
    let fmt = "6 qseqid sseqid staxids sallseqid";
    let out = blastn(
        s,
        "ours/db",
        &["-negative_taxids", "9606", "-no_taxid_expansion"],
        fmt,
    );
    assert!(
        !out.lines().any(|l| l.split('\t').nth(2) == Some("9606")),
        "{out}"
    );
    let decoy: Vec<&str> = out.lines().filter(|l| l.contains("_decoy_")).collect();
    assert!(!decoy.is_empty(), "{out}");
    for l in decoy {
        assert!(l.ends_with("\tgnl|mouse_a|decoy"), "{l}");
    }
}

#[test]
fn staxids_lists_every_taxon() {
    let s = setup();
    let out = blastn(s, "ours/db", &[], "6 qseqid sseqid staxids");
    let tax = |pat: &str| -> Vec<String> {
        out.lines()
            .filter(|l| l.starts_with('q') && l.contains(pat))
            .map(|l| l.split('\t').nth(2).unwrap().to_string())
            .collect()
    };
    let decoy = tax("_decoy_a\t");
    assert!(!decoy.is_empty());
    assert!(decoy.iter().all(|t| t == "9606;10090"), "{decoy:?}");
    let chr3 = tax("_chr3_a\t");
    assert!(!chr3.is_empty());
    assert!(chr3.iter().all(|t| t == "9606"), "{chr3:?}");
    let u1 = tax("_u1_a\t");
    assert!(u1.iter().all(|t| t == "0"), "{u1:?}");
}

#[test]
fn outfmt0_no_doubled_name() {
    let s = setup();
    let out = blastn(s, "ours/db", &[], "0");
    let mut seen_title = false;
    for l in out.lines() {
        let t = l.trim_start_matches('>').trim_start();
        if let Some((id, rest)) = t.split_once(' ')
            && let Some((_, name)) = id.split_once(':')
        {
            assert!(
                !rest.starts_with(&format!("{name} ")) && rest != name,
                "doubled name: {l}"
            );
        }
        if l.starts_with('>') && l.contains(" Homo sapiens chr1; 2 assemblies; length:300; md5:") {
            seen_title = true;
        }
    }
    assert!(seen_title, "{out}");
    assert!(
        out.contains(" human_a:chr1 Homo sapiens human_a (GCF_TEST.1)"),
        "{out}"
    );
    assert!(
        out.contains(" human_b:chr1 Homo sapiens human_b (GCA_TEST.1)"),
        "{out}"
    );
    assert!(
        out.contains(" mouse_a:decoy Mus musculus mouse_a (GCF_000002.1)"),
        "{out}"
    );
}

#[test]
fn sscinames_with_taxdb() {
    let Some(taxdb) = std::env::var_os("GTARS_RGBLAST_TAXDB") else {
        eprintln!(
            "skipped: set GTARS_RGBLAST_TAXDB to an unpacked NCBI taxdb (taxdb.btd, taxdb.bti, taxonomy4blast.sqlite3; https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz, 66 MB)"
        );
        return;
    };
    let s = setup();
    let dir = tempdir();
    let work = dir.path().to_path_buf();
    copy_dir(&s.work.join("ours"), &work.join("ours"), &|_| false);
    std::fs::copy(s.work.join("q.fa"), work.join("q.fa")).unwrap();
    for f in ["taxdb.btd", "taxdb.bti", "taxonomy4blast.sqlite3"] {
        let from = Path::new(&taxdb).join(f);
        if std::fs::hard_link(&from, work.join(f)).is_err() {
            std::fs::copy(&from, work.join(f)).unwrap();
        }
    }
    let run = |extra: &[&str]| {
        let mut a = vec![
            "-task",
            "megablast",
            "-db",
            "ours/db",
            "-query",
            "q.fa",
            "-outfmt",
            "6 qseqid sseqid sscinames staxids",
        ];
        a.extend(extra);
        s.blast.ok(&work, "blastn", &a)
    };
    let out = run(&[]);
    assert!(
        out.lines()
            .any(|l| l.contains("_decoy_a\t") && l.contains("\tHomo sapiens;Mus musculus\t")),
        "{out}"
    );
    // Mammalia (40674) expands to both species through taxonomy4blast.sqlite3
    let mam = run(&["-taxids", "40674"]);
    assert!(mam.lines().any(|l| l.contains("\t10090")), "{mam}");
    assert!(mam.lines().any(|l| l.contains("\t9606")), "{mam}");
    assert!(
        !mam.lines().any(|l| l.ends_with("\t0")),
        "nofhr (taxid 0) is not a mammal: {mam}"
    );
}

/// `blastdbcmd -info` sequence count and total bases.
fn info(s: &Setup, db: &str) -> (u64, u64) {
    let out = s.blast.ok(&s.work, "blastdbcmd", &["-db", db, "-info"]);
    let l = out
        .lines()
        .find(|l| l.contains(" sequences; "))
        .unwrap_or_else(|| panic!("{out}"));
    let nums: Vec<u64> = l
        .split(|c: char| !c.is_ascii_digit() && c != ',')
        .filter(|x| !x.is_empty())
        .map(|x| x.replace(',', "").parse().unwrap())
        .collect();
    (nums[0], nums[1])
}

fn nal_value(p: &Path, key: &str) -> u64 {
    std::fs::read_to_string(p)
        .unwrap()
        .lines()
        .find_map(|l| l.strip_prefix(&format!("{key} ")))
        .unwrap()
        .parse()
        .unwrap()
}

/// The blastn `-outfmt 0` footer's letters and sequences.
fn footer(out: &str) -> (u64, u64) {
    let num = |key: &str| -> u64 {
        out.lines()
            .find_map(|l| l.trim().strip_prefix(key))
            .unwrap_or_else(|| panic!("no '{key}' in output"))
            .trim()
            .replace(',', "")
            .parse()
            .unwrap()
    };
    (
        num("Number of sequences in database:"),
        num("Number of letters in database:"),
    )
}

#[test]
fn view_statistics() {
    let s = setup();
    let f = taxa_encoded();
    let idx = gtars_rgblast::views::read_index(&s.work.join("ours/db.views/views.tsv")).unwrap();
    assert_eq!(idx.len(), 5 + 2);
    let mut control_checked = false;
    for row in &idx {
        let (name, kind) = (&row[0], &row[1]);
        let view = format!("ours/db.views/{name}");
        let nal = s.work.join(format!("{view}.nal"));
        let (nseq, len) = (nal_value(&nal, "NSEQ"), nal_value(&nal, "LENGTH"));
        assert_eq!(info(s, &view), (nseq, len), "{name}: blastdbcmd -info");
        let out0 = blastn(s, &view, &[], "0");
        assert_eq!(footer(&out0), (nseq, len), "{name}: blastn footer");
        // the same selection built fresh gives the same hits, E-values and bit scores
        let sel = if kind == "per-collection" {
            Selection::explicit(&[row[2].as_str()])
        } else {
            Selection::attributes::<&str>(&[], &[row[3].as_str()])
        };
        let fresh = s.work.join(format!("fresh/{name}"));
        build(&f.spec(), sel, &fresh);
        let fdb = format!("fresh/{name}");
        assert_eq!(info(s, &fdb), (nseq, len), "{name}: fresh build");
        for task in ["megablast", "blastn"] {
            let a = s.blast.ok(
                &s.work,
                "blastn",
                &[
                    "-task", task, "-db", &view, "-query", "q.fa", "-outfmt", TAB,
                ],
            );
            let b = s.blast.ok(
                &s.work,
                "blastn",
                &["-task", task, "-db", &fdb, "-query", "q.fa", "-outfmt", TAB],
            );
            assert!(
                !a.is_empty() || name == "human_masked",
                "{name} {task}: no hits"
            );
            assert_eq!(by_digest(s, &a), by_digest(s, &b), "{name} {task}");
        }
        // negative control: blastdb_aliastool writes SeqDB's approximate LENGTH
        let ids: String = std::fs::read_to_string(s.work.join("ours/db.members.tsv"))
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .filter_map(|l| {
                let x: Vec<&str> = l.split('\t').collect();
                row[2]
                    .split(',')
                    .any(|c| c == x[2])
                    .then(|| format!("{}:{}\n", x[2], x[4]))
            })
            .collect();
        let list = format!("ctl_{name}.txt");
        std::fs::write(s.work.join(&list), ids).unwrap();
        s.blast.ok(
            &s.work,
            "blastdb_aliastool",
            &[
                "-db",
                "ours/db",
                "-dbtype",
                "nucl",
                "-seqidlist",
                &list,
                "-out",
                &format!("ctl_{name}"),
            ],
        );
        let approx = nal_value(&s.work.join(format!("ctl_{name}.nal")), "LENGTH");
        let oids = sidecars::read_oids(&s.work.join("ours/db.oids.tsv")).unwrap();
        let members = sidecars::read_members(&s.work.join("ours/db.members.tsv")).unwrap();
        let mut set: Vec<u32> = members
            .iter()
            .filter(|m| row[2].split(',').any(|c| c == m.collection))
            .map(|m| m.oid)
            .collect();
        set.sort();
        set.dedup();
        let formula: u64 = set
            .iter()
            .map(|&o| 4 * (oids[o as usize].length / 4) + (o as u64 % 4))
            .sum();
        assert_eq!(approx, formula, "{name}: SeqDB's approximation (M7)");
        if approx != len {
            control_checked = true;
        }
    }
    assert!(
        control_checked,
        "no view where the approximation differs from the exact length"
    );
}

#[test]
fn view_full_db_stats() {
    let s = setup();
    let view = gtars_rgblast::api::view(
        &ViewRequest {
            db: s.work.join("ours/db"),
            kind: ViewKind::Named,
            name: Some("mouse_full".into()),
            selection: Some(Selection::explicit(&["mouse_a"])),
            out_dir: Some(s.work.join("fullstats")),
            title: None,
            full_db_stats: true,
        },
        &mut |_| true,
    )
    .unwrap();
    let v = &view.views[0];
    assert_eq!(
        info(s, "fullstats/mouse_full"),
        (v.nseq, v.length),
        "-info shows the view itself"
    );
    for task in ["megablast", "blastn"] {
        let a = s.blast.ok(
            &s.work,
            "blastn",
            &[
                "-task",
                task,
                "-db",
                "fullstats/mouse_full",
                "-query",
                "q.fa",
                "-outfmt",
                TAB,
            ],
        );
        let full = s.blast.ok(
            &s.work,
            "blastn",
            &[
                "-task", task, "-db", "ours/db", "-query", "q.fa", "-outfmt", TAB,
            ],
        );
        let full: HashMap<(String, String), Vec<String>> =
            by_digest(s, &full)
                .into_iter()
                .fold(HashMap::new(), |mut m, l| {
                    let f: Vec<&str> = l.split('\t').collect();
                    m.entry((f[0].to_string(), f[1].to_string()))
                        .or_default()
                        .push(format!("{}\t{}", f[2], f[3]));
                    m
                });
        let a = by_digest(s, &a);
        assert!(!a.is_empty());
        for l in &a {
            let f: Vec<&str> = l.split('\t').collect();
            let want = &full[&(f[0].to_string(), f[1].to_string())];
            assert!(
                want.contains(&format!("{}\t{}", f[2], f[3])),
                "{task}: {l} vs {want:?}"
            );
        }
    }
}
