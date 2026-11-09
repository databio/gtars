//! The synthetic `origin_store`: two species and a collection with no FHR, every sequence
//! 1.2-3 kb of seeded random DNA.
//!
//! - species A (9606): A1 = chr1, chr2, chrM; A2 = chr1, chr2 with 5 SNPs, chrM;
//!   A3 = chr1, chr2, a decoy chrX9 only A3 holds, and vec1;
//! - species B (10090): B1 = chr1b, chrMb, chr2b; B2 = chr1b, chrMb, vec1 (also in A3);
//! - C1 (no FHR): chrU.
//!
//! The representatives are A1 and B1 (most-shared; B1 is larger) and C1 (only).

use gtars_refget::store::{FastaImportOptions, RefgetStore, StorageMode};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn dna(&mut self, n: usize) -> String {
        (0..n)
            .map(|_| b"ACGT"[(self.next() % 4) as usize] as char)
            .collect()
    }
}

pub struct Seqs {
    pub s1: String,
    pub s2: String,
    pub s2p: String,
    pub sm: String,
    pub d: String,
    pub v: String,
    pub u: String,
    pub b1: String,
    pub bm: String,
    pub b2: String,
    pub rand: String,
}

pub fn seqs() -> Seqs {
    let mut r = Rng(0x9e3779b97f4a7c15);
    let s2 = r.dna(2500);
    let mut s2p = s2.clone().into_bytes();
    for p in [300, 500, 700, 900, 1100] {
        s2p[p] = match s2p[p] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    Seqs {
        s1: r.dna(3000),
        s2,
        s2p: String::from_utf8(s2p).unwrap(),
        sm: r.dna(2000),
        d: r.dna(2200),
        v: r.dna(1200),
        u: r.dna(2000),
        b1: r.dna(3000),
        bm: r.dna(2000),
        b2: r.dna(2600),
        rand: r.dna(300),
    }
}

/// (alias, species taxid and name, records)
pub type Collection = (
    &'static str,
    Option<(u32, &'static str)>,
    Vec<(&'static str, String)>,
);

pub fn collections(s: &Seqs) -> Vec<Collection> {
    let h = Some((9606, "Homo sapiens"));
    let m = Some((10090, "Mus musculus"));
    vec![
        (
            "A1",
            h,
            vec![
                ("chr1", s.s1.clone()),
                ("chr2", s.s2.clone()),
                ("chrM", s.sm.clone()),
            ],
        ),
        (
            "A2",
            h,
            vec![
                ("chr1", s.s1.clone()),
                ("chr2", s.s2p.clone()),
                ("chrM", s.sm.clone()),
            ],
        ),
        (
            "A3",
            h,
            vec![
                ("chr1", s.s1.clone()),
                ("chr2", s.s2.clone()),
                ("chrX9", s.d.clone()),
                ("vec1", s.v.clone()),
            ],
        ),
        (
            "B1",
            m,
            vec![
                ("chr1b", s.b1.clone()),
                ("chrMb", s.bm.clone()),
                ("chr2b", s.b2.clone()),
            ],
        ),
        (
            "B2",
            m,
            vec![
                ("chr1b", s.b1.clone()),
                ("chrMb", s.bm.clone()),
                ("vec1", s.v.clone()),
            ],
        ),
        ("C1", None, vec![("chrU", s.u.clone())]),
    ]
}

/// The queries, in order.
pub fn queries(s: &Seqs) -> Vec<(&'static str, String)> {
    vec![
        ("q_s1", s.s1[100..1100].to_string()),
        ("q_snp", s.s2p[200..1200].to_string()),
        ("q_decoy", s.d[100..900].to_string()),
        ("q_vec", s.v[50..650].to_string()),
        ("q_bm", s.bm[100..1100].to_string()),
        (
            "q_chim",
            format!("{}{}", &s.s1[2000..2400], &s.v[700..1100]),
        ),
        ("q_u", s.u[100..700].to_string()),
        ("q_rand", s.rand.clone()),
    ]
}

fn fasta(records: &[(&str, String)]) -> String {
    records
        .iter()
        .map(|(n, s)| format!(">{n}\n{s}\n"))
        .collect()
}

/// Build the store at `dir` (which must not exist) and write the queries to `queries`.
pub fn build_origin_store(dir: &Path, queries_path: &Path) {
    let s = seqs();
    let mut st = RefgetStore::on_disk(dir).unwrap();
    st.set_quiet(true);
    st.set_encoding_mode(StorageMode::Encoded).unwrap();
    let src = dir.with_extension("src");
    std::fs::create_dir_all(&src).unwrap();
    for (alias, taxon, recs) in collections(&s) {
        let fa = src.join(format!("{alias}.fa"));
        std::fs::write(&fa, fasta(&recs)).unwrap();
        let (m, _) = st
            .add_sequence_collection_from_fasta(
                &fa,
                FastaImportOptions::new().collection_alias("name", alias),
            )
            .unwrap();
        if let Some((taxid, name)) = taxon {
            let j = src.join(format!("{alias}.fhr.json"));
            std::fs::write(
                &j,
                format!(
                    "{{\"schema\":\"https://raw.githubusercontent.com/FAIR-bioHeaders/FHR-Specification/main/fhr.json\",\"schemaVersion\":1.0,\"genome\":\"{alias}\",\"taxon\":{{\"name\":\"{name}\",\"uri\":\"https://identifiers.org/taxonomy:{taxid}\"}}}}"
                ),
            )
            .unwrap();
            let fhr = gtars_refget::store::load_from_json(j.to_str().unwrap()).unwrap();
            st.set_fhr_metadata(&m.digest, fhr).unwrap();
        }
    }
    st.write().unwrap();
    std::fs::write(queries_path, fasta(&queries(&s))).unwrap();
}

pub struct OriginFixture {
    _dir: tempfile::TempDir,
    pub store: PathBuf,
    pub queries: PathBuf,
}

/// The origin store, built once per test binary.
pub fn origin_store() -> &'static OriginFixture {
    static F: OnceLock<OriginFixture> = OnceLock::new();
    F.get_or_init(|| {
        let dir = super::store::tempdir();
        let store = dir.path().join("origin_store");
        let queries = dir.path().join("queries.fa");
        build_origin_store(&store, &queries);
        OriginFixture {
            _dir: dir,
            store,
            queries,
        }
    })
}

pub const HITS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/fixtures/origin/origin_store_hits"
);
pub const ORIGIN_FIXTURES: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/fixtures/origin");
