//! The `txome` fixture (plan 8 step 33): a deterministic transcriptome-and-microbiome store with
//! thousands of short subjects, generated at test time (seeded, no committed FASTA) and built
//! once per storage mode per test binary. `tests/data/fixtures/txome.expected.json` pins the
//! generated FASTAs (md5) and the plan totals; `txome_fixture_matches_expected` checks them.
//!
//! - `txome_a.fa`: 1,500 transcripts of 150–4,000 bp; 3 % carry an N run, 1 % IUPAC codes, 1 %
//!   repeat an earlier transcript's sequence under a new name, and 30 isoform groups of 3–25
//!   transcripts share one 300-bp exon per group.
//! - `txome_b.fa`: 800 transcripts; 300 share a sequence with `txome_a` (150 under the same name,
//!   150 under another).
//! - `bins/`: 40 collections of 10–30 contigs of 500–20,000 bp, each with its own `name:` alias
//!   and FHR taxid. Contig names carry the bin's alias, so makeblastdb's `-taxid_map` (which
//!   looks ids up by bare tag) can express the taxa for the canary; the many-bins spec repeats
//!   names across bins instead, as assembler output does.

use super::store::{Fixture, tempdir};
use gtars_refget::store::{FastaImportOptions, RefgetStore, StorageMode};
use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// How sequences are laid out on disk. Only the standalone layout (one `.seq` file per
/// sequence) exists today; the packed layout (plan 8 step 18) adds a variant here and to the
/// fixture cache key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Layout {
    Standalone,
}

/// What the generator writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TxomeSpec {
    /// Write `txome_a.fa` and `txome_b.fa`.
    pub transcripts: bool,
    /// Number of `bins/` collections.
    pub bins: usize,
    /// Contigs per bin (inclusive range).
    pub bin_contigs: (u32, u32),
    /// Contig lengths (inclusive range).
    pub bin_len: (u32, u32),
    /// Prefix contig names with the bin's alias (else names repeat across bins).
    pub bin_unique_names: bool,
}

impl TxomeSpec {
    /// The `txome` fixture.
    pub const STANDARD: TxomeSpec = TxomeSpec {
        transcripts: true,
        bins: 40,
        bin_contigs: (10, 30),
        bin_len: (500, 20_000),
        bin_unique_names: true,
    };

    /// `n` tiny bins and nothing else (plan 8 step 24: thousands of collections).
    pub fn many_bins(n: usize) -> TxomeSpec {
        TxomeSpec {
            transcripts: false,
            bins: n,
            bin_contigs: (1, 3),
            bin_len: (100, 400),
            bin_unique_names: false,
        }
    }
}

/// xorshift64*, seeded: the same spec always gives the same bytes.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    /// Uniform in [lo, hi].
    fn range(&mut self, lo: u32, hi: u32) -> u32 {
        lo + (self.next() % (hi - lo + 1) as u64) as u32
    }
    fn chance(&mut self, percent: u64) -> bool {
        self.next() % 100 < percent
    }
    fn acgt(&mut self, n: u32) -> Vec<u8> {
        (0..n)
            .map(|_| b"ACGT"[(self.next() >> 33) as usize & 3])
            .collect()
    }
}

/// One collection the generator wrote: its FASTA, `name:` alias and FHR file (bins only).
#[derive(Debug, Clone)]
pub struct Source {
    pub fasta: PathBuf,
    pub alias: String,
    pub fhr: Option<PathBuf>,
    /// FHR taxid of a bin.
    pub taxid: Option<i32>,
}

fn write_fasta(path: &Path, records: &[(String, Vec<u8>)]) {
    let mut s = String::new();
    for (name, seq) in records {
        writeln!(s, ">{name}").unwrap();
        for line in seq.chunks(60) {
            s.push_str(std::str::from_utf8(line).unwrap());
            s.push('\n');
        }
    }
    std::fs::write(path, s).unwrap();
}

/// A transcript length: mostly short, a long tail up to 4,000 bp.
fn tx_len(r: &mut Rng) -> u32 {
    let u = r.next() % 10_000;
    150 + (3_850 * u * u / 100_000_000) as u32
}

fn transcripts_a(r: &mut Rng) -> Vec<(String, Vec<u8>)> {
    let mut out: Vec<(String, Vec<u8>)> = Vec::with_capacity(1_500);
    // 30 isoform groups: a shared 300-bp exon between a random head and tail
    for g in 0..30u32 {
        let exon = r.acgt(300);
        let members = if g % 5 == 0 { 25 } else { r.range(3, 12) };
        for k in 0..members {
            let head = r.range(0, 1_500);
            let tail = r.range(0, 1_500);
            let mut s = r.acgt(head);
            s.extend_from_slice(&exon);
            s.extend(r.acgt(tail));
            out.push((format!("gene{g:02}-iso{k:02}"), s));
        }
    }
    let mut i = 0;
    while out.len() < 1_500 {
        i += 1;
        let name = format!("TXA{i:06}.1");
        if out.len() > 100 && r.chance(1) {
            // the same sequence as an earlier transcript, under a new name
            let j = (r.next() % out.len() as u64) as usize;
            let s = out[j].1.clone();
            out.push((name, s));
            continue;
        }
        let len = tx_len(r);
        let mut s = r.acgt(len);
        if r.chance(3) {
            let run = r.range(1, 60).min(s.len() as u32 - 1);
            let at = r.range(0, s.len() as u32 - run) as usize;
            s[at..at + run as usize].fill(b'N');
        } else if r.chance(1) {
            let codes = b"RYKMSWBDHV";
            for p in (r.range(0, 30) as usize..s.len()).step_by(37) {
                s[p] = codes[(r.next() % codes.len() as u64) as usize];
            }
        }
        out.push((name, s));
    }
    out
}

fn transcripts_b(r: &mut Rng, a: &[(String, Vec<u8>)]) -> Vec<(String, Vec<u8>)> {
    let mut out = Vec::with_capacity(800);
    for k in 0..300usize {
        let (name, s) = &a[(k * 5 + 3) % a.len()];
        let name = if k % 2 == 0 {
            name.clone()
        } else {
            format!("TXB_shared{k:04}.2")
        };
        out.push((name, s.clone()));
    }
    for k in 0..500u32 {
        let len = tx_len(r);
        out.push((format!("TXB{k:06}.1"), r.acgt(len)));
    }
    out
}

/// Write the fixture's FASTAs (and bin FHR files) under `dir`; returns one source per collection.
pub fn write_txome_sources(dir: &Path, spec: &TxomeSpec) -> Vec<Source> {
    std::fs::create_dir_all(dir.join("bins")).unwrap();
    let mut r = Rng(0x7478_6F6D_6520_2026);
    let mut out = Vec::new();
    if spec.transcripts {
        let a = transcripts_a(&mut r);
        let b = transcripts_b(&mut r, &a);
        for (name, recs) in [("txome_a", &a), ("txome_b", &b)] {
            let fasta = dir.join(format!("{name}.fa"));
            write_fasta(&fasta, recs);
            out.push(Source {
                fasta,
                alias: name.to_string(),
                fhr: None,
                taxid: None,
            });
        }
    }
    let width = spec.bins.to_string().len().max(2);
    for b in 0..spec.bins {
        let alias = format!("bin{:0width$}", b + 1);
        let n = r.range(spec.bin_contigs.0, spec.bin_contigs.1);
        let first = r.range(1, 500);
        let recs: Vec<(String, Vec<u8>)> = (0..n)
            .map(|k| {
                let len = r.range(spec.bin_len.0, spec.bin_len.1);
                let name = format!("k141_{}", first + 3 * k);
                let name = if spec.bin_unique_names {
                    format!("{alias}_{name}")
                } else {
                    name
                };
                (name, r.acgt(len))
            })
            .collect();
        let fasta = dir.join("bins").join(format!("{alias}.fa"));
        write_fasta(&fasta, &recs);
        let taxid = 100_000 + 7 * b as i32;
        let fhr = dir.join("bins").join(format!("{alias}.fhr.json"));
        std::fs::write(
            &fhr,
            format!(
                "{{\n  \"schema\": \"https://raw.githubusercontent.com/FAIR-bioHeaders/FHR-Specification/main/fhr.json\",\n  \"genome\": \"Bacterium fixturus {b}\",\n  \"taxon\": {{\"name\": \"Bacterium fixturus {b}\", \"uri\": \"https://identifiers.org/taxonomy:{taxid}\"}}\n}}\n"
            ),
        )
        .unwrap();
        out.push(Source {
            fasta,
            alias,
            fhr: Some(fhr),
            taxid: Some(taxid),
        });
    }
    out
}

/// Generate the sources into `<dir>.src` and build the store in `dir` (which must not exist
/// yet). Returns `name:` alias -> collection digest.
pub fn build_txome_store(
    dir: &Path,
    mode: StorageMode,
    layout: Layout,
    spec: &TxomeSpec,
) -> HashMap<String, String> {
    let mut s = match layout {
        Layout::Standalone => RefgetStore::on_disk(dir).unwrap(),
    };
    s.set_quiet(true);
    s.set_encoding_mode(mode).unwrap();
    let mut digests = HashMap::new();
    for src in write_txome_sources(&dir.with_extension("src"), spec) {
        let (m, _) = s
            .add_sequence_collection_from_fasta(
                &src.fasta,
                FastaImportOptions::new().collection_alias("name", &src.alias),
            )
            .unwrap();
        if let Some(f) = &src.fhr {
            s.load_fhr_metadata(&m.digest, f.to_str().unwrap()).unwrap();
        }
        digests.insert(src.alias, m.digest.clone());
    }
    s.write().unwrap();
    digests
}

const MODES: [StorageMode; 3] = [StorageMode::Encoded, StorageMode::Raw, StorageMode::Zstd];
const LAYOUTS: [Layout; 1] = [Layout::Standalone];

static TXOME: [[OnceLock<Fixture>; 1]; 3] =
    [[OnceLock::new()], [OnceLock::new()], [OnceLock::new()]];

/// The `txome` fixture store in a storage mode and layout, built once per test binary.
pub fn txome_in(mode: StorageMode, layout: Layout) -> &'static Fixture {
    let m = MODES.iter().position(|x| *x == mode).unwrap();
    let l = LAYOUTS.iter().position(|x| *x == layout).unwrap();
    TXOME[m][l].get_or_init(|| {
        let dir = tempdir();
        let path = dir.path().join("txome");
        let digests = build_txome_store(&path, mode, layout, &TxomeSpec::STANDARD);
        RefgetStore::open_local(&path).unwrap();
        Fixture {
            _dir: dir,
            path,
            digests,
            orphan: String::new(),
        }
    })
}

/// The `txome` fixture store (standalone) in a storage mode.
pub fn txome(mode: StorageMode) -> &'static Fixture {
    txome_in(mode, Layout::Standalone)
}
