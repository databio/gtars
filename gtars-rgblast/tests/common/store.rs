//! Fixture stores built from `tests/data/fixtures` with the gtars-refget API, once per
//! storage mode and layout per test binary, and helpers to build databases from them.

use gtars_refget::store::{FastaImportOptions, RefgetStore, StorageMode};
use gtars_rgblast::plan::{PlanOptions, SubjectPlan};
use gtars_rgblast::select::nstats::NstatsOptions;
use gtars_rgblast::select::{ResolveOptions, Selection};
use gtars_rgblast::{BuildOptions, BuildRequest, BuildSummary, OpenedStore, ProgressEvent};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

pub const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/fixtures");
pub const TAXA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/fixtures/taxa");
pub const EXPECTED: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/expected");
/// Fixed database date, so outputs are deterministic.
pub const DATE: &str = "Nov 8, 2026  12:51 PM";
pub const TITLE: &str = "fixture";

/// The collections of the fixture store, by their `name:` alias.
pub const COLLECTIONS: [&str; 5] = ["asmA", "asmB", "asmC", "asmD", "numeric"];
/// The collections of the `taxa` fixture store, by their `name:` alias.
pub const TAXA_COLLECTIONS: [&str; 5] = ["human_a", "human_b", "human_masked", "mouse_a", "nofhr"];

pub struct Fixture {
    pub(crate) _dir: tempfile::TempDir,
    /// The store directory.
    pub path: PathBuf,
    /// `name:` alias -> collection digest.
    pub digests: HashMap<String, String>,
    /// Digest of the orphan sequence.
    pub orphan: String,
}

impl Fixture {
    pub fn spec(&self) -> String {
        self.path.display().to_string()
    }
    pub fn digest(&self, alias: &str) -> String {
        self.digests[alias].clone()
    }
}

fn tmp_root() -> PathBuf {
    let root = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("rgblast-tests");
    std::fs::create_dir_all(&root).unwrap();
    // remove fixture copies left by earlier runs (statics are never dropped)
    if let Ok(rd) = std::fs::read_dir(&root) {
        for e in rd.flatten() {
            let old = e
                .metadata()
                .and_then(|m| m.modified())
                .map(|t| t.elapsed().map(|d| d.as_secs() > 3600).unwrap_or(false))
                .unwrap_or(false);
            if old {
                let _ = std::fs::remove_dir_all(e.path());
            }
        }
    }
    root
}

/// A fresh temporary directory under the cargo target tmp dir.
pub fn tempdir() -> tempfile::TempDir {
    tempfile::Builder::new()
        .prefix("t-")
        .tempdir_in(tmp_root())
        .unwrap()
}

fn fasta(name: &str) -> PathBuf {
    Path::new(FIXTURES).join(name)
}

/// Copy fixture files into `dir` and return their paths there. gtars caches a `<fasta>.rgsi`
/// next to each imported FASTA, so imports must never read from the source tree.
pub fn fixture_copies(dir: &Path, names: &[&str]) -> Vec<PathBuf> {
    std::fs::create_dir_all(dir).unwrap();
    names
        .iter()
        .map(|n| {
            let to = dir.join(n);
            std::fs::copy(fasta(n), &to).unwrap();
            to
        })
        .collect()
}

/// Read `aliases.tsv`: (kind, namespace, alias, collection, name).
fn alias_spec() -> Vec<[String; 5]> {
    alias_spec_in(&fasta("aliases.tsv"))
}

fn alias_spec_in(p: &Path) -> Vec<[String; 5]> {
    std::fs::read_to_string(p)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            [0, 1, 2, 3, 4].map(|i| f.get(i).unwrap_or(&"").to_string())
        })
        .collect()
}

/// How sequences are laid out on disk: one `.seq` file per sequence, or records concatenated
/// into pack files (`RefgetStore::on_disk_packed`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Layout {
    Standalone,
    Packed,
}

pub const LAYOUTS: [Layout; 2] = [Layout::Standalone, Layout::Packed];

/// The pack size cap of the small packed fixtures, so their records spread over several
/// packs and most sit at a nonzero, unaligned offset.
pub const PACK_CAP: u64 = 64 << 10;

/// A new store in `dir` (which must not exist yet) in `layout`, with gtars-refget's default
/// pack cap.
pub fn open_layout(dir: &Path, layout: Layout) -> RefgetStore {
    let mut s = match layout {
        Layout::Standalone => RefgetStore::on_disk(dir).unwrap(),
        Layout::Packed => RefgetStore::on_disk_packed(dir).unwrap(),
    };
    s.set_quiet(true);
    s
}

/// `open_layout` with small packs (`PACK_CAP`).
pub fn open_small(dir: &Path, layout: Layout) -> RefgetStore {
    let mut s = open_layout(dir, layout);
    if layout == Layout::Packed {
        s.set_pack_cap_bytes(PACK_CAP);
    }
    s
}

/// The store at `dir` is packed: no `.seq` file under `sequences/`, at least one pack, and the
/// manifest says so.
pub fn assert_packed(dir: &Path) {
    fn seq_files(d: &Path) -> usize {
        std::fs::read_dir(d).map_or(0, |rd| {
            rd.flatten()
                .map(|e| {
                    let p = e.path();
                    if p.is_dir() {
                        seq_files(&p)
                    } else {
                        p.extension().is_some_and(|x| x == "seq") as usize
                    }
                })
                .sum()
        })
    }
    assert_eq!(seq_files(&dir.join("sequences")), 0, "{}", dir.display());
    let packs = std::fs::read_dir(dir.join("packs"))
        .unwrap()
        .flatten()
        .filter(|e| e.path().extension().is_some_and(|x| x == "pack"))
        .count();
    assert!(packs > 0, "{}: no packs", dir.display());
    let manifest = std::fs::read_to_string(dir.join("rgstore.json")).unwrap();
    assert!(manifest.contains("\"layout\": \"packed\""), "{manifest}");
}

/// Build the fixture store in `dir` (which must not exist yet) in `mode` and `layout`.
pub fn build_fixture_store(
    dir: &Path,
    mode: StorageMode,
    layout: Layout,
) -> (HashMap<String, String>, String) {
    let mut s = open_small(dir, layout);
    s.set_encoding_mode(mode).unwrap();
    let src = dir.with_extension("src");
    fixture_copies(
        &src,
        &["asmA.fa", "asmB.fa", "asmC.fa", "asmD.fa", "numeric.fa"],
    );
    let mut digests = HashMap::new();
    for (file, alias) in [
        ("asmA.fa", "asmA"),
        ("asmB.fa", "asmB"),
        ("asmC.fa", "asmC"),
        ("asmD.fa", "asmD"),
        ("numeric.fa", "numeric"),
    ] {
        let (m, _) = s
            .add_sequence_collection_from_fasta(
                src.join(file),
                FastaImportOptions::new().collection_alias("name", alias),
            )
            .unwrap();
        digests.insert(alias.to_string(), m.digest.clone());
    }
    for [kind, ns, alias, coll, name] in alias_spec() {
        let cd = &digests[&coll];
        if kind == "collection" {
            s.add_collection_alias(&ns, &alias, cd).unwrap();
        } else {
            let sd = s
                .get_sequence_by_name(cd, &name)
                .unwrap()
                .metadata()
                .sha512t24u
                .clone();
            s.add_sequence_alias(&ns, &alias, &sd).unwrap();
        }
    }
    s.load_fhr_metadata(&digests["asmA"], fasta("asmA.fhr.json").to_str().unwrap())
        .unwrap();
    let orphan_fa = std::fs::read_to_string(fasta("orphan.fa")).unwrap();
    let seq: String = orphan_fa.lines().filter(|l| !l.starts_with('>')).collect();
    let orphan = s.ingest_sequence("orphan1", seq.as_bytes()).unwrap();
    s.write().unwrap();
    (digests, orphan)
}

/// Build a store in `dir` (which must not exist yet) holding one collection imported from
/// `fasta` with the `name:` alias `alias`; returns the collection digest.
pub fn build_store_from_fasta(
    dir: &Path,
    fasta: &Path,
    alias: &str,
    mode: StorageMode,
    layout: Layout,
) -> String {
    let mut s = open_layout(dir, layout);
    s.set_encoding_mode(mode).unwrap();
    let (m, _) = s
        .add_sequence_collection_from_fasta(
            fasta,
            FastaImportOptions::new().collection_alias("name", alias),
        )
        .unwrap();
    s.write().unwrap();
    m.digest.clone()
}

/// Build the `taxa` fixture store in `dir` (which must not exist yet) in `mode` and `layout`:
/// one collection per FASTA (no alias at import), the aliases of `taxa/aliases.tsv` in file
/// order, the FHR of `taxa/fhr/<collection>.json`, and the orphan.
pub fn build_taxa_store(
    dir: &Path,
    mode: StorageMode,
    layout: Layout,
) -> (HashMap<String, String>, String) {
    let mut s = open_small(dir, layout);
    s.set_encoding_mode(mode).unwrap();
    let src = dir.with_extension("src");
    std::fs::create_dir_all(&src).unwrap();
    let taxa = Path::new(TAXA);
    let mut digests = HashMap::new();
    for c in TAXA_COLLECTIONS {
        let to = src.join(format!("{c}.fa"));
        std::fs::copy(taxa.join(format!("{c}.fa")), &to).unwrap();
        let (m, _) = s
            .add_sequence_collection_from_fasta(&to, FastaImportOptions::new())
            .unwrap();
        digests.insert(c.to_string(), m.digest.clone());
    }
    for [kind, ns, alias, coll, name] in alias_spec_in(&taxa.join("aliases.tsv")) {
        let cd = &digests[&coll];
        if kind == "collection" {
            s.add_collection_alias(&ns, &alias, cd).unwrap();
        } else {
            let sd = s
                .get_sequence_by_name(cd, &name)
                .unwrap()
                .metadata()
                .sha512t24u
                .clone();
            s.add_sequence_alias(&ns, &alias, &sd).unwrap();
        }
    }
    for c in TAXA_COLLECTIONS {
        let f = taxa.join("fhr").join(format!("{c}.json"));
        if f.exists() {
            let m = gtars_refget::store::load_from_json(f.to_str().unwrap()).unwrap();
            s.set_fhr_metadata(&digests[c], m).unwrap();
        }
    }
    let orphan_fa = std::fs::read_to_string(taxa.join("orphan.fa")).unwrap();
    let seq: String = orphan_fa.lines().filter(|l| !l.starts_with('>')).collect();
    let orphan = s.ingest_sequence("o1", seq.as_bytes()).unwrap();
    s.write().unwrap();
    (digests, orphan)
}

/// Index of a (mode, layout) fixture cache slot.
fn slot(mode: StorageMode, layout: Layout) -> (usize, usize) {
    let m = match mode {
        StorageMode::Encoded => 0,
        StorageMode::Raw => 1,
        StorageMode::Zstd => 2,
    };
    (m, LAYOUTS.iter().position(|l| *l == layout).unwrap())
}

type Slots = [[OnceLock<Fixture>; 2]; 3];

const fn slots() -> Slots {
    [
        [OnceLock::new(), OnceLock::new()],
        [OnceLock::new(), OnceLock::new()],
        [OnceLock::new(), OnceLock::new()],
    ]
}

static TAXA_BY_MODE: Slots = slots();

/// The `taxa` fixture store (standalone) in a storage mode, built once per test binary.
pub fn taxa(mode: StorageMode) -> &'static Fixture {
    taxa_in(mode, Layout::Standalone)
}

/// The `taxa` fixture store in a storage mode and layout, built once per test binary. The
/// store directory is called `taxa` in both layouts.
pub fn taxa_in(mode: StorageMode, layout: Layout) -> &'static Fixture {
    let (m, l) = slot(mode, layout);
    TAXA_BY_MODE[m][l].get_or_init(|| {
        let dir = tempdir();
        let path = dir.path().join("taxa");
        let (digests, orphan) = build_taxa_store(&path, mode, layout);
        RefgetStore::open_local(&path).unwrap();
        Fixture {
            _dir: dir,
            path,
            digests,
            orphan,
        }
    })
}

pub fn taxa_encoded() -> &'static Fixture {
    taxa(StorageMode::Encoded)
}

/// A per-test-binary N-count cache, so tests never touch the user's cache.
pub fn nstats_cache() -> PathBuf {
    static C: OnceLock<tempfile::TempDir> = OnceLock::new();
    C.get_or_init(|| {
        tempfile::Builder::new()
            .prefix("nstats-")
            .tempdir_in(tmp_root())
            .unwrap()
    })
    .path()
    .to_path_buf()
}

pub fn nstats_opts() -> NstatsOptions {
    NstatsOptions {
        cache_dir: Some(nstats_cache()),
        ..Default::default()
    }
}

static FIXTURES_BY_MODE: Slots = slots();

/// The fixture store (standalone) in a storage mode, built once per test binary.
pub fn fixture(mode: StorageMode) -> &'static Fixture {
    fixture_in(mode, Layout::Standalone)
}

/// The fixture store in a storage mode and layout, built once per test binary. The store
/// directory is called `fixture` in both layouts, so labels and manifests match.
pub fn fixture_in(mode: StorageMode, layout: Layout) -> &'static Fixture {
    let (m, l) = slot(mode, layout);
    FIXTURES_BY_MODE[m][l].get_or_init(|| {
        let dir = tempdir();
        let path = dir.path().join("fixture");
        let (digests, orphan) = build_fixture_store(&path, mode, layout);
        // reopen, so reads go through disk
        RefgetStore::open_local(&path).unwrap();
        Fixture {
            _dir: dir,
            path,
            digests,
            orphan,
        }
    })
}

pub fn encoded() -> &'static Fixture {
    fixture(StorageMode::Encoded)
}

/// A store with one collection whose names collide after sanitizing (`a|b`, `a_b`).
pub fn build_collide_store(dir: &Path) {
    let src = fixture_copies(&dir.with_extension("src"), &["collide.fa"]);
    let mut s = RefgetStore::on_disk(dir).unwrap();
    s.set_quiet(true);
    s.add_sequence_collection_from_fasta(
        &src[0],
        FastaImportOptions::new().collection_alias("name", "collide"),
    )
    .unwrap();
    s.write().unwrap();
}

/// Default test options: fixed date, fixed title, 2 threads, no progress throttling.
pub fn opts(out: &Path) -> BuildOptions {
    BuildOptions {
        out: out.to_path_buf(),
        title: Some(TITLE.into()),
        threads: 2,
        date: Some(DATE.into()),
        progress_interval_ms: 0,
        nstats: nstats_opts(),
        ..Default::default()
    }
}

pub fn request(store: &str, selection: Selection, out: &Path) -> BuildRequest {
    BuildRequest {
        store: store.into(),
        cache_dir: None,
        selection,
        options: opts(out),
    }
}

pub fn build_req(req: &BuildRequest) -> BuildSummary {
    gtars_rgblast::build(req, &mut |_: &ProgressEvent| true).unwrap()
}

/// Build a database of `selection` from a store.
pub fn build(store: &str, selection: Selection, out: &Path) -> BuildSummary {
    build_req(&request(store, selection, out))
}

/// Resolve options for tests (the test N-count cache).
pub fn resolve_opts() -> ResolveOptions {
    ResolveOptions {
        nstats: nstats_opts(),
        threads: 2,
        ..Default::default()
    }
}

/// Open a store and build the plan for a selection (default plan options).
pub fn plan_for(store: &str, selection: &Selection) -> (OpenedStore, SubjectPlan) {
    plan_with(store, selection, &PlanOptions::default())
}

pub fn plan_with(
    store: &str,
    selection: &Selection,
    o: &PlanOptions,
) -> (OpenedStore, SubjectPlan) {
    let mut st = OpenedStore::open(store, None).unwrap();
    let (sel, _) =
        gtars_rgblast::api::resolve_in_store(&mut st, selection, &resolve_opts(), &mut |_| true)
            .unwrap();
    let plan = gtars_rgblast::api::plan_for_selection(&mut st, &sel, o).unwrap();
    (st, plan)
}

/// A database manifest without the store field (where the store is, and its digests).
pub fn manifest_sans_store(b: &[u8]) -> serde_json::Value {
    let mut v: serde_json::Value = serde_json::from_slice(b).unwrap();
    v.as_object_mut().unwrap().remove("store");
    v
}

/// Every database file identical; manifests equal apart from where the store is.
pub fn same_database(a: &Path, b: &Path) {
    let (a, b) = (files_in(a), files_in(b));
    assert_eq!(a.len(), b.len());
    for ((n1, b1), (n2, b2)) in a.iter().zip(&b) {
        assert_eq!(n1, n2);
        if n1.ends_with(".manifest.json") {
            let (ma, mb) = (manifest_sans_store(b1), manifest_sans_store(b2));
            assert_eq!(ma, mb, "manifest beyond the store location");
            continue;
        }
        assert!(b1 == b2, "{n1} differs between the two builds");
    }
}

/// Every file in a directory, by name.
pub fn files_in(dir: &Path) -> Vec<(String, Vec<u8>)> {
    let mut v: Vec<(String, Vec<u8>)> = std::fs::read_dir(dir)
        .unwrap()
        .flatten()
        .filter(|e| e.path().is_file())
        .map(|e| {
            (
                e.file_name().to_string_lossy().into_owned(),
                std::fs::read(e.path()).unwrap(),
            )
        })
        .collect();
    v.sort();
    v
}

/// Copy a directory tree.
pub fn copy_dir(from: &Path, to: &Path, skip: &dyn Fn(&Path) -> bool) {
    std::fs::create_dir_all(to).unwrap();
    for e in std::fs::read_dir(from).unwrap().flatten() {
        let p = e.path();
        if skip(&p) {
            continue;
        }
        let t = to.join(e.file_name());
        if p.is_dir() {
            copy_dir(&p, &t, skip);
        } else {
            std::fs::copy(&p, &t).unwrap();
        }
    }
}

pub fn md5_hex(b: &[u8]) -> String {
    use md5::{Digest, Md5};
    Md5::digest(b).iter().map(|x| format!("{x:02x}")).collect()
}
