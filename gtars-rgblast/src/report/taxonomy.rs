//! The optional NCBI taxdump (`nodes.dmp`, `names.dmp`, `merged.dmp`) for lineage, LCA calls
//! and rank folding. FHR has only a taxid and a species name, so lineage comes from here.
//! Nothing is downloaded.

use crate::error::Error;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

pub const TAXDUMP_HINT: &str = "for LCA calls: curl -O https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && tar xzf taxdump.tar.gz nodes.dmp names.dmp merged.dmp, then --taxdump .";

/// `--taxdump DIR`, else `$GTARS_RGBLAST_TAXDUMP`.
pub fn taxdump_dir(explicit: Option<&Path>) -> Option<PathBuf> {
    explicit.map(Path::to_path_buf).or_else(|| {
        std::env::var_os("GTARS_RGBLAST_TAXDUMP")
            .filter(|v| !v.is_empty())
            .map(PathBuf::from)
    })
}

/// One `\t|\t`-separated dump line, the trailing `\t|` removed.
fn fields(line: &str) -> Vec<&str> {
    line.trim_end_matches(['\n', '\r'])
        .trim_end_matches("\t|")
        .split("\t|\t")
        .collect()
}

pub struct Taxonomy {
    dir: PathBuf,
    /// Indexed by taxid; 0 = absent (root's parent is itself).
    parent: Vec<u32>,
    rank: Vec<u8>,
    ranks: Vec<String>,
    merged: HashMap<u32, u32>,
    names: HashMap<u32, String>,
    nodes_date: String,
    /// Unknown taxids met so far (each warned once, then a leaf under root).
    pub warnings: Vec<String>,
}

impl Taxonomy {
    /// Load `nodes.dmp` and (if present) `merged.dmp` from `dir`.
    pub fn load(dir: &Path) -> Result<Self, Error> {
        let nodes = dir.join("nodes.dmp");
        let f = std::fs::File::open(&nodes).map_err(|e| Error::io(&nodes, e))?;
        let mut t = Taxonomy {
            dir: dir.to_path_buf(),
            parent: Vec::new(),
            rank: Vec::new(),
            ranks: Vec::new(),
            merged: HashMap::new(),
            names: HashMap::new(),
            nodes_date: modified_date(&nodes),
            warnings: Vec::new(),
        };
        let mut rank_ix: HashMap<String, u8> = HashMap::new();
        for (n, line) in BufReader::new(f).lines().enumerate() {
            let line = line.map_err(|e| Error::io(&nodes, e))?;
            let fs = fields(&line);
            if fs.len() < 3 {
                continue;
            }
            let bad = || {
                Error::InvalidOptions(format!(
                    "{}: line {}: not a nodes.dmp line",
                    nodes.display(),
                    n + 1
                ))
            };
            let id: u32 = fs[0].trim().parse().map_err(|_| bad())?;
            let parent: u32 = fs[1].trim().parse().map_err(|_| bad())?;
            let r = fs[2].trim().to_string();
            let next = rank_ix.len() as u8;
            let ri = *rank_ix.entry(r.clone()).or_insert_with(|| {
                t.ranks.push(r);
                next
            });
            t.set(id, parent, ri);
        }
        let merged = dir.join("merged.dmp");
        if let Ok(f) = std::fs::File::open(&merged) {
            for line in BufReader::new(f).lines().map_while(Result::ok) {
                let fs = fields(&line);
                if let (Some(a), Some(b)) = (fs.first(), fs.get(1))
                    && let (Ok(a), Ok(b)) = (a.trim().parse(), b.trim().parse())
                {
                    t.merged.insert(a, b);
                }
            }
        }
        if !t.has(1) {
            let ri = t.rank_index("no rank");
            t.set(1, 1, ri);
        }
        Ok(t)
    }

    fn rank_index(&mut self, r: &str) -> u8 {
        if let Some(i) = self.ranks.iter().position(|x| x == r) {
            return i as u8;
        }
        self.ranks.push(r.to_string());
        (self.ranks.len() - 1) as u8
    }

    fn set(&mut self, id: u32, parent: u32, rank: u8) {
        let i = id as usize;
        if self.parent.len() <= i {
            self.parent.resize(i + 1, 0);
            self.rank.resize(i + 1, 0);
        }
        self.parent[i] = parent;
        self.rank[i] = rank;
    }

    pub fn has(&self, id: u32) -> bool {
        self.parent.get(id as usize).is_some_and(|p| *p != 0)
    }

    /// The current taxid for `id`: itself, or its `merged.dmp` replacement. An id that is still
    /// unknown is warned about once and becomes a leaf under root.
    pub fn resolve(&mut self, id: u32) -> u32 {
        if self.has(id) {
            return id;
        }
        if let Some(&m) = self.merged.get(&id)
            && self.has(m)
        {
            return m;
        }
        self.warnings.push(format!(
            "taxid {id} is not in {}; treated as a leaf under root",
            self.dir.join("nodes.dmp").display()
        ));
        let ri = self.rank_index("no rank");
        self.set(id, 1, ri);
        id
    }

    pub fn parent(&self, id: u32) -> u32 {
        self.parent.get(id as usize).copied().unwrap_or(1)
    }

    /// The rank as `nodes.dmp` gives it.
    pub fn rank(&self, id: u32) -> &str {
        self.rank
            .get(id as usize)
            .and_then(|r| self.ranks.get(*r as usize))
            .map(|s| s.as_str())
            .unwrap_or("no rank")
    }

    /// `id`, its parent, ..., root (1).
    pub fn lineage(&self, id: u32) -> Vec<u32> {
        let mut v = vec![id];
        let mut t = id;
        while t != 1 && self.has(t) && v.len() < 1000 {
            t = self.parent(t);
            v.push(t);
        }
        if *v.last().unwrap() != 1 {
            v.push(1);
        }
        v
    }

    /// The deepest node common to all lineages (root when there is nothing deeper).
    pub fn lca(&self, ids: &[u32]) -> u32 {
        let Some((&first, rest)) = ids.split_first() else {
            return 1;
        };
        let mut common = self.lineage(first);
        for &id in rest {
            let set: HashSet<u32> = self.lineage(id).into_iter().collect();
            common.retain(|t| set.contains(t));
        }
        common.first().copied().unwrap_or(1)
    }

    /// The ancestor of `id` (or itself) at `rank`.
    pub fn ancestor_at(&self, id: u32, rank: &str) -> Option<u32> {
        self.lineage(id).into_iter().find(|t| self.rank(*t) == rank)
    }

    /// Stream `names.dmp` once, keeping the scientific names of `wanted`.
    pub fn load_names(&mut self, wanted: &HashSet<u32>) -> Result<(), Error> {
        let missing: HashSet<u32> = wanted
            .iter()
            .filter(|t| !self.names.contains_key(t))
            .copied()
            .collect();
        if missing.is_empty() {
            return Ok(());
        }
        let p = self.dir.join("names.dmp");
        let f = std::fs::File::open(&p).map_err(|e| Error::io(&p, e))?;
        for line in BufReader::new(f).lines() {
            let line = line.map_err(|e| Error::io(&p, e))?;
            let fs = fields(&line);
            if fs.len() >= 4
                && fs[3].trim() == "scientific name"
                && let Ok(id) = fs[0].trim().parse::<u32>()
                && missing.contains(&id)
            {
                self.names.insert(id, fs[1].to_string());
            }
        }
        Ok(())
    }

    pub fn name(&self, id: u32) -> Option<&str> {
        self.names.get(&id).map(|s| s.as_str())
    }

    /// The number of names loaded (for tests: only the output's taxids).
    pub fn names_loaded(&self) -> usize {
        self.names.len()
    }

    /// "taxdump /data/taxdump (nodes.dmp 2026-11-01)"
    pub fn describe(&self) -> String {
        format!(
            "taxdump {} (nodes.dmp {})",
            self.dir.display(),
            self.nodes_date
        )
    }
}

fn modified_date(p: &Path) -> String {
    let Ok(t) = std::fs::metadata(p).and_then(|m| m.modified()) else {
        return "?".into();
    };
    let secs = t
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    crate::blastdb::BlastDate::from_epoch(secs)
        .iso()
        .chars()
        .take(10)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    pub fn mini(dir: &Path) {
        std::fs::write(
            dir.join("nodes.dmp"),
            "1\t|\t1\t|\tno rank\t|\n2759\t|\t1\t|\tsuperkingdom\t|\n314146\t|\t2759\t|\tsuperorder\t|\n9606\t|\t314146\t|\tspecies\t|\n10090\t|\t314146\t|\tspecies\t|\n10091\t|\t10090\t|\tsubspecies\t|\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("names.dmp"),
            "1\t|\troot\t|\t\t|\tscientific name\t|\n314146\t|\tEuarchontoglires\t|\t\t|\tscientific name\t|\n9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|\n9606\t|\thuman\t|\t\t|\tgenbank common name\t|\n10090\t|\tMus musculus\t|\t\t|\tscientific name\t|\n2759\t|\tEukaryota\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();
        std::fs::write(dir.join("merged.dmp"), "63221\t|\t9606\t|\n").unwrap();
    }

    #[test]
    fn lineage_lca_merged_unknown_names() {
        let d = tempfile::tempdir().unwrap();
        mini(d.path());
        let mut t = Taxonomy::load(d.path()).unwrap();
        assert_eq!(t.lca(&[9606, 10090]), 314146);
        assert_eq!(t.rank(314146), "superorder");
        assert_eq!(t.lca(&[9606]), 9606);
        assert_eq!(t.ancestor_at(10091, "species"), Some(10090));
        assert_eq!(t.resolve(63221), 9606);
        assert!(t.warnings.is_empty());
        assert_eq!(t.resolve(777), 777);
        assert_eq!(t.warnings.len(), 1);
        assert_eq!(t.parent(777), 1);
        assert_eq!(t.lca(&[777, 9606]), 1);
        t.load_names(&[9606, 314146].into_iter().collect()).unwrap();
        assert_eq!(t.name(9606), Some("Homo sapiens"));
        assert_eq!(t.names_loaded(), 2);
        assert!(t.describe().starts_with("taxdump "));
    }
}
