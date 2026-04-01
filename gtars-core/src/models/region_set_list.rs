use crate::models::region::Region;
use crate::models::region_set::RegionSet;
use crate::utils::read_bedset_file;

use anyhow::Result;
use md5::{Digest, Md5};
use std::fmt::Debug;
use std::path::{Path, PathBuf};

/// A collection of named region sets — the gtars equivalent of GRangesList.
///
/// `RegionSetList` holds multiple [`RegionSet`]s with optional names.
/// It supports iteration, indexing, and flattening all contained regions
/// into a single `RegionSet` via [`concat`](RegionSetList::concat).
#[derive(Clone, Debug)]
pub struct RegionSetList {
    pub region_sets: Vec<RegionSet>,
    pub names: Option<Vec<String>>,
    pub path: Option<PathBuf>,
}

pub struct RegionSetListIterator<'a> {
    list: &'a RegionSetList,
    index: usize,
}

impl RegionSetList {
    /// Create a new `RegionSetList` from a vector of region sets.
    pub fn new(sets: Vec<RegionSet>) -> Self {
        RegionSetList {
            region_sets: sets,
            names: None,
            path: None,
        }
    }

    /// Create a new `RegionSetList` with names for each region set.
    pub fn with_names(sets: Vec<RegionSet>, names: Vec<String>) -> Self {
        RegionSetList {
            region_sets: sets,
            names: Some(names),
            path: None,
        }
    }

    pub fn add(&mut self, region_set: RegionSet) {
        self.region_sets.push(region_set);
    }

    pub fn is_empty(&self) -> bool {
        self.region_sets.is_empty()
    }

    pub fn len(&self) -> usize {
        self.region_sets.len()
    }

    /// Get a reference to the region set at the given index.
    pub fn get(&self, index: usize) -> Option<&RegionSet> {
        self.region_sets.get(index)
    }

    /// Iterate over the contained region sets.
    pub fn iter(&self) -> impl Iterator<Item = &RegionSet> {
        self.region_sets.iter()
    }

    /// Flatten all region sets into a single `RegionSet`.
    ///
    /// Regions are concatenated without merging or deduplication.
    /// Apply `disjoin()` or `reduce()` on the result as needed.
    pub fn concat(&self) -> RegionSet {
        let regions: Vec<Region> = self
            .region_sets
            .iter()
            .flat_map(|rs| rs.regions.iter().cloned())
            .collect();
        RegionSet {
            regions,
            header: None,
            path: None,
        }
    }

    pub fn identifier(&self) -> String {
        let mut bedfile_ids = Vec::new();
        for rs in &self.region_sets {
            let id = rs.identifier();
            bedfile_ids.push(id);
        }
        bedfile_ids.sort();
        let mut hasher = Md5::new();
        let combined = bedfile_ids.join("");
        hasher.update(combined.as_bytes());

        let hash = hasher.finalize();
        format!("{:x}", hash)
    }
}

impl TryFrom<&Path> for RegionSetList {
    type Error = anyhow::Error;

    fn try_from(value: &Path) -> Result<Self> {
        let region_sets_paths = read_bedset_file(value)?;
        let mut region_sets = Vec::new();
        for region_set_path in region_sets_paths {
            let region_set = RegionSet::try_from(region_set_path)?;
            region_sets.push(region_set);
        }
        Ok(RegionSetList {
            region_sets,
            names: None,
            path: Some(value.to_owned()),
        })
    }
}

impl TryFrom<&str> for RegionSetList {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        RegionSetList::try_from(Path::new(value))
    }
}

impl TryFrom<String> for RegionSetList {
    type Error = anyhow::Error;

    fn try_from(value: String) -> Result<Self> {
        RegionSetList::try_from(Path::new(&value))
    }
}

impl TryFrom<PathBuf> for RegionSetList {
    type Error = anyhow::Error;

    fn try_from(value: PathBuf) -> Result<Self> {
        RegionSetList::try_from(value.as_path())
    }
}

impl TryFrom<Vec<&Path>> for RegionSetList {
    type Error = anyhow::Error;

    fn try_from(value: Vec<&Path>) -> Result<Self> {
        let mut region_sets = Vec::new();
        for region_set_path in value {
            let region_set = RegionSet::try_from(region_set_path)?;
            region_sets.push(region_set);
        }
        Ok(RegionSetList {
            region_sets,
            names: None,
            path: None,
        })
    }
}

impl TryFrom<Vec<&str>> for RegionSetList {
    type Error = anyhow::Error;

    fn try_from(value: Vec<&str>) -> Result<Self> {
        let mut region_sets = Vec::new();
        for region_set_path in value {
            let region_set = RegionSet::try_from(region_set_path)?;
            region_sets.push(region_set);
        }
        Ok(RegionSetList {
            region_sets,
            names: None,
            path: None,
        })
    }
}

impl TryFrom<Vec<String>> for RegionSetList {
    type Error = anyhow::Error;

    fn try_from(value: Vec<String>) -> Result<Self> {
        let mut region_sets = Vec::new();
        for region_set_path in value {
            let region_set = RegionSet::try_from(region_set_path)?;
            region_sets.push(region_set);
        }
        Ok(RegionSetList {
            region_sets,
            names: None,
            path: None,
        })
    }
}

impl TryFrom<Vec<PathBuf>> for RegionSetList {
    type Error = anyhow::Error;

    fn try_from(value: Vec<PathBuf>) -> Result<Self> {
        let mut region_sets = Vec::new();
        for region_set_path in value {
            let region_set = RegionSet::try_from(region_set_path)?;
            region_sets.push(region_set);
        }
        Ok(RegionSetList {
            region_sets,
            names: None,
            path: None,
        })
    }
}

impl From<Vec<RegionSet>> for RegionSetList {
    fn from(region_sets: Vec<RegionSet>) -> Self {
        RegionSetList {
            region_sets,
            names: None,
            path: None,
        }
    }
}

impl<'a> Iterator for RegionSetListIterator<'a> {
    type Item = &'a RegionSet;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.list.region_sets.len() {
            let region_set = &self.list.region_sets[self.index];
            self.index += 1;
            Some(region_set)
        } else {
            None
        }
    }
}

impl<'a> IntoIterator for &'a RegionSetList {
    type Item = &'a RegionSet;
    type IntoIter = RegionSetListIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        RegionSetListIterator {
            list: self,
            index: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::read_dir;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn get_test_bed_dir() -> PathBuf {
        std::env::current_dir()
            .unwrap()
            .join("../tests/data/bedset")
    }

    fn get_test_bed_paths() -> Vec<PathBuf> {
        let dir = get_test_bed_dir();
        let mut paths = vec![];
        for entry in read_dir(dir).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            paths.push(path);
        }
        paths.sort();
        paths
    }

    fn write_temp_bedset_list(paths: Vec<PathBuf>) -> NamedTempFile {
        let mut temp_file = NamedTempFile::new().expect("Failed to create temporary file");

        for path in paths {
            let abs_path = std::fs::canonicalize(&path)
                .unwrap_or_else(|_| panic!("Failed to canonicalize path: {}", path.display()));
            writeln!(temp_file, "{}", abs_path.display())
                .unwrap_or_else(|_| panic!("Failed to write to temp file: {}", abs_path.display()));
        }

        temp_file
    }

    fn sample_region_set(chr: &str, ranges: &[(u32, u32)]) -> RegionSet {
        RegionSet {
            regions: ranges
                .iter()
                .map(|(s, e)| Region {
                    chr: chr.to_string(),
                    start: *s,
                    end: *e,
                    rest: None,
                })
                .collect(),
            header: None,
            path: None,
        }
    }

    #[test]
    fn test_open_from_file_path() {
        let paths = get_test_bed_paths();
        let temp_file = write_temp_bedset_list(paths);
        assert!(RegionSetList::try_from(temp_file.path()).is_ok());
    }

    #[test]
    fn test_try_from_pathbuf_vec() {
        let paths = get_test_bed_paths();
        assert!(RegionSetList::try_from(paths).is_ok());
    }

    #[test]
    fn test_try_from_str_vec() {
        let paths = get_test_bed_paths();
        let path_strs: Vec<String> = paths.iter().map(|p| p.to_string_lossy().into()).collect();
        assert!(RegionSetList::try_from(path_strs).is_ok());
    }

    #[test]
    fn test_try_from_path_vec() {
        let paths = get_test_bed_paths();
        let refs: Vec<&Path> = paths.iter().map(|p| p.as_path()).collect();
        assert!(RegionSetList::try_from(refs).is_ok());
    }

    #[test]
    fn test_add_and_len() {
        let paths = get_test_bed_paths();
        let mut rsl = RegionSetList::try_from(paths[..1].to_vec()).unwrap();
        let len_before = rsl.len();

        let additional = RegionSet::try_from(paths[1].as_path()).unwrap();
        rsl.add(additional);
        assert_eq!(rsl.len(), len_before + 1);
    }

    #[test]
    fn test_is_empty() {
        let paths = get_test_bed_paths();
        let rsl = RegionSetList::try_from(paths).unwrap();
        assert!(!rsl.is_empty());
    }

    #[test]
    fn test_from_vec_regionset() {
        let paths = get_test_bed_paths();
        let region_sets: Vec<RegionSet> = paths
            .iter()
            .map(|p| RegionSet::try_from(p.as_path()).unwrap())
            .collect();
        let rsl = RegionSetList::from(region_sets);
        assert_eq!(rsl.len(), paths.len());
    }

    #[test]
    fn test_calculate_identifier() {
        let paths = get_test_bed_paths();
        let temp_file = write_temp_bedset_list(paths);

        let rsl = RegionSetList::try_from(temp_file.path()).unwrap();
        assert_eq!("17a10ce63638431b34e7d044c3eac186", rsl.identifier());
    }

    // --- New tests for RegionSetList-specific functionality ---

    #[test]
    fn test_new_empty() {
        let rsl = RegionSetList::new(vec![]);
        assert_eq!(rsl.len(), 0);
        assert!(rsl.is_empty());
        assert!(rsl.names.is_none());
    }

    #[test]
    fn test_new_from_sets() {
        let rs1 = sample_region_set("chr1", &[(0, 100), (200, 300)]);
        let rs2 = sample_region_set("chr1", &[(50, 150)]);
        let rsl = RegionSetList::new(vec![rs1, rs2]);
        assert_eq!(rsl.len(), 2);
        assert!(!rsl.is_empty());
    }

    #[test]
    fn test_with_names() {
        let rs1 = sample_region_set("chr1", &[(0, 100)]);
        let rs2 = sample_region_set("chr2", &[(0, 50)]);
        let rsl = RegionSetList::with_names(
            vec![rs1, rs2],
            vec!["file1.bed".into(), "file2.bed".into()],
        );
        assert_eq!(rsl.names.as_ref().unwrap(), &["file1.bed", "file2.bed"]);
    }

    #[test]
    fn test_get() {
        let rs1 = sample_region_set("chr1", &[(0, 100)]);
        let rs2 = sample_region_set("chr2", &[(0, 50)]);
        let rsl = RegionSetList::new(vec![rs1, rs2]);

        let first = rsl.get(0).unwrap();
        assert_eq!(first.regions.len(), 1);
        assert_eq!(first.regions[0].chr, "chr1");

        assert!(rsl.get(5).is_none());
    }

    #[test]
    fn test_iter() {
        let rs1 = sample_region_set("chr1", &[(0, 100)]);
        let rs2 = sample_region_set("chr2", &[(0, 50)]);
        let rsl = RegionSetList::new(vec![rs1, rs2]);

        let total_regions: usize = rsl.iter().map(|rs| rs.regions.len()).sum();
        assert_eq!(total_regions, 2);
    }

    #[test]
    fn test_concat_flattens() {
        let rs1 = sample_region_set("chr1", &[(0, 100), (200, 300)]);
        let rs2 = sample_region_set("chr1", &[(50, 150)]);
        let rs3 = sample_region_set("chr2", &[(0, 500)]);
        let rsl = RegionSetList::new(vec![rs1, rs2, rs3]);

        let combined = rsl.concat();
        assert_eq!(combined.regions.len(), 4);
    }

    #[test]
    fn test_concat_empty() {
        let rsl = RegionSetList::new(vec![]);
        let combined = rsl.concat();
        assert_eq!(combined.regions.len(), 0);
    }

    #[test]
    fn test_concat_single() {
        let rs = sample_region_set("chr1", &[(0, 100), (200, 300)]);
        let rsl = RegionSetList::new(vec![rs]);
        let combined = rsl.concat();
        assert_eq!(combined.regions.len(), 2);
    }
}
