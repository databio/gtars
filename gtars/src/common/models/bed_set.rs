use crate::common::models::region_set::RegionSet;
use crate::common::utils::read_bedset_file;
use anyhow::Result;
use md5::{Digest, Md5};
use std::fmt::Debug;
use std::path::{Path, PathBuf};
use std::io::Error;
// use std::io::Error;

#[derive(Clone, Debug)]
pub struct BedSet {
    pub region_sets: Vec<RegionSet>,
    pub path: Option<PathBuf>,
}

// pub struct BedSetIterator<'a> {
//     bed_set: &'a BedSet,
//     index: usize,
// }

impl TryFrom<&Path> for BedSet {
    type Error = anyhow::Error;

    fn try_from(value: &Path) -> Result<Self> {
        let region_sets_paths = read_bedset_file(value)?;
        let mut region_sets = Vec::new();
        for region_set_path in region_sets_paths {
            let region_set = RegionSet::try_from(region_set_path)?;
            region_sets.push(region_set);
        }
        Ok(BedSet {
            region_sets,
            path: Some(value.to_owned()),
        })
    }
}

impl TryFrom<&str> for BedSet {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        BedSet::try_from(Path::new(value))
    }
}

impl TryFrom<String> for BedSet {
    type Error = anyhow::Error;

    fn try_from(value: String) -> Result<Self> {
        BedSet::try_from(Path::new(&value))
    }
}

impl TryFrom<PathBuf> for BedSet {
    type Error = anyhow::Error;

    fn try_from(value: PathBuf) -> Result<Self> {
        BedSet::try_from(value.as_path())
    }
}

impl TryFrom<Vec<&Path>> for BedSet {
    type Error = anyhow::Error;

    fn try_from(value: Vec<&Path>) -> Result<Self> {
        let mut region_sets = Vec::new();
        for region_set_path in value {
            let region_set = RegionSet::try_from(region_set_path)?;
            region_sets.push(region_set);
        }
        Ok(BedSet {
            region_sets: region_sets,
            path: None,
        })
    }
}

impl TryFrom<Vec<&str>> for BedSet {
    type Error = anyhow::Error;

    fn try_from(value: Vec<&str>) -> Result<Self> {
        let mut region_sets = Vec::new();
        for region_set_path in value {
            let region_set = RegionSet::try_from(region_set_path)?;
            region_sets.push(region_set);
        }
        Ok(BedSet {
            region_sets: region_sets,
            path: None,
        })
    }
}

impl TryFrom<Vec<String>> for BedSet {
    type Error = anyhow::Error;

    fn try_from(value: Vec<String>) -> Result<Self> {
        let mut region_sets = Vec::new();
        for region_set_path in value {
            let region_set = RegionSet::try_from(region_set_path)?;
            region_sets.push(region_set);
        }
        Ok(BedSet {
            region_sets: region_sets,
            path: None,
        })
    }
}

impl TryFrom<Vec<PathBuf>> for BedSet {
    type Error = anyhow::Error;

    fn try_from(value: Vec<PathBuf>) -> Result<Self> {
        let mut region_sets = Vec::new();
        for region_set_path in value {
            let region_set = RegionSet::try_from(region_set_path)?;
            region_sets.push(region_set);
        }
        Ok(BedSet {
            region_sets: region_sets,
            path: None,
        })
    }
}

impl From<Vec<RegionSet>> for BedSet {
    fn from(region_sets: Vec<RegionSet>) -> Self {
        BedSet {
            region_sets: region_sets,
            path: None,
        }
    }
}

impl BedSet {
    pub fn add(&mut self, region_set: RegionSet) {
        self.region_sets.push(region_set);
    }

    pub fn is_empty(&self) -> bool {
        if self.region_sets.len() == 0 {
            return true;
        }
        false
    }

    pub fn len(&self) -> usize {
        self.region_sets.len()
    }

    // pub fn identifier(&self) -> String {
    //     let mut bedfile_ids = Vec::new();
    //     for rs in &self.region_sets {
    //         if let id = rs.identifier() {
    //             bedfile_ids.push(id);
    //         }
    //     }
    //     let mut hasher = Md5::new();
    //     let combined = bedfile_ids.join("");
    //     hasher.update(combined.as_bytes());

    //     let hash = hasher.finalize();
    //     let bedset_digest: String = format!("{:x}", hash);

    //     bedset_digest
    // }
}



#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::read_dir;
    use tempfile::NamedTempFile;
    use std::io::Write;
 
    fn get_test_bed_dir() ->  Result<PathBuf, Error> {
        let folder_path = std::env::current_dir()
            .unwrap()
            .join("tests/data/bedset");
        
        Ok(folder_path)
    }



    fn get_test_bed_paths() -> Vec<PathBuf> {
        let dir = get_test_bed_dir();
        let mut paths = vec![];
        for entry in read_dir(dir.unwrap()).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.extension().map(|ext| ext == "bed" || ext == "narrowPeak").unwrap_or(false) {
                paths.push(path);
            }
        }
        paths.sort(); // for consistent identifier
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

    #[test]
    fn test_open_from_file_path() {
        let paths = get_test_bed_paths(); // Returns Vec<PathBuf>
        let temp_file = write_temp_bedset_list(paths);
        
        // let bs = BedSet::try_from(file_path.unwrap().as_path());
        assert!(BedSet::try_from(temp_file.path()).is_ok());
    }

    #[test]
    fn test_try_from_pathbuf_vec() {
        let paths = get_test_bed_paths();
        let bedset = BedSet::try_from(paths.clone());
        assert!(bedset.is_ok());
        let bedset = bedset.unwrap();
        assert_eq!(bedset.len(), paths.len());
    }

    #[test]
    fn test_try_from_str_vec() {
        let paths = get_test_bed_paths();
        let path_strs: Vec<String> = paths.iter().map(|p| p.to_string_lossy().into()).collect();
        let bedset = BedSet::try_from(path_strs);
        assert!(bedset.is_ok());
    }

    #[test]
    fn test_try_from_path_vec() {
        let paths = get_test_bed_paths();
        let refs: Vec<&Path> = paths.iter().map(|p| p.as_path()).collect();
        let bedset = BedSet::try_from(refs);
        assert!(bedset.is_ok());
    }

    #[test]
    fn test_bedset_add_and_len() {
        let paths = get_test_bed_paths();
        let mut bedset = BedSet::try_from(paths[..1].to_vec()).unwrap();
        let len_before = bedset.len();

        let additional = RegionSet::try_from(paths[1].as_path()).unwrap();
        bedset.add(additional);
        assert_eq!(bedset.len(), len_before + 1);
    }

    #[test]
    fn test_bedset_is_empty() {
        let paths = get_test_bed_paths();
        let bedset = BedSet::try_from(paths).unwrap();
        assert!(!bedset.is_empty());
    }

    // #[test]
    // fn test_bedset_identifier_consistency() {
    //     let paths = get_test_bed_paths();
    //     let bedset1 = BedSet::try_from(paths.clone()).unwrap();
    //     let bedset2 = BedSet::try_from(paths).unwrap();
    //     assert_eq!(bedset1.identifier(), bedset2.identifier());
    // }

    #[test]
    fn test_bedset_from_vec_regionset() {
        let paths = get_test_bed_paths();
        let region_sets: Vec<RegionSet> = paths
            .iter()
            .map(|p| RegionSet::try_from(p.as_path()).unwrap())
            .collect();
        let bedset = BedSet::from(region_sets);
        assert_eq!(bedset.len(), paths.len());
    }
}
