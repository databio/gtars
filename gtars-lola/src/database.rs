//! RegionDB: LOLA database loading, wrapping IGD + annotations + original region sets.

use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use gtars_core::models::RegionSet;
use gtars_igd::igd::Igd;

use crate::errors::LolaError;

/// Collection-level annotation from collection.txt.
#[derive(Debug, Clone, Default)]
pub struct CollectionAnno {
    pub collector: String,
    pub date: String,
    pub source: String,
    pub description: String,
    pub collection_name: String,
}

/// Per-file annotation from index.txt.
#[derive(Debug, Clone, Default)]
pub struct RegionSetAnno {
    pub filename: String,
    pub cell_type: String,
    pub description: String,
    pub tissue: String,
    pub data_source: String,
    pub antibody: String,
    pub treatment: String,
    /// Which collection this file belongs to.
    pub collection: String,
}

/// A LOLA region database: IGD index + original region sets + annotations.
pub struct RegionDB {
    /// The overlap index (all region sets indexed together).
    pub igd: Igd,
    /// Original regions per file, in the same order as IGD file indices.
    pub region_sets: Vec<RegionSet>,
    /// Path to the LOLA database folder (None if built from API/memory).
    pub db_location: Option<PathBuf>,
    /// Collection-level annotations.
    pub collection_anno: Vec<CollectionAnno>,
    /// Per-file annotations (same length and order as region_sets / IGD file_info).
    pub region_anno: Vec<RegionSetAnno>,
}

impl RegionDB {
    /// Load a RegionDB from a LOLA-format folder structure.
    ///
    /// Expects:
    /// ```text
    /// db_path/
    /// ├── collection1/
    /// │   ├── collection.txt
    /// │   ├── index.txt
    /// │   └── regions/
    /// │       ├── file1.bed
    /// │       └── file2.bed
    /// ├── collection2/
    /// │   └── ...
    /// ```
    pub fn from_lola_folder(
        db_path: &Path,
        collections_filter: Option<&[&str]>,
        limit: Option<usize>,
    ) -> Result<Self, LolaError> {
        let mut all_region_sets: Vec<RegionSet> = Vec::new();
        let mut all_region_anno: Vec<RegionSetAnno> = Vec::new();
        let mut all_collection_anno: Vec<CollectionAnno> = Vec::new();
        let mut igd_inputs: Vec<(String, Vec<(String, i32, i32)>)> = Vec::new();

        // Discover collections (subdirectories with a regions/ folder)
        let mut collections: Vec<PathBuf> = Vec::new();
        for entry in fs::read_dir(db_path).map_err(|e| LolaError::Other(e.into()))? {
            let entry = entry.map_err(|e| LolaError::Other(e.into()))?;
            let path = entry.path();
            if path.is_dir() && path.join("regions").is_dir() {
                let name = path
                    .file_name()
                    .map(|n| n.to_string_lossy().into_owned())
                    .unwrap_or_default();

                // Apply collection filter
                if let Some(filter) = collections_filter {
                    if !filter.iter().any(|f| *f == name) {
                        continue;
                    }
                }
                collections.push(path);
            }
        }
        collections.sort();

        for coll_path in &collections {
            let coll_name = coll_path
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default();

            // Parse collection.txt
            let coll_anno = parse_collection_txt(&coll_path.join("collection.txt"), &coll_name);
            all_collection_anno.push(coll_anno);

            // Parse index.txt for annotations (keyed by filename)
            let index_annos = parse_index_txt(&coll_path.join("index.txt"), &coll_name);
            let anno_map: HashMap<String, RegionSetAnno> = index_annos
                .into_iter()
                .map(|a| (a.filename.clone(), a))
                .collect();

            // Discover all BED files in regions/ directory (matching R LOLA behavior)
            let regions_dir = coll_path.join("regions");
            let mut bed_files: Vec<String> = Vec::new();
            if let Ok(entries) = fs::read_dir(&regions_dir) {
                for entry in entries.flatten() {
                    let fname = entry.file_name().to_string_lossy().into_owned();
                    if entry.path().is_file() {
                        bed_files.push(fname);
                    }
                }
            }
            bed_files.sort();

            let mut files_loaded = 0;

            for fname in &bed_files {
                if let Some(max) = limit {
                    if files_loaded >= max {
                        break;
                    }
                }

                let bed_path = regions_dir.join(fname);

                match RegionSet::try_from(bed_path.as_path()) {
                    Ok(region_set) => {
                        // Build IGD input tuples
                        let regions: Vec<(String, i32, i32)> = region_set
                            .regions
                            .iter()
                            .map(|r| (r.chr.clone(), r.start as i32, r.end as i32))
                            .collect();

                        igd_inputs.push((fname.clone(), regions));
                        all_region_sets.push(region_set);

                        // Use annotation from index.txt if available, otherwise default
                        let anno = anno_map.get(fname).cloned().unwrap_or(RegionSetAnno {
                            filename: fname.clone(),
                            collection: coll_name.clone(),
                            ..Default::default()
                        });
                        all_region_anno.push(anno);
                        files_loaded += 1;
                    }
                    Err(_) => continue,
                }
            }
        }

        // Build IGD from all loaded region sets
        let igd = Igd::from_region_sets(igd_inputs);

        Ok(RegionDB {
            igd,
            region_sets: all_region_sets,
            db_location: Some(db_path.to_path_buf()),
            collection_anno: all_collection_anno,
            region_anno: all_region_anno,
        })
    }

    /// Build a RegionDB from pre-loaded IGD, region sets, and annotations.
    /// For BEDbase integration where data comes from an API.
    pub fn from_igd_with_regions(
        igd: Igd,
        region_sets: Vec<RegionSet>,
        region_anno: Vec<RegionSetAnno>,
    ) -> Self {
        RegionDB {
            igd,
            region_sets,
            db_location: None,
            collection_anno: Vec::new(),
            region_anno,
        }
    }

    /// Merge two RegionDBs.
    pub fn merge(a: RegionDB, b: RegionDB) -> Self {
        let mut all_sets: Vec<(String, Vec<(String, i32, i32)>)> = Vec::new();
        let mut region_sets = Vec::new();
        let mut region_anno = Vec::new();

        // Add all from a
        for (i, rs) in a.region_sets.into_iter().enumerate() {
            let filename = if i < a.region_anno.len() {
                a.region_anno[i].filename.clone()
            } else if i < a.igd.file_info.len() {
                a.igd.file_info[i].filename.clone()
            } else {
                format!("set_{}", i)
            };

            let regions: Vec<(String, i32, i32)> = rs
                .regions
                .iter()
                .map(|r| (r.chr.clone(), r.start as i32, r.end as i32))
                .collect();
            all_sets.push((filename, regions));
            region_sets.push(rs);
        }
        region_anno.extend(a.region_anno);

        // Add all from b
        for (i, rs) in b.region_sets.into_iter().enumerate() {
            let filename = if i < b.region_anno.len() {
                b.region_anno[i].filename.clone()
            } else if i < b.igd.file_info.len() {
                b.igd.file_info[i].filename.clone()
            } else {
                format!("set_{}", i)
            };

            let regions: Vec<(String, i32, i32)> = rs
                .regions
                .iter()
                .map(|r| (r.chr.clone(), r.start as i32, r.end as i32))
                .collect();
            all_sets.push((filename, regions));
            region_sets.push(rs);
        }
        region_anno.extend(b.region_anno);

        let igd = Igd::from_region_sets(all_sets);

        let mut collection_anno = a.collection_anno;
        collection_anno.extend(b.collection_anno);

        RegionDB {
            igd,
            region_sets,
            db_location: None,
            collection_anno,
            region_anno,
        }
    }

    /// List all region set filenames, optionally filtered by collection.
    pub fn list_region_sets(&self, collections: Option<&[&str]>) -> Vec<String> {
        self.region_anno
            .iter()
            .filter(|a| {
                if let Some(filter) = collections {
                    filter.iter().any(|f| *f == a.collection)
                } else {
                    true
                }
            })
            .map(|a| a.filename.clone())
            .collect()
    }

    /// Get region sets by filename, optionally filtered by collection.
    pub fn get_region_set(
        &self,
        filenames: &[&str],
        collections: Option<&[&str]>,
    ) -> Vec<&RegionSet> {
        self.region_anno
            .iter()
            .enumerate()
            .filter(|(_, a)| {
                let name_match = filenames.iter().any(|f| *f == a.filename);
                let coll_match = if let Some(filter) = collections {
                    filter.iter().any(|f| *f == a.collection)
                } else {
                    true
                };
                name_match && coll_match
            })
            .filter_map(|(i, _)| self.region_sets.get(i))
            .collect()
    }

    /// Number of region sets in this database.
    pub fn num_region_sets(&self) -> usize {
        self.region_sets.len()
    }
}

// ---------------------------------------------------------------------------
// TSV/metadata parsers
// ---------------------------------------------------------------------------

/// Parse collection.txt (TSV with header).
fn parse_collection_txt(path: &Path, collection_name: &str) -> CollectionAnno {
    if !path.exists() {
        return CollectionAnno {
            collection_name: collection_name.to_string(),
            ..Default::default()
        };
    }

    let file = match fs::File::open(path) {
        Ok(f) => f,
        Err(_) => {
            return CollectionAnno {
                collection_name: collection_name.to_string(),
                ..Default::default()
            }
        }
    };

    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Read header to determine column order
    let header = match lines.next() {
        Some(Ok(h)) => h,
        _ => {
            return CollectionAnno {
                collection_name: collection_name.to_string(),
                ..Default::default()
            }
        }
    };

    let col_names: Vec<&str> = header.split('\t').collect();
    let col_map: HashMap<&str, usize> = col_names
        .iter()
        .enumerate()
        .map(|(i, &name)| (name.trim(), i))
        .collect();

    // Read first data line
    let data_line = match lines.next() {
        Some(Ok(l)) => l,
        _ => {
            return CollectionAnno {
                collection_name: collection_name.to_string(),
                ..Default::default()
            }
        }
    };

    let fields: Vec<&str> = data_line.split('\t').collect();

    let get = |key: &str| -> String {
        col_map
            .get(key)
            .and_then(|&i| fields.get(i))
            .map(|s| s.trim().to_string())
            .unwrap_or_default()
    };

    CollectionAnno {
        collector: get("collector"),
        date: get("date"),
        source: get("source"),
        description: get("description"),
        collection_name: collection_name.to_string(),
    }
}

/// Parse index.txt (TSV with header), returning per-file annotations.
fn parse_index_txt(path: &Path, collection_name: &str) -> Vec<RegionSetAnno> {
    if !path.exists() {
        return Vec::new();
    }

    let file = match fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return Vec::new(),
    };

    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Read header
    let header = match lines.next() {
        Some(Ok(h)) => h,
        _ => return Vec::new(),
    };

    let col_names: Vec<&str> = header.split('\t').collect();
    let col_map: HashMap<&str, usize> = col_names
        .iter()
        .enumerate()
        .map(|(i, &name)| (name.trim(), i))
        .collect();

    let mut annos = Vec::new();

    for line in lines {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };

        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        let get = |key: &str| -> String {
            col_map
                .get(key)
                .and_then(|&i| fields.get(i))
                .map(|s| s.trim().to_string())
                .unwrap_or_default()
        };

        annos.push(RegionSetAnno {
            filename: get("filename"),
            cell_type: get("cellType"),
            description: get("description"),
            tissue: get("tissue"),
            data_source: get("dataSource"),
            antibody: get("antibody"),
            treatment: get("treatment"),
            collection: collection_name.to_string(),
        });
    }

    annos
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn write_file(path: &Path, content: &str) {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).unwrap();
        }
        fs::write(path, content).unwrap();
    }

    fn create_test_lola_db(tmpdir: &Path) {
        // Create collection "coll1" with 2 BED files
        let coll1 = tmpdir.join("coll1");
        let regions = coll1.join("regions");
        fs::create_dir_all(&regions).unwrap();

        write_file(
            &coll1.join("collection.txt"),
            "collector\tdate\tsource\tdescription\nJohn\t2024\tENCODE\tTest collection\n",
        );

        write_file(
            &coll1.join("index.txt"),
            "filename\tcellType\tdescription\ttissue\tdataSource\tantibody\ttreatment\n\
             file1.bed\tK562\tHistone marks\tblood\tENCODE\tH3K27ac\tNone\n\
             file2.bed\tHeLa\tTF binding\tcervix\tENCODE\tCTCF\tNone\n",
        );

        write_file(
            &regions.join("file1.bed"),
            "chr1\t100\t200\nchr1\t300\t400\nchr1\t500\t600\n",
        );
        write_file(
            &regions.join("file2.bed"),
            "chr1\t150\t250\nchr1\t350\t450\n",
        );
    }

    #[test]
    fn test_load_from_lola_folder() {
        let tmpdir = tempfile::tempdir().unwrap();
        create_test_lola_db(tmpdir.path());

        let db = RegionDB::from_lola_folder(tmpdir.path(), None, None).unwrap();

        assert_eq!(db.num_region_sets(), 2);
        assert_eq!(db.igd.num_files(), 2);
        assert_eq!(db.collection_anno.len(), 1);
        assert_eq!(db.collection_anno[0].collector, "John");
        assert_eq!(db.region_anno.len(), 2);
        assert_eq!(db.region_anno[0].cell_type, "K562");
        assert_eq!(db.region_anno[1].cell_type, "HeLa");

        // Region sets should have correct counts
        assert_eq!(db.region_sets[0].regions.len(), 3); // file1: 3 regions
        assert_eq!(db.region_sets[1].regions.len(), 2); // file2: 2 regions
    }

    #[test]
    fn test_load_with_collection_filter() {
        let tmpdir = tempfile::tempdir().unwrap();
        create_test_lola_db(tmpdir.path());

        // Also create coll2
        let coll2 = tmpdir.path().join("coll2");
        let regions2 = coll2.join("regions");
        fs::create_dir_all(&regions2).unwrap();
        write_file(
            &coll2.join("index.txt"),
            "filename\tcellType\n\
             file3.bed\tGM12878\n",
        );
        write_file(&regions2.join("file3.bed"), "chr1\t100\t200\n");

        // Filter to only coll1
        let db =
            RegionDB::from_lola_folder(tmpdir.path(), Some(&["coll1"]), None).unwrap();
        assert_eq!(db.num_region_sets(), 2); // only coll1's files
    }

    #[test]
    fn test_load_with_limit() {
        let tmpdir = tempfile::tempdir().unwrap();
        create_test_lola_db(tmpdir.path());

        let db = RegionDB::from_lola_folder(tmpdir.path(), None, Some(1)).unwrap();
        assert_eq!(db.num_region_sets(), 1); // limited to 1 file per collection
    }

    #[test]
    fn test_list_region_sets() {
        let tmpdir = tempfile::tempdir().unwrap();
        create_test_lola_db(tmpdir.path());

        let db = RegionDB::from_lola_folder(tmpdir.path(), None, None).unwrap();
        let all = db.list_region_sets(None);
        assert_eq!(all.len(), 2);

        let filtered = db.list_region_sets(Some(&["coll1"]));
        assert_eq!(filtered.len(), 2);

        let empty = db.list_region_sets(Some(&["nonexistent"]));
        assert!(empty.is_empty());
    }

    #[test]
    fn test_get_region_set() {
        let tmpdir = tempfile::tempdir().unwrap();
        create_test_lola_db(tmpdir.path());

        let db = RegionDB::from_lola_folder(tmpdir.path(), None, None).unwrap();
        let sets = db.get_region_set(&["file1.bed"], None);
        assert_eq!(sets.len(), 1);
        assert_eq!(sets[0].regions.len(), 3);
    }

    #[test]
    fn test_merge_region_dbs() {
        let tmpdir = tempfile::tempdir().unwrap();
        create_test_lola_db(tmpdir.path());

        let db1 = RegionDB::from_lola_folder(tmpdir.path(), None, Some(1)).unwrap();
        let db2 = RegionDB::from_lola_folder(tmpdir.path(), None, Some(1)).unwrap();

        let merged = RegionDB::merge(db1, db2);
        assert_eq!(merged.num_region_sets(), 2);
        assert_eq!(merged.igd.num_files(), 2);
    }

    #[test]
    fn test_from_igd_with_regions() {
        let igd = Igd::from_region_sets(vec![(
            "test.bed".to_string(),
            vec![("chr1".to_string(), 100, 200)],
        )]);

        let region_set = RegionSet::from(vec![gtars_core::models::Region {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
            rest: None,
        }]);

        let anno = RegionSetAnno {
            filename: "test.bed".to_string(),
            cell_type: "K562".to_string(),
            ..Default::default()
        };

        let db = RegionDB::from_igd_with_regions(igd, vec![region_set], vec![anno]);
        assert_eq!(db.num_region_sets(), 1);
        assert!(db.db_location.is_none());
    }

    #[test]
    fn test_query_loaded_db() {
        // End-to-end: load DB, query with run_lola
        let tmpdir = tempfile::tempdir().unwrap();
        create_test_lola_db(tmpdir.path());

        let db = RegionDB::from_lola_folder(tmpdir.path(), None, None).unwrap();

        let user = RegionSet::from(vec![gtars_core::models::Region {
            chr: "chr1".to_string(),
            start: 150,
            end: 180,
            rest: None,
        }]);

        let universe = RegionSet::from(vec![
            gtars_core::models::Region {
                chr: "chr1".to_string(),
                start: 50,
                end: 650,
                rest: None,
            },
            gtars_core::models::Region {
                chr: "chr1".to_string(),
                start: 700,
                end: 800,
                rest: None,
            },
        ]);

        let config = crate::models::LolaConfig::default();
        let results =
            crate::enrichment::run_lola(&db.igd, &[user], &universe, &config).unwrap();

        assert_eq!(results.len(), 2); // 2 DB sets
        // User [150,180) overlaps file1's [100,200) and file2's [150,250)
        let r0 = results.iter().find(|r| r.db_set == 0).unwrap();
        let r1 = results.iter().find(|r| r.db_set == 1).unwrap();
        assert_eq!(r0.support, 1);
        assert_eq!(r1.support, 1);
    }
}
