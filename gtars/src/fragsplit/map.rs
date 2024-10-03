use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};

pub struct BarcodeToClusterMap {
    map: HashMap<String, char>,
    cluster_labels: HashSet<char>,
}

pub trait ClusterLookup {
    fn get_cluster_from_barcode(&self, barcode: &str) -> Option<char>;
}

pub trait ClusterCount {
    fn n_clusters(&self) -> u16;
}

impl ClusterLookup for BarcodeToClusterMap {
    fn get_cluster_from_barcode(&self, barcode: &str) -> Option<char> {
        self.map.get(barcode).copied()
    }
}

impl ClusterCount for BarcodeToClusterMap {
    fn n_clusters(&self) -> u16 {
        self.cluster_labels.len() as u16
    }
}

impl BarcodeToClusterMap {
    pub fn from_file(file: &Path) -> Result<Self> {
        let file = File::open(file).with_context(|| format!("Couldn't open file: {:?}", file))?;

        let mut map: HashMap<String, char> = HashMap::new();
        let mut cluster_labels: HashSet<char> = HashSet::new();

        let reader = BufReader::new(file);

        for (index, line) in reader.lines().enumerate() {
            let line =
                line.with_context(|| format!("There was an error reading line {}", index + 1))?;

            let mut parts = line.split('\t');

            let barcode = parts.next();
            let cluster_id = parts.next();

            if barcode.is_none() || cluster_id.is_none() {
                anyhow::bail!(
                    "Invalid line format: Expected two tab-separated values, found: {:?}",
                    line
                );
            }

            if let (Some(barcode), Some(cluster_id)) = (barcode, cluster_id) {
                if cluster_id.len() > 1 {
                    anyhow::bail!(
                        "Invalid cluster id: Must be coercible to a char type. Found: {:?}",
                        cluster_id
                    );
                }
                let cluster_id = cluster_id.chars().next().unwrap();
                map.insert(barcode.to_string(), cluster_id);
                if !cluster_labels.contains(&cluster_id) {
                    cluster_labels.insert(cluster_id);
                }
            } else {
                anyhow::bail!(
                    "There was an error parsing the cluster map file for the following line: {:?}",
                    line
                )
            }
        }

        Ok(BarcodeToClusterMap {
            map,
            cluster_labels,
        })
    }

    pub fn get_cluster_labels(&self) -> HashSet<char> {
        self.cluster_labels.clone()
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn barcode_cluster_map_file() -> &'static str {
        "tests/data/barcode_cluster_map.tsv"
    }

    #[rstest]
    fn make_map_from_file(barcode_cluster_map_file: &str) {
        let path = Path::new(barcode_cluster_map_file);
        let mapping = BarcodeToClusterMap::from_file(path);

        assert_eq!(mapping.is_ok(), true);
        assert_eq!(mapping.unwrap().get_cluster_labels().len(), 3);
    }

    #[rstest]
    fn test_get_cluster_label(barcode_cluster_map_file: &str) {
        let path = Path::new(barcode_cluster_map_file);
        let mapping = BarcodeToClusterMap::from_file(path).unwrap();

        let cluster_id_none = mapping.get_cluster_from_barcode("AAACGCAAGCAAAGGATCGGCT");
        let cluster_id_some = mapping.get_cluster_from_barcode("AAACGCAAGCAACTGCGTCTTT");

        assert_eq!(cluster_id_none.is_none(), true);
        assert_eq!(cluster_id_some.is_some(), true);
    }
}
