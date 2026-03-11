//! FAIR Headers Reference genome (FHR) metadata for sequence collections.
//!
//! This module contains the FHR data types, sidecar JSON I/O functions,
//! and RefgetStore bridge methods for managing FHR metadata.
//!
//! See: https://github.com/FAIR-bioHeaders/FHR-Specification

use super::*;
use super::readonly::ReadonlyRefgetStore;
use super::core::RefgetStore;

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::hashkeyable::{HashKeyable, key_to_digest_string};

// ============================================================================
// Types
// ============================================================================

/// FAIR Headers Reference genome (FHR) metadata for a sequence collection.
///
/// All fields are optional to allow partial metadata. RefgetStore does not
/// enforce FHR schema compliance -- that's the user's responsibility.
#[derive(Clone, Debug, Serialize, Deserialize, Default)]
#[serde(rename_all = "camelCase")]
pub struct FhrMetadata {
    /// URL to the FHR JSON schema
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub schema: Option<String>,

    /// FHR schema version (numeric per spec, e.g. 1 or 1.0)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub schema_version: Option<serde_json::Number>,

    /// Genome name (e.g., "Homo sapiens")
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome: Option<String>,

    /// Taxonomy information
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub taxon: Option<FhrTaxon>,

    /// Genome version (e.g., "GRCh38.p14")
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,

    /// Who created the metadata (ORCID URIs)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub metadata_author: Option<Vec<FhrAuthor>>,

    /// Who assembled the genome
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assembly_author: Option<Vec<FhrAuthor>>,

    /// Assembly creation date (ISO 8601)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub date_created: Option<String>,

    /// Description of the physical sample
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub voucher_specimen: Option<String>,

    /// Masking type
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub masking: Option<String>,

    /// File-level checksum (SHA2-512/256 per FHR spec)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub checksum: Option<String>,

    /// Alternative common names for this genome
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_synonym: Option<Vec<String>>,

    /// Database accession identifier (single object per spec)
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        rename = "accessionID"
    )]
    pub accession_id: Option<FhrIdentifier>,

    /// Sequencing instruments used
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub instrument: Option<Vec<String>>,

    /// DOI or scholarly article reference (single string per spec)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scholarly_article: Option<String>,

    /// Documentation about the genome
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub documentation: Option<String>,

    /// Identifiers of the genome (namespace:value format)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub identifier: Option<Vec<String>>,

    /// License information
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub license: Option<String>,

    /// Related URLs
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub related_link: Option<Vec<String>>,

    /// Funding information (single string per spec)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub funding: Option<String>,

    /// General statistics about the genome assembly
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vital_stats: Option<FhrVitalStats>,

    /// Seqcol digest (added by RefgetStore, not part of FHR 1.0)
    #[serde(skip)]
    pub seqcol_digest: Option<String>,

    /// Catch-all for any other FHR fields or custom extensions
    #[serde(flatten)]
    pub extra: HashMap<String, serde_json::Value>,
}

/// General statistics about a genome assembly.
#[derive(Clone, Debug, Serialize, Deserialize, Default)]
#[serde(rename_all = "camelCase")]
pub struct FhrVitalStats {
    #[serde(default, skip_serializing_if = "Option::is_none", rename = "L50")]
    pub l50: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none", rename = "N50")]
    pub n50: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none", rename = "L90")]
    pub l90: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub total_base_pairs: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub number_contigs: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub number_scaffolds: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub read_technology: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FhrTaxon {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uri: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FhrAuthor {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uri: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FhrIdentifier {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
}

// ============================================================================
// Disk I/O helpers -- called by RefgetStore, no dependency on it
// ============================================================================

pub(crate) const SIDECAR_EXTENSION: &str = ".fhr.json";

/// Load all FHR sidecar files from a collections directory.
///
/// Scans for `*.fhr.json` files, parses each one, and returns a map
/// keyed by collection digest. Malformed files are skipped with a warning to stderr.
pub fn load_sidecars(collections_dir: &Path) -> HashMap<[u8; 32], FhrMetadata> {
    let mut map = HashMap::new();
    if !collections_dir.exists() {
        return map;
    }
    let entries = match fs::read_dir(collections_dir) {
        Ok(e) => e,
        Err(e) => {
            eprintln!(
                "Warning: could not read FHR sidecar directory {}: {}",
                collections_dir.display(),
                e
            );
            return map;
        }
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if name.ends_with(SIDECAR_EXTENSION) {
                let digest_str = &name[..name.len() - SIDECAR_EXTENSION.len()];
                let key = digest_str.to_key();
                match fs::read_to_string(&path) {
                    Ok(json) => {
                        match serde_json::from_str::<FhrMetadata>(&json) {
                            Ok(fhr) => {
                                map.insert(key, fhr);
                            }
                            Err(e) => {
                                eprintln!(
                                    "Warning: skipping malformed FHR sidecar {}: {}",
                                    path.display(),
                                    e
                                );
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!(
                            "Warning: could not read FHR sidecar {}: {}",
                            path.display(),
                            e
                        );
                    }
                }
            }
        }
    }
    map
}

/// Write all FHR sidecar files to a collections directory.
pub fn write_sidecars(
    collections_dir: &Path,
    metadata: &HashMap<[u8; 32], FhrMetadata>,
) -> Result<()> {
    for (key, fhr) in metadata {
        let digest_str = key_to_digest_string(key);
        let path = collections_dir.join(format!("{}{}", digest_str, SIDECAR_EXTENSION));
        write_sidecar(&path, fhr)?;
    }
    Ok(())
}

/// Write a single FHR sidecar JSON file.
pub fn write_sidecar(path: &Path, metadata: &FhrMetadata) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let json = serde_json::to_string_pretty(metadata)?;
    fs::write(path, json)?;
    Ok(())
}

/// Remove a single FHR sidecar file (if it exists). Returns quietly if missing.
pub fn remove_sidecar(collections_dir: &Path, digest_str: &str) {
    let path = collections_dir.join(format!("{}{}", digest_str, SIDECAR_EXTENSION));
    let _ = fs::remove_file(path);
}

/// Build the sidecar file path for a given digest.
pub fn sidecar_path(collections_dir: &Path, digest_str: &str) -> std::path::PathBuf {
    collections_dir.join(format!("{}{}", digest_str, SIDECAR_EXTENSION))
}

/// Load FHR metadata from a standalone JSON file.
pub fn load_from_json(path: &str) -> Result<FhrMetadata> {
    let json = fs::read_to_string(path)
        .context(format!("Failed to read FHR metadata from {}", path))?;
    serde_json::from_str(&json).context("Failed to parse FHR JSON")
}

// ============================================================================
// ReadonlyRefgetStore FHR bridge methods
// ============================================================================

impl ReadonlyRefgetStore {
    /// Set FHR metadata for a collection.
    pub fn set_fhr_metadata(
        &mut self,
        collection_digest: &str,
        metadata: FhrMetadata,
    ) -> Result<()> {
        let key = collection_digest.to_key();
        if !self.collections.contains_key(&key) {
            return Err(anyhow::anyhow!("Collection not found: {}", collection_digest));
        }
        if self.persist_to_disk {
            if let Some(ref local_path) = self.local_path {
                let path = sidecar_path(
                    &local_path.join("collections"),
                    collection_digest,
                );
                write_sidecar(&path, &metadata)?;
            }
        }
        self.fhr_metadata.insert(key, metadata);
        Ok(())
    }

    /// Get FHR metadata for a collection. Returns None if missing.
    pub fn get_fhr_metadata(&self, collection_digest: &str) -> Option<&FhrMetadata> {
        let key = collection_digest.to_key();
        self.fhr_metadata.get(&key)
    }

    /// Remove FHR metadata for a collection.
    pub fn remove_fhr_metadata(&mut self, collection_digest: &str) -> bool {
        let key = collection_digest.to_key();
        if self.persist_to_disk {
            if let Some(ref local_path) = self.local_path {
                remove_sidecar(
                    &local_path.join("collections"),
                    collection_digest,
                );
            }
        }
        self.fhr_metadata.remove(&key).is_some()
    }

    /// List all collection digests that have FHR metadata.
    pub fn list_fhr_metadata(&self) -> Vec<String> {
        self.fhr_metadata
            .keys()
            .map(|key| key_to_digest_string(key))
            .collect()
    }

    /// Load FHR metadata from a JSON file and attach it to a collection.
    pub fn load_fhr_metadata(&mut self, collection_digest: &str, path: &str) -> Result<()> {
        let metadata = load_from_json(path)?;
        self.set_fhr_metadata(collection_digest, metadata)
    }
}

// ============================================================================
// RefgetStore FHR wrapper delegates
// ============================================================================

impl RefgetStore {
    /// Set FHR metadata for a collection.
    pub fn set_fhr_metadata(&mut self, collection_digest: &str, metadata: FhrMetadata) -> Result<()> {
        self.inner.set_fhr_metadata(collection_digest, metadata)
    }

    /// Remove FHR metadata for a collection.
    pub fn remove_fhr_metadata(&mut self, collection_digest: &str) -> bool {
        self.inner.remove_fhr_metadata(collection_digest)
    }

    /// Load FHR metadata from a JSON file.
    pub fn load_fhr_metadata(&mut self, collection_digest: &str, path: &str) -> Result<()> {
        self.inner.load_fhr_metadata(collection_digest, path)
    }

    /// Pull FHR metadata sidecars from the remote store.
    ///
    /// If `digest` is Some, pulls only that collection's FHR.
    /// If `digest` is None, pulls FHR for all known collections.
    pub fn pull_fhr(
        &mut self,
        digest: Option<&str>,
        strategy: SyncStrategy,
    ) -> Result<PullResult> {
        let mut result = PullResult::default();

        let digests: Vec<String> = match digest {
            Some(d) => vec![d.to_string()],
            None => self
                .inner
                .collections
                .values()
                .map(|r| r.metadata().digest.to_string())
                .collect(),
        };

        for digest_str in &digests {
            let relative_path = format!("collections/{}.fhr.json", digest_str);

            match strategy {
                SyncStrategy::KeepOurs => {
                    let was_local = self
                        .inner
                        .local_path
                        .as_ref()
                        .map(|p| p.join(&relative_path).exists())
                        .unwrap_or(false);
                    match ReadonlyRefgetStore::fetch_file(
                        &self.inner.local_path,
                        &self.inner.remote_source,
                        &relative_path,
                        self.inner.persist_to_disk,
                        false,
                    ) {
                        Ok(data) => {
                            if was_local {
                                result.skipped += 1;
                            } else {
                                if let Ok(fhr) = serde_json::from_slice::<FhrMetadata>(&data) {
                                    let key = digest_str.to_key();
                                    self.inner.fhr_metadata.insert(key, fhr);
                                }
                                result.pulled += 1;
                            }
                        }
                        Err(_) => {
                            result.not_found += 1;
                        }
                    }
                }
                SyncStrategy::KeepTheirs => {
                    match ReadonlyRefgetStore::fetch_file(
                        &self.inner.local_path,
                        &self.inner.remote_source,
                        &relative_path,
                        self.inner.persist_to_disk,
                        true,
                    ) {
                        Ok(data) => {
                            if let Ok(fhr) = serde_json::from_slice::<FhrMetadata>(&data) {
                                let key = digest_str.to_key();
                                self.inner.fhr_metadata.insert(key, fhr);
                            }
                            result.pulled += 1;
                        }
                        Err(_) => {
                            result.not_found += 1;
                        }
                    }
                }
                SyncStrategy::Notify => {
                    let local_exists = self
                        .inner
                        .local_path
                        .as_ref()
                        .map(|p| p.join(&relative_path).exists())
                        .unwrap_or(false);

                    if local_exists {
                        match ReadonlyRefgetStore::fetch_file(
                            &None,
                            &self.inner.remote_source,
                            &relative_path,
                            false,
                            false,
                        ) {
                            Ok(remote_data) => {
                                let local_path = self
                                    .inner
                                    .local_path
                                    .as_ref()
                                    .unwrap()
                                    .join(&relative_path);
                                let local_data = fs::read(&local_path)?;
                                if local_data != remote_data {
                                    result.conflicts.push(relative_path);
                                } else {
                                    result.skipped += 1;
                                }
                            }
                            Err(_) => {
                                result.not_found += 1;
                            }
                        }
                    } else {
                        match ReadonlyRefgetStore::fetch_file(
                            &None,
                            &self.inner.remote_source,
                            &relative_path,
                            false,
                            false,
                        ) {
                            Ok(_) => {
                                result.conflicts.push(relative_path);
                            }
                            Err(_) => {
                                result.not_found += 1;
                            }
                        }
                    }
                }
            }
        }

        Ok(result)
    }
}

// ============================================================================
// Tests -- serialization, roundtripping, disk I/O, and store integration
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_json_roundtrip() {
        let fhr = FhrMetadata {
            schema: Some("https://raw.githubusercontent.com/FAIR-bioHeaders/FHR-Specification/main/fhr.json".to_string()),
            schema_version: Some(serde_json::Number::from_f64(1.0).unwrap()),
            genome: Some("Homo sapiens".to_string()),
            taxon: Some(FhrTaxon {
                name: Some("Homo sapiens".to_string()),
                uri: Some("https://identifiers.org/taxonomy:9606".to_string()),
            }),
            version: Some("GRCh38.p14".to_string()),
            masking: Some("soft-masked".to_string()),
            genome_synonym: Some(vec!["hg38".to_string()]),
            scholarly_article: Some("10.1371/journal.pntd.0008755".to_string()),
            funding: Some("NIH R01".to_string()),
            accession_id: Some(FhrIdentifier {
                name: Some("GCA_000001405.29".to_string()),
                url: Some("https://www.ncbi.nlm.nih.gov/assembly/GCA_000001405.29".to_string()),
            }),
            ..Default::default()
        };

        let json = serde_json::to_string_pretty(&fhr).unwrap();
        let roundtripped: FhrMetadata = serde_json::from_str(&json).unwrap();
        assert_eq!(roundtripped.genome, fhr.genome);
        assert_eq!(roundtripped.taxon.as_ref().unwrap().name, fhr.taxon.as_ref().unwrap().name);
        assert_eq!(roundtripped.genome_synonym, fhr.genome_synonym);
        assert_eq!(roundtripped.scholarly_article, Some("10.1371/journal.pntd.0008755".to_string()));
        assert_eq!(roundtripped.funding, Some("NIH R01".to_string()));
        assert!(roundtripped.accession_id.is_some());
    }

    #[test]
    fn test_extra_fields_preserved() {
        let json = r#"{
            "genome": "Test",
            "customField": "custom_value",
            "anotherCustom": [1, 2, 3]
        }"#;
        let fhr: FhrMetadata = serde_json::from_str(json).unwrap();
        assert_eq!(fhr.genome, Some("Test".to_string()));
        assert!(fhr.extra.contains_key("customField"));

        let json_out = serde_json::to_string(&fhr).unwrap();
        assert!(json_out.contains("customField"));
        assert!(json_out.contains("custom_value"));
    }

    #[test]
    fn test_camel_case_serialization() {
        let fhr = FhrMetadata {
            schema_version: Some(serde_json::Number::from_f64(1.0).unwrap()),
            genome_synonym: Some(vec!["hg38".to_string()]),
            date_created: Some("2024-01-01".to_string()),
            ..Default::default()
        };
        let json = serde_json::to_string(&fhr).unwrap();
        assert!(json.contains("schemaVersion"));
        assert!(json.contains("genomeSynonym"));
        assert!(json.contains("dateCreated"));
        assert!(!json.contains("schema_version"));
        assert!(!json.contains("genome_synonym"));
    }

    #[test]
    fn test_default_is_empty() {
        let fhr = FhrMetadata::default();
        let json = serde_json::to_string(&fhr).unwrap();
        assert_eq!(json, "{}");
    }

    #[test]
    fn test_write_and_load_sidecar() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.fhr.json");

        let fhr = FhrMetadata {
            genome: Some("Test".to_string()),
            version: Some("1.0".to_string()),
            ..Default::default()
        };

        write_sidecar(&path, &fhr).unwrap();
        assert!(path.exists());

        let loaded = load_from_json(path.to_str().unwrap()).unwrap();
        assert_eq!(loaded.genome, Some("Test".to_string()));
        assert_eq!(loaded.version, Some("1.0".to_string()));
    }

    #[test]
    fn test_load_sidecars_empty_dir() {
        let dir = tempdir().unwrap();
        let map = load_sidecars(dir.path());
        assert!(map.is_empty());
    }

    #[test]
    fn test_load_sidecars_nonexistent_dir() {
        let map = load_sidecars(Path::new("/nonexistent/path"));
        assert!(map.is_empty());
    }

    #[test]
    fn test_remove_sidecar_missing_is_ok() {
        let dir = tempdir().unwrap();
        remove_sidecar(dir.path(), "nonexistent_digest");
    }

    #[test]
    fn test_accession_id_casing() {
        let fhr = FhrMetadata {
            accession_id: Some(FhrIdentifier {
                name: Some("GCA_000001405.29".to_string()),
                url: Some("https://ncbi.nlm.nih.gov".to_string()),
            }),
            ..Default::default()
        };
        let json = serde_json::to_string(&fhr).unwrap();
        assert!(json.contains("accessionID"));
        assert!(!json.contains("accessionId"));
    }

    #[test]
    fn test_schema_version_as_number() {
        let json = r#"{"schemaVersion": 1}"#;
        let fhr: FhrMetadata = serde_json::from_str(json).unwrap();
        assert!(fhr.schema_version.is_some());
        let ver = fhr.schema_version.unwrap();
        assert_eq!(ver.to_string(), "1");

        let json = r#"{"schemaVersion": 1.0}"#;
        let fhr: FhrMetadata = serde_json::from_str(json).unwrap();
        assert!(fhr.schema_version.is_some());
        let ver = fhr.schema_version.unwrap();
        assert_eq!(ver.to_string(), "1.0");
    }

    #[test]
    fn test_vital_stats_roundtrip() {
        let fhr = FhrMetadata {
            vital_stats: Some(FhrVitalStats {
                l50: Some(42),
                n50: Some(1_000_000),
                l90: Some(100),
                total_base_pairs: Some(3_000_000_000),
                number_contigs: Some(500),
                number_scaffolds: Some(24),
                read_technology: Some("hifi".to_string()),
            }),
            ..Default::default()
        };
        let json = serde_json::to_string_pretty(&fhr).unwrap();
        assert!(json.contains("\"L50\""));
        assert!(json.contains("\"N50\""));
        assert!(json.contains("\"L90\""));
        assert!(json.contains("\"totalBasePairs\""));
        assert!(json.contains("\"numberContigs\""));
        let roundtripped: FhrMetadata = serde_json::from_str(&json).unwrap();
        let stats = roundtripped.vital_stats.unwrap();
        assert_eq!(stats.l50, Some(42));
        assert_eq!(stats.n50, Some(1_000_000));
        assert_eq!(stats.read_technology, Some("hifi".to_string()));
    }

    #[test]
    fn test_spec_example_roundtrip() {
        let json = r#"{
            "schema":"https://raw.githubusercontent.com/FAIR-bioHeaders/FHR-Specification/main/fhr.jso",
            "schemaVersion": 1.0,
            "taxon": {"name":"Bombas huntii", "uri": "https://identifiers.org/taxonomy:9606"},
            "genome": "Bombas huntii",
            "genomeSynonym": ["B. huntii"],
            "version": "0.0.1",
            "metadataAuthor": [{"name":"Adam Wright", "uri":"https://orcid.org/0000-0002-5719-4024"}],
            "assemblyAuthor": [{"name":"David Molik", "url":"https://orcid.org/0000-0003-3192-6538"}],
            "dateCreated":"2022-03-21",
            "accessionID": {"name":"PBARC", "url":"https://www.ars.usda.gov/pacific-west-area/hilo-hi/daniel-k-inouye-us-pacific-basin-agricultural-research-center/"},
            "instrument": ["Sequel IIe", "Nanopore"],
            "voucherSpecimen":"Located in Freezer 33, Drawer 137",
            "scholarlyArticle":"10.1371/journal.pntd.0008755",
            "assemblySoftware":"HiFiASM",
            "funding":"funding",
            "reuseConditions":"public domain",
            "documentation":"Built assembly from... ",
            "masking":"soft-masked",
            "identifier": ["beetlebase:TC010103"],
            "relatedLink": ["http://wfleabase.org/genome/Daphnia_pulex/dpulex_jgi060905/fasta/"],
            "checksum":"md5:7582b26fcb0a9775b87c38f836e97c42"
        }"#;
        let fhr: FhrMetadata = serde_json::from_str(json).unwrap();
        assert_eq!(fhr.genome, Some("Bombas huntii".to_string()));
        assert_eq!(fhr.voucher_specimen, Some("Located in Freezer 33, Drawer 137".to_string()));
        assert_eq!(fhr.documentation, Some("Built assembly from... ".to_string()));
        assert_eq!(fhr.scholarly_article, Some("10.1371/journal.pntd.0008755".to_string()));
        assert_eq!(fhr.funding, Some("funding".to_string()));
        assert_eq!(fhr.identifier, Some(vec!["beetlebase:TC010103".to_string()]));
        assert!(fhr.accession_id.is_some());
        assert_eq!(fhr.accession_id.as_ref().unwrap().name, Some("PBARC".to_string()));
        assert!(fhr.extra.contains_key("assemblySoftware"));
        assert!(fhr.extra.contains_key("reuseConditions"));
    }

    #[test]
    fn test_seqcol_digest_skipped_in_json() {
        let mut fhr = FhrMetadata {
            genome: Some("Test".to_string()),
            ..Default::default()
        };
        fhr.seqcol_digest = Some("abc123".to_string());
        let json = serde_json::to_string(&fhr).unwrap();
        assert!(!json.contains("seqcolDigest"));
        assert!(!json.contains("seqcol_digest"));
        assert!(!json.contains("abc123"));
    }

    #[test]
    fn test_new_fields_present() {
        let fhr = FhrMetadata {
            voucher_specimen: Some("Freezer 33".to_string()),
            documentation: Some("Assembly notes".to_string()),
            identifier: Some(vec!["ncbi:GCA_000001405".to_string()]),
            ..Default::default()
        };
        let json = serde_json::to_string(&fhr).unwrap();
        assert!(json.contains("voucherSpecimen"));
        assert!(json.contains("documentation"));
        assert!(json.contains("identifier"));
    }

    #[test]
    fn test_load_sidecars_skips_malformed_json() {
        let dir = tempdir().unwrap();
        let bad_path = dir.path().join("baddigest.fhr.json");
        fs::write(&bad_path, "{ not valid json }").unwrap();
        let map = load_sidecars(dir.path());
        assert!(map.is_empty());
    }

    #[test]
    fn test_load_sidecars_loads_valid_skips_invalid() {
        let dir = tempdir().unwrap();

        let valid_fhr = FhrMetadata {
            genome: Some("ValidGenome".to_string()),
            ..Default::default()
        };
        write_sidecar(&dir.path().join("validdigest.fhr.json"), &valid_fhr).unwrap();

        fs::write(dir.path().join("baddigest.fhr.json"), "not json at all").unwrap();

        let map = load_sidecars(dir.path());
        assert_eq!(map.len(), 1);
    }

    // =========================================================================
    // Store-level FHR integration tests
    // =========================================================================

    #[test]
    fn test_fhr_metadata_empty_by_default() {
        let mut store = RefgetStore::in_memory();

        let (meta, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        assert!(store.get_fhr_metadata(&meta.digest).is_none());
        assert!(store.list_fhr_metadata().is_empty());
    }

    #[test]
    fn test_fhr_metadata_set_get() {
        let mut store = RefgetStore::in_memory();
        let (meta, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let mut fhr = FhrMetadata::default();
        fhr.genome = Some("Test genome".to_string());
        fhr.version = Some("1.0".to_string());
        fhr.masking = Some("not-masked".to_string());

        store.set_fhr_metadata(&meta.digest, fhr.clone()).unwrap();

        let retrieved = store.get_fhr_metadata(&meta.digest).unwrap();
        assert_eq!(retrieved.genome, Some("Test genome".to_string()));
        assert_eq!(retrieved.version, Some("1.0".to_string()));
    }

    #[test]
    fn test_fhr_metadata_nonexistent_collection() {
        let mut store = RefgetStore::in_memory();
        let fhr = FhrMetadata::default();
        assert!(store.set_fhr_metadata("nonexistent_digest", fhr).is_err());
    }

    #[test]
    fn test_fhr_metadata_remove() {
        let mut store = RefgetStore::in_memory();
        let (meta, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let fhr = FhrMetadata {
            genome: Some("Test".to_string()),
            ..Default::default()
        };
        store.set_fhr_metadata(&meta.digest, fhr).unwrap();

        assert!(store.get_fhr_metadata(&meta.digest).is_some());
        assert!(store.remove_fhr_metadata(&meta.digest));
        assert!(store.get_fhr_metadata(&meta.digest).is_none());
    }

    #[test]
    fn test_fhr_metadata_persistence() {
        let dir = tempdir().unwrap();
        let store_path = dir.path().join("store");
        let digest: String;

        {
            let mut store = RefgetStore::on_disk(&store_path).unwrap();
            let (meta, _) = store
                .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
                .unwrap();
            digest = meta.digest.clone();

            let fhr = FhrMetadata {
                genome: Some("Homo sapiens".to_string()),
                version: Some("GRCh38".to_string()),
                masking: Some("soft-masked".to_string()),
                ..Default::default()
            };
            store.set_fhr_metadata(&digest, fhr).unwrap();
        }

        {
            let store = RefgetStore::open_local(&store_path).unwrap();
            let fhr = store.get_fhr_metadata(&digest).unwrap();
            assert_eq!(fhr.genome, Some("Homo sapiens".to_string()));
            assert_eq!(fhr.version, Some("GRCh38".to_string()));
            assert_eq!(fhr.masking, Some("soft-masked".to_string()));
        }
    }

    #[test]
    fn test_fhr_list() {
        let mut store = RefgetStore::in_memory();
        assert!(store.list_fhr_metadata().is_empty());

        let (meta, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();
        let fhr = FhrMetadata {
            genome: Some("Test".to_string()),
            ..Default::default()
        };
        store.set_fhr_metadata(&meta.digest, fhr).unwrap();

        let list = store.list_fhr_metadata();
        assert_eq!(list.len(), 1);
        assert!(list.contains(&meta.digest));
    }

    #[test]
    fn test_remove_collection_cleans_up_fhr_metadata() {
        let dir = tempdir().unwrap();
        let fasta = dir.path().join("test.fa");
        std::fs::write(&fasta, ">chr1\nACGT\n").unwrap();

        let mut store = RefgetStore::in_memory();
        let (meta, _) = store
            .add_sequence_collection_from_fasta(&fasta, FastaImportOptions::new())
            .unwrap();
        let digest = meta.digest;

        let fhr = FhrMetadata::default();
        store.set_fhr_metadata(&digest, fhr).unwrap();
        assert!(store.get_fhr_metadata(&digest).is_some());

        store.remove_collection(&digest, false).unwrap();

        assert!(store.get_fhr_metadata(&digest).is_none());
    }

    // -----------------------------------------------------------------------
    // KeepOurs sync strategy tests (regression test for was_local ordering bug)
    // -----------------------------------------------------------------------

    /// Spin up a minimal HTTP server serving files from `serve_dir`.
    /// Returns `(base_url, shutdown_fn)`.
    fn start_file_server(serve_dir: std::path::PathBuf) -> (String, impl FnOnce()) {
        use std::io::{Read as _, Write as _};
        use std::net::TcpListener;
        use std::sync::{Arc, atomic::{AtomicBool, Ordering}};

        let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
        let port = listener.local_addr().unwrap().port();
        let base_url = format!("http://127.0.0.1:{}", port);
        let stop = Arc::new(AtomicBool::new(false));
        let stop_clone = Arc::clone(&stop);

        std::thread::spawn(move || {
            listener.set_nonblocking(false).ok();
            while !stop_clone.load(Ordering::Relaxed) {
                match listener.accept() {
                    Ok((mut stream, _)) => {
                        let mut buf = [0u8; 4096];
                        let n = stream.read(&mut buf).unwrap_or(0);
                        let request = std::str::from_utf8(&buf[..n]).unwrap_or("");
                        let path = request
                            .lines()
                            .next()
                            .and_then(|l| l.split_whitespace().nth(1))
                            .unwrap_or("/");
                        let rel = path.trim_start_matches('/');
                        let file_path = serve_dir.join(rel);
                        if file_path.exists() && file_path.is_file() {
                            let data = fs::read(&file_path).unwrap_or_default();
                            let header = format!(
                                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                                data.len()
                            );
                            let _ = stream.write_all(header.as_bytes());
                            let _ = stream.write_all(&data);
                        } else {
                            let body = b"Not Found";
                            let header = format!(
                                "HTTP/1.1 404 Not Found\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                                body.len()
                            );
                            let _ = stream.write_all(header.as_bytes());
                            let _ = stream.write_all(body);
                        }
                    }
                    Err(_) => break,
                }
            }
        });

        let shutdown = move || {
            stop.store(true, Ordering::Relaxed);
            let _ = std::net::TcpStream::connect(format!("127.0.0.1:{}", port));
        };

        (base_url, shutdown)
    }

    /// Pull an FHR sidecar that does NOT exist locally yet.
    /// KeepOurs: first pull should count as `pulled`, second pull as `skipped`.
    #[test]
    fn test_keep_ours_fhr_first_pull_counts_as_pulled() {
        // "Remote" store: directory with a pre-built FHR JSON sidecar.
        let remote_dir = tempdir().unwrap();
        let collections_dir = remote_dir.path().join("collections");
        fs::create_dir_all(&collections_dir).unwrap();

        // We need a fake digest string to use as the collection identity.
        let fake_digest = "SQ.aaaaaaaaaaaaaaaaaaaaaaaa";
        let sidecar_name = format!("{}.fhr.json", fake_digest);
        let fhr = FhrMetadata {
            genome: Some("TestGenome".to_string()),
            ..Default::default()
        };
        let sidecar_json = serde_json::to_string(&fhr).unwrap();
        fs::write(collections_dir.join(&sidecar_name), &sidecar_json).unwrap();

        // Start HTTP server.
        let (base_url, shutdown) = start_file_server(remote_dir.path().to_path_buf());

        // "Local" store: disk-backed, with a stub collection so pull_fhr has a digest.
        let local_dir = tempdir().unwrap();
        let local_store_path = local_dir.path().join("store");

        let mut store = RefgetStore::on_disk(&local_store_path).unwrap();
        store.inner.remote_source = Some(base_url);

        // Inject a minimal stub collection so pull_fhr iterates over it.
        use crate::hashkeyable::HashKeyable;
        use crate::digest::SequenceCollectionRecord;
        let key = fake_digest.to_key();
        let stub = crate::digest::SequenceCollectionMetadata {
            digest: fake_digest.to_string(),
            n_sequences: 0,
            names_digest: String::new(),
            sequences_digest: String::new(),
            lengths_digest: String::new(),
            name_length_pairs_digest: None,
            sorted_name_length_pairs_digest: None,
            sorted_sequences_digest: None,
            file_path: None,
        };
        store.inner.collections.insert(key, SequenceCollectionRecord::Stub(stub));

        // First pull: FHR sidecar not yet local → should be pulled.
        let result = store.pull_fhr(Some(fake_digest), SyncStrategy::KeepOurs).unwrap();
        assert_eq!(result.pulled, 1, "first pull should count as pulled, not skipped");
        assert_eq!(result.skipped, 0, "first pull should not be skipped");
        assert_eq!(result.not_found, 0);

        // Second pull: sidecar now on disk → should be skipped.
        let result2 = store.pull_fhr(Some(fake_digest), SyncStrategy::KeepOurs).unwrap();
        assert_eq!(result2.skipped, 1, "second pull should be skipped (file already local)");
        assert_eq!(result2.pulled, 0, "second pull should not count as pulled");

        shutdown();
    }
}
