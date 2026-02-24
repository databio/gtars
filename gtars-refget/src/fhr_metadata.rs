//! FAIR Headers Reference genome (FHR) metadata for sequence collections.
//!
//! This module is self-contained: it defines the FHR data types and provides
//! free functions for reading/writing FHR sidecar JSON files. RefgetStore
//! calls into these helpers but the module has no dependency on RefgetStore.
//!
//! See: https://github.com/FAIR-bioHeaders/FHR-Specification

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
/// enforce FHR schema compliance — that's the user's responsibility.
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
// Disk I/O helpers — called by RefgetStore, no dependency on it
// ============================================================================

const SIDECAR_EXTENSION: &str = ".fhr.json";

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
// Tests — serialization, roundtripping, disk I/O (no RefgetStore needed)
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

        // Round-trip preserves extra fields
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
        // snake_case should NOT appear
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
        // Should not panic
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
        // Must be uppercase "ID", not camelCase "Id"
        assert!(json.contains("accessionID"));
        assert!(!json.contains("accessionId"));
    }

    #[test]
    fn test_schema_version_as_number() {
        // Integer
        let json = r#"{"schemaVersion": 1}"#;
        let fhr: FhrMetadata = serde_json::from_str(json).unwrap();
        assert!(fhr.schema_version.is_some());
        let ver = fhr.schema_version.unwrap();
        assert_eq!(ver.to_string(), "1");

        // Float
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
        // Uppercase stat names
        assert!(json.contains("\"L50\""));
        assert!(json.contains("\"N50\""));
        assert!(json.contains("\"L90\""));
        assert!(json.contains("\"totalBasePairs\""));
        assert!(json.contains("\"numberContigs\""));
        // Roundtrip
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
        // Non-spec fields captured in extra
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
        // seqcol_digest should NOT appear in serialized JSON
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

        // Write a valid sidecar
        let valid_fhr = FhrMetadata {
            genome: Some("ValidGenome".to_string()),
            ..Default::default()
        };
        write_sidecar(&dir.path().join("validdigest.fhr.json"), &valid_fhr).unwrap();

        // Write an invalid sidecar
        fs::write(dir.path().join("baddigest.fhr.json"), "not json at all").unwrap();

        let map = load_sidecars(dir.path());
        assert_eq!(map.len(), 1);
    }
}
