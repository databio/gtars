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

use crate::hashkeyable::HashKeyable;

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

    /// FHR schema version (e.g., "1.0")
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub schema_version: Option<String>,

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

    /// Masking type
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub masking: Option<String>,

    /// File-level checksum (SHA2-512/256 per FHR spec)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub checksum: Option<String>,

    /// Alternative common names for this genome
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_synonym: Option<Vec<String>>,

    /// Database accession identifiers
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub accession_id: Option<Vec<FhrIdentifier>>,

    /// Sequencing instruments used
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub instrument: Option<Vec<String>>,

    /// DOI or scholarly article references
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scholarly_article: Option<Vec<String>>,

    /// License information
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub license: Option<String>,

    /// Related URLs
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub related_link: Option<Vec<String>>,

    /// Funding information
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub funding: Option<Vec<String>>,

    /// Seqcol digest (added by RefgetStore, not part of FHR 1.0)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seqcol_digest: Option<String>,

    /// Catch-all for any other FHR fields or custom extensions
    #[serde(flatten)]
    pub extra: HashMap<String, serde_json::Value>,
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

/// Convert a [u8; 32] key back to a digest string.
fn key_to_digest_string(key: &[u8; 32]) -> String {
    let len = key.iter().position(|&b| b == 0).unwrap_or(32);
    String::from_utf8_lossy(&key[..len]).to_string()
}

/// Load all FHR sidecar files from a collections directory.
///
/// Scans for `*.fhr.json` files, parses each one, and returns a map
/// keyed by collection digest. Malformed files are silently skipped.
pub fn load_sidecars(collections_dir: &Path) -> HashMap<[u8; 32], FhrMetadata> {
    let mut map = HashMap::new();
    if !collections_dir.exists() {
        return map;
    }
    let entries = match fs::read_dir(collections_dir) {
        Ok(e) => e,
        Err(_) => return map,
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if name.ends_with(SIDECAR_EXTENSION) {
                let digest_str = &name[..name.len() - SIDECAR_EXTENSION.len()];
                let key = digest_str.to_key();
                if let Ok(json) = fs::read_to_string(&path) {
                    if let Ok(fhr) = serde_json::from_str::<FhrMetadata>(&json) {
                        map.insert(key, fhr);
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
            schema_version: Some("1.0".to_string()),
            genome: Some("Homo sapiens".to_string()),
            taxon: Some(FhrTaxon {
                name: Some("Homo sapiens".to_string()),
                uri: Some("https://identifiers.org/taxonomy:9606".to_string()),
            }),
            version: Some("GRCh38.p14".to_string()),
            masking: Some("soft-masked".to_string()),
            genome_synonym: Some(vec!["hg38".to_string()]),
            ..Default::default()
        };

        let json = serde_json::to_string_pretty(&fhr).unwrap();
        let roundtripped: FhrMetadata = serde_json::from_str(&json).unwrap();
        assert_eq!(roundtripped.genome, fhr.genome);
        assert_eq!(roundtripped.taxon.as_ref().unwrap().name, fhr.taxon.as_ref().unwrap().name);
        assert_eq!(roundtripped.genome_synonym, fhr.genome_synonym);
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
            schema_version: Some("1.0".to_string()),
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
}
