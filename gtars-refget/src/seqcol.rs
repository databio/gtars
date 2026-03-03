//! Seqcol spec operations: comparison, level-based retrieval, attribute search.
//!
//! These methods are defined on `RefgetStore` but live in this file for
//! organizational clarity. They implement the GA4GH Sequence Collections
//! specification endpoints (level 1/2 retrieval, comparison, attribute search).

use crate::digest::{CollectionLevel1, CollectionLevel2, SeqColComparison};
use crate::hashkeyable::HashKeyable;
use crate::store::RefgetStore;
use anyhow::{anyhow, Result};

/// Warn users when brute-force scanning more than this many collections.
const ATTRIBUTE_SEARCH_WARN_THRESHOLD: usize = 10_000;

/// Error if brute-force scanning would exceed this many collections.
const ATTRIBUTE_SEARCH_ERROR_THRESHOLD: usize = 100_000;

impl RefgetStore {
    /// Enable computation and storage of ancillary digests (nlp, snlp, sorted_sequences).
    pub fn enable_ancillary_digests(&mut self) {
        self.ancillary_digests = true;
    }

    /// Disable computation and storage of ancillary digests.
    pub fn disable_ancillary_digests(&mut self) {
        self.ancillary_digests = false;
    }

    /// Returns whether ancillary digests are enabled.
    pub fn has_ancillary_digests(&self) -> bool {
        self.ancillary_digests
    }

    /// Returns whether the on-disk attribute index is enabled.
    pub fn has_attribute_index(&self) -> bool {
        self.attribute_index
    }

    /// Get collection at level 1 representation (attribute digests with spec field names).
    /// This is a lightweight operation that only reads metadata, no loading needed.
    pub fn get_collection_level1(&self, digest: &str) -> Result<CollectionLevel1> {
        let key = digest.to_key();
        let record = self
            .collections
            .get(&key)
            .ok_or_else(|| anyhow!("Collection not found: {}", digest))?;
        Ok(record.metadata().to_level1())
    }

    /// Get collection at level 2 representation (full arrays, spec format).
    /// May need to load the collection from disk/remote.
    pub fn get_collection_level2(&mut self, digest: &str) -> Result<CollectionLevel2> {
        let collection = self.get_collection(digest)?;
        Ok(collection.to_level2())
    }

    /// Compare two collections by digest. Loads both if needed.
    pub fn compare(&mut self, digest_a: &str, digest_b: &str) -> Result<SeqColComparison> {
        let coll_a = self.get_collection(digest_a)?;
        let coll_b = self.get_collection(digest_b)?;
        Ok(coll_a.compare(&coll_b))
    }

    /// Find all collections with a specific attribute digest.
    ///
    /// Dispatches to indexed lookup (if attribute_index enabled) or
    /// brute-force metadata scan (default).
    ///
    /// Supported attr_name values: "names", "lengths", "sequences",
    /// "name_length_pairs", "sorted_name_length_pairs", "sorted_sequences"
    pub fn find_collections_by_attribute(
        &self,
        attr_name: &str,
        attr_digest: &str,
    ) -> Result<Vec<String>> {
        if self.attribute_index {
            self.find_collections_by_attribute_indexed(attr_name, attr_digest)
        } else {
            self.find_collections_by_attribute_scan(attr_name, attr_digest)
        }
    }

    /// Brute-force scan of collection metadata.
    /// Warns at 10k collections, errors at 100k collections.
    fn find_collections_by_attribute_scan(
        &self,
        attr_name: &str,
        attr_digest: &str,
    ) -> Result<Vec<String>> {
        let count = self.collections.len();

        if count > ATTRIBUTE_SEARCH_ERROR_THRESHOLD {
            return Err(anyhow!(
                "Brute-force attribute search is limited to {} collections ({} in store). \
                 Indexed attribute lookup is planned for a future release.",
                ATTRIBUTE_SEARCH_ERROR_THRESHOLD,
                count
            ));
        }

        if count > ATTRIBUTE_SEARCH_WARN_THRESHOLD {
            eprintln!(
                "Warning: brute-force attribute search scanning {} collections. \
                 This may be slow.",
                count
            );
        }

        let mut results = Vec::new();
        for record in self.collections.values() {
            let meta = record.metadata();
            let matches = match attr_name {
                "names" => meta.names_digest == attr_digest,
                "lengths" => meta.lengths_digest == attr_digest,
                "sequences" => meta.sequences_digest == attr_digest,
                "name_length_pairs" => meta
                    .name_length_pairs_digest
                    .as_deref()
                    .map_or(false, |d| d == attr_digest),
                "sorted_name_length_pairs" => meta
                    .sorted_name_length_pairs_digest
                    .as_deref()
                    .map_or(false, |d| d == attr_digest),
                "sorted_sequences" => meta
                    .sorted_sequences_digest
                    .as_deref()
                    .map_or(false, |d| d == attr_digest),
                _ => {
                    return Err(anyhow!(
                        "Unknown attribute: '{}'. Supported: names, lengths, sequences, \
                         name_length_pairs, sorted_name_length_pairs, sorted_sequences",
                        attr_name
                    ))
                }
            };
            if matches {
                results.push(meta.digest.clone());
            }
        }
        Ok(results)
    }

    /// Get the raw attribute array for a given attribute digest.
    /// Finds a collection with this attribute (via search), loads it,
    /// and extracts the array.
    ///
    /// Returns the array as a serde_json::Value (array of strings or numbers).
    /// Returns Ok(None) if no collection has this attribute digest.
    pub fn get_attribute(
        &mut self,
        attr_name: &str,
        attr_digest: &str,
    ) -> Result<Option<serde_json::Value>> {
        let collections = self.find_collections_by_attribute(attr_name, attr_digest)?;
        if collections.is_empty() {
            return Ok(None);
        }

        // Load the first matching collection
        let collection = self.get_collection(&collections[0])?;
        let lvl2 = collection.to_level2();

        let value = match attr_name {
            "names" => serde_json::Value::Array(
                lvl2.names
                    .iter()
                    .map(|s| serde_json::Value::String(s.clone()))
                    .collect(),
            ),
            "lengths" => serde_json::Value::Array(
                lvl2.lengths
                    .iter()
                    .map(|l| serde_json::Value::Number(serde_json::Number::from(*l)))
                    .collect(),
            ),
            "sequences" => serde_json::Value::Array(
                lvl2.sequences
                    .iter()
                    .map(|s| serde_json::Value::String(s.clone()))
                    .collect(),
            ),
            _ => {
                return Err(anyhow!(
                    "Cannot retrieve attribute array for '{}'. \
                     Only 'names', 'lengths', and 'sequences' have raw arrays.",
                    attr_name
                ))
            }
        };

        Ok(Some(value))
    }

    /// Enable indexed attribute lookup (not yet implemented).
    ///
    /// Note: The indexed lookup feature is planned for a future release.
    /// Enabling this will cause `find_collections_by_attribute()` to return
    /// a "not implemented" error until the feature is complete.
    pub fn enable_attribute_index(&mut self) {
        self.attribute_index = true;
    }

    /// Disable indexed attribute lookup, using brute-force scan instead.
    pub fn disable_attribute_index(&mut self) {
        self.attribute_index = false;
    }

    /// Indexed lookup from on-disk reverse index.
    /// Stub: not yet implemented. Planned for a future release.
    fn find_collections_by_attribute_indexed(
        &self,
        _attr_name: &str,
        _attr_digest: &str,
    ) -> Result<Vec<String>> {
        Err(anyhow!(
            "Indexed attribute lookup is not yet implemented. \
             This feature is planned for a future release. \
             For now, use the brute-force scan by keeping attribute_index disabled."
        ))
    }
}

#[cfg(test)]
mod tests {
    use crate::store::{FastaImportOptions, RefgetStore};
    use std::path::PathBuf;

    /// Copy a test FASTA to a temp directory to avoid writing RGSI cache files
    /// into the test data directory.
    fn copy_test_fasta(temp_dir: &std::path::Path, name: &str) -> PathBuf {
        let src = format!("../tests/data/fasta/{}", name);
        let dst = temp_dir.join(name);
        std::fs::copy(&src, &dst)
            .unwrap_or_else(|e| panic!("Failed to copy {} to tempdir: {}", src, e));
        dst
    }

    #[test]
    fn test_ancillary_digests_computed() {
        let mut store = RefgetStore::in_memory();
        assert!(store.has_ancillary_digests());

        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        // Ancillary digests should be present
        assert!(metadata.name_length_pairs_digest.is_some());
        assert!(metadata.sorted_name_length_pairs_digest.is_some());
        assert!(metadata.sorted_sequences_digest.is_some());

        // The stored collection should also have them
        let coll_meta = store.get_collection_metadata(&metadata.digest).unwrap();
        assert!(coll_meta.name_length_pairs_digest.is_some());
        assert!(coll_meta.sorted_name_length_pairs_digest.is_some());
        assert!(coll_meta.sorted_sequences_digest.is_some());
    }

    #[test]
    fn test_ancillary_digests_disabled() {
        let mut store = RefgetStore::in_memory();
        store.disable_ancillary_digests();
        assert!(!store.has_ancillary_digests());

        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        // Ancillary digests should NOT be present
        assert!(metadata.name_length_pairs_digest.is_none());
        assert!(metadata.sorted_name_length_pairs_digest.is_none());
        assert!(metadata.sorted_sequences_digest.is_none());
    }

    #[test]
    fn test_collection_level1() {
        let mut store = RefgetStore::in_memory();
        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let lvl1 = store.get_collection_level1(&metadata.digest).unwrap();
        assert_eq!(lvl1.names, metadata.names_digest);
        assert_eq!(lvl1.lengths, metadata.lengths_digest);
        assert_eq!(lvl1.sequences, metadata.sequences_digest);
        assert!(lvl1.name_length_pairs.is_some());
        assert!(lvl1.sorted_name_length_pairs.is_some());
        assert!(lvl1.sorted_sequences.is_some());
    }

    #[test]
    fn test_collection_level2() {
        let mut store = RefgetStore::in_memory();
        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let lvl2 = store.get_collection_level2(&metadata.digest).unwrap();
        assert_eq!(lvl2.names.len(), 3); // chrX, chr1, chr2
        assert_eq!(lvl2.lengths.len(), 3);
        assert_eq!(lvl2.sequences.len(), 3);

        // Check that sequences have SQ. prefix
        for seq in &lvl2.sequences {
            assert!(seq.starts_with("SQ."), "Expected SQ. prefix, got: {}", seq);
        }

        // Check lengths match the FASTA data
        assert!(lvl2.lengths.contains(&8)); // chrX = TTGGGGAA
        assert!(lvl2.lengths.contains(&4)); // chr1 = GGAA, chr2 = GCGC
    }

    #[test]
    fn test_compare_collections() {
        let mut store = RefgetStore::in_memory();
        let (meta_a, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();
        let (meta_b, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa", FastaImportOptions::new())
            .unwrap();

        // Self-compare: identical digests, all elements same order
        let self_result = store.compare(&meta_a.digest, &meta_a.digest).unwrap();
        assert_eq!(self_result.digests.a, self_result.digests.b);
        assert_eq!(self_result.attributes.a_and_b.len(), 6); // 3 core + 3 ancillary
        for attr in &self_result.attributes.a_and_b {
            assert_eq!(self_result.array_elements.a_and_b_same_order[attr], Some(true));
        }

        // Cross-compare: different digests, all 6 attributes shared
        let cross_result = store.compare(&meta_a.digest, &meta_b.digest).unwrap();
        assert_ne!(cross_result.digests.a, cross_result.digests.b);
        assert_eq!(cross_result.attributes.a_and_b.len(), 6);
        assert!(cross_result.attributes.a_only.is_empty());
        assert!(cross_result.attributes.b_only.is_empty());
    }

    #[test]
    fn test_compare_mixed_ancillary() {
        let mut store = RefgetStore::in_memory();
        let (meta_a, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();
        store.disable_ancillary_digests();
        let (meta_b, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa", FastaImportOptions::new())
            .unwrap();

        let result = store.compare(&meta_a.digest, &meta_b.digest).unwrap();
        assert_eq!(result.attributes.a_and_b.len(), 3);
        assert_eq!(result.attributes.a_only.len(), 3);
        assert!(result.attributes.b_only.is_empty());
    }

    #[test]
    fn test_find_collections_by_attribute() {
        let mut store = RefgetStore::in_memory();
        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        // Search by names digest should find our collection
        let results = store
            .find_collections_by_attribute("names", &metadata.names_digest)
            .unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0], metadata.digest);

        // Search by lengths digest
        let results = store
            .find_collections_by_attribute("lengths", &metadata.lengths_digest)
            .unwrap();
        assert_eq!(results.len(), 1);

        // Search by sequences digest
        let results = store
            .find_collections_by_attribute("sequences", &metadata.sequences_digest)
            .unwrap();
        assert_eq!(results.len(), 1);

        // Search by ancillary digest
        let nlp = metadata.name_length_pairs_digest.as_ref().unwrap();
        let results = store
            .find_collections_by_attribute("name_length_pairs", nlp)
            .unwrap();
        assert_eq!(results.len(), 1);

        // Search with nonexistent digest
        let results = store
            .find_collections_by_attribute("names", "nonexistent")
            .unwrap();
        assert!(results.is_empty());

        // Unknown attribute should error
        assert!(store
            .find_collections_by_attribute("unknown", "digest")
            .is_err());
    }

    #[test]
    fn test_get_attribute() {
        let mut store = RefgetStore::in_memory();
        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        // Get names array
        let result = store
            .get_attribute("names", &metadata.names_digest)
            .unwrap();
        assert!(result.is_some());
        let names = result.unwrap();
        assert!(names.is_array());
        assert_eq!(names.as_array().unwrap().len(), 3);

        // Get lengths array
        let result = store
            .get_attribute("lengths", &metadata.lengths_digest)
            .unwrap();
        assert!(result.is_some());

        // Nonexistent digest returns None
        let result = store.get_attribute("names", "nonexistent").unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_rgci_roundtrip_with_ancillary() {
        let dir = tempfile::tempdir().unwrap();
        let dir_path = dir.path();
        let temp_fasta = copy_test_fasta(dir_path, "base.fa");

        // Create and save a store with ancillary digests
        {
            let mut store = RefgetStore::on_disk(dir_path).unwrap();
            store
                .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
                .unwrap();
            store.write().unwrap();
        }

        // Reload and verify ancillary digests survived
        {
            let store = RefgetStore::open_local(dir_path).unwrap();
            let collections = store.list_collections();
            assert_eq!(collections.len(), 1);

            let meta = &collections[0];
            assert!(meta.name_length_pairs_digest.is_some());
            assert!(meta.sorted_name_length_pairs_digest.is_some());
            assert!(meta.sorted_sequences_digest.is_some());
        }
    }

    // ================================================================
    // Compliance tests against Python refget test_fasta_digests.json
    // To add new test cases, edit tests/data/fasta/test_fasta_digests.json
    // and add corresponding .fa files to tests/data/fasta/.
    // ================================================================

    #[test]
    fn test_compliance_digests_from_fixture() {
        let fixture_path = "../tests/data/fasta/test_fasta_digests.json";
        let fixture_str = std::fs::read_to_string(fixture_path)
            .unwrap_or_else(|e| panic!("Failed to read {}: {}", fixture_path, e));
        let fixture: serde_json::Value = serde_json::from_str(&fixture_str)
            .unwrap_or_else(|e| panic!("Failed to parse {}: {}", fixture_path, e));

        let mut store = RefgetStore::in_memory();
        store.enable_ancillary_digests();

        for (fa_name, bundle) in fixture.as_object().unwrap() {
            let fasta_path = format!("../tests/data/fasta/{}", fa_name);
            let (meta, _) = store
                .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
                .unwrap_or_else(|e| panic!("{}: {}", fa_name, e));

            let lvl1 = bundle["level1"].as_object().unwrap();
            let expected_digest = bundle["top_level_digest"].as_str().unwrap();

            assert_eq!(meta.digest, expected_digest, "{}: top_level_digest", fa_name);
            assert_eq!(meta.names_digest, lvl1["names"].as_str().unwrap(), "{}: names", fa_name);
            assert_eq!(meta.lengths_digest, lvl1["lengths"].as_str().unwrap(), "{}: lengths", fa_name);
            assert_eq!(meta.sequences_digest, lvl1["sequences"].as_str().unwrap(), "{}: sequences", fa_name);
            assert_eq!(
                meta.sorted_sequences_digest.as_deref(),
                Some(lvl1["sorted_sequences"].as_str().unwrap()),
                "{}: sorted_sequences", fa_name
            );
            assert_eq!(
                meta.name_length_pairs_digest.as_deref(),
                Some(lvl1["name_length_pairs"].as_str().unwrap()),
                "{}: name_length_pairs", fa_name
            );
            assert_eq!(
                meta.sorted_name_length_pairs_digest.as_deref(),
                Some(lvl1["sorted_name_length_pairs"].as_str().unwrap()),
                "{}: sorted_name_length_pairs", fa_name
            );
        }
    }

    #[test]
    fn test_store_config_persisted() {
        let dir = tempfile::tempdir().unwrap();
        let dir_path = dir.path();
        let temp_fasta = copy_test_fasta(dir_path, "base.fa");

        // Create store with ancillary enabled (default)
        {
            let mut store = RefgetStore::on_disk(dir_path).unwrap();
            assert!(store.has_ancillary_digests());
            assert!(!store.has_attribute_index());
            store
                .add_sequence_collection_from_fasta(&temp_fasta, FastaImportOptions::new())
                .unwrap();
            store.write().unwrap();
        }

        // Reload and verify config
        {
            let store = RefgetStore::open_local(dir_path).unwrap();
            assert!(store.has_ancillary_digests());
            assert!(!store.has_attribute_index());
        }
    }
}
