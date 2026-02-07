//! Seqcol spec operations: comparison, level-based retrieval, attribute search.
//!
//! These methods are defined on `RefgetStore` but live in this file for
//! organizational clarity. They implement the GA4GH Sequence Collections
//! specification endpoints (level 1/2 retrieval, comparison, attribute search).

use crate::digest::{
    CollectionLevel1, CollectionLevel2, SeqColComparison, SequenceCollectionMetadata,
};
use crate::hashkeyable::HashKeyable;
use crate::store::RefgetStore;
use anyhow::{anyhow, Result};

/// Parse a single line from an RGCI (collection index) file.
///
/// RGCI format is tab-separated with 5+ columns:
/// digest, n_sequences, names_digest, sequences_digest, lengths_digest,
/// [name_length_pairs_digest, sorted_name_length_pairs_digest, sorted_sequences_digest]
///
/// Lines starting with '#' are treated as comments and return None.
/// Lines with fewer than 5 columns return None.
/// Columns 5-7 are optional ancillary digests (empty string = None).
pub(crate) fn parse_rgci_line(line: &str) -> Option<SequenceCollectionMetadata> {
    if line.starts_with('#') {
        return None;
    }
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 5 {
        return None;
    }
    // Parse optional ancillary digest columns (empty string -> None)
    let opt_col = |i: usize| -> Option<String> {
        parts.get(i).and_then(|s| {
            if s.is_empty() { None } else { Some(s.to_string()) }
        })
    };
    Some(SequenceCollectionMetadata {
        digest: parts[0].to_string(),
        n_sequences: parts[1].parse().ok()?,
        names_digest: parts[2].to_string(),
        sequences_digest: parts[3].to_string(),
        lengths_digest: parts[4].to_string(),
        name_length_pairs_digest: opt_col(5),
        sorted_name_length_pairs_digest: opt_col(6),
        sorted_sequences_digest: opt_col(7),
        file_path: None,
    })
}

/// Default limit for brute-force attribute search
pub const DEFAULT_ATTRIBUTE_SEARCH_LIMIT: usize = 10_000;

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

    /// Set the maximum number of collections for brute-force attribute search.
    /// Set to 0 for unlimited.
    pub fn set_attribute_search_limit(&mut self, limit: usize) {
        self.attribute_search_limit = limit;
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
    /// Returns error if collection count exceeds attribute_search_limit (unless 0).
    fn find_collections_by_attribute_scan(
        &self,
        attr_name: &str,
        attr_digest: &str,
    ) -> Result<Vec<String>> {
        let count = self.collections.len();
        if self.attribute_search_limit > 0 && count > self.attribute_search_limit {
            return Err(anyhow!(
                "Too many collections ({}) for brute-force search (limit: {}). \
                 Use set_attribute_search_limit(0) for unlimited or enable attribute index.",
                count,
                self.attribute_search_limit
            ));
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

    /// Indexed lookup from on-disk reverse index.
    /// Stub in Part 1; Part 2 replaces with real implementation.
    fn find_collections_by_attribute_indexed(
        &self,
        _attr_name: &str,
        _attr_digest: &str,
    ) -> Result<Vec<String>> {
        Err(anyhow!(
            "Attribute index not available. \
             Use enable_attribute_index() or fall back to brute-force search."
        ))
    }
}

#[cfg(test)]
mod tests {
    use crate::store::RefgetStore;
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
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
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
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
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
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
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
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
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
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
            .unwrap();
        let (meta_b, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa")
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
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
            .unwrap();
        store.disable_ancillary_digests();
        let (meta_b, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa")
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
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
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
    fn test_find_collections_by_attribute_limit() {
        let mut store = RefgetStore::in_memory();
        store.set_attribute_search_limit(0); // unlimited first

        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
            .unwrap();

        // Should work with limit=0 (unlimited)
        assert!(store
            .find_collections_by_attribute("names", &metadata.names_digest)
            .is_ok());

        // Now set limit very low - but we only have 1 collection so it should still work
        store.set_attribute_search_limit(1);
        assert!(store
            .find_collections_by_attribute("names", &metadata.names_digest)
            .is_ok());
    }

    #[test]
    fn test_get_attribute() {
        let mut store = RefgetStore::in_memory();
        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa")
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
                .add_sequence_collection_from_fasta(&temp_fasta)
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
                .add_sequence_collection_from_fasta(&fasta_path)
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
                .add_sequence_collection_from_fasta(&temp_fasta)
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
