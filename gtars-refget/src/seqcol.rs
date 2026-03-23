//! Seqcol spec operations: comparison, level-based retrieval, attribute search.
//!
//! These methods are defined on `RefgetStore` but live in this file for
//! organizational clarity. They implement the GA4GH Sequence Collections
//! specification endpoints (level 1/2 retrieval, comparison, attribute search).

use crate::digest::{CollectionLevel1, CollectionLevel2, SeqColComparison, SequenceCollectionMetadata};
use crate::digest::types::{compare_arrays, level2_to_comparison_arrays};
use crate::hashkeyable::HashKeyable;
use crate::store::{PagedResult, ReadonlyRefgetStore};
use anyhow::{anyhow, Result};

/// Trait defining the read-only operations a seqcol API server needs from its backend.
///
/// This is the canonical interface for seqcol server operations. Any backend
/// (filesystem-backed `ReadonlyRefgetStore`, Postgres, remote proxy, in-memory mock)
/// can implement this trait and be used via `Arc<dyn SeqColService>`.
///
/// All methods take `&self` and return concrete types, making the trait object-safe.
pub trait SeqColService {
    /// Level 1: attribute digests only.
    fn get_collection_level1(&self, digest: &str) -> Result<CollectionLevel1>;

    /// Level 2: full arrays with values.
    fn get_collection_level2(&self, digest: &str) -> Result<CollectionLevel2>;

    /// Compare two collections by digest.
    fn compare(&self, digest_a: &str, digest_b: &str) -> Result<SeqColComparison>;

    /// Compare a stored collection against an externally-provided level-2 body.
    fn compare_with_level2(
        &self,
        digest_a: &str,
        external: &CollectionLevel2,
    ) -> Result<SeqColComparison>;

    /// Find collections sharing an attribute digest.
    fn find_collections_by_attribute(
        &self,
        attr_name: &str,
        attr_digest: &str,
    ) -> Result<Vec<String>>;

    /// Get the raw attribute array for a given attribute digest.
    fn get_attribute(
        &self,
        attr_name: &str,
        attr_digest: &str,
    ) -> Result<Option<serde_json::Value>>;

    /// List collections with pagination and optional attribute filters.
    fn list_collections(
        &self,
        page: usize,
        page_size: usize,
        filters: &[(&str, &str)],
    ) -> Result<PagedResult<SequenceCollectionMetadata>>;

    /// Total number of collections in the store.
    fn collection_count(&self) -> usize;
}

/// Check if a collection's metadata matches a single attribute filter.
/// Returns Err if attr_name is unrecognized.
pub(crate) fn metadata_matches_attribute(
    meta: &SequenceCollectionMetadata,
    attr_name: &str,
    attr_digest: &str,
) -> Result<bool> {
    match attr_name {
        "names" => Ok(meta.names_digest == attr_digest),
        "lengths" => Ok(meta.lengths_digest == attr_digest),
        "sequences" => Ok(meta.sequences_digest == attr_digest),
        "name_length_pairs" => Ok(meta
            .name_length_pairs_digest
            .as_deref()
            .map_or(false, |d| d == attr_digest)),
        "sorted_name_length_pairs" => Ok(meta
            .sorted_name_length_pairs_digest
            .as_deref()
            .map_or(false, |d| d == attr_digest)),
        "sorted_sequences" => Ok(meta
            .sorted_sequences_digest
            .as_deref()
            .map_or(false, |d| d == attr_digest)),
        _ => Err(anyhow!(
            "Unknown attribute: '{}'. Supported: names, lengths, sequences, \
             name_length_pairs, sorted_name_length_pairs, sorted_sequences",
            attr_name
        )),
    }
}

/// Warn users when brute-force scanning more than this many collections.
const ATTRIBUTE_SEARCH_WARN_THRESHOLD: usize = 10_000;

/// Error if brute-force scanning would exceed this many collections.
const ATTRIBUTE_SEARCH_ERROR_THRESHOLD: usize = 100_000;

impl ReadonlyRefgetStore {
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
    pub fn get_collection_level2(&self, digest: &str) -> Result<CollectionLevel2> {
        let collection = self.get_collection(digest)?;
        Ok(collection.to_level2())
    }

    /// Compare two collections by digest. Both must be preloaded.
    pub fn compare(&self, digest_a: &str, digest_b: &str) -> Result<SeqColComparison> {
        let coll_a = self.get_collection(digest_a)?;
        let coll_b = self.get_collection(digest_b)?;
        Ok(coll_a.compare(&coll_b))
    }

    /// Compare a stored collection (by digest) against an externally-provided level-2 body.
    ///
    /// Used for the seqcol spec `POST /comparison/:digest1` endpoint where the client
    /// submits a local collection as JSON rather than referencing a stored digest.
    ///
    /// The returned `SeqColComparison` has `digests.a` set to the stored collection's
    /// digest and `digests.b` set to `None` because the external collection has no
    /// server-side digest.
    pub fn compare_with_level2(
        &self,
        digest_a: &str,
        external: &CollectionLevel2,
    ) -> Result<SeqColComparison> {
        let coll_a = self.get_collection(digest_a)?;
        let arrays_a = coll_a.to_comparison_arrays();
        let arrays_b = level2_to_comparison_arrays(external);
        Ok(compare_arrays(
            arrays_a,
            arrays_b,
            coll_a.metadata.digest.clone(),
            None,
        ))
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
            if metadata_matches_attribute(meta, attr_name, attr_digest)? {
                results.push(meta.digest.clone());
            }
        }
        Ok(results)
    }

    /// Get the raw attribute array for a given attribute digest.
    /// Finds a collection with this attribute (via search), loads it,
    /// and extracts the array.
    ///
    /// Supported attr_name values: "names", "lengths", "sequences",
    /// "name_length_pairs", "sorted_name_length_pairs", "sorted_sequences"
    ///
    /// Returns the array as a serde_json::Value (array of strings or numbers).
    /// Returns Ok(None) if no collection has this attribute digest.
    pub fn get_attribute(
        &self,
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
            "sorted_sequences" => serde_json::Value::Array(
                collection
                    .build_sorted_sequences()
                    .into_iter()
                    .map(serde_json::Value::String)
                    .collect(),
            ),
            "name_length_pairs" => {
                serde_json::Value::Array(collection.build_name_length_pairs())
            }
            "sorted_name_length_pairs" => {
                serde_json::Value::Array(collection.build_sorted_name_length_pairs())
            }
            _ => {
                return Err(anyhow!(
                    "Unknown attribute: '{}'. Supported: names, lengths, sequences, \
                     name_length_pairs, sorted_name_length_pairs, sorted_sequences",
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

    /// Total number of collections in the store.
    pub fn collection_count(&self) -> usize {
        self.collections.len()
    }
}

impl SeqColService for ReadonlyRefgetStore {
    fn get_collection_level1(&self, digest: &str) -> Result<CollectionLevel1> {
        ReadonlyRefgetStore::get_collection_level1(self, digest)
    }

    fn get_collection_level2(&self, digest: &str) -> Result<CollectionLevel2> {
        ReadonlyRefgetStore::get_collection_level2(self, digest)
    }

    fn compare(&self, digest_a: &str, digest_b: &str) -> Result<SeqColComparison> {
        ReadonlyRefgetStore::compare(self, digest_a, digest_b)
    }

    fn compare_with_level2(
        &self,
        digest_a: &str,
        external: &CollectionLevel2,
    ) -> Result<SeqColComparison> {
        ReadonlyRefgetStore::compare_with_level2(self, digest_a, external)
    }

    fn find_collections_by_attribute(
        &self,
        attr_name: &str,
        attr_digest: &str,
    ) -> Result<Vec<String>> {
        ReadonlyRefgetStore::find_collections_by_attribute(self, attr_name, attr_digest)
    }

    fn get_attribute(
        &self,
        attr_name: &str,
        attr_digest: &str,
    ) -> Result<Option<serde_json::Value>> {
        ReadonlyRefgetStore::get_attribute(self, attr_name, attr_digest)
    }

    fn list_collections(
        &self,
        page: usize,
        page_size: usize,
        filters: &[(&str, &str)],
    ) -> Result<PagedResult<SequenceCollectionMetadata>> {
        ReadonlyRefgetStore::list_collections(self, page, page_size, filters)
    }

    fn collection_count(&self) -> usize {
        ReadonlyRefgetStore::collection_count(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::store::{FastaImportOptions, ReadonlyRefgetStore, RefgetStore};
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
        assert_eq!(Some(self_result.digests.a.as_str()), self_result.digests.b.as_deref());
        assert_eq!(self_result.attributes.a_and_b.len(), 6); // 3 core + 3 ancillary
        for attr in &self_result.attributes.a_and_b {
            assert_eq!(self_result.array_elements.a_and_b_same_order[attr], Some(true));
        }

        // Cross-compare: different digests, all 6 attributes shared
        let cross_result = store.compare(&meta_a.digest, &meta_b.digest).unwrap();
        assert_ne!(Some(cross_result.digests.a.as_str()), cross_result.digests.b.as_deref());
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
    fn test_get_attribute_sorted_sequences() {
        let mut store = RefgetStore::in_memory();
        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let digest = metadata.sorted_sequences_digest.as_ref().unwrap();
        let result = store.get_attribute("sorted_sequences", digest).unwrap();
        assert!(result.is_some());
        let arr = result.unwrap();
        assert!(arr.is_array());
        let items = arr.as_array().unwrap();
        assert_eq!(items.len(), 3);

        // All items should be strings with SQ. prefix
        for item in items {
            let s = item.as_str().unwrap();
            assert!(s.starts_with("SQ."), "Expected SQ. prefix, got: {}", s);
        }

        // Should be sorted lexicographically
        let strings: Vec<&str> = items.iter().map(|v| v.as_str().unwrap()).collect();
        let mut sorted = strings.clone();
        sorted.sort();
        assert_eq!(strings, sorted, "sorted_sequences should be in sorted order");
    }

    #[test]
    fn test_get_attribute_name_length_pairs() {
        let mut store = RefgetStore::in_memory();
        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let digest = metadata.name_length_pairs_digest.as_ref().unwrap();
        let result = store.get_attribute("name_length_pairs", digest).unwrap();
        assert!(result.is_some());
        let arr = result.unwrap();
        assert!(arr.is_array());
        let items = arr.as_array().unwrap();
        assert_eq!(items.len(), 3);

        // Each item should be an object with "name" (string) and "length" (number) keys
        for item in items {
            let obj = item.as_object().unwrap();
            assert!(obj.contains_key("name"), "Expected 'name' key in object");
            assert!(obj.contains_key("length"), "Expected 'length' key in object");
            assert!(obj["name"].is_string(), "name should be a string");
            assert!(obj["length"].is_number(), "length should be a number");
        }
    }

    #[test]
    fn test_get_attribute_sorted_name_length_pairs() {
        use crate::digest::algorithms::{canonicalize_json, sha512t24u};

        let mut store = RefgetStore::in_memory();
        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let digest = metadata.sorted_name_length_pairs_digest.as_ref().unwrap();
        let result = store
            .get_attribute("sorted_name_length_pairs", digest)
            .unwrap();
        assert!(result.is_some());
        let arr = result.unwrap();
        assert!(arr.is_array());
        let items = arr.as_array().unwrap();
        assert_eq!(items.len(), 3);

        // Each item should be an object with "name" and "length" keys
        for item in items {
            let obj = item.as_object().unwrap();
            assert!(obj.contains_key("name"));
            assert!(obj.contains_key("length"));
        }

        // Verify the objects are sorted by their canonical JSON digest
        let digests: Vec<String> = items
            .iter()
            .map(|v| sha512t24u(canonicalize_json(v).as_bytes()))
            .collect();
        let mut sorted_digests = digests.clone();
        sorted_digests.sort();
        assert_eq!(
            digests, sorted_digests,
            "sorted_name_length_pairs objects should be in sorted digest order"
        );
    }

    #[test]
    fn test_get_attribute_ancillary_not_computed() {
        let mut store = RefgetStore::in_memory();
        store.disable_ancillary_digests();
        let (metadata, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        // Ancillary digests not present, so search returns empty, get_attribute returns None
        assert!(metadata.name_length_pairs_digest.is_none());
        let result = store
            .get_attribute("name_length_pairs", "some_digest")
            .unwrap();
        assert!(
            result.is_none(),
            "Expected None when no ancillary digests are computed"
        );
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
            let collections = store.list_collections(0, usize::MAX, &[]).unwrap();
            assert_eq!(collections.results.len(), 1);

            let meta = &collections.results[0];
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

    // ================================================================
    // Tests for compare_with_level2
    // ================================================================

    /// Test 1: compare_with_level2 produces same result as compare when inputs are identical.
    #[test]
    fn test_compare_with_level2_self_identical() {
        let mut store = RefgetStore::in_memory();
        let (meta, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        // Get the level-2 representation of the stored collection
        let level2 = store.get_collection_level2(&meta.digest).unwrap();

        // Compare the stored collection against its own level-2 body
        let result = store.compare_with_level2(&meta.digest, &level2).unwrap();

        // digests.a is the stored digest, digests.b is None (no server-side digest for external body)
        assert_eq!(result.digests.a, meta.digest);
        assert!(result.digests.b.is_none(), "digests.b should be None for external level-2 comparison");

        // All three core attributes should be shared (a_and_b)
        assert!(result.attributes.a_and_b.contains(&"names".to_string()));
        assert!(result.attributes.a_and_b.contains(&"lengths".to_string()));
        assert!(result.attributes.a_and_b.contains(&"sequences".to_string()));

        // All shared core attributes should be in same order
        for attr in &["names", "lengths", "sequences"] {
            assert_eq!(
                result.array_elements.a_and_b_same_order[*attr],
                Some(true),
                "{} should be in same order",
                attr
            );
        }
    }

    /// Test 2: compare_with_level2 with a different external collection produces equivalent
    /// result to compare(), except digests.b is None instead of Some(digest_b).
    #[test]
    fn test_compare_with_level2_cross_compare() {
        let mut store = RefgetStore::in_memory();
        let (meta_a, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();
        let (meta_b, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa", FastaImportOptions::new())
            .unwrap();

        // Get level-2 of B
        let level2_b = store.get_collection_level2(&meta_b.digest).unwrap();

        // compare via digest-digest
        let compare_result = store.compare(&meta_a.digest, &meta_b.digest).unwrap();
        // compare via digest + external level-2
        let with_level2_result = store.compare_with_level2(&meta_a.digest, &level2_b).unwrap();

        // digests.b is None in compare_with_level2 result, Some(digest_b) in compare result
        assert_eq!(compare_result.digests.b, Some(meta_b.digest.clone()));
        assert!(with_level2_result.digests.b.is_none());

        // Everything else should match
        assert_eq!(compare_result.digests.a, with_level2_result.digests.a);
        // The a_and_b core attributes should agree (both compare the same data for core attributes)
        for attr in &["names", "lengths", "sequences"] {
            assert!(
                with_level2_result.attributes.a_and_b.contains(&attr.to_string())
                    || with_level2_result.attributes.a_only.contains(&attr.to_string())
                    || with_level2_result.attributes.b_only.contains(&attr.to_string()),
                "attr {} must appear somewhere",
                attr
            );
        }
    }

    /// Test 3: compare_with_level2 with an unknown digest returns an error.
    #[test]
    fn test_compare_with_level2_unknown_digest_returns_error() {
        let store = RefgetStore::in_memory();
        // Build a minimal CollectionLevel2 body
        use crate::digest::CollectionLevel2;
        let level2 = CollectionLevel2 {
            names: vec!["chr1".to_string()],
            lengths: vec![100],
            sequences: vec!["SQ.aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa".to_string()],
        };

        let result = store.compare_with_level2("nonexistent_digest", &level2);
        assert!(result.is_err(), "Expected error for unknown digest");
    }

    /// Test 4: ancillary attributes from the stored collection appear in a_only when
    /// the external level-2 body does not contain them.
    #[test]
    fn test_compare_with_level2_ancillary_in_a_only() {
        let mut store = RefgetStore::in_memory();
        // Import with ancillary digests enabled (default)
        let (meta, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        // Verify ancillary digests were computed
        assert!(meta.name_length_pairs_digest.is_some());
        assert!(meta.sorted_name_length_pairs_digest.is_some());
        assert!(meta.sorted_sequences_digest.is_some());

        // The level-2 body has only the three core attributes — no ancillary
        let level2 = store.get_collection_level2(&meta.digest).unwrap();

        let result = store.compare_with_level2(&meta.digest, &level2).unwrap();

        // Ancillary attributes that are in stored collection but not in level-2 body
        // should appear in a_only
        assert!(
            result.attributes.a_only.contains(&"sorted_sequences".to_string()),
            "sorted_sequences should be in a_only"
        );
        assert!(
            result.attributes.a_only.contains(&"name_length_pairs".to_string()),
            "name_length_pairs should be in a_only"
        );
        assert!(
            result.attributes.a_only.contains(&"sorted_name_length_pairs".to_string()),
            "sorted_name_length_pairs should be in a_only"
        );

        // b_only should be empty (level-2 body only has core attributes)
        assert!(
            result.attributes.b_only.is_empty(),
            "b_only should be empty when level-2 has only core attributes"
        );
    }

    // =========================================================================
    // list_collections pagination and filtering tests
    // =========================================================================

    #[test]
    fn test_list_collections_paged_no_filters() {
        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa", FastaImportOptions::new())
            .unwrap();

        let result = store.list_collections(0, 2, &[]).unwrap();
        assert_eq!(result.results.len(), 2);
        assert_eq!(result.pagination.total, 2);
        assert_eq!(result.pagination.page, 0);
        assert_eq!(result.pagination.page_size, 2);
        // Results should be sorted by digest
        assert!(result.results[0].digest <= result.results[1].digest);
    }

    #[test]
    fn test_list_collections_paged_second_page() {
        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa", FastaImportOptions::new())
            .unwrap();

        // Page 0 with page_size 1
        let page0 = store.list_collections(0, 1, &[]).unwrap();
        assert_eq!(page0.results.len(), 1);
        assert_eq!(page0.pagination.total, 2);

        // Page 1 with page_size 1
        let page1 = store.list_collections(1, 1, &[]).unwrap();
        assert_eq!(page1.results.len(), 1);
        assert_eq!(page1.pagination.total, 2);

        // Different results on different pages
        assert_ne!(page0.results[0].digest, page1.results[0].digest);
    }

    #[test]
    fn test_list_collections_paged_beyond_end() {
        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let result = store.list_collections(10, 100, &[]).unwrap();
        assert!(result.results.is_empty());
        assert_eq!(result.pagination.total, 1);
    }

    #[test]
    fn test_list_collections_single_filter() {
        let mut store = RefgetStore::in_memory();
        let (meta, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa", FastaImportOptions::new())
            .unwrap();

        let result = store.list_collections(0, 100, &[("names", &meta.names_digest)]).unwrap();
        assert_eq!(result.results.len(), 1);
        assert_eq!(result.results[0].digest, meta.digest);
        assert_eq!(result.pagination.total, 1);
    }

    #[test]
    fn test_list_collections_multi_filter_and() {
        let mut store = RefgetStore::in_memory();
        let (meta, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa", FastaImportOptions::new())
            .unwrap();

        // Filter by both names AND lengths -- should match only base.fa
        let result = store.list_collections(0, 100, &[
            ("names", &meta.names_digest),
            ("lengths", &meta.lengths_digest),
        ]).unwrap();
        assert_eq!(result.results.len(), 1);
        assert_eq!(result.results[0].digest, meta.digest);
    }

    #[test]
    fn test_list_collections_filter_no_match() {
        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let result = store.list_collections(0, 100, &[("names", "nonexistent_digest")]).unwrap();
        assert!(result.results.is_empty());
        assert_eq!(result.pagination.total, 0);
    }

    #[test]
    fn test_list_collections_invalid_attribute() {
        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let result = store.list_collections(0, 100, &[("unknown_attr", "digest")]);
        assert!(result.is_err());
    }

    #[test]
    fn test_list_collections_filter_with_pagination() {
        let mut store = RefgetStore::in_memory();
        // base.fa and different_names.fa share the same lengths
        let (meta_a, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();
        let (_meta_b, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa", FastaImportOptions::new())
            .unwrap();

        // Both collections share lengths -- filter by lengths, page_size=1
        let page0 = store.list_collections(0, 1, &[("lengths", &meta_a.lengths_digest)]).unwrap();
        assert_eq!(page0.results.len(), 1);
        assert_eq!(page0.pagination.total, 2); // both match lengths filter

        let page1 = store.list_collections(1, 1, &[("lengths", &meta_a.lengths_digest)]).unwrap();
        assert_eq!(page1.results.len(), 1);
        assert_eq!(page1.pagination.total, 2);

        assert_ne!(page0.results[0].digest, page1.results[0].digest);
    }

    // =========================================================================
    // SeqColService trait tests
    // =========================================================================

    #[test]
    fn test_trait_object_safety() {
        use super::SeqColService;
        use std::sync::Arc;

        let mut store = RefgetStore::in_memory();
        let (meta, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let readonly = store.into_readonly();
        let service: Arc<dyn SeqColService + Send + Sync> = Arc::new(readonly);

        // Call each trait method through the trait object
        let lvl1 = service.get_collection_level1(&meta.digest).unwrap();
        assert_eq!(lvl1.names, meta.names_digest);

        let lvl2 = service.get_collection_level2(&meta.digest).unwrap();
        assert_eq!(lvl2.names.len(), 3);

        let cmp = service.compare(&meta.digest, &meta.digest).unwrap();
        assert_eq!(cmp.digests.a, meta.digest);

        let cmp2 = service.compare_with_level2(&meta.digest, &lvl2).unwrap();
        assert_eq!(cmp2.digests.a, meta.digest);

        let found = service
            .find_collections_by_attribute("names", &meta.names_digest)
            .unwrap();
        assert_eq!(found.len(), 1);

        let attr = service
            .get_attribute("names", &meta.names_digest)
            .unwrap();
        assert!(attr.is_some());

        let paged = service.list_collections(0, 10, &[]).unwrap();
        assert_eq!(paged.results.len(), 1);

        assert_eq!(service.collection_count(), 1);
    }

    #[test]
    fn test_collection_count() {
        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();
        store
            .add_sequence_collection_from_fasta("../tests/data/fasta/different_names.fa", FastaImportOptions::new())
            .unwrap();

        assert_eq!(store.collection_count(), 2);
    }

    #[test]
    fn test_trait_methods_match_concrete() {
        use super::SeqColService;

        let mut store = RefgetStore::in_memory();
        let (meta, _) = store
            .add_sequence_collection_from_fasta("../tests/data/fasta/base.fa", FastaImportOptions::new())
            .unwrap();

        let readonly = store.into_readonly();

        // Call via concrete method
        let concrete_lvl1 = ReadonlyRefgetStore::get_collection_level1(&readonly, &meta.digest).unwrap();
        // Call via trait
        let trait_ref: &dyn SeqColService = &readonly;
        let trait_lvl1 = trait_ref.get_collection_level1(&meta.digest).unwrap();

        assert_eq!(concrete_lvl1.names, trait_lvl1.names);
        assert_eq!(concrete_lvl1.lengths, trait_lvl1.lengths);
        assert_eq!(concrete_lvl1.sequences, trait_lvl1.sequences);

        // list_collections
        let concrete_list = ReadonlyRefgetStore::list_collections(&readonly, 0, 10, &[]).unwrap();
        let trait_list = trait_ref.list_collections(0, 10, &[]).unwrap();
        assert_eq!(concrete_list.results.len(), trait_list.results.len());
        assert_eq!(concrete_list.pagination.total, trait_list.pagination.total);

        // collection_count
        assert_eq!(
            ReadonlyRefgetStore::collection_count(&readonly),
            trait_ref.collection_count()
        );
    }
}
