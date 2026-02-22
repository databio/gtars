use extendr_api::prelude::*;
use gtars_refget::collection::{
    SeqColDigestLvl1, SequenceCollection, SequenceCollectionMetadata, SequenceMetadata,
    SequenceRecord,
};
use gtars_refget::digest::{md5, sha512t24u, AlphabetType};
use gtars_refget::fhr_metadata::FhrMetadata;
use gtars_refget::store::{RefgetStore, RetrievedSequence, StorageMode};

// =========================================================================
// Helper macro for extracting store from external pointer
// =========================================================================

macro_rules! with_store {
    ($store_ptr:expr, $store:ident, $body:block) => {{
        let mut ext_ptr = <ExternalPtr<RefgetStore>>::try_from($store_ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid store pointer".into()))?;
        let $store = &mut *ext_ptr;
        $body
    }};
}

macro_rules! with_store_ref {
    ($store_ptr:expr, $store:ident, $body:block) => {{
        let ext_ptr = <ExternalPtr<RefgetStore>>::try_from($store_ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid store pointer".into()))?;
        let $store = &*ext_ptr;
        $body
    }};
}

// =========================================================================
// Digest Functions
// =========================================================================

/// Create sha512t24u digest
/// @export
/// @param readable A readable string representing a sequence.
#[extendr]
pub fn sha512t24u_digest(readable: &str) -> String {
    sha512t24u(readable.as_bytes())
}

/// Create md5 digest
/// @export
/// @param readable A readable string representing a sequence.
#[extendr]
pub fn md5_digest(readable: &str) -> String {
    md5(readable.as_bytes())
}

/// Digest fasta file
/// @param fasta A filepath string to a fasta file.
#[extendr]
pub fn digest_fasta_raw(fasta: &str) -> extendr_api::Result<List> {
    match gtars_refget::fasta::digest_fasta(fasta) {
        Ok(sequence_collection) => Ok(sequence_collection_to_list(sequence_collection)),
        Err(e) => Err(format!("Error processing FASTA file: {}", e).into()),
    }
}

// =========================================================================
// Store Constructors
// =========================================================================

/// Create an in-memory RefgetStore
/// @param mode Storage mode: "raw" or "encoded"
#[extendr]
pub fn refget_store_raw(mode: &str) -> extendr_api::Result<Robj> {
    let storage_mode = match mode.to_lowercase().as_str() {
        "raw" => StorageMode::Raw,
        "encoded" => StorageMode::Encoded,
        _ => return Err(format!("Invalid mode: {}", mode).into()),
    };

    let mut store = RefgetStore::in_memory();
    store.set_encoding_mode(storage_mode);
    let ext_ptr = ExternalPtr::new(store);

    Ok(ext_ptr.into())
}

/// Create a disk-backed RefgetStore
/// @export
/// @param path Path for storing sequences and metadata
#[extendr]
pub fn on_disk_store(path: &str) -> extendr_api::Result<Robj> {
    match RefgetStore::on_disk(path) {
        Ok(store) => {
            let ext_ptr = ExternalPtr::new(store);
            Ok(ext_ptr.into())
        }
        Err(e) => Err(format!("Error creating disk-backed store: {}", e).into()),
    }
}

/// Open store from local directory
/// @export
/// @param root_path Path to read store from
#[extendr]
pub fn open_local_store(root_path: &str) -> extendr_api::Result<Robj> {
    match RefgetStore::open_local(root_path) {
        Ok(store) => {
            let ext_ptr = ExternalPtr::new(store);
            Ok(ext_ptr.into())
        }
        Err(e) => Err(format!("Error opening store from directory: {}", e).into()),
    }
}

/// Open a remote RefgetStore with local caching
/// @export
/// @param cache_path Local directory to cache downloaded metadata and sequences
/// @param remote_url Base URL of the remote refget store
#[extendr]
pub fn open_remote_store(cache_path: &str, remote_url: &str) -> extendr_api::Result<Robj> {
    match RefgetStore::open_remote(cache_path, remote_url.to_string()) {
        Ok(store) => {
            let ext_ptr = ExternalPtr::new(store);
            Ok(ext_ptr.into())
        }
        Err(e) => Err(format!("Error opening remote store: {}", e).into()),
    }
}

// =========================================================================
// Encoding Mode Controls
// =========================================================================

/// Enable 2-bit encoding for space efficiency
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn enable_encoding_store(store_ptr: Robj) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store.enable_encoding();
        Ok(())
    })
}

/// Disable encoding, use raw byte storage
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn disable_encoding_store(store_ptr: Robj) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store.disable_encoding();
        Ok(())
    })
}

/// Get current storage mode
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn get_storage_mode_store(store_ptr: Robj) -> extendr_api::Result<String> {
    with_store_ref!(store_ptr, store, {
        Ok(match store.storage_mode() {
            StorageMode::Raw => "raw".to_string(),
            StorageMode::Encoded => "encoded".to_string(),
        })
    })
}

// =========================================================================
// Quiet Mode
// =========================================================================

/// Set quiet mode
/// @param store_ptr External pointer to RefgetStore
/// @param quiet Whether to suppress output
#[extendr]
pub fn set_quiet_store(store_ptr: Robj, quiet: bool) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store.set_quiet(quiet);
        Ok(())
    })
}

/// Get quiet mode status
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn get_quiet_store(store_ptr: Robj) -> extendr_api::Result<bool> {
    with_store_ref!(store_ptr, store, { Ok(store.is_quiet()) })
}

// =========================================================================
// Persistence Controls
// =========================================================================

/// Check if store is persisting to disk
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn is_persisting_store(store_ptr: Robj) -> extendr_api::Result<bool> {
    with_store_ref!(store_ptr, store, { Ok(store.is_persisting()) })
}

/// Enable disk persistence
/// @param store_ptr External pointer to RefgetStore
/// @param path Directory for storing sequences and metadata
#[extendr]
pub fn enable_persistence_store(store_ptr: Robj, path: &str) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store
            .enable_persistence(path.to_string())
            .map_err(|e| format!("Error enabling persistence: {}", e).into())
    })
}

/// Disable disk persistence
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn disable_persistence_store(store_ptr: Robj) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store.disable_persistence();
        Ok(())
    })
}

// =========================================================================
// Store Path Accessors
// =========================================================================

/// Get cache path
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn get_cache_path_store(store_ptr: Robj) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        Ok(store
            .local_path()
            .map(|p| Robj::from(p.display().to_string()))
            .unwrap_or_else(|| Robj::from(())))
    })
}

/// Get remote URL
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn get_remote_url_store(store_ptr: Robj) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        Ok(store
            .remote_source()
            .map(|s| Robj::from(s.to_string()))
            .unwrap_or_else(|| Robj::from(())))
    })
}

// =========================================================================
// Adding Sequences and Collections
// =========================================================================

/// Import FASTA file into store
/// @param store_ptr External pointer to RefgetStore
/// @param file_path Path to FASTA file
/// @param force Whether to overwrite existing sequences/collections
#[extendr]
pub fn import_fasta_store(store_ptr: Robj, file_path: &str) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store
            .add_sequence_collection_from_fasta(file_path)
            .map(|_| ())
            .map_err(|e| format!("Error importing FASTA: {}", e).into())
    })
}

/// Add a FASTA file with force option
/// @param store_ptr External pointer to RefgetStore
/// @param file_path Path to FASTA file
/// @param force Whether to overwrite existing sequences/collections
#[extendr]
pub fn add_fasta_store(store_ptr: Robj, file_path: &str, force: bool) -> extendr_api::Result<List> {
    with_store!(store_ptr, store, {
        let result = if force {
            store.add_sequence_collection_from_fasta_force(file_path)
        } else {
            store.add_sequence_collection_from_fasta(file_path)
        };
        result
            .map(|(metadata, was_new)| {
                list!(
                    digest = metadata.digest,
                    n_sequences = metadata.n_sequences as i32,
                    names_digest = metadata.names_digest,
                    sequences_digest = metadata.sequences_digest,
                    lengths_digest = metadata.lengths_digest,
                    was_new = was_new
                )
            })
            .map_err(|e| format!("Error importing FASTA: {}", e).into())
    })
}

// =========================================================================
// Sequence Retrieval
// =========================================================================

/// Get sequence by digest from store (supports SHA512t24u or MD5)
/// @param store_ptr External pointer to RefgetStore
/// @param digest Sequence digest (SHA512t24u or MD5)
#[extendr]
pub fn get_sequence_store(store_ptr: Robj, digest: &str) -> extendr_api::Result<Robj> {
    let digest = strip_sq_prefix(digest);
    with_store!(store_ptr, store, {
        if let Ok(record) = store.get_sequence(digest) {
            Ok(record_to_list(record.clone()).into())
        } else {
            Ok(Robj::from(())) // NULL
        }
    })
}

/// Get sequence by collection and name
/// @param store_ptr External pointer to RefgetStore
/// @param collection_digest Sequence collection digest
/// @param sequence_name Sequence name
#[extendr]
pub fn get_sequence_by_name_store(
    store_ptr: Robj,
    collection_digest: &str,
    sequence_name: &str,
) -> extendr_api::Result<Robj> {
    let collection_digest = strip_sq_prefix(collection_digest);
    with_store!(store_ptr, store, {
        let result = store.get_sequence_by_name(collection_digest, sequence_name);
        if let Ok(record) = result {
            Ok(record_to_list(record.clone()).into())
        } else {
            Ok(Robj::from(())) // NULL
        }
    })
}

/// Get sequence metadata (no sequence data)
/// @param store_ptr External pointer to RefgetStore
/// @param digest Sequence digest
#[extendr]
pub fn get_sequence_metadata_store(store_ptr: Robj, digest: &str) -> extendr_api::Result<Robj> {
    let digest = strip_sq_prefix(digest);
    with_store_ref!(store_ptr, store, {
        Ok(store
            .get_sequence_metadata(digest.as_bytes())
            .map(|meta| metadata_to_list(meta.clone()).into())
            .unwrap_or_else(|| Robj::from(())))
    })
}

/// Get substring from sequence
/// @param store_ptr External pointer to RefgetStore
/// @param seq_digest Sequence digest
/// @param start Start position
/// @param end End position
#[extendr]
pub fn get_substring_store(
    store_ptr: Robj,
    seq_digest: &str,
    start: i32,
    end: i32,
) -> extendr_api::Result<Robj> {
    let seq_digest = strip_sq_prefix(seq_digest);
    with_store!(store_ptr, store, {
        if let Ok(substr) = store.get_substring(seq_digest, start as usize, end as usize) {
            Ok(Robj::from(substr))
        } else {
            Ok(Robj::from(())) // NULL
        }
    })
}

// =========================================================================
// Collection Retrieval
// =========================================================================

/// List all collections in the store
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn list_collections_store(store_ptr: Robj) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        let collections: Vec<Robj> = store
            .list_collections()
            .into_iter()
            .map(|meta| collection_metadata_to_list(meta).into())
            .collect();
        Ok(collections.into())
    })
}

/// List all sequences in the store
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn list_sequences_store(store_ptr: Robj) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        let sequences: Vec<Robj> = store
            .list_sequences()
            .into_iter()
            .map(|meta| metadata_to_list(meta).into())
            .collect();
        Ok(sequences.into())
    })
}

/// Get a collection with all sequences loaded
/// @param store_ptr External pointer to RefgetStore
/// @param digest Collection digest
#[extendr]
pub fn get_collection_store(store_ptr: Robj, digest: &str) -> extendr_api::Result<Robj> {
    with_store!(store_ptr, store, {
        store
            .get_collection(digest)
            .map(|coll| sequence_collection_to_list(coll).into())
            .map_err(|e| format!("Error loading collection: {}", e).into())
    })
}

/// Get collection metadata
/// @param store_ptr External pointer to RefgetStore
/// @param digest Collection digest
#[extendr]
pub fn get_collection_metadata_store(store_ptr: Robj, digest: &str) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        Ok(store
            .get_collection_metadata(digest)
            .map(|meta| collection_metadata_to_list(meta.clone()).into())
            .unwrap_or_else(|| Robj::from(())))
    })
}

/// Check if collection is fully loaded
/// @param store_ptr External pointer to RefgetStore
/// @param digest Collection digest
#[extendr]
pub fn is_collection_loaded_store(store_ptr: Robj, digest: &str) -> extendr_api::Result<bool> {
    with_store_ref!(store_ptr, store, { Ok(store.is_collection_loaded(digest)) })
}

/// Iterate over all collections with sequences loaded
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn iter_collections_store(store_ptr: Robj) -> extendr_api::Result<Robj> {
    with_store!(store_ptr, store, {
        let collections: Vec<Robj> = store
            .iter_collections()
            .map(|coll| sequence_collection_to_list(coll).into())
            .collect();
        Ok(collections.into())
    })
}

/// Iterate over all sequences with data loaded
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn iter_sequences_store(store_ptr: Robj) -> extendr_api::Result<Robj> {
    with_store!(store_ptr, store, {
        let sequences: Vec<Robj> = store
            .iter_sequences()
            .map(|rec| record_to_list(rec).into())
            .collect();
        Ok(sequences.into())
    })
}

/// Get store statistics
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn stats_store(store_ptr: Robj) -> extendr_api::Result<List> {
    with_store_ref!(store_ptr, store, {
        let stats = store.stats_extended();
        Ok(list!(
            n_sequences = stats.n_sequences as i32,
            n_sequences_loaded = stats.n_sequences_loaded as i32,
            n_collections = stats.n_collections as i32,
            n_collections_loaded = stats.n_collections_loaded as i32,
            storage_mode = stats.storage_mode,
            total_disk_size = stats.total_disk_size as i64
        ))
    })
}

// =========================================================================
// Seqcol Spec Operations
// =========================================================================

/// Get level 1 representation (attribute digests) for a collection
/// @param store_ptr External pointer to RefgetStore
/// @param digest Collection digest
#[extendr]
pub fn get_collection_level1_store(store_ptr: Robj, digest: &str) -> extendr_api::Result<List> {
    with_store_ref!(store_ptr, store, {
        store
            .get_collection_level1(digest)
            .map(|lvl1| {
                let mut result = list!(
                    names = lvl1.names,
                    lengths = lvl1.lengths,
                    sequences = lvl1.sequences
                );
                if let Some(ref v) = lvl1.name_length_pairs {
                    result.set_attrib("name_length_pairs", v.clone()).ok();
                }
                if let Some(ref v) = lvl1.sorted_name_length_pairs {
                    result
                        .set_attrib("sorted_name_length_pairs", v.clone())
                        .ok();
                }
                if let Some(ref v) = lvl1.sorted_sequences {
                    result.set_attrib("sorted_sequences", v.clone()).ok();
                }
                result
            })
            .map_err(|e| format!("Error getting level 1: {}", e).into())
    })
}

/// Get level 2 representation (full arrays) for a collection
/// @param store_ptr External pointer to RefgetStore
/// @param digest Collection digest
#[extendr]
pub fn get_collection_level2_store(store_ptr: Robj, digest: &str) -> extendr_api::Result<List> {
    with_store!(store_ptr, store, {
        store
            .get_collection_level2(digest)
            .map(|lvl2| {
                let lengths: Vec<i64> = lvl2.lengths.iter().map(|&l| l as i64).collect();
                list!(
                    names = lvl2.names,
                    lengths = lengths,
                    sequences = lvl2.sequences
                )
            })
            .map_err(|e| format!("Error getting level 2: {}", e).into())
    })
}

/// Compare two collections
/// @param store_ptr External pointer to RefgetStore
/// @param digest_a First collection digest
/// @param digest_b Second collection digest
#[extendr]
pub fn compare_store(
    store_ptr: Robj,
    digest_a: &str,
    digest_b: &str,
) -> extendr_api::Result<List> {
    with_store!(store_ptr, store, {
        store
            .compare(digest_a, digest_b)
            .map(|comparison| {
                let digests = list!(a = comparison.digests.a, b = comparison.digests.b);
                let attributes = list!(
                    a_only = comparison.attributes.a_only,
                    b_only = comparison.attributes.b_only,
                    a_and_b = comparison.attributes.a_and_b
                );

                // Convert HashMap to R list for array_elements
                let a_count: Vec<(String, i32)> = comparison
                    .array_elements
                    .a_count
                    .into_iter()
                    .map(|(k, v)| (k, v as i32))
                    .collect();
                let b_count: Vec<(String, i32)> = comparison
                    .array_elements
                    .b_count
                    .into_iter()
                    .map(|(k, v)| (k, v as i32))
                    .collect();
                let a_and_b_count: Vec<(String, i32)> = comparison
                    .array_elements
                    .a_and_b_count
                    .into_iter()
                    .map(|(k, v)| (k, v as i32))
                    .collect();

                list!(
                    digests = digests,
                    attributes = attributes,
                    a_count = hashmap_to_list(a_count),
                    b_count = hashmap_to_list(b_count),
                    a_and_b_count = hashmap_to_list(a_and_b_count)
                )
            })
            .map_err(|e| format!("Error comparing collections: {}", e).into())
    })
}

/// Find collections by attribute digest
/// @param store_ptr External pointer to RefgetStore
/// @param attr_name Attribute name (names, lengths, sequences, etc.)
/// @param attr_digest The digest to search for
#[extendr]
pub fn find_collections_by_attribute_store(
    store_ptr: Robj,
    attr_name: &str,
    attr_digest: &str,
) -> extendr_api::Result<Vec<String>> {
    with_store_ref!(store_ptr, store, {
        store
            .find_collections_by_attribute(attr_name, attr_digest)
            .map_err(|e| format!("Error finding collections: {}", e).into())
    })
}

/// Get attribute array by digest
/// @param store_ptr External pointer to RefgetStore
/// @param attr_name Attribute name
/// @param attr_digest Attribute digest
#[extendr]
pub fn get_attribute_store(
    store_ptr: Robj,
    attr_name: &str,
    attr_digest: &str,
) -> extendr_api::Result<Robj> {
    with_store!(store_ptr, store, {
        let result = store
            .get_attribute(attr_name, attr_digest)
            .map_err(|e| -> extendr_api::Error { format!("{}", e).into() })?;

        match result {
            None => Ok(Robj::from(())),
            Some(value) => {
                // Convert serde_json::Value to R
                json_value_to_robj(value)
            }
        }
    })
}

/// Enable ancillary digest computation
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn enable_ancillary_digests_store(store_ptr: Robj) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store.enable_ancillary_digests();
        Ok(())
    })
}

/// Disable ancillary digest computation
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn disable_ancillary_digests_store(store_ptr: Robj) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store.disable_ancillary_digests();
        Ok(())
    })
}

/// Check if ancillary digests are enabled
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn has_ancillary_digests_store(store_ptr: Robj) -> extendr_api::Result<bool> {
    with_store_ref!(store_ptr, store, { Ok(store.has_ancillary_digests()) })
}

/// Enable attribute index
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn enable_attribute_index_store(store_ptr: Robj) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store.enable_attribute_index();
        Ok(())
    })
}

/// Disable attribute index
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn disable_attribute_index_store(store_ptr: Robj) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store.disable_attribute_index();
        Ok(())
    })
}

/// Check if attribute index is enabled
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn has_attribute_index_store(store_ptr: Robj) -> extendr_api::Result<bool> {
    with_store_ref!(store_ptr, store, { Ok(store.has_attribute_index()) })
}

// =========================================================================
// Write/Export Operations
// =========================================================================

/// Write store to directory
/// @param store_ptr External pointer to RefgetStore
/// @param root_path Path to write store
/// @param seqdata_path_template Path template name (optional, pass empty string for default)
#[extendr]
pub fn write_store_to_directory_store(
    store_ptr: Robj,
    root_path: &str,
    seqdata_path_template: &str,
) -> extendr_api::Result<()> {
    with_store_ref!(store_ptr, store, {
        let template = if seqdata_path_template.is_empty() {
            None
        } else {
            Some(seqdata_path_template)
        };
        store
            .write_store_to_dir(root_path, template)
            .map_err(|e| format!("Error writing store to directory: {}", e).into())
    })
}

/// Write store using configured paths
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn write_store(store_ptr: Robj) -> extendr_api::Result<()> {
    with_store_ref!(store_ptr, store, {
        store
            .write()
            .map_err(|e| format!("Error writing store: {}", e).into())
    })
}

/// Export sequences from a collection to FASTA
/// @param store_ptr External pointer to RefgetStore
/// @param collection_digest Collection digest
/// @param output_path Output FASTA file path
/// @param sequence_names Optional vector of sequence names to export (NULL for all)
/// @param line_width Bases per line (default 80)
#[extendr]
pub fn export_fasta_store(
    store_ptr: Robj,
    collection_digest: &str,
    output_path: &str,
    sequence_names: Robj,
    line_width: i32,
) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        let names: Option<Vec<String>> = if sequence_names.is_null() {
            None
        } else {
            Some(sequence_names.as_str_vector().unwrap_or_default().iter().map(|s| s.to_string()).collect())
        };
        let names_refs = names
            .as_ref()
            .map(|v| v.iter().map(|s| s.as_str()).collect::<Vec<&str>>());
        let width = if line_width <= 0 {
            None
        } else {
            Some(line_width as usize)
        };
        store
            .export_fasta(collection_digest, output_path, names_refs, width)
            .map_err(|e| format!("Error exporting FASTA: {}", e).into())
    })
}

/// Export sequences by their digests to FASTA
/// @param store_ptr External pointer to RefgetStore
/// @param seq_digests Vector of sequence digests
/// @param output_path Output FASTA file path
/// @param line_width Bases per line (default 80)
#[extendr]
pub fn export_fasta_by_digests_store(
    store_ptr: Robj,
    seq_digests: Vec<String>,
    output_path: &str,
    line_width: i32,
) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        let digests_refs: Vec<&str> = seq_digests.iter().map(|s| s.as_str()).collect();
        let width = if line_width <= 0 {
            None
        } else {
            Some(line_width as usize)
        };
        store
            .export_fasta_by_digests(digests_refs, output_path, width)
            .map_err(|e| format!("Error exporting FASTA: {}", e).into())
    })
}

/// Extract BED file sequences from store as FASTA
/// @param store_ptr External pointer to RefgetStore
/// @param collection_digest Sequence collection digest
/// @param bed_file_path Path to BED file
/// @param output_file_path Path to output FASTA file
#[extendr]
pub fn get_seqs_bed_file_store(
    store_ptr: Robj,
    collection_digest: &str,
    bed_file_path: &str,
    output_file_path: &str,
) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store
            .export_fasta_from_regions(collection_digest, bed_file_path, output_file_path)
            .map_err(|e| format!("Error writing sequences to file: {}", e).into())
    })
}

/// Extract BED file sequences from store into memory
/// @param store_ptr External pointer to RefgetStore
/// @param collection_digest Sequence collection digest
/// @param bed_file_path Path to BED file
#[extendr]
pub fn get_seqs_bed_file_to_vec_store(
    store_ptr: Robj,
    collection_digest: &str,
    bed_file_path: &str,
) -> extendr_api::Result<Robj> {
    with_store!(store_ptr, store, {
        match store.substrings_from_regions(collection_digest, bed_file_path) {
            Ok(iter) => {
                let mut r_results: Vec<Robj> = Vec::new();
                for result in iter {
                    match result {
                        Ok(retrieved_seq) => {
                            r_results.push(retrieved_sequence_to_list(retrieved_seq).into());
                        }
                        Err(e) => {
                            eprintln!("Warning: {}", e);
                        }
                    }
                }
                Ok(r_results.into())
            }
            Err(e) => Err(format!("Error retrieving sequences from BED file: {}", e).into()),
        }
    })
}

// =========================================================================
// Sequence Alias Operations
// =========================================================================

/// Add a sequence alias
/// @param store_ptr External pointer to RefgetStore
/// @param namespace Alias namespace
/// @param alias Alias name
/// @param digest Sequence digest
#[extendr]
pub fn add_sequence_alias_store(
    store_ptr: Robj,
    namespace: &str,
    alias: &str,
    digest: &str,
) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store
            .add_sequence_alias(namespace, alias, digest)
            .map_err(|e| format!("Error adding sequence alias: {}", e).into())
    })
}

/// Get sequence by alias
/// @param store_ptr External pointer to RefgetStore
/// @param namespace Alias namespace
/// @param alias Alias name
#[extendr]
pub fn get_sequence_by_alias_store(
    store_ptr: Robj,
    namespace: &str,
    alias: &str,
) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        Ok(store
            .get_sequence_by_alias(namespace, alias)
            .map(|r| record_to_list(r.clone()).into())
            .unwrap_or_else(|| Robj::from(())))
    })
}

/// Get all aliases for a sequence
/// @param store_ptr External pointer to RefgetStore
/// @param digest Sequence digest
#[extendr]
pub fn get_aliases_for_sequence_store(store_ptr: Robj, digest: &str) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        let aliases: Vec<Robj> = store
            .get_aliases_for_sequence(digest)
            .into_iter()
            .map(|(ns, alias)| list!(namespace = ns, alias = alias).into())
            .collect();
        Ok(aliases.into())
    })
}

/// List all sequence alias namespaces
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn list_sequence_alias_namespaces_store(
    store_ptr: Robj,
) -> extendr_api::Result<Vec<String>> {
    with_store_ref!(store_ptr, store, {
        Ok(store.list_sequence_alias_namespaces())
    })
}

/// List all aliases in a sequence namespace
/// @param store_ptr External pointer to RefgetStore
/// @param namespace Namespace name
#[extendr]
pub fn list_sequence_aliases_store(
    store_ptr: Robj,
    namespace: &str,
) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        Ok(store
            .list_sequence_aliases(namespace)
            .map(|v| Robj::from(v))
            .unwrap_or_else(|| Robj::from(())))
    })
}

/// Remove a sequence alias
/// @param store_ptr External pointer to RefgetStore
/// @param namespace Namespace name
/// @param alias Alias name
#[extendr]
pub fn remove_sequence_alias_store(
    store_ptr: Robj,
    namespace: &str,
    alias: &str,
) -> extendr_api::Result<bool> {
    with_store!(store_ptr, store, {
        store
            .remove_sequence_alias(namespace, alias)
            .map_err(|e| format!("Error removing sequence alias: {}", e).into())
    })
}

/// Load sequence aliases from TSV file
/// @param store_ptr External pointer to RefgetStore
/// @param namespace Namespace to load into
/// @param path Path to TSV file
#[extendr]
pub fn load_sequence_aliases_store(
    store_ptr: Robj,
    namespace: &str,
    path: &str,
) -> extendr_api::Result<i32> {
    with_store!(store_ptr, store, {
        store
            .load_sequence_aliases(namespace, path)
            .map(|n| n as i32)
            .map_err(|e| format!("Error loading sequence aliases: {}", e).into())
    })
}

// =========================================================================
// Collection Alias Operations
// =========================================================================

/// Add a collection alias
/// @param store_ptr External pointer to RefgetStore
/// @param namespace Alias namespace
/// @param alias Alias name
/// @param digest Collection digest
#[extendr]
pub fn add_collection_alias_store(
    store_ptr: Robj,
    namespace: &str,
    alias: &str,
    digest: &str,
) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store
            .add_collection_alias(namespace, alias, digest)
            .map_err(|e| format!("Error adding collection alias: {}", e).into())
    })
}

/// Get collection by alias
/// @param store_ptr External pointer to RefgetStore
/// @param namespace Alias namespace
/// @param alias Alias name
#[extendr]
pub fn get_collection_by_alias_store(
    store_ptr: Robj,
    namespace: &str,
    alias: &str,
) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        Ok(store
            .get_collection_by_alias(namespace, alias)
            .map(|r| collection_metadata_to_list(r.metadata().clone()).into())
            .unwrap_or_else(|| Robj::from(())))
    })
}

/// Get all aliases for a collection
/// @param store_ptr External pointer to RefgetStore
/// @param digest Collection digest
#[extendr]
pub fn get_aliases_for_collection_store(
    store_ptr: Robj,
    digest: &str,
) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        let aliases: Vec<Robj> = store
            .get_aliases_for_collection(digest)
            .into_iter()
            .map(|(ns, alias)| list!(namespace = ns, alias = alias).into())
            .collect();
        Ok(aliases.into())
    })
}

/// List all collection alias namespaces
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn list_collection_alias_namespaces_store(
    store_ptr: Robj,
) -> extendr_api::Result<Vec<String>> {
    with_store_ref!(store_ptr, store, {
        Ok(store.list_collection_alias_namespaces())
    })
}

/// List all aliases in a collection namespace
/// @param store_ptr External pointer to RefgetStore
/// @param namespace Namespace name
#[extendr]
pub fn list_collection_aliases_store(
    store_ptr: Robj,
    namespace: &str,
) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        Ok(store
            .list_collection_aliases(namespace)
            .map(|v| Robj::from(v))
            .unwrap_or_else(|| Robj::from(())))
    })
}

/// Remove a collection alias
/// @param store_ptr External pointer to RefgetStore
/// @param namespace Namespace name
/// @param alias Alias name
#[extendr]
pub fn remove_collection_alias_store(
    store_ptr: Robj,
    namespace: &str,
    alias: &str,
) -> extendr_api::Result<bool> {
    with_store!(store_ptr, store, {
        store
            .remove_collection_alias(namespace, alias)
            .map_err(|e| format!("Error removing collection alias: {}", e).into())
    })
}

/// Load collection aliases from TSV file
/// @param store_ptr External pointer to RefgetStore
/// @param namespace Namespace to load into
/// @param path Path to TSV file
#[extendr]
pub fn load_collection_aliases_store(
    store_ptr: Robj,
    namespace: &str,
    path: &str,
) -> extendr_api::Result<i32> {
    with_store!(store_ptr, store, {
        store
            .load_collection_aliases(namespace, path)
            .map(|n| n as i32)
            .map_err(|e| format!("Error loading collection aliases: {}", e).into())
    })
}

// =========================================================================
// FHR Metadata Operations
// =========================================================================

/// Set FHR metadata for a collection
/// @param store_ptr External pointer to RefgetStore
/// @param collection_digest Collection digest
/// @param metadata_json JSON string with FHR metadata
#[extendr]
pub fn set_fhr_metadata_store(
    store_ptr: Robj,
    collection_digest: &str,
    metadata_json: &str,
) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        let metadata: FhrMetadata = serde_json::from_str(metadata_json)
            .map_err(|e| -> extendr_api::Error { format!("Invalid JSON: {}", e).into() })?;
        store
            .set_fhr_metadata(collection_digest, metadata)
            .map_err(|e| format!("Error setting FHR metadata: {}", e).into())
    })
}

/// Get FHR metadata for a collection
/// @param store_ptr External pointer to RefgetStore
/// @param collection_digest Collection digest
#[extendr]
pub fn get_fhr_metadata_store(
    store_ptr: Robj,
    collection_digest: &str,
) -> extendr_api::Result<Robj> {
    with_store_ref!(store_ptr, store, {
        Ok(store
            .get_fhr_metadata(collection_digest)
            .map(|fhr| {
                serde_json::to_string(fhr)
                    .map(|s| Robj::from(s))
                    .unwrap_or_else(|_| Robj::from(()))
            })
            .unwrap_or_else(|| Robj::from(())))
    })
}

/// Remove FHR metadata for a collection
/// @param store_ptr External pointer to RefgetStore
/// @param collection_digest Collection digest
#[extendr]
pub fn remove_fhr_metadata_store(store_ptr: Robj, collection_digest: &str) -> extendr_api::Result<bool> {
    with_store!(store_ptr, store, {
        Ok(store.remove_fhr_metadata(collection_digest))
    })
}

/// List all collection digests with FHR metadata
/// @param store_ptr External pointer to RefgetStore
#[extendr]
pub fn list_fhr_metadata_store(store_ptr: Robj) -> extendr_api::Result<Vec<String>> {
    with_store_ref!(store_ptr, store, { Ok(store.list_fhr_metadata()) })
}

/// Load FHR metadata from JSON file
/// @param store_ptr External pointer to RefgetStore
/// @param collection_digest Collection digest
/// @param path Path to JSON file
#[extendr]
pub fn load_fhr_metadata_store(
    store_ptr: Robj,
    collection_digest: &str,
    path: &str,
) -> extendr_api::Result<()> {
    with_store!(store_ptr, store, {
        store
            .load_fhr_metadata(collection_digest, path)
            .map_err(|e| format!("Error loading FHR metadata: {}", e).into())
    })
}

// =========================================================================
// Helper Functions
// =========================================================================

/// Strip "SQ." prefix from digest if present
fn strip_sq_prefix(digest: &str) -> &str {
    if digest.len() > 3 {
        let prefix = &digest[..3];
        if prefix.eq_ignore_ascii_case("SQ.") {
            return &digest[3..];
        }
    }
    digest
}

fn alphabet_to_string(alphabet: AlphabetType) -> &'static str {
    match alphabet {
        AlphabetType::Dna2bit => "dna2bit",
        AlphabetType::Dna3bit => "dna3bit",
        AlphabetType::DnaIupac => "dnaio",
        AlphabetType::Protein => "protein",
        AlphabetType::Ascii => "ASCII",
        AlphabetType::Unknown => "Unknown",
    }
}

fn metadata_to_list(metadata: SequenceMetadata) -> List {
    list!(
        name = metadata.name,
        length = metadata.length as i32,
        sha512t24u = metadata.sha512t24u,
        md5 = metadata.md5,
        alphabet = alphabet_to_string(metadata.alphabet),
        description = metadata.description.unwrap_or_default()
    )
}

fn record_to_list(record: SequenceRecord) -> List {
    list!(
        metadata = metadata_to_list(record.metadata().clone()),
        data = record.decode().unwrap_or_default()
    )
}

fn lvl1_to_list(lvl1: SeqColDigestLvl1) -> List {
    list!(
        sequences_digest = lvl1.sequences_digest,
        names_digest = lvl1.names_digest,
        lengths_digest = lvl1.lengths_digest
    )
}

fn collection_metadata_to_list(metadata: SequenceCollectionMetadata) -> List {
    list!(
        digest = metadata.digest,
        n_sequences = metadata.n_sequences as i32,
        names_digest = metadata.names_digest,
        sequences_digest = metadata.sequences_digest,
        lengths_digest = metadata.lengths_digest,
        name_length_pairs_digest = metadata.name_length_pairs_digest.unwrap_or_default(),
        sorted_name_length_pairs_digest = metadata.sorted_name_length_pairs_digest.unwrap_or_default(),
        sorted_sequences_digest = metadata.sorted_sequences_digest.unwrap_or_default()
    )
}

fn sequence_collection_to_list(collection: SequenceCollection) -> List {
    let sequences: Vec<Robj> = collection
        .sequences
        .into_iter()
        .map(|seq_record| record_to_list(seq_record).into())
        .collect();

    let lvl1 = lvl1_to_list(collection.metadata.to_lvl1());
    list!(
        sequences = sequences,
        digest = collection.metadata.digest,
        lvl1 = lvl1,
        file_path = collection
            .metadata
            .file_path
            .map(|p| p.to_string_lossy().to_string())
            .unwrap_or_default()
    )
}

fn retrieved_sequence_to_list(retrieved_sequence: RetrievedSequence) -> List {
    list!(
        sequence = retrieved_sequence.sequence,
        chrom_name = retrieved_sequence.chrom_name,
        start = retrieved_sequence.start as i32,
        end = retrieved_sequence.end as i32
    )
}

fn hashmap_to_list(items: Vec<(String, i32)>) -> List {
    let keys: Vec<&str> = items.iter().map(|(k, _)| k.as_str()).collect();
    let mut named_list = List::new(items.len());
    for (i, (_, val)) in items.iter().enumerate() {
        named_list.set_elt(i, Robj::from(*val)).ok();
    }
    let _ = named_list.set_names(keys);
    named_list
}

fn json_value_to_robj(value: serde_json::Value) -> extendr_api::Result<Robj> {
    match value {
        serde_json::Value::Null => Ok(Robj::from(())),
        serde_json::Value::Bool(b) => Ok(Robj::from(b)),
        serde_json::Value::Number(n) => {
            if let Some(i) = n.as_i64() {
                Ok(Robj::from(i as i32))
            } else if let Some(f) = n.as_f64() {
                Ok(Robj::from(f))
            } else {
                Ok(Robj::from(n.to_string()))
            }
        }
        serde_json::Value::String(s) => Ok(Robj::from(s)),
        serde_json::Value::Array(arr) => {
            // Try to create a typed vector
            let first = arr.first();
            match first {
                Some(serde_json::Value::Number(_)) => {
                    // Try integer first
                    let ints: Option<Vec<i32>> = arr
                        .iter()
                        .map(|v| v.as_i64().map(|i| i as i32))
                        .collect();
                    if let Some(v) = ints {
                        return Ok(Robj::from(v));
                    }
                    // Fall back to floats
                    let floats: Option<Vec<f64>> = arr.iter().map(|v| v.as_f64()).collect();
                    if let Some(v) = floats {
                        return Ok(Robj::from(v));
                    }
                }
                Some(serde_json::Value::String(_)) => {
                    let strings: Option<Vec<String>> =
                        arr.iter().map(|v| v.as_str().map(|s| s.to_string())).collect();
                    if let Some(v) = strings {
                        return Ok(Robj::from(v));
                    }
                }
                _ => {}
            }
            // Fall back to list
            let items: extendr_api::Result<Vec<Robj>> =
                arr.into_iter().map(json_value_to_robj).collect();
            Ok(items?.into())
        }
        serde_json::Value::Object(obj) => {
            let mut list = List::new(obj.len());
            let keys: Vec<String> = obj.keys().cloned().collect();
            for (i, (_, v)) in obj.into_iter().enumerate() {
                list.set_elt(i, json_value_to_robj(v)?).ok();
            }
            let key_refs: Vec<&str> = keys.iter().map(|s| s.as_str()).collect();
            let _ = list.set_names(key_refs);
            Ok(list.into())
        }
    }
}

// =========================================================================
// Module Registration
// =========================================================================

extendr_module! {
    mod refget;

    // Digest functions
    fn sha512t24u_digest;
    fn md5_digest;
    fn digest_fasta_raw;

    // Store constructors
    fn refget_store_raw;
    fn on_disk_store;
    fn open_local_store;
    fn open_remote_store;

    // Encoding mode
    fn enable_encoding_store;
    fn disable_encoding_store;
    fn get_storage_mode_store;

    // Quiet mode
    fn set_quiet_store;
    fn get_quiet_store;

    // Persistence
    fn is_persisting_store;
    fn enable_persistence_store;
    fn disable_persistence_store;

    // Path accessors
    fn get_cache_path_store;
    fn get_remote_url_store;

    // Adding sequences/collections
    fn import_fasta_store;
    fn add_fasta_store;

    // Sequence retrieval
    fn get_sequence_store;
    fn get_sequence_by_name_store;
    fn get_sequence_metadata_store;
    fn get_substring_store;

    // Collection retrieval
    fn list_collections_store;
    fn list_sequences_store;
    fn get_collection_store;
    fn get_collection_metadata_store;
    fn is_collection_loaded_store;
    fn iter_collections_store;
    fn iter_sequences_store;
    fn stats_store;

    // Seqcol spec operations
    fn get_collection_level1_store;
    fn get_collection_level2_store;
    fn compare_store;
    fn find_collections_by_attribute_store;
    fn get_attribute_store;
    fn enable_ancillary_digests_store;
    fn disable_ancillary_digests_store;
    fn has_ancillary_digests_store;
    fn enable_attribute_index_store;
    fn disable_attribute_index_store;
    fn has_attribute_index_store;

    // Write/export
    fn write_store_to_directory_store;
    fn write_store;
    fn export_fasta_store;
    fn export_fasta_by_digests_store;
    fn get_seqs_bed_file_store;
    fn get_seqs_bed_file_to_vec_store;

    // Sequence aliases
    fn add_sequence_alias_store;
    fn get_sequence_by_alias_store;
    fn get_aliases_for_sequence_store;
    fn list_sequence_alias_namespaces_store;
    fn list_sequence_aliases_store;
    fn remove_sequence_alias_store;
    fn load_sequence_aliases_store;

    // Collection aliases
    fn add_collection_alias_store;
    fn get_collection_by_alias_store;
    fn get_aliases_for_collection_store;
    fn list_collection_alias_namespaces_store;
    fn list_collection_aliases_store;
    fn remove_collection_alias_store;
    fn load_collection_aliases_store;

    // FHR metadata
    fn set_fhr_metadata_store;
    fn get_fhr_metadata_store;
    fn remove_fhr_metadata_store;
    fn list_fhr_metadata_store;
    fn load_fhr_metadata_store;
}
