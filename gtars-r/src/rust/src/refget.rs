use extendr_api::prelude::*;
use gtars_refget::alphabet::AlphabetType;
use gtars_refget::collection::{
    SeqColDigestLvl1, SequenceCollection, SequenceMetadata, SequenceRecord,
};
use gtars_refget::digest::{md5, sha512t24u};
use gtars_refget::store::{RefgetStore, RetrievedSequence, StorageMode};

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

/// Create a RefgetStore
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
    let store = Box::new(store);
    let ptr = unsafe { Robj::make_external_ptr(Box::into_raw(store), Robj::from(())) };

    Ok(ptr)
}

/// Import FASTA file into store
/// @param store_ptr External pointer to RefgetStore
/// @param file_path Path to FASTA file
#[extendr]
pub fn import_fasta_store(store_ptr: Robj, file_path: &str) -> extendr_api::Result<()> {
    let store_raw_ptr = unsafe { store_ptr.external_ptr_addr::<RefgetStore>() };
    if store_raw_ptr.is_null() {
        return Err("Invalid store pointer".into());
    }

    let store = unsafe { &mut *store_raw_ptr };

    store
        .add_sequence_collection_from_fasta(file_path)
        .map(|_| ())
        .map_err(|e| format!("Error importing FASTA: {}", e).into())
}

/// Get sequence by digest from store (supports SHA512t24u or MD5)
/// @param store_ptr External pointer to RefgetStore
/// @param digest Sequence digest (SHA512t24u or MD5)
#[extendr]
pub fn get_sequence_store(store_ptr: Robj, digest: &str) -> extendr_api::Result<Robj> {
    let store_raw_ptr = unsafe { store_ptr.external_ptr_addr::<RefgetStore>() };
    if store_raw_ptr.is_null() {
        return Err("Invalid store pointer".into());
    }

    let store = unsafe { &mut *store_raw_ptr };

    // get_sequence handles both SHA512t24u and MD5 lookups
    if let Ok(record) = store.get_sequence(digest) {
        Ok(record_to_list(record.clone()).into())
    } else {
        Ok(Robj::from(())) // NULL
    }
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
    let store_raw_ptr = unsafe { store_ptr.external_ptr_addr::<RefgetStore>() };
    if store_raw_ptr.is_null() {
        return Err("Invalid store pointer".into());
    }

    let store = unsafe { &mut *store_raw_ptr };

    let result = store.get_sequence_by_name(collection_digest, sequence_name);

    if let Ok(record) = result {
        Ok(record_to_list(record.clone()).into())
    } else {
        Ok(Robj::from(())) // NULL
    }
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
    let store_raw_ptr = unsafe { store_ptr.external_ptr_addr::<RefgetStore>() };
    if store_raw_ptr.is_null() {
        return Err("Invalid store pointer".into());
    }

    let store = unsafe { &mut *store_raw_ptr };

    if let Ok(substr) = store.get_substring(seq_digest, start as usize, end as usize) {
        Ok(Robj::from(substr))
    } else {
        Ok(Robj::from(())) // NULL
    }
}

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
    let store_raw_ptr = unsafe { store_ptr.external_ptr_addr::<RefgetStore>() };
    if store_raw_ptr.is_null() {
        return Err("Invalid store pointer".into());
    }

    let store = unsafe { &*store_raw_ptr };

    let template = if seqdata_path_template.is_empty() {
        None
    } else {
        Some(seqdata_path_template)
    };

    store
        .write_store_to_dir(root_path, template)
        .map_err(|e| format!("Error writing store to directory: {}", e).into())
}

/// Open store from directory
/// @export
/// @param root_path Path to read store from
#[extendr]
pub fn open_local_store(root_path: &str) -> extendr_api::Result<Robj> {
    match RefgetStore::open_local(root_path) {
        Ok(store) => {
            let boxed_store = Box::new(store);
            let ptr =
                unsafe { Robj::make_external_ptr(Box::into_raw(boxed_store), Robj::from(())) };
            Ok(ptr)
        }
        Err(e) => Err(format!("Error opening store from directory: {}", e).into()),
    }
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
    let store_raw_ptr = unsafe { store_ptr.external_ptr_addr::<RefgetStore>() };
    if store_raw_ptr.is_null() {
        return Err("Invalid store pointer".into());
    }

    let store = unsafe { &mut *store_raw_ptr };

    store
        .export_fasta_from_regions(collection_digest, bed_file_path, output_file_path)
        .map_err(|e| format!("Error writing sequences to file: {}", e).into())
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
    let store_raw_ptr = unsafe { store_ptr.external_ptr_addr::<RefgetStore>() };
    if store_raw_ptr.is_null() {
        return Err("Invalid store pointer".into());
    }

    let store = unsafe { &mut *store_raw_ptr };

    match store.substrings_from_regions(collection_digest, bed_file_path) {
        Ok(iter) => {
            let mut r_results: Vec<Robj> = Vec::new();
            for result in iter {
                match result {
                    Ok(retrieved_seq) => {
                        r_results.push(retrieved_sequence_to_list(retrieved_seq).into());
                    }
                    Err(e) => {
                        return Err(format!("Error retrieving sequence: {}", e).into());
                    }
                }
            }
            Ok(r_results.into())
        }
        Err(e) => Err(format!("Error retrieving sequences from BED file: {}", e).into()),
    }
}

fn alphabet_to_string(alphabet: AlphabetType) -> &'static str {
    match alphabet {
        AlphabetType::Dna2bit => "dna2bit", // 2-bit DNA encoding
        AlphabetType::Dna3bit => "dna3bit", // 3-bit DNA encoding
        AlphabetType::DnaIupac => "dnaio",  // IUPAC DNA (includes ambiguous bases)
        AlphabetType::Protein => "protein", // Amino acid sequences
        AlphabetType::Ascii => "ASCII",     // Plain ASCII text
        AlphabetType::Unknown => "Unknown", // Unrecognized format
    }
}

fn metadata_to_list(metadata: SequenceMetadata) -> List {
    list!(
        name = metadata.name,
        length = metadata.length as i32,
        sha512t24u = metadata.sha512t24u,
        md5 = metadata.md5,
        alphabet = alphabet_to_string(metadata.alphabet) // Convert enum to string
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

fn sequence_collection_to_list(collection: SequenceCollection) -> List {
    let sequences: Vec<Robj> = collection
        .sequences
        .into_iter()
        .map(|seq_record| record_to_list(seq_record).into())
        .collect();

    list!(
        sequences = sequences,
        digest = collection.metadata.digest,
        lvl1 = lvl1_to_list(collection.metadata.to_lvl1()),
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

extendr_module! {
    mod refget;
    fn sha512t24u_digest;
    fn md5_digest;
    fn digest_fasta_raw;
    fn refget_store_raw;
    fn import_fasta_store;
    fn get_sequence_store;
    fn get_sequence_by_name_store;
    fn get_substring_store;
    fn write_store_to_directory_store;
    fn open_local_store;
    fn get_seqs_bed_file_store;
    fn get_seqs_bed_file_to_vec_store;
}
