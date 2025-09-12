use extendr_api::prelude::*;
use refget::digest::{md5, sha512t24u};
use refget::collection::{SequenceCollection, SequenceMetadata, SequenceRecord, SeqColDigestLvl1};
use refget::alphabet::AlphabetType;

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
    match refget::fasta::digest_fasta(fasta) {
        Ok(sequence_collection) => {
            Ok(sequence_collection_to_list(sequence_collection))
        }
        Err(e) => Err(format!("Error processing FASTA file: {}", e).into())
    }
}

fn alphabet_to_string(alphabet: AlphabetType) -> &'static str {
    match alphabet {
        AlphabetType::Dna2bit => "dna2bit",     // 2-bit DNA encoding
        AlphabetType::Dna3bit => "dna3bit",     // 3-bit DNA encoding  
        AlphabetType::DnaIupac => "dnaio",      // IUPAC DNA (includes ambiguous bases)
        AlphabetType::Protein => "protein",     // Amino acid sequences
        AlphabetType::Ascii => "ASCII",         // Plain ASCII text
        AlphabetType::Unknown => "Unknown",     // Unrecognized format
    }
}

fn metadata_to_list(metadata: SequenceMetadata) -> List {
    list!(
        name = metadata.name,
        length = metadata.length as i32,
        sha512t24u = metadata.sha512t24u,
        md5 = metadata.md5,
        alphabet = alphabet_to_string(metadata.alphabet)  // Convert enum to string
    )
}

fn record_to_list(record: SequenceRecord) -> List {
    list!(
        metadata = metadata_to_list(record.metadata),
        data = record.data.unwrap_or(Vec::new())
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
    let sequences: Vec<Robj> = collection.sequences
        .into_iter()
        .map(|seq_record| record_to_list(seq_record).into())
        .collect();

    list!(
        sequences = sequences,
        digest = collection.digest,
        lvl1 = lvl1_to_list(collection.lvl1),
        file_path = collection.file_path
            .map(|p| p.to_string_lossy().to_string())
            .unwrap_or_else(|| String::new()),
        has_data = collection.has_data
    )
}

extendr_module! {
    mod refget;
    fn sha512t24u_digest;
    fn md5_digest;
    fn digest_fasta_raw;
}
