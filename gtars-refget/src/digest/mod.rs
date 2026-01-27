//! Core digest and encoding functionality - WASM-safe.
//!
//! This module contains all the WASM-compatible code for computing
//! sequence digests. No filesystem dependencies.
//!
//! # Submodules
//!
//! - `algorithms` - Core hash functions (sha512t24u, md5, canonicalize_json)
//! - `alphabet` - Alphabet types and encoding/decoding arrays
//! - `auto_decompress` - Auto-decompressing writer (push-based gzip handling)
//! - `encoder` - Sequence bit-packing and decoding
//! - `types` - Data structures (SequenceRecord, SequenceCollection, etc.)
//! - `fasta` - Bytes-based FASTA parsing (WASM-compatible)
//! - `stream` - Streaming FASTA hasher for chunk-by-chunk processing

pub mod algorithms;
pub mod alphabet;
pub mod auto_decompress;
pub mod encoder;
pub mod fasta;
pub mod stream;
pub mod types;

// Re-export commonly used items at the module level
pub use algorithms::{canonicalize_json, md5, sha512t24u};
pub use alphabet::{
    guess_alphabet, lookup_alphabet, Alphabet, AlphabetGuesser, AlphabetType,
    ASCII_ALPHABET, DNA_2BIT_ALPHABET, DNA_3BIT_ALPHABET, DNA_IUPAC_ALPHABET, PROTEIN_ALPHABET,
};
pub use encoder::{
    decode_string_from_bytes, decode_substring_from_bytes, encode_sequence, SequenceEncoder,
};
pub use fasta::{digest_fasta_bytes, load_fasta_bytes, parse_fasta_header, ParseOptions};
pub use auto_decompress::AutoDecompressWriter;
pub use stream::FastaStreamHasher;
pub use types::{
    digest_sequence, digest_sequence_with_description, parse_rgsi_line,
    FaiMetadata, SeqColDigestLvl1, SequenceCollection, SequenceCollectionMetadata,
    SequenceCollectionRecord, SequenceMetadata, SequenceRecord,
};
