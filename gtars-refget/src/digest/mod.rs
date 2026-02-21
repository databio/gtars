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
    ASCII_ALPHABET, Alphabet, AlphabetGuesser, AlphabetType, DNA_2BIT_ALPHABET, DNA_3BIT_ALPHABET,
    DNA_IUPAC_ALPHABET, PROTEIN_ALPHABET, guess_alphabet, lookup_alphabet,
};
pub use auto_decompress::AutoDecompressWriter;
pub use encoder::{
    SequenceEncoder, decode_string_from_bytes, decode_substring_from_bytes, encode_sequence,
};
pub use fasta::{ParseOptions, digest_fasta_bytes, load_fasta_bytes, parse_fasta_header};
pub use stream::FastaStreamHasher;
pub use types::{
    ArrayElementComparison, AttributeComparison, CollectionLevel1, CollectionLevel2,
    ComparisonDigests, FaiMetadata, SeqColComparison, SeqColDigestLvl1, SequenceCollection,
    SequenceCollectionMetadata, SequenceCollectionRecord, SequenceMetadata, SequenceRecord,
    digest_sequence, digest_sequence_with_description, parse_rgci_line, parse_rgsi_line,
};
