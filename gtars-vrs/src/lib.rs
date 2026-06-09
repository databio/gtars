//! # GA4GH VRS Allele Digest Computation
//!
//! This crate computes GA4GH VRS (Variation Representation Specification) identifiers
//! for genomic variants. It provides:
//!
//! - VRS data models (Allele, SequenceLocation, etc.)
//! - VRS digest computation (canonical JSON serialization + SHA-512/24u)
//! - Allele normalization (port of bioutils `normalize()`)
//! - VCF batch processing (parse VCF → compute VRS IDs)

pub mod digest;
pub mod hgvs;
pub mod models;
pub mod normalize;
pub mod provider;
// WASM-safe VCF→VRS primitives (parsing + the reader-generic readonly streaming
// core). No filesystem/threading deps, so it compiles for the wasm build too.
pub mod vcf_core;
// The full VCF pipeline uses BGZF/threads/filesystem and the mutable
// `RefgetStore`; it is not part of the WASM surface.
#[cfg(feature = "filesystem")]
pub mod vcf;

pub use models::{Allele, AlleleState, SequenceLocation, SequenceReference};
pub use digest::{allele_digest, allele_identifier, sequence_location_digest};
pub use normalize::normalize;
// Readonly HGVS->VRS entries are WASM-safe; the mutable `&mut RefgetStore`
// entries depend on the filesystem-backed VCF helpers.
pub use hgvs::bridge::{BridgeError, hgvs_str_to_vrs_id_readonly, hgvs_to_allele_readonly};
#[cfg(feature = "filesystem")]
pub use hgvs::bridge::{hgvs_str_to_vrs_id, hgvs_to_allele};
pub use provider::{NoTranscriptProvider, ProviderError, TranscriptProvider};
// WASM-safe VCF surface: parsing + the reader-generic readonly streaming core.
pub use vcf_core::{
    ParsedRecord, VrsResult, compute_vrs_ids_streaming_readonly_from_reader, is_real_alt,
    parse_vcf_record,
};
#[cfg(feature = "filesystem")]
pub use vcf::{
    compute_vrs_ids_from_vcf, compute_vrs_ids_from_vcf_readonly,
    compute_vrs_ids_parallel_encoded, compute_vrs_ids_streaming,
    compute_vrs_ids_streaming_readonly,
};
