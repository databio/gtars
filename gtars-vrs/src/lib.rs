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
// The VCF pipeline uses BGZF/threads/filesystem and the mutable `RefgetStore`;
// it is not part of the WASM surface.
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
#[cfg(feature = "filesystem")]
pub use vcf::{
    VrsResult, compute_vrs_ids_from_vcf, compute_vrs_ids_from_vcf_readonly,
    compute_vrs_ids_parallel_encoded, compute_vrs_ids_streaming,
    compute_vrs_ids_streaming_readonly,
};
