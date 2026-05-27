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
pub mod vcf;

pub use models::{Allele, AlleleState, SequenceLocation, SequenceReference};
pub use digest::{allele_digest, allele_identifier, sequence_location_digest};
pub use normalize::normalize;
pub use hgvs::bridge::{
    BridgeError, hgvs_str_to_vrs_id, hgvs_str_to_vrs_id_readonly, hgvs_to_allele,
    hgvs_to_allele_readonly,
};
pub use provider::{NoTranscriptProvider, ProviderError, TranscriptProvider};
pub use vcf::{
    VrsResult, compute_vrs_ids_from_vcf, compute_vrs_ids_from_vcf_readonly,
    compute_vrs_ids_streaming, compute_vrs_ids_streaming_readonly,
};
