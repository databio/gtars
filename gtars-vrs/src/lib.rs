//! # GA4GH VRS Allele Digest Computation
//!
//! This crate computes GA4GH VRS (Variation Representation Specification) identifiers
//! for genomic variants. It provides:
//!
//! - VRS data models (Allele, SequenceLocation, etc.)
//! - VRS digest computation (canonical JSON serialization + SHA-512/24u)
//! - Allele normalization (port of bioutils `normalize()`)
//! - VCF batch processing (parse VCF → compute VRS IDs)

pub mod models;
pub mod digest;
pub mod normalize;
pub mod vcf;

pub use models::{Allele, AlleleState, SequenceLocation, SequenceReference};
pub use digest::{allele_digest, allele_identifier, sequence_location_digest};
pub use normalize::normalize;
pub use vcf::{VrsResult, compute_vrs_ids_from_vcf};
