//! HGVS variant parsing and VRS bridge.
//!
//! Supports g., c., n., and m. reference types with substitution, deletion,
//! duplication, insertion, delins, and identity edits. Single positions and
//! ranges are both supported, including 5' UTR (`c.-N`), 3' UTR (`c.*N`),
//! and intronic offsets (`c.93+5`, `c.94-3`).
//!
//! ```ignore
//! use gtars_vrs::hgvs::{parse, ReferenceType};
//!
//! let v = parse("NM_004333.6:c.1799T>A").unwrap();
//! assert_eq!(v.accession, "NM_004333.6");
//! assert!(matches!(v.reference_type, ReferenceType::C));
//! ```

pub mod ast;
pub mod bridge;
pub mod error;
pub mod parser;

pub use ast::{
    Datum, Edit, HgvsVariant, HgvsVariantOwned, LocationRange, PosEdit, Position, ReferenceType,
};
pub use bridge::{BridgeError, hgvs_str_to_vrs_id_readonly, hgvs_to_allele_readonly};
#[cfg(feature = "filesystem")]
pub use bridge::{hgvs_str_to_vrs_id, hgvs_to_allele};
pub use error::HgvsError;
pub use parser::parse;
