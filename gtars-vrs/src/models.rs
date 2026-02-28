//! VRS data models for digest computation.
//!
//! Minimal structs representing GA4GH VRS objects — just enough
//! for canonical JSON serialization and digest computation.

/// A reference to a specific sequence identified by its refget accession.
#[derive(Debug, Clone)]
pub struct SequenceReference {
    /// GA4GH refget accession, e.g. "SQ.F-LrL..."
    pub refget_accession: String,
}

/// A location on a sequence defined by start/end coordinates (interbase, 0-based).
#[derive(Debug, Clone)]
pub struct SequenceLocation {
    pub sequence_reference: SequenceReference,
    pub start: u64,
    pub end: u64,
}

/// The state (alternate allele) of a VRS Allele.
#[derive(Debug, Clone)]
pub enum AlleleState {
    /// A literal sequence expression (SNV, indel, MNV).
    LiteralSequenceExpression { sequence: String },
    /// A reference-length expression (for CNVs/repeats).
    ReferenceLengthExpression {
        length: u64,
        repeat_subunit_length: u64,
        sequence: Option<String>,
    },
}

/// A VRS Allele: a specific sequence state at a specific genomic location.
#[derive(Debug, Clone)]
pub struct Allele {
    pub location: SequenceLocation,
    pub state: AlleleState,
}
