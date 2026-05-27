//! HGVS abstract syntax tree.

/// Reference sequence type of an HGVS expression.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceType {
    /// Genomic (`g.`)
    G,
    /// Coding transcript (`c.`)
    C,
    /// Non-coding transcript (`n.`)
    N,
    /// Mitochondrial (`m.`) — handled as genomic.
    M,
    /// RNA (`r.`) — handled similarly to `n.`; `U` is treated as `T`.
    R,
    /// Protein (`p.`) — parsed but not bridged in v1.
    P,
}

/// A parsed HGVS variant (borrowing the input string for zero-copy fields).
#[derive(Debug, Clone, PartialEq)]
pub struct HgvsVariant<'a> {
    /// Reference accession (e.g. `NM_004333.6`) or gene symbol (e.g. `BRAF`).
    pub accession: &'a str,
    /// Optional gene symbol in parentheses, e.g. `NM_004333.6(BRAF):c.1799T>A`.
    pub gene: Option<&'a str>,
    /// Reference type prefix (`g.`, `c.`, ...).
    pub reference_type: ReferenceType,
    /// Position + edit payload.
    pub posedit: PosEdit<'a>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PosEdit<'a> {
    pub pos: LocationRange,
    pub edit: Edit<'a>,
    /// True if wrapped in parentheses, i.e. uncertain.
    pub uncertain: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum LocationRange {
    Single(Position),
    Range { start: Position, end: Position },
    /// Whole sequence (no position specified, e.g., `g.=` or `g.?`)
    WholeSequence,
    /// Range with uncertain start: `(low_high)_end` e.g., `c.(4_6)_246del`
    UncertainStart {
        start_low: Option<Position>,  // None = `?`
        start_high: Option<Position>, // None = `?`
        end: Position,
    },
    /// Range with uncertain end: `start_(low_high)` e.g., `c.4_(245_246)del`
    UncertainEnd {
        start: Position,
        end_low: Option<Position>,  // None = `?`
        end_high: Option<Position>, // None = `?`
    },
    /// Range with both uncertain: `(low_high)_(low_high)`
    UncertainBoth {
        start_low: Option<Position>,
        start_high: Option<Position>,
        end_low: Option<Position>,
        end_high: Option<Position>,
    },
}

/// A single position with intronic offset and CDS-anchored datum.
///
/// HGVS examples:
/// - `c.100` → base=100, offset=0, datum=CdsStart
/// - `c.100+5` → base=100, offset=5, datum=CdsStart
/// - `c.100-3` → base=100, offset=-3, datum=CdsStart
/// - `c.-14` → base=-14, offset=0, datum=CdsStart
/// - `c.*37` → base=37, offset=0, datum=CdsEnd
/// - `c.*37+1` → base=37, offset=1, datum=CdsEnd
/// - `n.500` / `g.12345` → base=N, offset=0, datum=SeqStart
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Position {
    pub base: i64,
    pub offset: i64,
    pub datum: Datum,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Datum {
    /// Relative to sequence start (`g.`, `n.`, `m.`, `r.`).
    SeqStart,
    /// Relative to CDS start (`c.`).
    CdsStart,
    /// Relative to CDS end (`c.*N`).
    CdsEnd,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Edit<'a> {
    /// `A>T`
    Sub {
        reference: &'a str,
        alternate: &'a str,
    },
    /// `del` / `delAGT`
    Del { reference: Option<&'a str> },
    /// `dup` / `dupA`
    Dup { reference: Option<&'a str> },
    /// `insATG`
    Ins { alternate: &'a str },
    /// `delinsCT` / `delAinsCT`
    DelIns {
        reference: Option<&'a str>,
        alternate: &'a str,
    },
    /// `inv` — parsed but bridge rejects.
    Inv { reference: Option<&'a str> },
    /// `=` no change.
    Identity,
    /// `?` unknown edit.
    Unknown,
    /// `copyN` — Invitae extension for copy number (parsed but bridge rejects).
    Copy { count: u32 },
    /// `SEQ[N]` — repeat notation (e.g., `CA[4]` means 4 copies of CA).
    Repeat { sequence: &'a str, count: u32 },
}

// ─── Owned variants (for FFI / Python) ──────────────────────────────────────

#[derive(Debug, Clone, PartialEq)]
pub struct HgvsVariantOwned {
    pub accession: String,
    pub gene: Option<String>,
    pub reference_type: ReferenceType,
    pub posedit: PosEditOwned,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PosEditOwned {
    pub pos: LocationRange,
    pub edit: EditOwned,
    pub uncertain: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum EditOwned {
    Sub {
        reference: String,
        alternate: String,
    },
    Del {
        reference: Option<String>,
    },
    Dup {
        reference: Option<String>,
    },
    Ins {
        alternate: String,
    },
    DelIns {
        reference: Option<String>,
        alternate: String,
    },
    Inv {
        reference: Option<String>,
    },
    Identity,
    Unknown,
    Copy { count: u32 },
    Repeat { sequence: String, count: u32 },
}

impl<'a> From<HgvsVariant<'a>> for HgvsVariantOwned {
    fn from(v: HgvsVariant<'a>) -> Self {
        Self {
            accession: v.accession.to_string(),
            gene: v.gene.map(|s| s.to_string()),
            reference_type: v.reference_type,
            posedit: PosEditOwned {
                pos: v.posedit.pos,
                edit: match v.posedit.edit {
                    Edit::Sub { reference, alternate } => EditOwned::Sub {
                        reference: reference.to_string(),
                        alternate: alternate.to_string(),
                    },
                    Edit::Del { reference } => EditOwned::Del {
                        reference: reference.map(|s| s.to_string()),
                    },
                    Edit::Dup { reference } => EditOwned::Dup {
                        reference: reference.map(|s| s.to_string()),
                    },
                    Edit::Ins { alternate } => EditOwned::Ins {
                        alternate: alternate.to_string(),
                    },
                    Edit::DelIns {
                        reference,
                        alternate,
                    } => EditOwned::DelIns {
                        reference: reference.map(|s| s.to_string()),
                        alternate: alternate.to_string(),
                    },
                    Edit::Inv { reference } => EditOwned::Inv {
                        reference: reference.map(|s| s.to_string()),
                    },
                    Edit::Identity => EditOwned::Identity,
                    Edit::Unknown => EditOwned::Unknown,
                    Edit::Copy { count } => EditOwned::Copy { count },
                    Edit::Repeat { sequence, count } => EditOwned::Repeat {
                        sequence: sequence.to_string(),
                        count,
                    },
                },
                uncertain: v.posedit.uncertain,
            },
        }
    }
}
