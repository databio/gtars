//! Python bindings for `gtars-vrs::hgvs`: AST + parser, exposed under
//! `gtars.vrs.hgvs`.

use std::fmt::Display;

use pyo3::create_exception;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use gtars_vrs::hgvs::ast::{
    Datum as RsDatum, EditOwned, HgvsVariantOwned, LocationRange, PosEditOwned, Position,
    ReferenceType as RsReferenceType,
};
use gtars_vrs::hgvs::parser::parse as rs_parse;

create_exception!(gtars.vrs.hgvs, HgvsError, pyo3::exceptions::PyException);

#[inline]
fn map_hgvs_err<E: Display>(e: E) -> PyErr {
    HgvsError::new_err(format!("{}", e))
}

// ---------------------------------------------------------------------------
// Reference type
// ---------------------------------------------------------------------------

#[pyclass(name = "ReferenceType", module = "gtars.vrs.hgvs", eq, eq_int)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceTypePy {
    Coding,
    NonCoding,
    Genomic,
    Mitochondrial,
    Protein,
    Rna,
}

#[pymethods]
impl ReferenceTypePy {
    fn __repr__(&self) -> String {
        format!("ReferenceType.{:?}", self)
    }
}

impl ReferenceTypePy {
    fn from_rs(r: RsReferenceType) -> Self {
        match r {
            RsReferenceType::G => ReferenceTypePy::Genomic,
            RsReferenceType::C => ReferenceTypePy::Coding,
            RsReferenceType::N => ReferenceTypePy::NonCoding,
            RsReferenceType::M => ReferenceTypePy::Mitochondrial,
            RsReferenceType::R => ReferenceTypePy::Rna,
            RsReferenceType::P => ReferenceTypePy::Protein,
        }
    }
}

// ---------------------------------------------------------------------------
// Datum
// ---------------------------------------------------------------------------

#[pyclass(name = "Datum", module = "gtars.vrs.hgvs", eq, eq_int)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DatumPy {
    SeqStart,
    Cds,
    CdsStop,
}

#[pymethods]
impl DatumPy {
    fn __repr__(&self) -> String {
        format!("Datum.{:?}", self)
    }
}

impl DatumPy {
    fn from_rs(d: RsDatum) -> Self {
        match d {
            RsDatum::SeqStart => DatumPy::SeqStart,
            RsDatum::CdsStart => DatumPy::Cds,
            RsDatum::CdsEnd => DatumPy::CdsStop,
        }
    }

    fn to_str(self) -> &'static str {
        match self {
            DatumPy::SeqStart => "seq_start",
            DatumPy::Cds => "cds",
            DatumPy::CdsStop => "cds_stop",
        }
    }
}

// ---------------------------------------------------------------------------
// Position
// ---------------------------------------------------------------------------

#[pyclass(name = "Position", module = "gtars.vrs.hgvs")]
#[derive(Debug, Clone, Copy)]
pub struct PositionPy {
    #[pyo3(get)]
    pub base: i64,
    #[pyo3(get)]
    pub offset: i64,
    #[pyo3(get)]
    pub datum: DatumPy,
}

#[pymethods]
impl PositionPy {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("base", self.base)?;
        d.set_item("offset", self.offset)?;
        d.set_item("datum", self.datum.to_str())?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "Position(base={}, offset={}, datum={:?})",
            self.base, self.offset, self.datum
        )
    }
}

impl PositionPy {
    fn from_rs(p: Position) -> Self {
        Self {
            base: p.base,
            offset: p.offset,
            datum: DatumPy::from_rs(p.datum),
        }
    }
}

// ---------------------------------------------------------------------------
// Edit
// ---------------------------------------------------------------------------

/// Edit payload. The `kind` discriminator is one of:
/// `"substitution"`, `"deletion"`, `"insertion"`, `"delins"`, `"duplication"`,
/// `"inversion"`, `"identity"`, `"unknown"`.
#[pyclass(name = "Edit", module = "gtars.vrs.hgvs")]
#[derive(Debug, Clone)]
pub struct EditPy {
    #[pyo3(get)]
    pub kind: String,
    /// Reference allele (substitution / del / dup / delins / inv); may be None.
    #[pyo3(get, name = "ref")]
    pub reference: Option<String>,
    /// Alternate allele (substitution / ins / delins).
    #[pyo3(get)]
    pub alt: Option<String>,
}

#[pymethods]
impl EditPy {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("kind", &self.kind)?;
        d.set_item("ref", &self.reference)?;
        d.set_item("alt", &self.alt)?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "Edit(kind={:?}, ref={:?}, alt={:?})",
            self.kind, self.reference, self.alt
        )
    }
}

impl EditPy {
    fn from_rs(e: EditOwned) -> Self {
        match e {
            EditOwned::Sub {
                reference,
                alternate,
            } => Self {
                kind: "substitution".to_string(),
                reference: Some(reference),
                alt: Some(alternate),
            },
            EditOwned::Del { reference } => Self {
                kind: "deletion".to_string(),
                reference,
                alt: None,
            },
            EditOwned::Dup { reference } => Self {
                kind: "duplication".to_string(),
                reference,
                alt: None,
            },
            EditOwned::Ins { alternate } => Self {
                kind: "insertion".to_string(),
                reference: None,
                alt: Some(alternate),
            },
            EditOwned::DelIns {
                reference,
                alternate,
            } => Self {
                kind: "delins".to_string(),
                reference,
                alt: Some(alternate),
            },
            EditOwned::Inv { reference } => Self {
                kind: "inversion".to_string(),
                reference,
                alt: None,
            },
            EditOwned::Identity => Self {
                kind: "identity".to_string(),
                reference: None,
                alt: None,
            },
            EditOwned::Unknown => Self {
                kind: "unknown".to_string(),
                reference: None,
                alt: None,
            },
            EditOwned::Copy { count } => Self {
                kind: "copy".to_string(),
                reference: None,
                alt: Some(format!("[{}]", count)),
            },
            EditOwned::Repeat { sequence, count } => Self {
                kind: "repeat".to_string(),
                reference: None,
                alt: Some(format!("{}[{}]", sequence, count)),
            },
        }
    }
}

// ---------------------------------------------------------------------------
// PositionBound
// ---------------------------------------------------------------------------

/// One side (start or end) of a location range.
///
/// `kind` is one of:
/// - `"certain"`: a concrete point, available in `position` (`low`/`high` None).
/// - `"uncertain"`: an interval bounded by `low` and `high`, where a `None`
///   bound represents the HGVS `?` token (unbounded on that side).
#[pyclass(name = "PositionBound", module = "gtars.vrs.hgvs")]
#[derive(Debug, Clone)]
pub struct PositionBoundPy {
    /// "certain" | "uncertain"
    #[pyo3(get)]
    pub kind: String,
    /// For kind="certain": the position. None for uncertain.
    #[pyo3(get)]
    pub position: Option<PositionPy>,
    /// For kind="uncertain": low bound; None means `?` (unbounded-low).
    #[pyo3(get)]
    pub low: Option<PositionPy>,
    /// For kind="uncertain": high bound; None means `?` (unbounded-high).
    #[pyo3(get)]
    pub high: Option<PositionPy>,
}

#[pymethods]
impl PositionBoundPy {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("kind", &self.kind)?;
        d.set_item("position", self.position.map(|p| p.to_dict(py)).transpose()?)?;
        d.set_item("low", self.low.map(|p| p.to_dict(py)).transpose()?)?;
        d.set_item("high", self.high.map(|p| p.to_dict(py)).transpose()?)?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "PositionBound(kind={:?}, position={:?}, low={:?}, high={:?})",
            self.kind, self.position, self.low, self.high
        )
    }
}

impl PositionBoundPy {
    fn certain(p: Position) -> Self {
        Self {
            kind: "certain".to_string(),
            position: Some(PositionPy::from_rs(p)),
            low: None,
            high: None,
        }
    }

    fn uncertain(low: Option<Position>, high: Option<Position>) -> Self {
        Self {
            kind: "uncertain".to_string(),
            position: None,
            low: low.map(PositionPy::from_rs),
            high: high.map(PositionPy::from_rs),
        }
    }
}

// ---------------------------------------------------------------------------
// PosEdit
// ---------------------------------------------------------------------------

#[pyclass(name = "PosEdit", module = "gtars.vrs.hgvs")]
#[derive(Debug, Clone)]
pub struct PosEditPy {
    /// One of: "single", "range", "whole_sequence".
    #[pyo3(get)]
    pub location_kind: String,
    /// Start endpoint. None only for whole_sequence.
    #[pyo3(get)]
    pub start: Option<PositionBoundPy>,
    /// End endpoint. None for single and whole_sequence.
    #[pyo3(get)]
    pub end: Option<PositionBoundPy>,
    #[pyo3(get)]
    pub edit: EditPy,
    #[pyo3(get)]
    pub uncertain: bool,
}

#[pymethods]
impl PosEditPy {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("location_kind", &self.location_kind)?;
        d.set_item(
            "start",
            self.start.as_ref().map(|b| b.to_dict(py)).transpose()?,
        )?;
        d.set_item(
            "end",
            self.end.as_ref().map(|b| b.to_dict(py)).transpose()?,
        )?;
        d.set_item("edit", self.edit.to_dict(py)?)?;
        d.set_item("uncertain", self.uncertain)?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "PosEdit(location_kind={:?}, start={:?}, end={:?}, edit={:?}, uncertain={})",
            self.location_kind, self.start, self.end, self.edit, self.uncertain
        )
    }
}

impl PosEditPy {
    fn from_rs(pe: PosEditOwned) -> Self {
        let (location_kind, start, end) = match pe.pos {
            LocationRange::Single(p) => {
                ("single", Some(PositionBoundPy::certain(p)), None)
            }
            LocationRange::Range { start, end } => (
                "range",
                Some(PositionBoundPy::certain(start)),
                Some(PositionBoundPy::certain(end)),
            ),
            LocationRange::WholeSequence => ("whole_sequence", None, None),
            LocationRange::UncertainStart {
                start_low,
                start_high,
                end,
            } => (
                "range",
                Some(PositionBoundPy::uncertain(start_low, start_high)),
                Some(PositionBoundPy::certain(end)),
            ),
            LocationRange::UncertainEnd {
                start,
                end_low,
                end_high,
            } => (
                "range",
                Some(PositionBoundPy::certain(start)),
                Some(PositionBoundPy::uncertain(end_low, end_high)),
            ),
            LocationRange::UncertainBoth {
                start_low,
                start_high,
                end_low,
                end_high,
            } => (
                "range",
                Some(PositionBoundPy::uncertain(start_low, start_high)),
                Some(PositionBoundPy::uncertain(end_low, end_high)),
            ),
        };
        Self {
            location_kind: location_kind.to_string(),
            start,
            end,
            edit: EditPy::from_rs(pe.edit),
            uncertain: pe.uncertain,
        }
    }
}

// ---------------------------------------------------------------------------
// HgvsVariant
// ---------------------------------------------------------------------------

#[pyclass(name = "HgvsVariant", module = "gtars.vrs.hgvs")]
#[derive(Debug, Clone)]
pub struct HgvsVariantPy {
    #[pyo3(get)]
    pub accession: String,
    #[pyo3(get)]
    pub gene: Option<String>,
    #[pyo3(get)]
    pub reference_type: ReferenceTypePy,
    #[pyo3(get)]
    pub pos_edit: PosEditPy,
}

#[pymethods]
impl HgvsVariantPy {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("accession", &self.accession)?;
        d.set_item("gene", &self.gene)?;
        d.set_item(
            "reference_type",
            match self.reference_type {
                ReferenceTypePy::Coding => "c",
                ReferenceTypePy::NonCoding => "n",
                ReferenceTypePy::Genomic => "g",
                ReferenceTypePy::Mitochondrial => "m",
                ReferenceTypePy::Protein => "p",
                ReferenceTypePy::Rna => "r",
            },
        )?;
        d.set_item("pos_edit", self.pos_edit.to_dict(py)?)?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "HgvsVariant(accession={:?}, reference_type={:?}, pos_edit={:?})",
            self.accession, self.reference_type, self.pos_edit
        )
    }
}

impl From<HgvsVariantOwned> for HgvsVariantPy {
    fn from(v: HgvsVariantOwned) -> Self {
        Self {
            accession: v.accession,
            gene: v.gene,
            reference_type: ReferenceTypePy::from_rs(v.reference_type),
            pos_edit: PosEditPy::from_rs(v.posedit),
        }
    }
}

// ---------------------------------------------------------------------------
// parse_hgvs
// ---------------------------------------------------------------------------

/// Parse an HGVS string into a structured `HgvsVariant`.
///
/// Raises `HgvsError` on invalid input.
#[pyfunction]
#[pyo3(signature = (s))]
pub fn parse_hgvs(py: Python<'_>, s: &str) -> PyResult<HgvsVariantPy> {
    let owned = py
        .allow_threads(|| {
            rs_parse(s).map(HgvsVariantOwned::from)
        })
        .map_err(map_hgvs_err)?;
    Ok(HgvsVariantPy::from(owned))
}

// ---------------------------------------------------------------------------
// Bridge: hgvs_to_vrs_id
// ---------------------------------------------------------------------------

/// Parse + bridge + normalize + digest in one call.
///
/// Returns the canonical `ga4gh:VA.<digest>` identifier matching what the
/// equivalent VCF row would produce.
///
/// Args:
///     hgvs_str (str): HGVS expression (e.g. `"NM_004333.6:c.1799T>A"`).
///     provider (gtars.reftx.ReftxProvider): Transcript provider.
///     refget (gtars.refget.RefgetStore): Mutable RefgetStore.
///     collection_digest (str): Sequence collection digest (level-0).
///
/// Returns:
///     str: VRS allele identifier (`ga4gh:VA.<digest>`).
///
/// Raises:
///     HgvsError: parse, mapping, or refget failure.
#[pyfunction]
#[pyo3(signature = (hgvs_str, provider, refget, collection_digest))]
pub fn hgvs_to_vrs_id(
    py: Python<'_>,
    hgvs_str: &str,
    provider: &crate::reftx::ReftxProviderPy,
    refget: &mut crate::refget::PyRefgetStore,
    collection_digest: &str,
) -> PyResult<String> {
    let s = hgvs_str.to_string();
    let coll = collection_digest.to_string();
    let prov = provider.inner.clone();
    // NOTE: warnings are currently dropped at the FFI boundary; surfacing them
    // to Python is out of scope for this change.
    py.allow_threads(|| {
        gtars_vrs::hgvs::bridge::hgvs_str_to_vrs_id(
            &s,
            &prov,
            &mut refget.inner,
            &coll,
        )
    })
    .map(|b| b.value)
    .map_err(map_hgvs_err)
}

// ---------------------------------------------------------------------------
// Module registration
// ---------------------------------------------------------------------------

pub fn register(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ReferenceTypePy>()?;
    m.add_class::<DatumPy>()?;
    m.add_class::<PositionPy>()?;
    m.add_class::<EditPy>()?;
    m.add_class::<PositionBoundPy>()?;
    m.add_class::<PosEditPy>()?;
    m.add_class::<HgvsVariantPy>()?;
    m.add_function(wrap_pyfunction!(parse_hgvs, m)?)?;
    m.add_function(wrap_pyfunction!(hgvs_to_vrs_id, m)?)?;
    m.add("HgvsError", py.get_type::<HgvsError>())?;
    Ok(())
}
