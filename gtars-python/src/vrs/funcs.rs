//! Free functions exposed at the `gtars.vrs` top level.

use pyo3::prelude::*;
use pyo3::types::PyDict;

use gtars_vrs::digest::{allele_digest, allele_identifier, sequence_location_digest};
use gtars_vrs::models::{Allele, AlleleState, SequenceLocation, SequenceReference};
use gtars_vrs::normalize::normalize;

/// Build an Allele with a LiteralSequenceExpression from primitive parameters.
fn build_literal_allele(seq_digest: &str, start: u64, end: u64, alt: &str) -> Allele {
    Allele {
        location: SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: seq_digest.to_string(),
            },
            start,
            end,
        },
        state: AlleleState::LiteralSequenceExpression {
            sequence: alt.to_string(),
        },
    }
}

/// Compute the VRS digest for a single allele (without "ga4gh:VA." prefix).
#[pyfunction]
#[pyo3(signature = (seq_digest, start, end, alt))]
pub fn vrs_digest(seq_digest: &str, start: u64, end: u64, alt: &str) -> PyResult<String> {
    Ok(allele_digest(&build_literal_allele(
        seq_digest, start, end, alt,
    )))
}

/// Compute the full VRS identifier for a single allele ("ga4gh:VA.<digest>").
#[pyfunction]
#[pyo3(signature = (seq_digest, start, end, alt))]
pub fn vrs_id(seq_digest: &str, start: u64, end: u64, alt: &str) -> PyResult<String> {
    Ok(allele_identifier(&build_literal_allele(
        seq_digest, start, end, alt,
    )))
}

/// Normalize an allele against a reference sequence.
#[pyfunction]
#[pyo3(signature = (sequence, start, ref_allele, alt_allele))]
pub fn normalize_allele<'py>(
    py: Python<'py>,
    sequence: &str,
    start: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let result = normalize(
        sequence.as_bytes(),
        start,
        ref_allele.as_bytes(),
        alt_allele.as_bytes(),
    )
    .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("{}", e)))?;
    let dict = PyDict::new(py);
    dict.set_item("start", result.start)?;
    dict.set_item("end", result.end)?;
    dict.set_item(
        "allele",
        String::from_utf8(result.allele).map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("Invalid UTF-8 in allele: {}", e))
        })?,
    )?;
    Ok(dict)
}

/// Compute the VRS SequenceLocation digest.
#[pyfunction]
#[pyo3(signature = (seq_digest, start, end))]
pub fn location_digest(seq_digest: &str, start: u64, end: u64) -> PyResult<String> {
    let loc = SequenceLocation {
        sequence_reference: SequenceReference {
            refget_accession: seq_digest.to_string(),
        },
        start,
        end,
    };
    Ok(sequence_location_digest(&loc))
}
