//! Python bindings for VRS allele digest computation.

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
///
/// Args:
///     seq_digest (str): Refget sequence accession (e.g. "SQ.xxx")
///     start (int): 0-based interbase start coordinate
///     end (int): 0-based interbase end coordinate
///     alt (str): Alternate allele sequence
///
/// Returns:
///     str: The VRS allele digest (32-char base64url)
#[pyfunction]
#[pyo3(signature = (seq_digest, start, end, alt))]
pub fn vrs_digest(seq_digest: &str, start: u64, end: u64, alt: &str) -> PyResult<String> {
    Ok(allele_digest(&build_literal_allele(seq_digest, start, end, alt)))
}

/// Compute the full VRS identifier for a single allele ("ga4gh:VA.<digest>").
///
/// Args:
///     seq_digest (str): Refget sequence accession (e.g. "SQ.xxx")
///     start (int): 0-based interbase start coordinate
///     end (int): 0-based interbase end coordinate
///     alt (str): Alternate allele sequence
///
/// Returns:
///     str: The VRS allele identifier (e.g. "ga4gh:VA.xxx")
#[pyfunction]
#[pyo3(signature = (seq_digest, start, end, alt))]
pub fn vrs_id(seq_digest: &str, start: u64, end: u64, alt: &str) -> PyResult<String> {
    Ok(allele_identifier(&build_literal_allele(seq_digest, start, end, alt)))
}

/// Normalize an allele against a reference sequence.
///
/// Args:
///     sequence (str): Full reference chromosome sequence
///     start (int): 0-based interbase start position
///     ref_allele (str): Reference allele
///     alt_allele (str): Alternate allele
///
/// Returns:
///     dict: {"start": int, "end": int, "allele": str}
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
///
/// Args:
///     seq_digest (str): Refget sequence accession (e.g. "SQ.xxx")
///     start (int): 0-based interbase start coordinate
///     end (int): 0-based interbase end coordinate
///
/// Returns:
///     str: The SequenceLocation digest (32-char base64url)
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

#[pymodule]
pub fn vrs(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(vrs_digest, m)?)?;
    m.add_function(wrap_pyfunction!(vrs_id, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_allele, m)?)?;
    m.add_function(wrap_pyfunction!(location_digest, m)?)?;
    Ok(())
}
