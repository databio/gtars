//! Python bindings for `gtars-vrs`.
//!
//! Top-level free functions live in `funcs.rs`; the HGVS submodule lives in
//! `hgvs.rs` and is registered as `gtars.vrs.hgvs` from `lib.rs`.

use pyo3::prelude::*;

pub mod funcs;
pub mod hgvs;

#[pymodule]
pub fn vrs(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(funcs::vrs_digest, m)?)?;
    m.add_function(wrap_pyfunction!(funcs::vrs_id, m)?)?;
    m.add_function(wrap_pyfunction!(funcs::normalize_allele, m)?)?;
    m.add_function(wrap_pyfunction!(funcs::location_digest, m)?)?;
    Ok(())
}
