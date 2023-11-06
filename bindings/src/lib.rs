use pyo3::prelude::*;

use crate::functions::prune_universe;

mod functions;

#[pymodule]
fn genimtools(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(prune_universe, m)?)?;
    Ok(())
}