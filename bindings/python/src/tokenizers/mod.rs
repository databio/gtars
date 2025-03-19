mod tokenizers;
mod universe;

use pyo3::prelude::*;

#[pymodule]
#[pyo3(name = "tokenizers")]
pub fn tokenizers_m(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // will place stuff here -- I nuked the module lol
    Ok(())
}
