mod py_tokenizers;
mod tokens;
mod universe;

use pyo3::prelude::*;

#[pymodule]
pub fn tokenizers(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // will place stuff here -- I nuked the module lol
    Ok(())
}
