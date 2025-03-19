mod py_tokenizers;
mod tokens;
mod universe;

use pyo3::prelude::*;

use crate::tokenizers::py_tokenizers::PyTokenizer;

#[pymodule]
pub fn tokenizers(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTokenizer>()?;
    Ok(())
}
