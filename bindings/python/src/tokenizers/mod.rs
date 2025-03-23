mod py_tokenizers;
mod universe;

use pyo3::prelude::*;

use crate::tokenizers::py_tokenizers::PyTokenizer;
// use crate::tokenizers::universe::PyUniverse;
// use crate::tokenizers::encoding::{PyBatchEncoding, PyEncoding};

#[pymodule]
pub fn tokenizers(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTokenizer>()?;
    Ok(())
}
