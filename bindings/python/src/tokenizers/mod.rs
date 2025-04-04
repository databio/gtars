mod encoding;
mod py_tokenizers;
mod universe;
mod utils;

use pyo3::prelude::*;

use crate::tokenizers::py_tokenizers::PyTokenizer;
use crate::tokenizers::utils::py_create_instances;
// use crate::tokenizers::universe::PyUniverse;
// use crate::tokenizers::encoding::{PyBatchEncoding, PyEncoding};

#[pymodule]
pub fn tokenizers(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTokenizer>()?;
    m.add_wrapped(wrap_pyfunction!(py_create_instances))?;
    Ok(())
}
