mod py_tokenizers;
mod tokens;
mod universe;

use pyo3::prelude::*;

use crate::tokenizers::py_tokenizers::PyTokenizer;
use crate::tokenizers::universe::PyUniverse;
use crate::tokenizers::tokens::PyTokenizedRegionSet;

#[pymodule]
pub fn tokenizers(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTokenizer>()?;
    m.add_class::<PyTokenizedRegionSet>()?;
    m.add_class::<PyUniverse>()?;
    Ok(())
}
