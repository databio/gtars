mod tree_tokenizer;

use pyo3::prelude::*;

pub use self::tree_tokenizer::PyTreeTokenizer;
pub use crate::models::{PyRegion, PyTokenizedRegion, PyTokenizedRegionSet, PyUniverse};

#[pymodule]
pub fn tokenizers(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyTreeTokenizer>()?;
    m.add_class::<PyRegion>()?;
    m.add_class::<PyTokenizedRegionSet>()?;
    m.add_class::<PyTokenizedRegion>()?;
    m.add_class::<PyUniverse>()?;
    Ok(())
}