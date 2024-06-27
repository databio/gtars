mod fragments_tokenizer;
mod meta_tokenizer;
mod tree_tokenizer;
// mod builder;

use pyo3::prelude::*;

pub use self::fragments_tokenizer::PyFragmentTokenizer;
pub use self::meta_tokenizer::PyMetaTokenizer;
pub use self::tree_tokenizer::PyTreeTokenizer;
pub use crate::models::{
    PyRegion, PyRegionSet, PyTokenizedRegion, PyTokenizedRegionSet, PyUniverse,
};

#[pymodule]
pub fn tokenizers(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTreeTokenizer>()?;
    m.add_class::<PyMetaTokenizer>()?;
    m.add_class::<PyFragmentTokenizer>()?;
    m.add_class::<PyRegion>()?;
    m.add_class::<PyTokenizedRegionSet>()?;
    m.add_class::<PyTokenizedRegion>()?;
    m.add_class::<PyUniverse>()?;
    m.add_class::<PyRegionSet>()?;
    Ok(())
}
