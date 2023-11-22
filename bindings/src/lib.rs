use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::vocab::prune_universe;
use crate::tokenizers::PyTreeTokenizer;
use crate::models::{PyRegion, PyTokenizedRegion, PyTokenizedRegionSet, PyUniverse};

mod vocab;
mod tokenizers;
mod models;

#[pymodule]
fn vocab_module(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(prune_universe))?;
    Ok(())
}

#[pymodule]
fn tokenize_module(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyTreeTokenizer>()?;
    m.add_class::<PyRegion>()?;
    m.add_class::<PyTokenizedRegionSet>()?;
    m.add_class::<PyTokenizedRegion>()?;
    m.add_class::<PyUniverse>()?;
    Ok(())
}

#[pymodule]
fn genimtools(py: Python, m: &PyModule) -> PyResult<()> {    
    let vocab_module = pyo3::wrap_pymodule!(vocab_module);
    let tokenize_module = pyo3::wrap_pymodule!(tokenize_module);

    m.add_wrapped(vocab_module)?;
    m.add_wrapped(tokenize_module)?;

    let sys = PyModule::import(py, "sys")?;
    let sys_modules: &PyDict = sys.getattr("modules")?.downcast()?;

    sys_modules.set_item("genimtools.vocab", m.getattr("vocab_module")?)?;
    sys_modules.set_item("genimtools.tokenize", m.getattr("tokenize_module")?)?;

    Ok(())
}

