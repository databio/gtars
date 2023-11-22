use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::vocab::prune_universe;
use crate::tokenizers::PyTreeTokenizer;
use crate::models::{PyRegion, PyTokenizedRegion, PyTokenizedRegionSet, PyUniverse};

mod vocab;
mod tokenizers;
mod models;
mod consts;

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
fn consts_module(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add("PAD_CHR", consts::PAD_CHR)?;
    m.add("PAD_START", consts::PAD_START)?;
    m.add("PAD_END", consts::PAD_END)?;
    m.add("UNKNOWN_CHR", consts::UNKNOWN_CHR)?;
    m.add("UNKNOWN_START", consts::UNKNOWN_START)?;
    m.add("UNKNOWN_END", consts::UNKNOWN_END)?;
    Ok(())
}

#[pymodule]
fn genimtools(py: Python, m: &PyModule) -> PyResult<()> {    
    let vocab_module = pyo3::wrap_pymodule!(vocab_module);
    let tokenize_module = pyo3::wrap_pymodule!(tokenize_module);
    let consts_module = pyo3::wrap_pymodule!(consts_module);

    m.add_wrapped(vocab_module)?;
    m.add_wrapped(tokenize_module)?;
    m.add_wrapped(consts_module)?;

    let sys = PyModule::import(py, "sys")?;
    let sys_modules: &PyDict = sys.getattr("modules")?.downcast()?;

    sys_modules.set_item("genimtools.vocab", m.getattr("vocab_module")?)?;
    sys_modules.set_item("genimtools.tokenize", m.getattr("tokenize_module")?)?;
    sys_modules.set_item("genimtools.consts", m.getattr("consts_module")?)?;

    Ok(())
}

