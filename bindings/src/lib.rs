use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::vocab::prune_universe;

mod vocab;
mod tokenizers;
mod models;

#[pymodule]
fn vocab_module(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(prune_universe))?;
    Ok(())
}

#[pymodule]
fn genimtools(py: Python, m: &PyModule) -> PyResult<()> {    
    let vocab_module = pyo3::wrap_pymodule!(vocab_module);
    // create more here

    m.add_wrapped(vocab_module)?;
    // add more here

    let sys = PyModule::import(py, "sys")?;
    let sys_modules: &PyDict = sys.getattr("modules")?.downcast()?;
    sys_modules.set_item("genimtools.vocab", m.getattr("vocab")?)?;
    sys_modules.set_item("genimtools.uniwig", m.getattr("uniwig")?)?;
    Ok(())
}

