use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::functions::{prune_universe, print_uniwig};

mod functions;

#[pymodule]
fn vocab(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(prune_universe))?;
    Ok(())
}

#[pymodule]
fn uniwig(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(print_uniwig))?;
    Ok(())
}

#[pymodule]
fn genimtools(py: Python, m: &PyModule) -> PyResult<()> {    
    let vocab_module = pyo3::wrap_pymodule!(vocab);
    let uniwig_module = pyo3::wrap_pymodule!(uniwig);

    m.add_wrapped(vocab_module)?;
    m.add_wrapped(uniwig_module)?;

    let sys = PyModule::import(py, "sys")?;
    let sys_modules: &PyDict = sys.getattr("modules")?.downcast()?;
    sys_modules.set_item("genimtools.vocab", m.getattr("vocab")?)?;
    sys_modules.set_item("genimtools.uniwig", m.getattr("uniwig")?)?;
    Ok(())
}

