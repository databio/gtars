use pyo3::prelude::*;
use pyo3::types::PyDict;

mod ailist;
mod models;
mod tokenizers;
mod utils;
mod digests;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[pymodule]
fn gtars(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    let tokenize_module = pyo3::wrap_pymodule!(tokenizers::tokenizers);
    let ailist_module = pyo3::wrap_pymodule!(ailist::ailist);
    let utils_module = pyo3::wrap_pymodule!(utils::utils);
    let models_module = pyo3::wrap_pymodule!(models::models);
    let digests_module = pyo3::wrap_pymodule!(digests::digests);

    m.add_wrapped(tokenize_module)?;
    m.add_wrapped(ailist_module)?;
    m.add_wrapped(utils_module)?;
    m.add_wrapped(models_module)?;
    m.add_wrapped(digests_module)?;

    let sys = PyModule::import_bound(py, "sys")?;
    let binding = sys.getattr("modules")?;
    let sys_modules: &Bound<'_, PyDict> = binding.downcast()?;

    // set names of submodules
    sys_modules.set_item("gtars.tokenizers", m.getattr("tokenizers")?)?;
    sys_modules.set_item("gtars.ailist", m.getattr("ailist")?)?;
    sys_modules.set_item("gtars.utils", m.getattr("utils")?)?;
    sys_modules.set_item("gtars.models", m.getattr("models")?)?;
    sys_modules.set_item("gtars.digests", m.getattr("digests")?)?;

    // add constants
    m.add("__version__", VERSION)?;

    Ok(())
}
