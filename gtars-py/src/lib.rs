use pyo3::prelude::*;
use pyo3::types::PyDict;

mod models;
mod refget;
mod tokenizers;
mod utils;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[pymodule]
fn gtars(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    let refget_module = pyo3::wrap_pymodule!(refget::refget);
    let tokenize_module = pyo3::wrap_pymodule!(tokenizers::tokenizers);
    let models_module = pyo3::wrap_pymodule!(models::models);
    let utils_module = pyo3::wrap_pymodule!(utils::utils);

    m.add_wrapped(refget_module)?;
    m.add_wrapped(tokenize_module)?;
    m.add_wrapped(models_module)?;
    m.add_wrapped(utils_module)?;

    let sys = PyModule::import_bound(py, "sys")?;
    let binding = sys.getattr("modules")?;
    let sys_modules: &Bound<'_, PyDict> = binding.downcast()?;

    // set names of submodules
    sys_modules.set_item("gtars.refget", m.getattr("refget")?)?;
    sys_modules.set_item("gtars.tokenizers", m.getattr("tokenizers")?)?;
    sys_modules.set_item("gtars.models", m.getattr("models")?)?;
    sys_modules.set_item("gtars.utils", m.getattr("utils")?)?;

    // add constants
    m.add("__version__", VERSION)?;

    Ok(())
}
