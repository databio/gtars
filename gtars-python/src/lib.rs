use pyo3::prelude::*;
use pyo3::types::PyDict;

mod genomic_distributions;
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
    let gd_module = pyo3::wrap_pymodule!(genomic_distributions::genomic_distributions);

    m.add_wrapped(refget_module)?;
    m.add_wrapped(tokenize_module)?;
    m.add_wrapped(models_module)?;
    m.add_wrapped(utils_module)?;
    m.add_wrapped(gd_module)?;

    let sys = PyModule::import(py, "sys")?;
    let binding = sys.getattr("modules")?;
    let sys_modules: &Bound<'_, PyDict> = binding.cast()?;

    // set names of submodules
    sys_modules.set_item("gtars.refget", m.getattr("refget")?)?;
    sys_modules.set_item("gtars.tokenizers", m.getattr("tokenizers")?)?;
    sys_modules.set_item("gtars.models", m.getattr("models")?)?;
    sys_modules.set_item("gtars.utils", m.getattr("utils")?)?;
    sys_modules.set_item(
        "gtars.genomic_distributions",
        m.getattr("genomic_distributions")?,
    )?;

    // add constants
    m.add("__version__", VERSION)?;

    Ok(())
}
