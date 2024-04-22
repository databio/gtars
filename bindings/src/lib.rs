use pyo3::prelude::*;
use pyo3::types::PyDict;

mod ailist;
mod models;
mod tokenizers;
mod utils;
mod vocab;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[pymodule]
fn genimtools(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    let vocab_module = pyo3::wrap_pymodule!(vocab::vocab);
    let tokenize_module = pyo3::wrap_pymodule!(tokenizers::tokenizers);
    let ailist_module = pyo3::wrap_pymodule!(ailist::ailist);
    let utils_module = pyo3::wrap_pymodule!(utils::utils);
    let models_modeule = pyo3::wrap_pymodule!(models::models);

    m.add_wrapped(vocab_module)?;
    m.add_wrapped(tokenize_module)?;
    m.add_wrapped(ailist_module)?;
    m.add_wrapped(utils_module)?;
    m.add_wrapped(models_modeule)?;

    let sys = PyModule::import_bound(py, "sys")?;
    let binding = sys.getattr("modules")?;
    let sys_modules: &Bound<'_, PyDict> = binding.downcast()?;

    // set names of submodules
    sys_modules.set_item("genimtools.vocab", m.getattr("vocab")?)?;
    sys_modules.set_item("genimtools.tokenizers", m.getattr("tokenizers")?)?;
    sys_modules.set_item("genimtools.ailist", m.getattr("ailist")?)?;
    sys_modules.set_item("genimtools.utils", m.getattr("utils")?)?;
    sys_modules.set_item("genimtools.models", m.getattr("models")?)?;

    // add constants
    m.add("__version__", VERSION)?;

    Ok(())
}
