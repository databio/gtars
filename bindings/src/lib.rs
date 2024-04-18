use pyo3::prelude::*;
use pyo3::types::PyDict;

mod ailist;
mod consts;
mod models;
mod tokenizers;
mod utils;
mod vocab;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[pymodule]
fn genimtools(py: Python, m: &PyModule) -> PyResult<()> {
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

    let sys = PyModule::import(py, "sys")?;
    let sys_modules: &PyDict = sys.getattr("modules")?.downcast()?;

    // set names of submodules
    sys_modules.set_item("genimtools.vocab", m.getattr("vocab")?)?;
    sys_modules.set_item("genimtools.tokenizers", m.getattr("tokenizers")?)?;
    sys_modules.set_item("genimtools.ailist", m.getattr("ailist")?)?;
    sys_modules.set_item("genimtools.utils", m.getattr("utils")?)?;
    sys_modules.set_item("genimtools.models", m.getattr("models")?)?;

    // add constants
    m.add("PAD_CHR", consts::special_tokens::PAD_CHR)?;
    m.add("PAD_START", consts::special_tokens::PAD_START)?;
    m.add("PAD_END", consts::special_tokens::PAD_END)?;
    m.add("UNKNOWN_CHR", consts::special_tokens::UNKNOWN_CHR)?;
    m.add("UNKNOWN_START", consts::special_tokens::UNKNOWN_START)?;
    m.add("UNKNOWN_END", consts::special_tokens::UNKNOWN_END)?;
    m.add("__version__", VERSION)?;

    Ok(())
}
