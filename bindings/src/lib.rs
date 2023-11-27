use pyo3::prelude::*;
use pyo3::types::PyDict;

mod vocab;
mod tokenizers;
mod models;
mod consts;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[pymodule]
fn genimtools(py: Python, m: &PyModule) -> PyResult<()> {

    let vocab_module = pyo3::wrap_pymodule!(vocab::vocab);
    let tokenize_module = pyo3::wrap_pymodule!(tokenizers::tokenizers);

    m.add_wrapped(vocab_module)?;
    m.add_wrapped(tokenize_module)?;

    let sys = PyModule::import(py, "sys")?;
    let sys_modules: &PyDict = sys.getattr("modules")?.downcast()?;

    // set names of submodules
    sys_modules.set_item("genimtools.vocab", m.getattr("vocab")?)?;
    sys_modules.set_item("genimtools.tokenizers", m.getattr("tokenizers")?)?;

    // add constants
    m.add("PAD_CHR", consts::PAD_CHR)?;
    m.add("PAD_START", consts::PAD_START)?;
    m.add("PAD_END", consts::PAD_END)?;
    m.add("UNKNOWN_CHR", consts::UNKNOWN_CHR)?;
    m.add("UNKNOWN_START", consts::UNKNOWN_START)?;
    m.add("UNKNOWN_END", consts::UNKNOWN_END)?;
    m.add("__version__", VERSION)?;

    Ok(())
}

