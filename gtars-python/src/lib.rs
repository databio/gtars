#[cfg(target_os = "linux")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use pyo3::prelude::*;
use pyo3::types::PyDict;

#[cfg(feature = "bm25")]
mod bm25;
#[cfg(feature = "genomic_distributions")]
mod genomic_distributions;
#[cfg(feature = "models")]
mod models;
#[cfg(feature = "refget")]
mod refget;
#[cfg(feature = "tokenizers")]
mod tokenizers;
#[cfg(any(feature = "utils", feature = "tokenizers"))]
mod utils;
#[cfg(feature = "lola")]
mod lola;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[pymodule]
fn gtars(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    let sys = PyModule::import(py, "sys")?;
    let binding = sys.getattr("modules")?;
    let sys_modules: &Bound<'_, PyDict> = binding.cast()?;

    #[cfg(feature = "bm25")]
    {
        let bm25_module = pyo3::wrap_pymodule!(bm25::bm25);
        m.add_wrapped(bm25_module)?;
        sys_modules.set_item("gtars.bm25", m.getattr("bm25")?)?;
    }

    #[cfg(feature = "refget")]
    {
        let refget_module = pyo3::wrap_pymodule!(refget::refget);
        m.add_wrapped(refget_module)?;
        sys_modules.set_item("gtars.refget", m.getattr("refget")?)?;
    }

    #[cfg(feature = "tokenizers")]
    {
        let tokenize_module = pyo3::wrap_pymodule!(tokenizers::tokenizers);
        m.add_wrapped(tokenize_module)?;
        sys_modules.set_item("gtars.tokenizers", m.getattr("tokenizers")?)?;
    }

    #[cfg(feature = "models")]
    {
        let models_module = pyo3::wrap_pymodule!(models::models);
        m.add_wrapped(models_module)?;
        sys_modules.set_item("gtars.models", m.getattr("models")?)?;
    }

    #[cfg(feature = "utils")]
    {
        let utils_module = pyo3::wrap_pymodule!(utils::utils);
        m.add_wrapped(utils_module)?;
        sys_modules.set_item("gtars.utils", m.getattr("utils")?)?;
    }

    #[cfg(feature = "genomic_distributions")]
    {
        let gd_module = pyo3::wrap_pymodule!(genomic_distributions::genomic_distributions);
        m.add_wrapped(gd_module)?;
        sys_modules.set_item(
            "gtars.genomic_distributions",
            m.getattr("genomic_distributions")?,
        )?;
    }

    #[cfg(feature = "lola")]
    {
        let lola_module = pyo3::wrap_pymodule!(lola::lola);
        m.add_wrapped(lola_module)?;
        sys_modules.set_item("gtars.lola", m.getattr("lola")?)?;
    }

    // add constants
    m.add("__version__", VERSION)?;

    Ok(())
}
