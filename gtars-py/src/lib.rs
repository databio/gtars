use pyo3::prelude::*;
use pyo3::types::PyDict;

mod refget;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[pymodule]
fn gtars(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {

    let refget_module = pyo3::wrap_pymodule!(refget::refget);

    m.add_wrapped(refget_module)?;

    let sys = PyModule::import_bound(py, "sys")?;
    let binding = sys.getattr("modules")?;
    let sys_modules: &Bound<'_, PyDict> = binding.downcast()?;

    // set names of submodules
    sys_modules.set_item("gtars.refget", m.getattr("refget")?)?;

    // add constants
    m.add("__version__", VERSION)?;

    Ok(())
}
