use gtars::tokenizers::config::TokenizerConfig;
use pyo3::prelude::*;
use pyo3::types::PyType;

use anyhow::Result;

use crate::tokenizers::universe::PyUniverse;
use gtars::tokenizers::Tokenizer;

#[pyclass(name = "Tokenizer", module = "gtars.tokenizers")]
pub struct PyTokenizer {
    tokenizer: Tokenizer,
    universe: Py<PyUniverse>, // this is a Py-wrapped version self.tokenizer.universe for performance reasons
}

#[pymethods]
impl PyTokenizer {
    #[classmethod]
    pub fn from_config(_cls: &Bound<'_, PyType>, cfg: &str) -> Result<Self> {
        Python::with_gil(|py| {
            let tokenizer = Tokenizer::from_config(cfg)?;
            let universe = PyUniverse::from(tokenizer.get_universe().clone());
            Ok(PyTokenizer {
                tokenizer,
                universe: Py::new(py, universe)?,
            })
        })
    }

    #[classmethod]
    pub fn from_bed(_cls: &Bound<'_, PyType>, path: &str) -> Result<Self> {
        Python::with_gil(|py| {
            let tokenizer = Tokenizer::from_bed(path)?;
            let universe = PyUniverse::from(tokenizer.get_universe().clone());
            Ok(PyTokenizer {
                tokenizer,
                universe: Py::new(py, universe)?,
            })
        })
    }
}
