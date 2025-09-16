use pyo3::prelude::*;

use gtars_tokenizers::universe::Universe;

#[pyclass(name = "Universe", module = "gtars.tokenizers")]
#[derive(Clone, Debug)]
pub struct PyUniverse {
    pub universe: Universe,
}

impl From<Universe> for PyUniverse {
    fn from(value: Universe) -> Self {
        PyUniverse { universe: value }
    }
}

#[pymethods]
impl PyUniverse {
    pub fn len(&self) -> usize {
        self.universe.len()
    }

    pub fn is_empty(&self) -> bool {
        self.universe.is_empty()
    }

    pub fn __len__(&self) -> usize {
        self.len()
    }

    pub fn __repr__(&self) -> String {
        format!("Universe with {} regions", self.len())
    }
}
