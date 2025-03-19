use pyo3::prelude::*;

use anyhow::Result;

use crate::models::PyRegion;
use gtars::tokenizers::universe::Universe;

#[pyclass(name = "Universe", module = "gtars.tokenizers")]
#[derive(Clone, Debug)]
pub struct PyUniverse {
    pub universe: Universe,
}

impl From<Universe> for PyUniverse {
    fn from(value: Universe) -> Self {
        PyUniverse {
            universe: value,
        }
    }
}

#[pymethods]
impl PyUniverse {
    pub fn add_token_to_universe(&mut self, region: &PyRegion) {
        let region = region.to_region();
        self.universe.add_token_to_universe(&region);
    }

    pub fn convert_region_to_id(&self, region: &PyRegion) -> Option<u32> {
        let region = region.to_region();
        self.universe.convert_region_to_id(&region)
    }

    pub fn convert_id_to_region(&self, id: u32) -> Option<PyRegion> {
        self.universe.convert_id_to_region(id).map(|region| PyRegion::from(region))
    }

    pub fn len(&self) -> usize {
        self.universe.len()
    }

    pub fn is_empty(&self) -> bool {
        self.universe.is_empty()
    }

    pub fn is_ordered(&self) -> bool {
        self.universe.is_ordered()
    }

    #[setter]
    pub fn ordered(&mut self, ordered: bool) {
        self.universe.ordered = ordered;
    }

    pub fn __len__(&self) -> usize {
        self.len()
    }

    pub fn __repr__(&self) -> String {
        format!("Universe with {} regions", self.len())
    }
}
