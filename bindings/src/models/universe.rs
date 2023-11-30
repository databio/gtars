use std::collections::HashMap;

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

use crate::models::PyRegion;

#[pyclass(name = "Universe")]
#[derive(Clone, Debug)]
pub struct PyUniverse {
    pub regions: Vec<PyRegion>,
    pub region_to_id: HashMap<PyRegion, u32>,
    pub length: u32,
}

#[pymethods]
impl PyUniverse {
    #[getter]
    pub fn regions(&self) -> PyResult<Vec<PyRegion>> {
        Ok(self.regions.to_owned())
    }
    pub fn region_to_id(&self, region: &PyAny) -> PyResult<u32> {
        
        // get the chr, start, end (use can provide anything that has these attributes)
        let chr: String = region.getattr("chr")?.extract()?;
        let start: u32 = region.getattr("start")?.extract()?;
        let end: u32 = region.getattr("end")?.extract()?;

        // use to create PyRegion
        let region = PyRegion { chr, start, end };

        let id = self.region_to_id.get(&region);
        match id {
            Some(id) => Ok(id.to_owned()),
            None => Err(PyTypeError::new_err("Region not found in universe")),
        }
    }
    pub fn __len__(&self) -> usize {
        self.length as usize
    }
}