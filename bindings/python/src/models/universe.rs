use std::collections::HashMap;

use pyo3::prelude::*;

use anyhow::Result;

use crate::models::PyRegion;
use gtars::common::models::Universe;

#[pyclass(name = "Universe", module="gtars.models")]
#[derive(Clone, Debug)]
pub struct PyUniverse {
    pub regions: Vec<PyRegion>,
    pub region_to_id: HashMap<PyRegion, u32>,
    pub id_to_region: HashMap<u32, PyRegion>,
}

impl From<Universe> for PyUniverse {
    fn from(value: Universe) -> Self {
        let regions = value.regions.into_iter().map(|r| r.into()).collect();

        let region_to_id = value
            .region_to_id
            .into_iter()
            .map(|(k, v)| (k.into(), v))
            .collect();

        let id_to_region = value
            .id_to_region
            .into_iter()
            .map(|(k, v)| (k, v.into()))
            .collect();

        PyUniverse {
            regions,
            region_to_id,
            id_to_region,
        }
    }
}

#[pymethods]
impl PyUniverse {
    pub fn insert_token(&mut self, region: &PyRegion) {
        let new_id = self.region_to_id.len() + 1;
        self.region_to_id.insert(region.to_owned(), new_id as u32);
        self.id_to_region.insert(new_id as u32, region.to_owned());
    }

    pub fn convert_region_to_id(&self, region: &PyRegion) -> Option<u32> {
        let id = self.region_to_id.get(region);
        id.map(|id| id.to_owned())
    }

    pub fn convert_chr_start_end_to_id(&self, chr: &str, start: u32, end: u32) -> Option<u32> {
        let region = PyRegion {
            chr: chr.to_string(),
            start,
            end,
        };
        self.convert_region_to_id(&region)
    }

    pub fn convert_id_to_region(&self, id: u32) -> Option<PyRegion> {
        self.id_to_region.get(&id).cloned()
    }

    pub fn len(&self) -> usize {
        self.region_to_id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.region_to_id.len() == 0
    }

    #[getter]
    pub fn regions(&self) -> Result<Vec<PyRegion>> {
        Ok(self.regions.clone())
    }

    pub fn __len__(&self) -> usize {
        self.len()
    }

    pub fn __repr__(&self) -> String {
        format!("Universe with {} regions", self.len())
    }
}
