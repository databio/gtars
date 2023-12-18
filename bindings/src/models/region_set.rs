use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

use genimtools::common::consts::{PAD_CHR, PAD_START, PAD_END};

use crate::models::{PyRegion, PyTokenizedRegion};

#[pyclass(name = "TokenizedRegionSet")]
#[derive(Clone, Debug)]
pub struct PyTokenizedRegionSet {
    pub regions: Vec<PyRegion>,
    pub ids: Vec<u32>,
    curr: usize,
}

#[pymethods]
impl PyTokenizedRegionSet {
    #[new]
    pub fn new(regions: Vec<PyRegion>, ids: Vec<u32>) -> Self {
        PyTokenizedRegionSet {
            regions,
            ids,
            curr: 0,
        }
    }

    #[getter]
    pub fn regions(&self) -> PyResult<Vec<PyRegion>> {
        Ok(self.regions.to_owned())
    }

    #[getter]
    pub fn ids(&self) -> PyResult<Vec<u32>> {
        Ok(self.ids.clone())
    }

    // this is wrong: the padding token might not be in the universe
    pub fn pad(&mut self, len: usize) {
        let pad_region = PyRegion {
            chr: PAD_CHR.to_string(),
            start: PAD_START as u32,
            end: PAD_END as u32,
        };
        let pad_id = self.ids[0];
        let pad_region_set = PyTokenizedRegionSet {
            regions: vec![pad_region; len],
            ids: vec![pad_id; len],
            curr: 0,
        };
        self.regions.extend(pad_region_set.regions);
        self.ids.extend(pad_region_set.ids);
    }

    pub fn __repr__(&self) -> String {
        format!("TokenizedRegionSet({} regions)", self.regions.len())
    }

    pub fn __len__(&self) -> usize {
        self.regions.len()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(&mut self) -> Option<PyTokenizedRegion> {
        if self.curr < self.regions.len() {
            let region = self.regions[self.curr].clone();
            let id = self.ids[self.curr];

            self.curr += 1;
            Some(PyTokenizedRegion::new(region, id))
        } else {
            None
        }
    }

    pub fn __getitem__(&self, indx: isize) -> PyResult<PyTokenizedRegion> {
        let indx = if indx < 0 {
            self.regions.len() as isize + indx
        } else {
            indx
        };
        if indx < 0 || indx >= self.regions.len() as isize {
            Err(PyIndexError::new_err("Index out of bounds"))
        } else {
            let region = self.regions[indx as usize].clone();
            let id = self.ids[indx as usize];
            Ok(PyTokenizedRegion::new(region, id))
        }
    }
}
