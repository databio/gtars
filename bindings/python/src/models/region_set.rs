use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use std::{path::Path, vec};

use numpy::ndarray::Array;
use numpy::{IntoPyArray, PyArray1};

use anyhow::Result;
use gtars::common::utils::extract_regions_from_bed_file;

use crate::models::{PyRegion, PyTokenizedRegion, PyUniverse};

#[pyclass(name = "RegionSet", module="gtars.models")]
#[derive(Clone, Debug)]
pub struct PyRegionSet {
    pub regions: Vec<PyRegion>,
    curr: usize,
}

impl From<Vec<PyRegion>> for PyRegionSet {
    fn from(value: Vec<PyRegion>) -> Self {
        PyRegionSet {
            regions: value,
            curr: 0,
        }
    }
}

#[pymethods]
impl PyRegionSet {
    #[new]
    pub fn new(path: String) -> Result<Self> {
        let path = Path::new(&path);
        let regions = extract_regions_from_bed_file(path)?;

        Ok(PyRegionSet {
            regions: regions.into_iter().map(|region| region.into()).collect(),
            curr: 0,
        })
    }

    pub fn __repr__(&self) -> String {
        format!("RegionSet({} regions)", self.regions.len())
    }

    pub fn __len__(&self) -> usize {
        self.regions.len()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(&mut self) -> Option<PyRegion> {
        if self.curr < self.regions.len() {
            let region = self.regions[self.curr].clone();
            self.curr += 1;

            Some(PyRegion {
                chr: region.chr,
                start: region.start,
                end: region.end,
            })
        } else {
            None
        }
    }

    pub fn __getitem__(&self, indx: isize) -> Result<PyRegion> {
        let indx = if indx < 0 {
            self.regions.len() as isize + indx
        } else {
            indx
        };
        if indx < 0 || indx >= self.regions.len() as isize {
            anyhow::bail!(PyIndexError::new_err("Index out of bounds"));
        } else {
            let r = self.regions[indx as usize].clone();
            Ok(PyRegion {
                chr: r.chr,
                start: r.start,
                end: r.end,
            })
        }
    }
}

#[pyclass(name = "TokenizedRegionSet", module="gtars.models")]
#[derive(Clone, Debug)]
pub struct PyTokenizedRegionSet {
    pub ids: Vec<u32>,
    pub universe: Py<PyUniverse>,
    pub curr: usize,
}

#[pymethods]
impl PyTokenizedRegionSet {
    #[getter]
    pub fn ids(&self) -> Result<Vec<u32>> {
        Ok(self.ids.clone())
    }

    fn to_numpy<'py>(&self, py: Python<'py>) -> Result<Bound<'py, PyArray1<u32>>> {
        let array = Array::from_vec(self.ids.clone()).into_pyarray_bound(py);

        Ok(array)
    }

    pub fn to_bit_vector(&self) -> Result<Vec<u8>> {
        Python::with_gil(|py| {
            let mut bit_vector = vec![0; self.universe.borrow(py).id_to_region.len()];

            for id in &self.ids {
                bit_vector[*id as usize] = 1;
            }

            Ok(bit_vector)
        })
    }

    pub fn to_regions(&self) -> Result<Vec<PyRegion>> {
        Python::with_gil(|py| {
            Ok(self
                .ids
                .iter()
                .map(|id| self.universe.borrow(py).id_to_region[id].clone())
                .collect())
        })
    }

    pub fn to_ids(&self) -> Result<Vec<u32>> {
        Ok(self.ids.clone())
    }

    // gensim needs strings as input, so to speed up
    // iterating over datasets, lets provide a rust
    // interface to directly convert to strings
    #[getter]
    pub fn ids_as_strs(&self) -> Result<Vec<String>> {
        Ok(self
            .ids
            .to_owned()
            .iter()
            .map(|id| id.to_string())
            .collect())
    }

    pub fn __repr__(&self) -> String {
        format!("TokenizedRegionSet({:?})", self.ids)
    }

    pub fn __len__(&self) -> usize {
        self.ids.len()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(&mut self) -> Option<PyTokenizedRegion> {
        Python::with_gil(|py| {
            if self.curr < self.ids.len() {
                let id = self.ids[self.curr];
                self.curr += 1;

                Some(PyTokenizedRegion {
                    universe: self.universe.clone_ref(py),
                    id,
                })
            } else {
                None
            }
        })
    }

    pub fn __getitem__(&self, indx: isize) -> Result<PyTokenizedRegion> {
        let indx = if indx < 0 {
            self.ids.len() as isize + indx
        } else {
            indx
        };
        if indx < 0 || indx >= self.ids.len() as isize {
            anyhow::bail!(PyIndexError::new_err("Index out of bounds"));
        } else {
            Ok(PyTokenizedRegion {
                universe: self.universe.to_owned(),
                id: self.ids[indx as usize],
            })
        }
    }
}
