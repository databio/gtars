use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

use anyhow::Result;
use genimtools::common::models::TokenizedRegionSet;

use crate::models::{PyRegion, PyTokenizedRegion, PyUniverse};

#[pyclass(name = "TokenizedRegionSet")]
#[derive(Clone, Debug)]
pub struct PyTokenizedRegionSet {
    pub ids: Vec<u32>,
    pub universe: PyUniverse,
    curr: usize,
}

impl From<TokenizedRegionSet<'_>> for PyTokenizedRegionSet {
    fn from(value: TokenizedRegionSet) -> Self {
        PyTokenizedRegionSet {
            ids: value.ids,
            universe: value.universe.to_owned().into(),
            curr: 0,
        }
    }
}

#[pymethods]
impl PyTokenizedRegionSet {
    #[getter]
    pub fn ids(&self) -> Result<Vec<u32>> {
        Ok(self.ids.clone())
    }

    #[getter]
    pub fn to_bit_vector(&self) -> Result<Vec<u8>> {
        let mut bit_vector = Vec::with_capacity(self.universe.len());

        for id in &self.ids {
            bit_vector[*id as usize] = 1;
        }

        Ok(bit_vector)
    }

    #[getter]
    pub fn to_regions(&self) -> Result<Vec<PyRegion>> {
        Ok(self
            .ids
            .iter()
            .map(|id| self.universe.id_to_region[&id].clone())
            .collect())
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
        format!("TokenizedRegionSet({} regions)", self.ids.len())
    }

    pub fn __len__(&self) -> usize {
        self.ids.len()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(&mut self) -> Option<PyTokenizedRegion> {
        if self.curr < self.ids.len() {
            let id = self.ids[self.curr];
            self.curr += 1;

            Some(PyTokenizedRegion {
                universe: self.universe.to_owned(),
                id,
            })
        } else {
            None
        }
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
