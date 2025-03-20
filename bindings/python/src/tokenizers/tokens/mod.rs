use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

use crate::models::PyRegion;
use crate::tokenizers::universe::PyUniverse;

use anyhow::Result;

#[pyclass(name = "TokenizedRegionSet", module = "gtars.models")]
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

    pub fn to_bit_vector(&self) -> Result<Vec<u8>> {
        Python::with_gil(|py| {
            let mut bit_vector = vec![0; self.universe.borrow(py).len()];

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
                .map(|id| {
                    self.universe
                        .borrow(py)
                        .convert_id_to_region(*id)
                        .unwrap()
                        .clone()
                })
                .collect())
        })
    }

    pub fn to_ids(&self) -> Result<Vec<u32>> {
        Ok(self.ids.clone())
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

#[pyclass(name = "TokenizedRegion", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyTokenizedRegion {
    pub id: u32,
    pub universe: Py<PyUniverse>,
}

impl From<PyTokenizedRegion> for PyRegion {
    fn from(value: PyTokenizedRegion) -> Self {
        Python::with_gil(|py| {
            value
                .universe
                .borrow(py)
                .convert_id_to_region(value.id)
                .unwrap()
        })
    }
}

#[pymethods]
impl PyTokenizedRegion {
    #[getter]
    pub fn chr(&self) -> Result<String> {
        Python::with_gil(|py| {
            Ok(self
                .universe
                .borrow(py)
                .convert_id_to_region(self.id)
                .unwrap()
                .chr)
        })
    }

    #[getter]
    pub fn start(&self) -> Result<u32> {
        Python::with_gil(|py| {
            Ok(self
                .universe
                .borrow(py)
                .convert_id_to_region(self.id)
                .unwrap()
                .start)
        })
    }

    #[getter]
    pub fn end(&self) -> Result<u32> {
        Python::with_gil(|py| {
            Ok(self
                .universe
                .borrow(py)
                .convert_id_to_region(self.id)
                .unwrap()
                .end)
        })
    }

    pub fn to_region(&self) -> Result<PyRegion> {
        Python::with_gil(|py| {
            Ok(self
                .universe
                .borrow(py)
                .convert_id_to_region(self.id)
                .unwrap())
        })
    }
    #[getter]
    pub fn id(&self) -> Result<u32> {
        Ok(self.id)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "TokenizedRegion({}, {}, {})",
            self.chr().unwrap(),
            self.start().unwrap(),
            self.end().unwrap()
        )
    }
}
