use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use std::{path::Path, vec};

use anyhow::Result;
use gtars::common::utils::extract_regions_from_bed_file;

use crate::models::{PyRegion, PyTokenizedRegion, PyUniverse};

use super::region::PyTokenizedRegionPointer;

#[pyclass(name = "RegionSet", module = "gtars.models")]
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

#[pyclass(name = "TokenizedRegionSet", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyTokenizedRegionSet {
    pub pointers: Vec<PyTokenizedRegionPointer>,
    pub universe: Py<PyUniverse>,
    pub curr: usize,
}

#[pymethods]
impl PyTokenizedRegionSet {
    #[getter]
    pub fn ids(&self) -> Result<Vec<u32>> {
        Ok(self.pointers.iter().map(|p| p.id).collect::<Vec<u32>>())
    }

    #[getter]
    pub fn source_chrom_ids(&self) -> Result<Vec<u16>> {
        Ok(self.pointers.iter().map(|p| p.source_chrom_id).collect::<Vec<u16>>())
    }

    #[getter]
    pub fn source_starts(&self) -> Result<Vec<u32>> {
        Ok(self.pointers.iter().map(|p| p.source_start).collect::<Vec<u32>>())
    }

    #[getter]
    pub fn source_ends(&self) -> Result<Vec<u32>> {
        Ok(self.pointers.iter().map(|p| p.source_end).collect::<Vec<u32>>())
    }

    #[getter]
    pub fn pointers(&self) -> Result<Vec<PyTokenizedRegionPointer>> {
        Ok(self.pointers.clone())
    }

    pub fn to_bit_vector(&self) -> Result<Vec<u8>> {
        Python::with_gil(|py| {
            let mut bit_vector = vec![0; self.universe.borrow(py).id_to_region.len()];
            let ids = self.ids()?;

            for id in ids {
                bit_vector[id as usize] = 1;
            }

            Ok(bit_vector)
        })
    }

    pub fn to_regions(&self) -> Result<Vec<PyRegion>> {
        Python::with_gil(|py| {
            Ok(self
                .ids()?
                .iter()
                .map(|id| self.universe.borrow(py).id_to_region[id].clone())
                .collect())
        })
    }

    pub fn to_ids(&self) -> Result<Vec<u32>> {
        self.ids()
    }

    pub fn __repr__(&self) -> String {
        format!(
            "TokenizedRegionSet({:?})",
            self.pointers.iter().map(|p| p.id).collect::<Vec<u32>>()
        )
    }

    pub fn __len__(&self) -> usize {
        self.pointers.len()
    }

    pub fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(&mut self) -> Option<PyTokenizedRegion> {
        Python::with_gil(|py| {
            if self.curr < self.pointers.len() {
                let pointer = self.pointers[self.curr];
                self.curr += 1;

                Some(PyTokenizedRegion {
                    universe: self.universe.clone_ref(py),
                    pointer,
                })
            } else {
                None
            }
        })
    }

    pub fn __getitem__(&self, indx: isize) -> Result<PyTokenizedRegion> {
        let indx = if indx < 0 {
            self.pointers.len() as isize + indx
        } else {
            indx
        };
        if indx < 0 || indx >= self.pointers.len() as isize {
            anyhow::bail!(PyIndexError::new_err("Index out of bounds"));
        } else {
            Ok(PyTokenizedRegion {
                universe: self.universe.to_owned(),
                pointer: self.pointers[indx as usize],
            })
        }
    }
}
