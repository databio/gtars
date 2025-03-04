use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use std::{path::Path, vec};

use numpy::ndarray::Array;
use numpy::{IntoPyArray, PyArray1};

use anyhow::Result;

use crate::models::{PyRegion, PyTokenizedRegion, PyUniverse};
use gtars::common::models::{Region, RegionSet};

#[pyclass(name = "RegionSet", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyRegionSet {
    pub regionset: RegionSet,
    curr: usize,
}

impl From<Vec<PyRegion>> for PyRegionSet {
    fn from(value: Vec<PyRegion>) -> Self {
        let mut rust_regions: Vec<Region> = Vec::new();
        for region in value {
            rust_regions.push(Region {
                chr: region.chr,
                start: region.start,
                end: region.end,
                rest: region.rest,
            })
        }
        PyRegionSet {
            regionset: RegionSet::from(rust_regions),
            curr: 0,
        }
    }
}

#[pymethods]
impl PyRegionSet {
    #[new]
    /// Create a new RegionSet object
    ///
    /// Args:
    ///     path: path to the bed, or bed.gz file
    ///
    /// Returns:
    ///     RegionSet object
    fn py_new(path: &str) -> PyResult<Self> {
        Ok(Self {
            regionset: RegionSet::try_from(path)?,
            curr: 0,
        })
    }

    #[getter]
    fn get_identifier(&self) -> PyResult<String> {
        Ok(self.regionset.identifier())
    }

    #[getter]
    fn get_path(&self) -> PyResult<String> {
        Ok(self
            .regionset
            .path
            .clone()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string())
    }

    #[getter]
    fn get_header(&self) -> PyResult<Option<String>> {
        Ok(self.regionset.header.clone())
    }

    fn __repr__(&self) -> String {
        self.regionset.to_string()
    }

    fn __str__(&self) -> String {
        self.regionset.to_string()
    }

    fn __len__(&self) -> usize {
        self.regionset.len()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(&mut self) -> Option<PyRegion> {
        if self.curr < self.regionset.regions.len() {
            let region = self.regionset.regions[self.curr].clone();
            self.curr += 1;

            Some(PyRegion {
                chr: region.chr,
                start: region.start,
                end: region.end,
                rest: region.rest,
            })
        } else {
            None
        }
    }

    pub fn __getitem__(&self, indx: isize) -> PyResult<PyRegion> {
        let len = self.regionset.regions.len() as isize;

        let indx = if indx < 0 {
            len + indx // Convert negative index to positive
        } else {
            indx
        };

        if indx < 0 || indx >= len {
            Err(PyIndexError::new_err("Index out of bounds"))
        } else {
            let r = &self.regionset.regions[indx as usize];
            Ok(PyRegion {
                chr: r.chr.clone(),
                start: r.start,
                end: r.end,
                rest: r.rest.clone(),
            })
        }
    }

    fn to_bigbed(&self, out_path: &str, chrom_size: &str) -> PyResult<()> {
        let chrom_size_path = Path::new(chrom_size);
        let out_path_path = Path::new(out_path);

        self.regionset.to_bigbed(out_path_path, chrom_size_path)?;
        Ok(())
    }

    fn to_bed(&self, path: &str) -> PyResult<()> {
        let path = Path::new(path);
        self.regionset.to_bed(path)?;
        Ok(())
    }

    fn to_bed_gz(&self, path: &str) -> PyResult<()> {
        let path = Path::new(path);
        self.regionset.to_bed_gz(path)?;
        Ok(())
    }

    fn sort(&mut self) -> PyResult<()> {
        self.regionset.sort();
        Ok(())
    }
}

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
