use std::hash::Hash;

use pyo3::class::basic::CompareOp;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

use anyhow::Result;
use gtars::common::models::region::Region;

use crate::models::PyUniverse;

#[pyclass(name = "Region", module="gtars.models")]
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub struct PyRegion {
    pub chr: String,
    pub start: u32,
    pub end: u32,
}

impl From<Region> for PyRegion {
    fn from(region: Region) -> Self {
        PyRegion {
            chr: region.chr,
            start: region.start,
            end: region.end,
        }
    }
}

impl PyRegion {
    pub fn to_region(&self) -> Region {
        Region {
            chr: self.chr.clone(),
            start: self.start,
            end: self.end,
        }
    }
}

#[pymethods]
impl PyRegion {
    #[new]
    pub fn new(chr: String, start: u32, end: u32) -> Self {
        PyRegion { chr, start, end }
    }

    #[getter]
    pub fn chr(&self) -> Result<&str> {
        Ok(&self.chr)
    }

    #[getter]
    pub fn start(&self) -> Result<u32> {
        Ok(self.start)
    }

    #[getter]
    pub fn end(&self) -> Result<u32> {
        Ok(self.end)
    }
    pub fn __repr__(&self) -> String {
        format!("Region({}, {}, {})", self.chr, self.start, self.end)
    }

    pub fn __richcmp__(&self, other: PyRef<PyRegion>, op: CompareOp) -> Result<bool> {
        match op {
            CompareOp::Eq => {
                Ok(self.chr == other.chr && self.start == other.start && self.end == other.end)
            }
            CompareOp::Ne => {
                Ok(self.chr != other.chr || self.start != other.start || self.end != other.end)
            }
            _ => anyhow::bail!(PyTypeError::new_err("Unsupported comparison operator")),
        }
    }
}

#[pyclass(name = "TokenizedRegion", module="gtars.models")]
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
