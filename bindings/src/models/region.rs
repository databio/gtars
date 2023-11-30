use std::hash::Hash;

use pyo3::class::basic::CompareOp;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

use genimtools::common::models::region::Region;

#[pyclass(name = "Region")]
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub struct PyRegion {
    pub chr: String,
    pub start: u32,
    pub end: u32,
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
    pub fn chr(&self) -> PyResult<&str> {
        Ok(&self.chr)
    }

    #[getter]
    pub fn start(&self) -> PyResult<u32> {
        Ok(self.start)
    }

    #[getter]
    pub fn end(&self) -> PyResult<u32> {
        Ok(self.end)
    }
    pub fn __repr__(&self) -> String {
        format!("Region({}, {}, {})", self.chr, self.start, self.end)
    }

    pub fn __richcmp__(&self, other: PyRef<PyRegion>, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => {
                Ok(self.chr == other.chr && self.start == other.start && self.end == other.end)
            }
            CompareOp::Ne => {
                Ok(self.chr != other.chr || self.start != other.start || self.end != other.end)
            }
            _ => Err(PyTypeError::new_err("Unsupported comparison operator")),
        }
    }
}

#[pyclass(name = "TokenizedRegion")]
#[derive(Clone, Debug)]
pub struct PyTokenizedRegion {
    pub region: PyRegion,
    pub id: u32,
}

#[pymethods]
impl PyTokenizedRegion {
    #[new]
    pub fn new(region: PyRegion, id: u32) -> Self {
        PyTokenizedRegion {
            region,
            id,
        }
    }

    #[getter]
    pub fn chr(&self) -> PyResult<&str> {
        Ok(&self.region.chr)
    }

    #[getter]
    pub fn start(&self) -> PyResult<u32> {
        Ok(self.region.start)
    }

    #[getter]
    pub fn end(&self) -> PyResult<u32> {
        Ok(self.region.end)
    }

    #[getter]
    pub fn region(&self) -> PyResult<PyRegion> {
        Ok(self.region.clone())
    }
    #[getter]
    pub fn id(&self) -> PyResult<u32> {
        Ok(self.id)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "TokenizedRegion({}, {}, {})",
            self.region.chr, self.region.start, self.region.end
        )
    }
}