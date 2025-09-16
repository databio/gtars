use std::hash::Hash;

use pyo3::class::basic::CompareOp;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

use anyhow::Result;
use gtars_core::models::region::Region;

#[pyclass(name = "Region", module = "gtars.models")]
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub struct PyRegion {
    #[pyo3(get)]
    pub chr: String,

    #[pyo3(get)]
    pub start: u32,

    #[pyo3(get)]
    pub end: u32,

    #[pyo3(get)]
    pub rest: Option<String>,
}

impl From<Region> for PyRegion {
    fn from(region: Region) -> Self {
        PyRegion {
            chr: region.chr,
            start: region.start,
            end: region.end,
            rest: region.rest,
        }
    }
}

impl PyRegion {
    pub fn to_region(&self) -> Region {
        Region {
            chr: self.chr.clone(),
            start: self.start,
            end: self.end,
            rest: self.rest.clone(),
        }
    }

    pub fn into_region(self) -> Region {
        Region {
            chr: self.chr,
            start: self.start,
            end: self.end,
            rest: self.rest,
        }
    }
}

#[pymethods]
impl PyRegion {
    #[new]
    pub fn new(chr: String, start: u32, end: u32, rest: Option<String>) -> Self {
        PyRegion {
            chr,
            start,
            end,
            rest,
        }
    }

    fn __repr__(&self) -> String {
        format!("Region -> {} {} {}", self.chr, self.start, self.end)
    }

    fn __str__(&self) -> String {
        format!(
            "{}\t{}\t{}{}",
            self.chr,
            self.start,
            self.end,
            self.rest
                .as_deref()
                .map_or(String::new(), |s| format!("\t{}", s)),
        )
    }

    fn __len__(&self) -> usize {
        (self.end - self.start) as usize
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
