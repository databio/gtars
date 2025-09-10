use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use std::collections::HashMap;

use crate::models::PyRegion;
use gtars_core::models::{Region, RegionSet};

#[pyclass(name = "RegionSet", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyRegionSet {
    pub regionset: RegionSet,
    curr: usize,
    identifier: Option<String>,
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
            identifier: None,
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
    fn py_new(path: &Bound<'_, PyAny>) -> PyResult<Self> {
        let path = path.to_string();
        Ok(Self {
            regionset: RegionSet::try_from(path)?,
            curr: 0,
            identifier: None,
        })
    }

    /// Alternate constructor from a list of PyRegion
    /// Args:
    ///     regions: a list/vec of PyRegion objects
    ///
    /// Returns:
    ///     RegionSet object
    #[staticmethod]
    fn from_regions(regions: Vec<PyRegion>) -> PyResult<Self> {
        let rust_regions: Vec<Region> = regions.into_iter().map(PyRegion::into_region).collect();

        Ok(Self {
            regionset: RegionSet::from(rust_regions),
            curr: 0,
            identifier: None,
        })
    }

    #[getter]
    fn get_identifier(&mut self) -> PyResult<String> {
        if self.identifier.is_none() {
            self.identifier = Some(self.regionset.identifier());
        };
        Ok(self.identifier.clone().unwrap())
    }

    #[getter]
    fn get_file_digest(&self) -> PyResult<String> {
        Ok(self.regionset.file_digest())
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

    fn is_empty(&self) -> PyResult<bool> {
        Ok(self.regionset.is_empty())
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
            self.curr = 0;
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

    fn to_bigbed(
        &self,
        out_path: &Bound<'_, PyAny>,
        chrom_size: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        self.regionset
            .to_bigbed(out_path.to_string(), chrom_size.to_string())?;
        Ok(())
    }

    fn to_bed(&self, path: &Bound<'_, PyAny>) -> PyResult<()> {
        self.regionset.to_bed(path.to_string())?;
        Ok(())
    }

    fn to_bed_gz(&self, path: &Bound<'_, PyAny>) -> PyResult<()> {
        self.regionset.to_bed_gz(path.to_string())?;
        Ok(())
    }

    fn sort(&mut self) -> PyResult<()> {
        self.regionset.sort();
        Ok(())
    }

    fn mean_region_width(&self) -> f64 {
        self.regionset.mean_region_width()
    }

    fn get_max_end_per_chr(&self) -> HashMap<String, u32> {
        self.regionset.get_max_end_per_chr()
    }

    fn get_nucleotide_length(&self) -> u32 {
        self.regionset.nucleotides_length()
    }
}
