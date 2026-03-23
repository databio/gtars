use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

use crate::models::PyRegionSet;
use gtars_core::models::{RegionSet, RegionSetList};
use gtars_genomicdist::pairwise_jaccard;

/// A collection of RegionSets — the gtars equivalent of GRangesList.
#[pyclass(name = "RegionSetList", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyRegionSetList {
    pub inner: RegionSetList,
}

impl PyRegionSetList {
    pub fn from_inner(inner: RegionSetList) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyRegionSetList {
    /// Create a RegionSetList from a list of RegionSet objects.
    #[new]
    fn py_new(sets: Vec<PyRef<PyRegionSet>>) -> Self {
        let region_sets: Vec<RegionSet> = sets.iter().map(|s| s.regionset.clone()).collect();
        Self {
            inner: RegionSetList::new(region_sets),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__(&self, index: isize) -> PyResult<PyRegionSet> {
        let len = self.inner.len() as isize;
        let idx = if index < 0 { len + index } else { index } as usize;
        match self.inner.get(idx) {
            Some(rs) => Ok(PyRegionSet::from_regionset(rs.clone())),
            None => Err(PyIndexError::new_err(format!(
                "index {} out of range for RegionSetList of length {}",
                index,
                self.inner.len()
            ))),
        }
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRegionSetListIterator {
        PyRegionSetListIterator {
            inner: slf.inner.clone(),
            index: 0,
        }
    }

    /// Flatten all region sets into a single RegionSet (no merging).
    fn concat(&self) -> PyRegionSet {
        PyRegionSet::from_regionset(self.inner.concat())
    }

    /// Get the names of the region sets, or None.
    #[getter]
    fn names(&self) -> Option<Vec<String>> {
        self.inner.names.clone()
    }

    fn __repr__(&self) -> String {
        format!("RegionSetList({} region sets)", self.inner.len())
    }

    /// Compute pairwise Jaccard similarity for all pairs of region sets.
    ///
    /// Returns an N x N list of lists (symmetric matrix) with 1.0 on the diagonal.
    fn pairwise_jaccard(&self) -> Vec<Vec<f64>> {
        let sets: Vec<RegionSet> = (0..self.inner.len())
            .filter_map(|i| self.inner.get(i).cloned())
            .collect();
        let flat = pairwise_jaccard(&sets);
        let n = sets.len();
        (0..n).map(|i| flat[i * n..(i + 1) * n].to_vec()).collect()
    }
}

#[pyclass]
struct PyRegionSetListIterator {
    inner: RegionSetList,
    index: usize,
}

#[pymethods]
impl PyRegionSetListIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<PyRegionSet> {
        if self.index < self.inner.len() {
            let rs = self.inner.get(self.index).unwrap().clone();
            self.index += 1;
            Some(PyRegionSet::from_regionset(rs))
        } else {
            None
        }
    }
}
