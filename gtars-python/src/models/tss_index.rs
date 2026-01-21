use crate::models::PyRegionSet;
use gtars_genomicdist::models::TssIndex;
use pyo3::prelude::*;

#[pyclass(name = "TssIndex")]
pub struct PyTssIndex {
    pub tss_index: TssIndex,
}

#[pymethods]
impl PyTssIndex {
    #[new]
    /// Create a new TSS index
    ///
    /// Args:
    ///     path: path to the bed, or bed.gz file that contains tss regions
    ///
    /// Returns:
    ///     TSS index object
    fn py_new(path: &Bound<'_, PyAny>) -> PyResult<Self> {
        let path = path.to_string();
        Ok(Self {
            tss_index: TssIndex::try_from(path)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?,
        })
    }

    pub fn calc_tss_distances(&self, rs: &PyRegionSet) -> PyResult<Vec<u32>> {
        let distances: Vec<u32> = self
            .tss_index
            .calc_tss_distances(&rs.regionset)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(distances)
    }

    pub fn __repr__(&self) -> String {
        self.tss_index.region_set.to_string()
    }

    pub fn __str__(&self) -> String {
        self.tss_index.region_set.to_string()
    }

    fn __len__(&self) -> usize {
        self.tss_index.region_set.len()
    }
}
