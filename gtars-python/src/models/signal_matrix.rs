use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use gtars_genomicdist::SignalMatrix;

#[pyclass(name = "SignalMatrix", module = "gtars.models")]
pub struct PySignalMatrix {
    pub signal_matrix: SignalMatrix,
}

#[pymethods]
impl PySignalMatrix {
    #[staticmethod]
    pub fn load_bin(path: &str) -> PyResult<Self> {
        let sm = SignalMatrix::load_bin(path)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { signal_matrix: sm })
    }

    #[staticmethod]
    pub fn from_tsv(path: &str) -> PyResult<Self> {
        let sm = SignalMatrix::from_tsv(path)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { signal_matrix: sm })
    }

    #[getter]
    fn condition_names(&self) -> Vec<String> {
        self.signal_matrix.condition_names.clone()
    }

    #[getter]
    fn n_conditions(&self) -> usize {
        self.signal_matrix.n_conditions
    }

    #[getter]
    fn n_regions(&self) -> usize {
        self.signal_matrix.regions.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "SignalMatrix(n_regions={}, n_conditions={})",
            self.signal_matrix.regions.len(),
            self.signal_matrix.n_conditions
        )
    }

    fn __len__(&self) -> usize {
        self.signal_matrix.regions.len()
    }
}
