use gtars_genomicdist::models::GenomeAssembly;
use pyo3::prelude::*;

#[pyclass(name = "GenomeAssembly")]
pub struct PyGenomeAssembly {
    pub genome_assembly: GenomeAssembly,
}

#[pymethods]
impl PyGenomeAssembly {
    #[new]
    /// Create a new GenomeAssembly object
    ///
    /// Args:
    ///     path: path to the fasta file
    ///
    /// Returns:
    ///     GenomeAssembly object
    pub fn new(path: &Bound<'_, PyAny>) -> PyResult<Self> {
        let genome_assembly = GenomeAssembly::try_from(path.to_string())
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(PyGenomeAssembly { genome_assembly })
    }
}
