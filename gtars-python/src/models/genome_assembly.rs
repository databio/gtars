use gtars_genomicdist::models::{BinaryGenomeAssembly, GenomeAssembly};
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

/// Memory-mapped genome assembly backed by a .fab binary file.
///
/// Near-instant construction via mmap with zero-copy sequence access.
/// Create .fab files with ``gtars prep --fasta <file>``.
#[pyclass(name = "BinaryGenomeAssembly")]
pub struct PyBinaryGenomeAssembly {
    pub assembly: BinaryGenomeAssembly,
}

#[pymethods]
impl PyBinaryGenomeAssembly {
    #[new]
    /// Open a .fab binary FASTA file.
    ///
    /// Args:
    ///     path: path to the .fab file
    ///
    /// Returns:
    ///     BinaryGenomeAssembly object
    pub fn new(path: &Bound<'_, PyAny>) -> PyResult<Self> {
        let assembly = BinaryGenomeAssembly::try_from(path.to_string())
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(PyBinaryGenomeAssembly { assembly })
    }
}
