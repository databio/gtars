use gtars_genomicdist::models::GenomeAssembly;
use pyo3::prelude::*;
use std::path::Path;

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
    ///     GenomeAssebly object
    pub fn new(path: &Bound<'_, PyAny>) -> Self {
        // let path = Path::new(&path.to_string());
        let genome_assembly = GenomeAssembly::try_from(path.to_string()).unwrap();
        PyGenomeAssembly { genome_assembly }
    }
}
