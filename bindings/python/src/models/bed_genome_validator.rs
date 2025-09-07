use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::convert::TryInto;
use crate::models::PyRegion;
use gtars::common::models::{ReferenceValidator, CompatibilityConcise, ReferenceGenomeMetadata};
use std::path::{Path, PathBuf};
use crate::models::PyRegionSet;

#[pyclass(name = "ReferenceGenomeMetadata", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyReferenceGenomeMetadata {
    pub genome_reference: ReferenceGenomeMetadata,
}



#[pymethods]
impl PyReferenceGenomeMetadata {
    #[new]
    /// Create a new Reference genome metadata model
    ///
    /// Args:
    ///     path: path to the folder with chrom sizes
    ///
    /// Returns:
    ///     ReferenceValidator object
    fn py_new(
        genome: &Bound<'_, PyAny>,
        digest: &Bound<'_, PyAny>,
        description: &Bound<'_, PyAny>,
        collection: HashMap<String, u32>
    ) -> PyResult<Self> {
        let genome = genome.to_string();
        let digest = digest.to_string();
        let description = description.to_string();

        Ok(PyReferenceGenomeMetadata {
            genome_reference: ReferenceGenomeMetadata {
                genome,
                digest,
                description,
                collection,
            }
        })
    }

    /// Alternate constructor from a path, or string
    /// Args:
    ///     path: path to the chrom sizes file
    ///
    /// Returns:
    ///     ReferenceGenomeMetadata
    #[staticmethod]
    fn from_path(path: &Bound<'_, PyAny>, name: &Bound<'_, PyAny>) -> PyResult<Self> {
        let path = path.to_string();
        let name = name.to_string();
        let path = Path::new(&path);

        Ok(PyReferenceGenomeMetadata {
            genome_reference: ReferenceGenomeMetadata::try_from(path, &name).unwrap()
        })
    }

    #[getter]
    fn get_name(&self) -> PyResult<String> {
        Ok(self.genome_reference.genome.clone())
    }

    #[getter]
    fn get_digest(&self) -> PyResult<String> {
        Ok(self.genome_reference.digest.clone())
    }

    #[getter]
    fn get_description(&self) -> PyResult<String> {
        Ok(self.genome_reference.description.clone())
    }

    #[getter]
    fn get_collection(&self) -> PyResult<HashMap<String, u32>> {
        Ok(self.genome_reference.collection.clone())
    }

}


#[pyclass(name = "ReferenceValidator", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyReferenceValidator {
    pub reference_genomes: ReferenceValidator,
}

#[pymethods]
impl PyReferenceValidator {
    #[new]
    /// Create a new ReferenceValidator
    ///
    /// Args:
    ///     ref_genome:
    ///
    /// Returns:
    ///     ReferenceValidator object
    // fn py_new(path: &Bound<'_, PyAny>) -> PyResult<Self> {
    fn py_new(ref_genome: Vec<PyReferenceGenomeMetadata>) -> PyResult<Self> {
        Ok(Self{reference_genomes: ReferenceValidator::new(ref_genome.into_iter().map(|r| r.genome_reference).collect())})
    }

    #[staticmethod]
    /// Create a new ReferenceValidator
    ///
    /// Args:
    ///     path: path to the folder with chrom sizes
    ///
    /// Returns:
    ///     ReferenceValidator object
    // fn py_new(path: &Bound<'_, PyAny>) -> PyResult<Self> {
    pub fn from_path() -> PyResult<Self> {
        // let path = path.to_string();
        // let path = PathBuf::from(path);
        let folder_path = Path::new("/home/bnt4me/virginia/repos/bedboss/bedboss/refgenome_validator/chrom_sizes");

        // Ok(Self{reference_genomes: ReferenceValidator::try_from(path)?})
        Ok(Self{reference_genomes: ReferenceValidator::try_from(folder_path)})
    }

    pub fn determine_compatibility(&self, rs: PyRegionSet) -> PyResult<HashMap<String,PyCompatibilityConcise>> {

        let start_time1 = std::time::Instant::now();
        let results: HashMap<String, CompatibilityConcise> = self.reference_genomes.determine_compatibility(&rs.regionset);

        let duration = start_time1.elapsed().as_secs_f64();
        println!("Rust Time taken for compatibility determination: {:?}", duration);

        let start_time = std::time::Instant::now();
        let mut return_results: HashMap<String, PyCompatibilityConcise> = HashMap::new();

        for (key, result) in results {
            return_results.insert(key, PyCompatibilityConcise{
                compatibility: result
                    }
                );
        }
        let elapsed_time = start_time.elapsed().as_secs_f64();
        println!("Rust Conversion into PyCompatibilityConcise: {:.6} seconds", elapsed_time);
        let total_time = duration + elapsed_time;
        println!("Rust Total time taken: {:.6} seconds", total_time);
        Ok(return_results)
    }

}


#[pyclass(name = "CompatibilityConcise", module = "gtars.models")]
#[derive(Debug)]
pub struct PyCompatibilityConcise {
    pub compatibility: CompatibilityConcise,
}


#[pymethods]
impl PyCompatibilityConcise {
    #[getter]
    fn get_xs(&self) -> PyResult<f64> {
        Ok(self.compatibility.xs)
    }

    #[getter]
    fn get_oobr(&self) -> PyResult<Option<f64>> {
        Ok(self.compatibility.oobr)
    }

    #[getter]
    fn get_sequence_fit(&self) -> PyResult<Option<f64>> {
        Ok(self.compatibility.sequence_fit)
    }

    #[getter]
    fn get_assigned_points(&self) -> PyResult<i32> {
        Ok(self.compatibility.assigned_points)
    }
    #[getter]
    fn get_tier_ranking(&self) -> PyResult<i32> {
        Ok(self.compatibility.tier_ranking)
    }
}
