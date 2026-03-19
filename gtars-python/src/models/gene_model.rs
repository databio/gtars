use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use gtars_genomicdist::GeneModel;

#[pyclass(name = "GeneModel", module = "gtars.models")]
pub struct PyGeneModel {
    pub gene_model: GeneModel,
}

#[pymethods]
impl PyGeneModel {
    #[staticmethod]
    #[pyo3(signature = (path, filter_protein_coding = true, convert_ensembl_ucsc = true))]
    pub fn from_gtf(
        path: &str,
        filter_protein_coding: bool,
        convert_ensembl_ucsc: bool,
    ) -> PyResult<Self> {
        let gene_model = GeneModel::from_gtf(path, filter_protein_coding, convert_ensembl_ucsc)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { gene_model })
    }

    #[getter]
    fn n_genes(&self) -> usize {
        self.gene_model.genes.len()
    }

    #[getter]
    fn n_exons(&self) -> usize {
        self.gene_model.exons.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "GeneModel(n_genes={}, n_exons={})",
            self.gene_model.genes.len(),
            self.gene_model.exons.len()
        )
    }
}
