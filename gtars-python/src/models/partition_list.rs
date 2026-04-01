use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::HashMap;

use crate::models::PyGeneModel;
use gtars_genomicdist::{genome_partition_list, GeneModel, PartitionList};

#[pyclass(name = "PartitionList", module = "gtars.models")]
pub struct PyPartitionList {
    pub partition_list: PartitionList,
}

#[pymethods]
impl PyPartitionList {
    #[staticmethod]
    #[pyo3(signature = (gene_model, core_prom, prox_prom, chrom_sizes = None))]
    pub fn from_gene_model(
        gene_model: &PyGeneModel,
        core_prom: u32,
        prox_prom: u32,
        chrom_sizes: Option<HashMap<String, u32>>,
    ) -> PyResult<Self> {
        let pl = genome_partition_list(
            &gene_model.gene_model,
            core_prom,
            prox_prom,
            chrom_sizes.as_ref(),
        );
        Ok(Self { partition_list: pl })
    }

    #[staticmethod]
    #[pyo3(signature = (path, core_prom, prox_prom, filter_protein_coding = true, convert_ensembl_ucsc = true, chrom_sizes = None))]
    pub fn from_gtf(
        path: &str,
        core_prom: u32,
        prox_prom: u32,
        filter_protein_coding: bool,
        convert_ensembl_ucsc: bool,
        chrom_sizes: Option<HashMap<String, u32>>,
    ) -> PyResult<Self> {
        let gene_model = GeneModel::from_gtf(path, filter_protein_coding, convert_ensembl_ucsc)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let pl = genome_partition_list(&gene_model, core_prom, prox_prom, chrom_sizes.as_ref());
        Ok(Self { partition_list: pl })
    }

    fn partition_names(&self) -> Vec<String> {
        self.partition_list
            .partitions
            .iter()
            .map(|(name, _)| name.clone())
            .collect()
    }

    fn __repr__(&self) -> String {
        let names = self.partition_names();
        format!("PartitionList(partitions={:?})", names)
    }

    fn __len__(&self) -> usize {
        self.partition_list.partitions.len()
    }
}
