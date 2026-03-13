use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::HashMap;

use crate::models::gene_model::PyGeneModel;
use crate::models::partition_list::PyPartitionList;
use crate::models::tss_index::PyTssIndex;
use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::models::{Strand, TssIndex};
use gtars_genomicdist::{genome_partition_list, GenomicDistAnnotation};

#[pyclass(name = "GenomicDistAnnotation", module = "gtars.models")]
pub struct PyGenomicDistAnnotation {
    inner: GenomicDistAnnotation,
}

#[pymethods]
impl PyGenomicDistAnnotation {
    #[staticmethod]
    pub fn load_bin(path: &str) -> PyResult<Self> {
        let gda = GenomicDistAnnotation::load_bin(path)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner: gda })
    }

    #[staticmethod]
    #[pyo3(signature = (path, filter_protein_coding = true, convert_ensembl_ucsc = true))]
    pub fn from_gtf(
        path: &str,
        filter_protein_coding: bool,
        convert_ensembl_ucsc: bool,
    ) -> PyResult<Self> {
        let gene_model = gtars_genomicdist::GeneModel::from_gtf(
            path,
            filter_protein_coding,
            convert_ensembl_ucsc,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self {
            inner: GenomicDistAnnotation { gene_model },
        })
    }

    fn gene_model(&self) -> PyGeneModel {
        PyGeneModel {
            gene_model: self.inner.gene_model.clone(),
        }
    }

    #[pyo3(signature = (core_prom, prox_prom, chrom_sizes = None))]
    fn partition_list(
        &self,
        core_prom: u32,
        prox_prom: u32,
        chrom_sizes: Option<HashMap<String, u32>>,
    ) -> PyPartitionList {
        let pl = genome_partition_list(
            &self.inner.gene_model,
            core_prom,
            prox_prom,
            chrom_sizes.as_ref(),
        );
        PyPartitionList { partition_list: pl }
    }

    fn tss_index(&self) -> PyResult<PyTssIndex> {
        let model = &self.inner.gene_model;
        let tss_regions: Vec<Region> = model
            .genes
            .inner
            .regions
            .iter()
            .zip(model.genes.strands.iter())
            .map(|(r, strand)| {
                let tss_pos = match strand {
                    Strand::Minus => r.end.saturating_sub(1),
                    _ => r.start,
                };
                Region {
                    chr: r.chr.clone(),
                    start: tss_pos,
                    end: tss_pos + 1,
                    rest: None,
                }
            })
            .collect();
        let tss_rs = RegionSet {
            regions: tss_regions,
            header: None,
            path: None,
        };
        let tss_index = TssIndex::try_from(tss_rs)
            .map_err(|e| PyValueError::new_err(format!("Building TssIndex from GDA: {}", e)))?;
        Ok(PyTssIndex { tss_index })
    }

    fn __repr__(&self) -> String {
        format!(
            "GenomicDistAnnotation(n_genes={}, n_exons={})",
            self.inner.gene_model.genes.len(),
            self.inner.gene_model.exons.len()
        )
    }
}
