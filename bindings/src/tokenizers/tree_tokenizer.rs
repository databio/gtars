use genimtools::tokenizers::traits::SpecialTokens;
use pyo3::prelude::*;
use pyo3::types::PyList;

use anyhow::Result;

use std::path::Path;

use genimtools::common::models::{Region, RegionSet};
use genimtools::tokenizers::{Tokenizer, TreeTokenizer};

use crate::models::{PyRegion, PyTokenizedRegionSet};

#[pyclass(name = "TreeTokenizer")]
pub struct PyTreeTokenizer {
    pub tokenizer: TreeTokenizer,
}

#[pymethods]
impl PyTreeTokenizer {
    #[new]
    pub fn new(path: String) -> Result<Self> {
        let path = Path::new(&path);
        let tokenizer = TreeTokenizer::try_from(path)?;

        Ok(PyTreeTokenizer { tokenizer })
    }

    #[getter]
    pub fn unknown_token(&self) -> Result<PyRegion> {
        Ok(self.tokenizer.unknown_token().into())
    }

    #[getter]
    pub fn padding_token(&self) -> Result<PyRegion> {
        Ok(self.tokenizer.padding_token().into())
    }

    #[getter]
    pub fn mask_token(&self) -> Result<PyRegion> {
        Ok(self.tokenizer.mask_token().into())
    }

    #[getter]
    pub fn cls_token(&self) -> Result<PyRegion> {
        Ok(self.tokenizer.cls_token().into())
    }

    #[getter]
    pub fn bos_token(&self) -> Result<PyRegion> {
        Ok(self.tokenizer.bos_token().into())
    }

    #[getter]
    pub fn eos_token(&self) -> Result<PyRegion> {
        Ok(self.tokenizer.eos_token().into())
    }

    #[getter]
    pub fn sep_token(&self) -> Result<PyRegion> {
        Ok(self.tokenizer.sep_token().into())
    }

    #[getter]
    pub fn vocab_size(&self) -> usize {
        self.tokenizer.vocab_size()
    }

    // tokenize just returns a list of regions
    pub fn tokenize(&self, regions: &PyList) -> Result<Vec<PyRegion>> {
        // attempt to map the list to a vector of regions
        let regions = regions
            .iter()
            .map(|x| {
                // extract chr, start, end
                // this lets us interface any python object with chr, start, end attributes
                let chr = x.getattr("chr").unwrap().extract::<String>().unwrap();
                let start = x.getattr("start").unwrap().extract::<u32>().unwrap();
                let end = x.getattr("end").unwrap().extract::<u32>().unwrap();

                Region { chr, start, end }
            })
            .collect::<Vec<_>>();

        // create RegionSet
        let rs = RegionSet::from(regions);

        // tokenize the RegionSet
        let tokenized = self.tokenizer.tokenize_region_set(&rs);

        let regions = tokenized.into_region_vec();

        Ok(regions.into_iter().map(|r| r.into()).collect())
    }

    // __call__ returns a TokenizedRegionSet
    pub fn __call__(&self, regions: &PyList) -> Result<PyTokenizedRegionSet> {
        // attempt to map the list to a vector of regions
        let regions = regions
            .iter()
            .map(|x| {
                // extract chr, start, end
                // this lets us interface any python object with chr, start, end attributes
                let chr = x.getattr("chr").unwrap().extract::<String>().unwrap();
                let start = x.getattr("start").unwrap().extract::<u32>().unwrap();
                let end = x.getattr("end").unwrap().extract::<u32>().unwrap();

                Region { chr, start, end }
            })
            .collect::<Vec<_>>();

        // create RegionSet
        let rs = RegionSet::from(regions);

        // tokenize the RegionSet
        let tokenized = self.tokenizer.tokenize_region_set(&rs);

        Ok(tokenized.into())
    }

    // encode returns a list of ids
    pub fn encode(&self, regions: &PyList) -> Result<Vec<u32>> {
        // attempt to map the list to a vector of regions
        let regions = regions
            .iter()
            .map(|x| {
                // extract chr, start, end
                // this lets us interface any python object with chr, start, end attributes
                let chr = x.getattr("chr").unwrap().extract::<String>().unwrap();
                let start = x.getattr("start").unwrap().extract::<u32>().unwrap();
                let end = x.getattr("end").unwrap().extract::<u32>().unwrap();

                Region { chr, start, end }
            })
            .collect::<Vec<_>>();

        // create RegionSet
        let rs = RegionSet::from(regions);

        // tokenize the RegionSet
        let tokenized = self.tokenizer.tokenize_region_set(&rs);

        Ok(tokenized.ids)
    }

    pub fn decode(&self, ids: Vec<u32>) -> Result<Vec<PyRegion>> {
        let regions = ids
            .iter()
            .map(|id| self.tokenizer.universe.id_to_region[id].clone().into())
            .collect();

        Ok(regions)
    }

    pub fn vocab(&self) -> Vec<(PyRegion, u32)> {
        self.tokenizer
            .universe
            .regions
            .iter()
            .map(|r| (r.clone().into(), self.tokenizer.universe.region_to_id[r]))
            .collect()
    }

    pub fn __len__(&self) -> usize {
        self.tokenizer.universe.len()
    }

    pub fn __repr__(&self) -> String {
        format!(
            "TreeTokenizer({} total regions)",
            self.tokenizer.universe.len()
        )
    }
}
