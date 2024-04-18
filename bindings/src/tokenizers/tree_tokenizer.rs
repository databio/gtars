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

    pub fn token_to_id(&self, region: &PyRegion) -> usize {
        let region = Region {
            chr: region.chr.to_string(),
            start: region.start,
            end: region.end,
        };
        self.tokenizer
            .universe
            .convert_region_to_id(&region)
            .unwrap() as usize
    }

    pub fn tokens_to_id(&self, regions: &PyList) -> Vec<usize> {
        let ids = regions
            .iter()
            .map(|x| {
                // extract chr, start, end
                // this lets us interface any python object with chr, start, end attributes
                let chr = x.getattr("chr").unwrap().extract::<String>().unwrap();
                let start = x.getattr("start").unwrap().extract::<u32>().unwrap();
                let end = x.getattr("end").unwrap().extract::<u32>().unwrap();

                let region = Region { chr, start, end };

                self.tokenizer
                    .universe
                    .convert_region_to_id(&region)
                    .unwrap() as usize
            })
            .collect::<Vec<_>>();

        ids
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

    ///
    /// Tokenize a list of regions
    ///
    /// # Arguments
    /// - `regions` - a list of regions
    ///
    /// # Returns
    /// A PyTokenizedRegionSet that contains regions, and ids
    pub fn tokenize(&self, regions: &PyList) -> Result<PyTokenizedRegionSet> {
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

        // tokenize
        let tokenized_regions = self.tokenizer.tokenize_region_set(&rs);

        // create pytokenizedregionset
        let regions = tokenized_regions
            .into_iter()
            .map(|x| PyRegion {
                chr: x.chr,
                start: x.start,
                end: x.end,
            })
            .collect::<Vec<_>>();

        let ids = tokenized_regions.to_region_ids();

        Ok(PyTokenizedRegionSet::new(regions, ids))
    }

    pub fn tokenize_bed_file(&self, path: String) -> Result<PyTokenizedRegionSet> {
        let bed_file = Path::new(&path);
        let tokens = self.tokenizer.tokenize_bed_file(bed_file)?;

        let regions = tokens
            .into_iter()
            .map(|x| PyRegion {
                chr: x.chr,
                start: x.start,
                end: x.end,
            })
            .collect::<Vec<_>>();

        let ids = tokens.to_region_ids();

        Ok(PyTokenizedRegionSet::new(regions, ids))
    }
}
