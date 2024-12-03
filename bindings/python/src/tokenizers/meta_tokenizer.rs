use gtars::tokenizers::traits::SpecialTokens;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use anyhow::Result;

use std::path::Path;

use gtars::common::models::RegionSet;
use gtars::tokenizers::{MetaTokenizer, Tokenizer};

use crate::models::{PyRegion, PyTokenizedRegionSet, PyUniverse};
use crate::utils::extract_regions_from_py_any;

#[pyclass(name = "MetaTokenizer", module="gtars.tokenizers")]
pub struct PyMetaTokenizer {
    pub tokenizer: MetaTokenizer,
    pub universe: Py<PyUniverse>, // this is a Py-wrapped version self.tokenizer.universe for performance reasons
}

#[pymethods]
impl PyMetaTokenizer {
    #[new]
    pub fn new(path: String) -> Result<Self> {
        Python::with_gil(|py| {
            let path = Path::new(&path);
            let tokenizer = MetaTokenizer::try_from(path)?;
            let py_universe: PyUniverse = tokenizer.universe.to_owned().into();
            let py_universe_bound = Py::new(py, py_universe)?;

            Ok(PyMetaTokenizer {
                tokenizer,
                universe: py_universe_bound,
            })
        })
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
    pub fn padding_token_id(&self) -> u32 {
        self.tokenizer.padding_token_id()
    }

    #[getter]
    pub fn mask_token_id(&self) -> u32 {
        self.tokenizer.mask_token_id()
    }

    #[getter]
    pub fn cls_token_id(&self) -> u32 {
        self.tokenizer.cls_token_id()
    }

    #[getter]
    pub fn bos_token_id(&self) -> u32 {
        self.tokenizer.bos_token_id()
    }

    #[getter]
    pub fn eos_token_id(&self) -> u32 {
        self.tokenizer.eos_token_id()
    }

    #[getter]
    pub fn sep_token_id(&self) -> u32 {
        self.tokenizer.sep_token_id()
    }

    #[getter]
    pub fn unknown_token_id(&self) -> u32 {
        self.tokenizer.unknown_token_id()
    }

    #[getter]
    pub fn vocab_size(&self) -> usize {
        self.tokenizer.vocab_size()
    }

    #[getter]
    pub fn universe(&self) -> PyUniverse {
        self.tokenizer.universe.clone().into()
    }

    // tokenize just returns a list of regions
    pub fn tokenize(&self, regions: &Bound<'_, PyAny>) -> Result<Vec<PyRegion>> {
        let rs = extract_regions_from_py_any(regions)?;

        // tokenize the RegionSet
        let tokenized = self.tokenizer.tokenize_region_set(&rs);

        let regions = tokenized.into_region_vec();

        Ok(regions.into_iter().map(|r| r.into()).collect())
    }

    pub fn tokenize_bed_file(&self, path: String) -> Result<Vec<PyRegion>> {
        let path = Path::new(&path);
        let regions = RegionSet::try_from(path)?;

        let tokenized = self.tokenizer.tokenize_region_set(&regions);

        let regions = tokenized.into_region_vec();

        Ok(regions.into_iter().map(|r| r.into()).collect())
    }

    pub fn export(&self, path: String) -> Result<()> {
        let path = Path::new(&path);
        self.tokenizer.export(path)
    }

    // __call__ returns a TokenizedRegionSet
    pub fn __call__(&self, regions: &Bound<'_, PyAny>) -> Result<PyTokenizedRegionSet> {
        // attempt to map the list to a vector of regions
        let rs = extract_regions_from_py_any(regions)?;

        // tokenize the RegionSet
        let tokenized = self.tokenizer.tokenize_region_set(&rs);

        Python::with_gil(|py| {
            let py_tokenized_region_set = PyTokenizedRegionSet {
                ids: tokenized.ids,
                curr: 0,
                universe: self.universe.clone_ref(py),
            };

            Ok(py_tokenized_region_set)
        })
    }

    // encode returns a list of ids
    pub fn encode(&self, regions: &Bound<'_, PyAny>) -> Result<Vec<u32>> {
        // attempt to map the list to a vector of regions
        let rs = extract_regions_from_py_any(regions)?;

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
            "MetaTokenizer({} total regions)",
            self.tokenizer.universe.len()
        )
    }
}
