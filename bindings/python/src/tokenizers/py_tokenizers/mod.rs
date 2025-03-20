use pyo3::prelude::*;
use pyo3::types::PyType;

use anyhow::Result;

use crate::tokenizers::universe::PyUniverse;
use crate::tokenizers::tokens::PyTokenizedRegionSet;
use crate::models::PyRegion;
use crate::utils::extract_regions_from_py_any;
use gtars::tokenizers::Tokenizer;

#[pyclass(name = "Tokenizer", module = "gtars.tokenizers")]
pub struct PyTokenizer {
    tokenizer: Tokenizer,
    universe: Py<PyUniverse>, // this is a Py-wrapped version self.tokenizer.universe for performance reasons
}

#[pymethods]
impl PyTokenizer {

    #[new]
    pub fn new(path: &str) -> Result<Self> {
        Python::with_gil(|py| {
            let tokenizer = Tokenizer::from_auto(path)?;
            let universe = PyUniverse::from(tokenizer.get_universe().clone());
            Ok(PyTokenizer {
                tokenizer,
                universe: Py::new(py, universe).unwrap(),
            })
        })
    }

    #[classmethod]
    pub fn from_config(_cls: &Bound<'_, PyType>, cfg: &str) -> Result<Self> {
        Python::with_gil(|py| {
            let tokenizer = Tokenizer::from_config(cfg)?;
            let universe = PyUniverse::from(tokenizer.get_universe().clone());
            Ok(PyTokenizer {
                tokenizer,
                universe: Py::new(py, universe)?,
            })
        })
    }

    #[classmethod]
    pub fn from_bed(_cls: &Bound<'_, PyType>, path: &str) -> Result<Self> {
        Python::with_gil(|py| {
            let tokenizer = Tokenizer::from_bed(path)?;
            let universe = PyUniverse::from(tokenizer.get_universe().clone());
            Ok(PyTokenizer {
                tokenizer,
                universe: Py::new(py, universe)?,
            })
        })
    }
    pub fn tokenize(&self, regions: &Bound<'_, PyAny>) -> Result<PyTokenizedRegionSet> {
        // attempt to map the list to a vector of regions
        let rs = extract_regions_from_py_any(regions)?;
        let tokenized = self.tokenizer.tokenize(&rs.regions)?;

        Python::with_gil(|py| {
            let py_tokenized_region_set = PyTokenizedRegionSet {
                ids: tokenized.ids,
                curr: 0,
                universe: self.universe.clone_ref(py),
            };

            Ok(py_tokenized_region_set)
        })

    }

    pub fn get_unk_token(&self) -> PyRegion {
        self.tokenizer.get_unk_token().to_owned().into()
    }

    pub fn get_pad_token(&self) -> PyRegion {
        self.tokenizer.get_pad_token().to_owned().into()
    }

    pub fn get_mask_token(&self) -> PyRegion {
        self.tokenizer.get_mask_token().to_owned().into()
    }

    pub fn get_cls_token(&self) -> PyRegion {
        self.tokenizer.get_cls_token().to_owned().into()
    }

    pub fn get_bos_token(&self) -> PyRegion {
        self.tokenizer.get_bos_token().to_owned().into()
    }

    pub fn get_eos_token(&self) -> PyRegion {
        self.tokenizer.get_eos_token().to_owned().into()
    }

    pub fn get_sep_token(&self) -> PyRegion {
        self.tokenizer.get_sep_token().to_owned().into()
    }

    pub fn get_pad_token_id(&self) -> u32 {
        self.tokenizer.get_pad_token_id()
    }

    pub fn get_mask_token_id(&self) -> u32 {
        self.tokenizer.get_mask_token_id()
    }

    pub fn get_cls_token_id(&self) -> u32 {
        self.tokenizer.get_cls_token_id()
    }

    pub fn get_bos_token_id(&self) -> u32 {
        self.tokenizer.get_bos_token_id()
    }

    pub fn get_eos_token_id(&self) -> u32 {
        self.tokenizer.get_eos_token_id()
    }

    pub fn get_sep_token_id(&self) -> u32 {
        self.tokenizer.get_sep_token_id()
    }

    pub fn get_unk_token_id(&self) -> u32 {
        self.tokenizer.get_unk_token_id()
    }

    pub fn get_vocab_size(&self) -> usize {
        self.tokenizer.get_vocab_size()
    }

    pub fn convert_id_to_region(&self, id: u32) -> Option<PyRegion> {
        Python::with_gil(|py| {
            let region = self.universe.borrow(py).convert_id_to_region(id);
            region
        })
    }

    pub fn convert_region_to_id(&self, region: &PyRegion) -> Option<u32> {
        Python::with_gil(|py| {
            let id = self.universe.borrow(py).convert_region_to_id(region);
            id
        })
    }

    // encode returns a list of ids
    pub fn encode(&self, regions: &Bound<'_, PyAny>) -> Result<Vec<u32>> {
        // attempt to map the list to a vector of regions
        let rs = extract_regions_from_py_any(regions)?;

        // tokenize the RegionSet
        let tokenized = self.tokenizer.tokenize(&rs.regions)?;

        Ok(tokenized.ids)
    }
    
    pub fn decode(&self, ids: Vec<u32>) -> Result<Vec<PyRegion>> {
        Python::with_gil(|py| {
            let regions = ids
                .iter()
                .map(|id| self.universe.borrow(py).convert_id_to_region(*id).unwrap_or(
                    PyRegion {
                        chr: "chrUNK".to_string(),
                        start: 0,
                        end: 0,
                        rest: None,
                    }
                ))
                .collect::<Vec<PyRegion>>();
            Ok(regions)
        })
    }

    pub fn __len__(&self) -> usize {
        self.tokenizer.get_vocab_size()
    }

    pub fn __repr__(&self) -> String {
        format!(
            "Tokenizer({} total regions)",
            self.tokenizer.get_vocab_size()
        )
    }

    pub fn __call__(&self, regions: &Bound<'_, PyAny>) -> Result<PyTokenizedRegionSet> {
        self.tokenize(regions)
    }
}
