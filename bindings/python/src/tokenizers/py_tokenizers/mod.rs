use pyo3::prelude::*;
use pyo3::types::PyType;

use anyhow::Result;

use crate::tokenizers::universe::PyUniverse;
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
    fn new(path: &str) -> Result<Self> {
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
    fn from_config(_cls: &Bound<'_, PyType>, cfg: &str) -> Result<Self> {
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
    fn from_bed(_cls: &Bound<'_, PyType>, path: &str) -> Result<Self> {
        Python::with_gil(|py| {
            let tokenizer = Tokenizer::from_bed(path)?;
            let universe = PyUniverse::from(tokenizer.get_universe().clone());
            Ok(PyTokenizer {
                tokenizer,
                universe: Py::new(py, universe)?,
            })
        })
    }

    fn tokenize(&self, regions: &Bound<'_, PyAny>) -> Result<Vec<String>> {
        // attempt to map the list to a vector of regions
        let rs = extract_regions_from_py_any(regions)?;
        let tokenized = self.tokenizer.tokenize(&rs.regions)?;
        Ok(tokenized)
    }

    // encode returns a list of ids
    fn encode(&self, regions: &Bound<'_, PyAny>) -> Result<Vec<u32>> {
        let rs = extract_regions_from_py_any(regions)?;
        let tokenized = self.tokenizer.encode(&rs.regions)?;

        Ok(tokenized)
        
    }

    fn decode(&self, ids: Vec<u32>) -> Result<Vec<String>> {
        
    }

    fn convert_ids_to_token(&self, id: u32) -> Option<String> {
        
    }

    fn convert_token_to_ids(&self, region: &String) -> Option<u32> {
        
    }

    #[getter]
    fn get_unk_token(&self) -> String {
        self.tokenizer.get_unk_token().to_owned().into()
    }

    #[getter]
    fn get_pad_token(&self) -> String {
        self.tokenizer.get_pad_token().to_owned().into()
    }

    #[getter]
    fn get_mask_token(&self) -> String {
        self.tokenizer.get_mask_token().to_owned().into()
    }

    #[getter]
    fn get_cls_token(&self) -> String {
        self.tokenizer.get_cls_token().to_owned().into()
    }

    #[getter]
    fn get_bos_token(&self) -> String {
        self.tokenizer.get_bos_token().to_owned().into()
    }

    #[getter]
    fn get_eos_token(&self) -> String {
        self.tokenizer.get_eos_token().to_owned().into()
    }

    #[getter]
    fn get_sep_token(&self) -> String {
        self.tokenizer.get_sep_token().to_owned().into()
    }

    #[getter]
    fn get_pad_token_id(&self) -> u32 {
        self.tokenizer.get_pad_token_id()
    }

    #[getter]
    fn get_mask_token_id(&self) -> u32 {
        self.tokenizer.get_mask_token_id()
    }

    #[getter]
    fn get_cls_token_id(&self) -> u32 {
        self.tokenizer.get_cls_token_id()
    }

    #[getter]
    fn get_bos_token_id(&self) -> u32 {
        self.tokenizer.get_bos_token_id()
    }

    #[getter]
    fn get_eos_token_id(&self) -> u32 {
        self.tokenizer.get_eos_token_id()
    }

    #[getter]
    fn get_sep_token_id(&self) -> u32 {
        self.tokenizer.get_sep_token_id()
    }

    #[getter]
    fn get_unk_token_id(&self) -> u32 {
        self.tokenizer.get_unk_token_id()
    }

    #[getter]
    fn get_vocab_size(&self) -> usize {
        self.tokenizer.get_vocab_size()
    }

    fn __len__(&self) -> usize {
        self.tokenizer.get_vocab_size()
    }

    fn __repr__(&self) -> String {
        format!(
            "Tokenizer({} total regions)",
            self.tokenizer.get_vocab_size()
        )
    }

    fn __call__(&self, regions: &Bound<'_, PyAny>) -> Result<> {
        self.tokenize(regions)
    }
}
