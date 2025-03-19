use gtars::tokenizers::config::TokenizerConfig;
use pyo3::prelude::*;
use pyo3::types::PyType;

use anyhow::Result;

use crate::tokenizers::universe::PyUniverse;
use crate::tokenizers::tokens::PyTokenizedRegionSet;
use crate::models::PyRegion;
use gtars::tokenizers::Tokenizer;

#[pyclass(name = "Tokenizer", module = "gtars.tokenizers")]
pub struct PyTokenizer {
    tokenizer: Tokenizer,
    universe: Py<PyUniverse>, // this is a Py-wrapped version self.tokenizer.universe for performance reasons
}

#[pymethods]
impl PyTokenizer {
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
}
