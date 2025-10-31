use std::collections::HashMap;

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::types::{PyDict, PyType};
use pyo3::{prelude::*, IntoPyObjectExt};

use anyhow::Result;

use crate::utils::extract_regions_from_py_any;
use gtars_tokenizers::Tokenizer;

use super::encoding::{PyBatchEncoding, PyEncoding};

#[pyclass(name = "Tokenizer", module = "gtars.tokenizers", subclass)]
pub struct PyTokenizer {
    tokenizer: Tokenizer,
}

#[pymethods]
impl PyTokenizer {
    #[new]
    fn new(path: &str) -> Result<Self> {
        let tokenizer = Tokenizer::from_auto(path)?;
        Ok(PyTokenizer { tokenizer })
    }

    #[classmethod]
    fn from_config(_cls: &Bound<'_, PyType>, cfg: &str) -> Result<Self> {
        let tokenizer = Tokenizer::from_config(cfg)?;
        Ok(PyTokenizer { tokenizer })
    }

    #[classmethod]
    fn from_bed(_cls: &Bound<'_, PyType>, path: &str) -> Result<Self> {
        let tokenizer = Tokenizer::from_bed(path)?;
        Ok(PyTokenizer { tokenizer })
    }

    #[classmethod]
    fn from_pretrained(_cls: &Bound<'_, PyType>, path: &str) -> Result<Self> {
        let tokenizer = Tokenizer::from_pretrained(path)?;
        Ok(PyTokenizer { tokenizer })
    }

    fn tokenize(&self, regions: &Bound<'_, PyAny>) -> Result<Vec<String>> {
        // attempt to map the list to a vector of regions
        let rs = extract_regions_from_py_any(regions)?;
        let tokenized = self.tokenizer.tokenize(&rs.regions)?;
        Ok(tokenized)
    }

    // encode returns a list of ids
    fn encode(&self, tokens: &Bound<'_, PyAny>) -> Result<Py<PyAny>, PyErr> {
        Python::attach(|py| {
            // if a single token is passed
            if let Ok(token) = tokens.extract::<String>() {
                Ok(vec![self
                    .tokenizer
                    .convert_token_to_id(&token)
                    .unwrap_or(self.get_unk_token_id())]
                .into_py_any(py)?)
            }
            // if a list of tokens is passed
            else if let Ok(tokens) = tokens.extract::<Vec<String>>() {
                let ids: Vec<u32> = tokens
                    .iter()
                    .map(|token| {
                        self.tokenizer
                            .convert_token_to_id(token)
                            .unwrap_or(self.get_unk_token_id())
                    })
                    .collect();
                Ok(ids.into_py_any(py)?)
            } else {
                Err(PyValueError::new_err(
                    "Invalid input type for convert_ids_to_token",
                ))
            }
        })
    }

    fn decode(&self, ids: &Bound<'_, PyAny>) -> Result<Py<PyAny>, PyErr> {
        Python::attach(|py| {
            // if a single id is passed
            if let Ok(id) = ids.extract::<u32>() {
                Ok(vec![self
                    .tokenizer
                    .convert_id_to_token(id)
                    .unwrap_or(self.get_unk_token())]
                .into_py_any(py)?)
            }
            // if a list of ids is passed
            else if let Ok(ids) = ids.extract::<Vec<u32>>() {
                let tokens: Vec<String> = ids
                    .iter()
                    .map(|&id| {
                        self.tokenizer
                            .convert_id_to_token(id)
                            .unwrap_or(self.get_unk_token())
                    })
                    .collect();
                Ok(tokens.into_py_any(py)?)
            } else {
                Err(PyValueError::new_err(
                    "Invalid input type for convert_ids_to_token",
                ))
            }
        })
    }

    fn convert_ids_to_tokens(&self, id: &Bound<'_, PyAny>) -> Result<Py<PyAny>, PyErr> {
        Python::attach(|py| {
            // if a single id is passed
            if let Ok(id) = id.extract::<u32>() {
                Ok(self
                    .tokenizer
                    .convert_id_to_token(id)
                    .unwrap_or(self.get_unk_token())
                    .into_py_any(py)?)
            }
            // if a list of ids is passed
            else if let Ok(ids) = id.extract::<Vec<u32>>() {
                let tokens: Vec<String> = ids
                    .iter()
                    .map(|&id| {
                        self.tokenizer
                            .convert_id_to_token(id)
                            .unwrap_or(self.get_unk_token())
                    })
                    .collect();
                Ok(tokens.into_py_any(py)?)
            } else {
                Err(PyValueError::new_err(
                    "Invalid input type for convert_ids_to_token",
                ))
            }
        })
    }

    fn convert_tokens_to_ids(&self, region: &Bound<'_, PyAny>) -> Result<Py<PyAny>, PyErr> {
        Python::attach(|py| {
            // if a single token is passed
            if let Ok(token) = region.extract::<String>() {
                let id = self
                    .tokenizer
                    .convert_token_to_id(&token)
                    .unwrap_or(self.get_unk_token_id());
                Ok(id.into_py_any(py)?)
            }
            // if a list of tokens is passed
            else if let Ok(tokens) = region.extract::<Vec<String>>() {
                let ids: Vec<u32> = tokens
                    .iter()
                    .map(|token| {
                        self.tokenizer
                            .convert_token_to_id(token)
                            .unwrap_or(self.get_unk_token_id())
                    })
                    .collect();
                Ok(ids.into_py_any(py)?)
            } else {
                Err(PyValueError::new_err(
                    "Invalid input type for convert_token_to_ids",
                ))
            }
        })
    }

    #[getter]
    fn get_unk_token(&self) -> String {
        self.tokenizer.get_unk_token().to_owned()
    }

    #[getter]
    fn get_pad_token(&self) -> String {
        self.tokenizer.get_pad_token().to_owned()
    }

    #[getter]
    fn get_mask_token(&self) -> String {
        self.tokenizer.get_mask_token().to_owned()
    }

    #[getter]
    fn get_cls_token(&self) -> String {
        self.tokenizer.get_cls_token().to_owned()
    }

    #[getter]
    fn get_bos_token(&self) -> String {
        self.tokenizer.get_bos_token().to_owned()
    }

    #[getter]
    fn get_eos_token(&self) -> String {
        self.tokenizer.get_eos_token().to_owned()
    }

    #[getter]
    fn get_sep_token(&self) -> String {
        self.tokenizer.get_sep_token().to_owned()
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

    #[getter]
    fn get_special_tokens_map(&self) -> PyResult<Py<PyDict>> {
        let special_tokens = self.tokenizer.get_special_tokens();
        Python::attach(|py| {
            let dict = PyDict::new(py);
            dict.set_item("unk_token", special_tokens.unk.clone())?;
            dict.set_item("pad_token", special_tokens.pad.clone())?;
            dict.set_item("mask_token", special_tokens.mask.clone())?;
            dict.set_item("cls_token", special_tokens.cls.clone())?;
            dict.set_item("eos_token", special_tokens.eos.clone())?;
            dict.set_item("bos_token", special_tokens.bos.clone())?;
            dict.set_item("sep_token", special_tokens.sep.clone())?;
            Ok(dict.into())
        })
    }

    fn get_special_tokens_mask(&self, tokens: Vec<String>) -> Vec<bool> {
        self.inner().get_special_tokens_mask(&tokens)
    }

    fn get_vocab(&self) -> HashMap<String, u32> {
        self.tokenizer.get_vocab()
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

    fn __call__(&self, regions: &Bound<'_, PyAny>) -> Result<Py<PyAny>, PyErr> {
        Python::attach(|py| {
            let rs = extract_regions_from_py_any(regions)?;
            let encoded = self
                .tokenizer
                .encode(&rs.regions)
                .map_err(|err| PyRuntimeError::new_err(err.to_string()))?;
            let attention_mask = encoded
                .iter()
                .map(|id| if *id == self.get_pad_token_id() { 0 } else { 1 })
                .collect();

            let encoding = PyEncoding {
                ids: encoded,
                attention_mask,
            };

            PyBatchEncoding {
                input_ids: encoding.clone().ids.into_py_any(py)?,
                attention_mask: encoding.clone().attention_mask.into_py_any(py)?,
                encodings: vec![encoding],
            }
            .into_py_any(py)
        })
    }
}

impl PyTokenizer {
    pub fn inner(&self) -> &Tokenizer {
        &self.tokenizer
    }
}
