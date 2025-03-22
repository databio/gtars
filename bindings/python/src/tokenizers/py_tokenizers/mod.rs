use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
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

    fn convert_ids_to_token(&self, id: &Bound<'_, PyAny>) -> Result<PyObject, PyErr> {
        Python::with_gil(|py| {
            // if a single id is passed
            if let Ok(id) = id.extract::<u32>() {
                return Ok(self.tokenizer.convert_id_to_token(id).unwrap_or(self.get_unk_token()).into_py(py))
            }
            // if a list of ids is passed
            else if let Ok(ids) = id.extract::<Vec<u32>>() {
                let tokens: Vec<String> = ids.iter()
                    .map(|&id| self.tokenizer.convert_id_to_token(id).unwrap_or(self.get_unk_token()))
                    .collect();
                return Ok(tokens.into_py(py))
            } else {
                return Err(PyValueError::new_err("Invalid input type for convert_ids_to_token"))
            }
        })
    }

    fn convert_token_to_ids(&self, region: &Bound<'_, PyAny>) -> Result<PyObject, PyErr> {
        Python::with_gil(|py| {
            // if a single token is passed
            if let Ok(token) = region.extract::<String>() {
                let id = self.tokenizer.convert_token_to_id(&token).unwrap_or(self.get_unk_token_id());
                return Ok(vec![id].into_py(py))
            }
            // if a list of tokens is passed
            else if let Ok(tokens) = region.extract::<Vec<String>>() {
                let ids: Vec<u32> = tokens.iter()
                    .map(|token| self.tokenizer.convert_token_to_id(token).unwrap_or(self.get_unk_token_id()))
                    .collect();
                return Ok(ids.into_py(py))
            } else {
                return Err(PyValueError::new_err("Invalid input type for convert_token_to_ids"))
            }
        })
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
