use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyclass(name="Encoding", module = "gtars.tokenizers")]
#[derive(Clone)]
pub struct PyEncoding {
    pub ids: Vec<u32>,
    pub attention_mask: Vec<u32>,
}

#[pyclass(extends=PyDict, name="BatchEncoding", module = "gtars.tokenizers")]
#[derive(Clone)]
pub struct PyBatchEncoding {
    #[pyo3(get)]
    pub encodings: Vec<PyEncoding>,
}