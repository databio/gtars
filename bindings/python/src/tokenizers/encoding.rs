use pyo3::prelude::*;

#[pyclass(name = "Encoding", module = "gtars.tokenizers")]
pub struct PyEncoding {
    pub input_ids: Vec<u32>,
    pub attention_mask: Vec<u8>,
}

#[pyclass(name = "BatchEncoding", module = "gtars.tokenizers")]
pub struct PyBatchEncoding {
    pub encodings: Vec<PyEncoding>,
}