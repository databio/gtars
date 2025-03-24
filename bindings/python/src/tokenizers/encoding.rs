use pyo3::prelude::*;

#[pyclass(name="Encoding", module = "gtars.tokenizers")]
#[derive(Clone)]
pub struct PyEncoding {
    pub ids: Vec<u32>,
    pub attention_mask: Vec<u32>,
}

#[pyclass(name="BatchEncoding", module = "gtars.tokenizers")]
#[derive(Clone)]
pub struct PyBatchEncoding {
    pub input_ids: PyObject,
    pub attention_mask: PyObject,
    #[pyo3(get)]
    pub encodings: Vec<PyEncoding>,
}

#[pymethods]
impl PyBatchEncoding {
    fn __getitem__(&self, key: &str) -> PyResult<PyObject> {
        match key {
            "input_ids" => Ok(self.input_ids.clone()),
            "attention_mask" => Ok(self.attention_mask.clone()),
            _ => Err(pyo3::exceptions::PyKeyError::new_err(format!("Invalid key: {}", key))),
        }
    }
}