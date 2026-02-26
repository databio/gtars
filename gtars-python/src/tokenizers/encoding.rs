use pyo3::prelude::*;

#[pyclass(name = "Encoding", module = "gtars.tokenizers")]
#[derive(Clone)]
pub struct PyEncoding {
    pub ids: Vec<u32>,
    pub attention_mask: Vec<u32>,
}

#[pyclass(name = "BatchEncoding", module = "gtars.tokenizers")]
pub struct PyBatchEncoding {
    pub input_ids: Py<PyAny>,
    pub attention_mask: Py<PyAny>,
    #[pyo3(get)]
    pub encodings: Vec<PyEncoding>,
}

#[pymethods]
impl PyBatchEncoding {
    fn __getitem__(&self, key: &str) -> PyResult<Py<PyAny>> {
        Python::attach(|py| match key {
            "input_ids" => Ok(self.input_ids.clone_ref(py)),
            "attention_mask" => Ok(self.attention_mask.clone_ref(py)),
            _ => Err(pyo3::exceptions::PyKeyError::new_err(format!(
                "Invalid key: {}",
                key
            ))),
        })
    }

    fn __len__(&self) -> PyResult<usize> {
        Ok(self.encodings.len())
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "BatchEncoding(input_ids={:?}, attention_mask={:?})",
            self.input_ids, self.attention_mask
        ))
    }
}
