use pyo3::prelude::*;

#[pyclass(name = "Encoding", module = "gtars.tokenizers")]
#[derive(Clone)]
pub struct PyEncoding {
    pub ids: Vec<u32>,
    pub attention_mask: Vec<u32>,
}

#[pyclass(name = "BatchEncoding", module = "gtars.tokenizers")]
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
            "input_ids" => Ok(self.input_ids),
            "attention_mask" => Ok(self.attention_mask),
            _ => Err(pyo3::exceptions::PyKeyError::new_err(format!(
                "Invalid key: {}",
                key
            ))),
        }
    }
    fn __copy__<'py>(slf: PyRef<'py, Self>, py: Python<'py>) -> PyResult<Py<Self>> {
        // Create a new PyBatchEncoding instance by "cloning" the fields.
        // For PyObject fields, use `.clone_ref(py)`.
        // For Vec<PyEncoding>, use `.clone()` assuming PyEncoding implements Clone.
        Ok(Py::new(py, PyBatchEncoding {
            input_ids: slf.input_ids.clone_ref(py),       // Use clone_ref for PyObject
            attention_mask: slf.attention_mask.clone_ref(py), // Use clone_ref for PyObject
            encodings: slf.encodings.clone(),              // This works if PyEncoding is Clone
        })?)
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
