use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyDict};

use super::PyTokenizer;
use gtars::tokenizers::utils::fragments::tokenize_fragment_file;

#[pyfunction(name = "tokenize_fragment_file")]
pub fn py_tokenize_fragment_file(file: String, tokenizer: &PyTokenizer) -> PyResult<Py<PyDict>> {
    let res = tokenize_fragment_file(&file, tokenizer.inner());
    match res {
        Ok(res) => Python::with_gil(|py| {
            let py_dict = res.into_py_dict_bound(py);
            Ok(py_dict.into())
        }),
        Err(res) => Err(PyRuntimeError::new_err(res.to_string())),
    }
}
