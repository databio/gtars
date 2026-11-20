//! JSON values to Python objects.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

/// `dict`, `list`, `str`, `int`, `float`, `bool` or `None`.
pub fn json_to_py(py: Python<'_>, v: &serde_json::Value) -> PyResult<Py<PyAny>> {
    use serde_json::Value as V;
    Ok(match v {
        V::Null => py.None(),
        V::Bool(b) => b.into_pyobject(py)?.to_owned().into_any().unbind(),
        V::Number(n) => {
            if let Some(i) = n.as_i64() {
                i.into_pyobject(py)?.into_any().unbind()
            } else if let Some(u) = n.as_u64() {
                u.into_pyobject(py)?.into_any().unbind()
            } else {
                n.as_f64()
                    .unwrap_or(f64::NAN)
                    .into_pyobject(py)?
                    .into_any()
                    .unbind()
            }
        }
        V::String(s) => s.into_pyobject(py)?.into_any().unbind(),
        V::Array(a) => {
            let l = PyList::empty(py);
            for x in a {
                l.append(json_to_py(py, x)?)?;
            }
            l.into_any().unbind()
        }
        V::Object(o) => {
            let d = PyDict::new(py);
            for (k, x) in o {
                d.set_item(k, json_to_py(py, x)?)?;
            }
            d.into_any().unbind()
        }
    })
}
