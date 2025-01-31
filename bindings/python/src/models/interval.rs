use pyo3::prelude::*;

#[pyclass(name = "Interval", module="gtars.models")]
pub struct PyInterval {
    #[pyo3(get, set)]
    pub start: u32,
    #[pyo3(get, set)]
    pub end: u32,
    #[pyo3(get, set)]
    pub data: PyObject,
}

#[pymethods]
impl PyInterval {
    #[new]
    pub fn new(start: u32, end: u32, data: PyObject) -> PyInterval {
        PyInterval { start, end, data }
    }
    pub fn __repr__(&self) -> String {
        format!("Interval({}, {})", self.start, self.end)
    }
}
