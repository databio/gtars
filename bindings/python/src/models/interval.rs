use pyo3::prelude::*;

#[pyclass(name = "Interval", module="gtars.models")]
pub struct PyInterval {
    #[pyo3(get, set)]
    pub start: u32,
    #[pyo3(get, set)]
    pub end: u32,
}

#[pymethods]
impl PyInterval {
    #[new]
    pub fn new(start: u32, end: u32) -> PyInterval {
        PyInterval { start, end }
    }
    pub fn __repr__(&self) -> String {
        format!("Interval({}, {})", self.start, self.end)
    }
}
