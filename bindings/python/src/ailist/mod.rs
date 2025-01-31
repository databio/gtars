use gtars::ailist::{AIList, Interval};
use pyo3::{prelude::*, pyclass};

use crate::models::PyInterval;

#[pyclass(name = "AIList", module="gtars.ailist")]
struct PyAIList {
    ailist: AIList<PyObject>,
}

#[pymethods]
impl PyAIList {
    #[new]
    fn new(
        py_interval_list: Vec<PyRef<PyInterval>>,
        minimum_coverage_length: Option<usize>,
    ) -> PyAIList {
        let interval_list = py_interval_list
            .into_iter()
            .map(|x| Interval {
                start: x.start,
                end: x.end,
                data: x.data.clone()
            })
            .collect();
        let ailist = AIList::new(interval_list, minimum_coverage_length.unwrap_or(3));
        PyAIList { ailist }
    }

    fn query(&self, py_interval: &PyInterval) -> Vec<PyInterval> {
        self.ailist
            .query(py_interval.start, py_interval.end)
            .into_iter()
            .map(|x| PyInterval {
                start: x.start,
                end: x.end,
                data: x.data
            })
            .collect()
    }
}

/// A Python module implemented in Rust.
#[pymodule]
pub fn ailist(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyInterval>()?;
    m.add_class::<PyAIList>()?;
    Ok(())
}
