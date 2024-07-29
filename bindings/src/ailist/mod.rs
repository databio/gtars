use gtars::ailist::{AIList, Interval};
use pyo3::{prelude::*, pyclass};

use crate::models::PyInterval;

#[pyclass(name = "AIList", module="gtars.ailist")]
struct PyAIList {
    ailist: AIList,
}

#[pymethods]
impl PyAIList {
    #[new]
    fn new(
        py_interval_list: Vec<PyRef<PyInterval>>,
        minimum_coverage_length: Option<usize>,
    ) -> PyAIList {
        let mut interval_list: Vec<Interval> = py_interval_list
            .into_iter()
            .map(|x| Interval {
                start: x.start,
                end: x.end,
            })
            .collect();
        let ailist = AIList::new(&mut interval_list, minimum_coverage_length.unwrap_or(3));
        PyAIList { ailist }
    }
    fn query(&self, py_interval: &PyInterval) -> Vec<PyInterval> {
        let interval: Interval = Interval {
            start: py_interval.start,
            end: py_interval.end,
        };
        self.ailist
            .query(&interval)
            .into_iter()
            .map(|x| PyInterval {
                start: x.start,
                end: x.end,
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
