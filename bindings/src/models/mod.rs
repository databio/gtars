use pyo3::prelude::*;

mod interval;
mod region;
mod region_set;
mod universe;

pub use self::interval::PyInterval;
pub use self::region::{PyRegion, PyTokenizedRegion};
pub use self::region_set::{PyRegionSet, PyTokenizedRegionSet};
pub use self::universe::PyUniverse;

#[pymodule]
pub fn models(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRegion>()?;
    m.add_class::<PyTokenizedRegion>()?;
    m.add_class::<PyTokenizedRegionSet>()?;
    m.add_class::<PyInterval>()?;
    m.add_class::<PyRegionSet>()?;
    Ok(())
}
