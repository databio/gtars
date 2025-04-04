use pyo3::prelude::*;

mod interval;
mod region;
mod region_set;

pub use self::interval::PyInterval;
pub use self::region::PyRegion;
pub use self::region_set::PyRegionSet;

#[pymodule]
pub fn models(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRegion>()?;
    m.add_class::<PyInterval>()?;
    m.add_class::<PyRegionSet>()?;
    Ok(())
}
