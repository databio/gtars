use pyo3::prelude::*;

mod genome_assembly;
mod interval;
mod region;
mod region_set;
mod tss_index;

pub use self::genome_assembly::PyGenomeAssembly;
pub use self::interval::PyInterval;
pub use self::region::PyRegion;
pub use self::region_set::PyChromosomeStatistics;
pub use self::region_set::PyRegionSet;
pub use self::tss_index::PyTssIndex;

#[pymodule]
pub fn models(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRegion>()?;
    m.add_class::<PyInterval>()?;
    m.add_class::<PyRegionSet>()?;
    m.add_class::<PyChromosomeStatistics>()?;
    m.add_class::<PyGenomeAssembly>()?;
    m.add_class::<PyTssIndex>()?;
    Ok(())
}
