use pyo3::prelude::*;

mod interval;
mod region;
mod region_set;
mod bed_genome_validator;

pub use self::interval::PyInterval;
pub use self::region::PyRegion;
pub use self::region_set::PyRegionSet;
pub use self::bed_genome_validator::PyReferenceValidator;
pub use self::bed_genome_validator::PyReferenceGenomeMetadata;
pub use self::bed_genome_validator::PyCompatibilityConcise;

#[pymodule]
pub fn models(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRegion>()?;
    m.add_class::<PyInterval>()?;
    m.add_class::<PyRegionSet>()?;
    m.add_class::<PyReferenceValidator>()?;
    m.add_class::<PyReferenceGenomeMetadata>()?;
    m.add_class::<PyCompatibilityConcise>()?;
    Ok(())
}
