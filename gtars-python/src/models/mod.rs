use pyo3::prelude::*;

mod gda;
mod gene_model;
mod genome_assembly;
mod interval;
mod partition_list;
mod region;
mod region_set;
mod region_set_list;
mod signal_matrix;
pub(crate) mod tss_index;

pub use self::gda::PyGenomicDistAnnotation;
pub use self::gene_model::PyGeneModel;
pub use self::genome_assembly::{PyBinaryGenomeAssembly, PyGenomeAssembly};
pub use self::interval::PyInterval;
pub use self::partition_list::PyPartitionList;
pub use self::region::PyRegion;
pub use self::region_set::{PyChromosomeStatistics, PyRegionSet};
pub use self::region_set_list::PyRegionSetList;
pub use self::signal_matrix::PySignalMatrix;
pub use self::tss_index::PyTssIndex;

#[pymodule]
pub fn models(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRegion>()?;
    m.add_class::<PyInterval>()?;
    m.add_class::<PyRegionSet>()?;
    m.add_class::<PyChromosomeStatistics>()?;
    m.add_class::<PyGenomeAssembly>()?;
    m.add_class::<PyBinaryGenomeAssembly>()?;
    m.add_class::<PyTssIndex>()?;
    m.add_class::<PyGeneModel>()?;
    m.add_class::<PyPartitionList>()?;
    m.add_class::<PySignalMatrix>()?;
    m.add_class::<PyGenomicDistAnnotation>()?;
    m.add_class::<PyRegionSetList>()?;
    Ok(())
}
