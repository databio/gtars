use pyo3::prelude::*;
use pyo3::types::{PyAny, PyIterator};

use anyhow::Result;
use genimtools::common::models::{Region, RegionSet};

pub fn extract_regions_from_py_any(regions: &Bound<'_, PyAny>) -> Result<RegionSet> {
    let regions = PyIterator::from_bound_object(regions)?;

    // attempt to map the list to a vector of regions
    let regions = regions
        .iter()?
        .map(|x| {
            let x = match x {
                Ok(x) => x,
                Err(_) => anyhow::bail!("Error iterating over regions"),
            };

            // extract chr, start, end
            // this lets us interface any python object with chr, start, end attributes
            let chr = x.getattr("chr").unwrap().extract::<String>().unwrap();
            let start = x.getattr("start").unwrap().extract::<u32>().unwrap();
            let end = x.getattr("end").unwrap().extract::<u32>().unwrap();

            Ok(Region { chr, start, end })
        })
        .collect::<Vec<_>>();

    let regions = regions.into_iter().collect::<Result<Vec<_>>>()?;

    Ok(RegionSet::from(regions))
}

#[pyfunction]
pub fn write_tokens_to_gtok(filename: &str, tokens: Vec<u32>) -> PyResult<()> {
    genimtools::io::write_tokens_to_gtok(filename, &tokens)?;
    Ok(())
}

#[pyfunction]
pub fn read_tokens_from_gtok(filename: &str) -> PyResult<Vec<u32>> {
    let tokens = genimtools::io::read_tokens_from_gtok(filename)?;
    Ok(tokens)
}

#[pymodule]
pub fn utils(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(write_tokens_to_gtok))?;
    m.add_wrapped(wrap_pyfunction!(read_tokens_from_gtok))?;
    Ok(())
}
