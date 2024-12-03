use std::path::Path;

use pyo3::prelude::*;
use pyo3::types::{PyAny, PyIterator};

use anyhow::Result;
use gtars::common::models::{Region, RegionSet};

// this is for internal use only
pub fn extract_regions_from_py_any(regions: &Bound<'_, PyAny>) -> Result<RegionSet> {
    // is a string?
    if let Ok(regions) = regions.extract::<String>() {
        let regions = Path::new(&regions);

        if !regions.exists() {
            return Err(pyo3::exceptions::PyFileNotFoundError::new_err(format!(
                "The file {} does not exist.",
                regions.display()
            ))
            .into());
        }

        let regions = gtars::common::utils::extract_regions_from_bed_file(regions);
        match regions {
            Ok(regions) => return Ok(RegionSet::from(regions)),
            Err(e) => return Err(pyo3::exceptions::PyValueError::new_err(e.to_string()).into()),
        }
    }

    let regions = PyIterator::from_bound_object(regions)?;

    // attempt to map the list to a vector of regions
    let regions = regions
        .iter()?
        .map(|x| {
            let x = match x {
                Ok(x) => x,
                Err(e) => anyhow::bail!("There was an error iterating over regions: {}", e),
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
    gtars::io::write_tokens_to_gtok(filename, &tokens)?;
    Ok(())
}

#[pyfunction]
pub fn read_tokens_from_gtok(filename: &str) -> PyResult<Vec<u32>> {
    let tokens = gtars::io::read_tokens_from_gtok(filename)?;
    Ok(tokens)
}

#[pyfunction]
pub fn read_tokens_from_gtok_as_strings(filename: &str) -> PyResult<Vec<String>> {
    let tokens = gtars::io::read_tokens_from_gtok(filename)?;
    let tokens = tokens.iter().map(|t| t.to_string()).collect();
    Ok(tokens)
}

#[pymodule]
pub fn utils(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(write_tokens_to_gtok))?;
    m.add_wrapped(wrap_pyfunction!(read_tokens_from_gtok))?;
    m.add_wrapped(wrap_pyfunction!(read_tokens_from_gtok_as_strings))?;
    Ok(())
}
