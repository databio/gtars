use std::path::Path;

use pyo3::prelude::*;
use pyo3::types::{PyAny, PyIterator};

use anyhow::Result;
use gtars_core::models::{Region, RegionSet};

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

        let regions = gtars_core::models::RegionSet::try_from(regions);
        match regions {
            Ok(regions) => return Ok(regions),
            Err(e) => return Err(pyo3::exceptions::PyValueError::new_err(e.to_string()).into()),
        }
    }

    let regions = PyIterator::from_object(regions)?;

    // attempt to map the list to a vector of regions
    let regions = regions
        .map(|x| {
            let x = match x {
                Ok(x) => x,
                Err(e) => anyhow::bail!("There was an error iterating over regions: {}", e),
            };

            // extract chr, start, end
            // this lets us interface any python object with chr, start, end attributes
            let chr = x.getattr("chr")
                .and_then(|v| v.extract::<String>())
                .map_err(|e| anyhow::anyhow!(
                    "Region object missing or invalid 'chr' attribute (expected str): {}", e
                ))?;
            let start = x.getattr("start")
                .and_then(|v| v.extract::<u32>())
                .map_err(|e| anyhow::anyhow!(
                    "Region object missing or invalid 'start' attribute (expected u32): {}", e
                ))?;
            let end = x.getattr("end")
                .and_then(|v| v.extract::<u32>())
                .map_err(|e| anyhow::anyhow!(
                    "Region object missing or invalid 'end' attribute (expected u32): {}", e
                ))?;

            Ok(Region {
                chr,
                start,
                end,
                rest: None,
            })
        })
        .collect::<Vec<_>>();

    let regions = regions.into_iter().collect::<Result<Vec<_>>>()?;

    Ok(RegionSet::from(regions))
}

#[cfg(feature = "utils")]
#[pyfunction]
pub fn write_tokens_to_gtok(filename: &str, tokens: Vec<u32>) -> PyResult<()> {
    gtars_io::write_tokens_to_gtok(filename, &tokens)?;
    Ok(())
}

#[cfg(feature = "utils")]
#[pyfunction]
pub fn read_tokens_from_gtok(filename: &str) -> PyResult<Vec<u32>> {
    let tokens = gtars_io::read_tokens_from_gtok(filename)?;
    Ok(tokens)
}

#[cfg(feature = "utils")]
#[pyfunction]
pub fn read_tokens_from_gtok_as_strings(filename: &str) -> PyResult<Vec<String>> {
    let tokens = gtars_io::read_tokens_from_gtok(filename)?;
    let tokens = tokens.iter().map(|t| t.to_string()).collect();
    Ok(tokens)
}

#[cfg(feature = "utils")]
#[pymodule]
pub fn utils(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(write_tokens_to_gtok))?;
    m.add_wrapped(wrap_pyfunction!(read_tokens_from_gtok))?;
    m.add_wrapped(wrap_pyfunction!(read_tokens_from_gtok_as_strings))?;
    Ok(())
}
