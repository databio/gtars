use pyo3::prelude::*;
use pyo3::types::PyList;

use anyhow::Result;

use std::collections::HashMap;
use std::path::Path;

use genimtools::common::consts::{
    MASK_CHR, MASK_END, MASK_START, PAD_CHR, PAD_END, PAD_START, UNKNOWN_CHR, UNKNOWN_END,
    UNKNOWN_START,
};
use genimtools::common::models::{Region, RegionSet};
use genimtools::tokenizers::{Tokenizer, TreeTokenizer};

use crate::models::{PyRegion, PyTokenizedRegion, PyTokenizedRegionSet, PyUniverse};

#[pyclass(name = "TreeTokenizer")]
pub struct PyTreeTokenizer {
    pub tokenizer: TreeTokenizer,
}

#[pymethods]
impl PyTreeTokenizer {
    #[new]
    pub fn new(path: String) -> Result<Self> {
        let path = Path::new(&path);
        let tokenizer = TreeTokenizer::try_from(path)?;

        Ok(PyTreeTokenizer { tokenizer })
    }

    #[getter]
    pub fn unknown_token(&self) -> Result<PyTokenizedRegion> {
        let region = PyRegion {
            chr: UNKNOWN_CHR.to_string(),
            start: UNKNOWN_START as u32,
            end: UNKNOWN_END as u32,
        };

        let id = self
            .tokenizer
            .universe
            .convert_region_to_id(&region.to_region());

        match id {
            Some(id) => Ok(PyTokenizedRegion { region, id }),
            None => {
                anyhow::bail!("Something went wrong -- could not find an id for the unknown token!")
            }
        }
    }

    #[getter]
    pub fn padding_token(&self) -> Result<PyTokenizedRegion> {
        let region = PyRegion {
            chr: PAD_CHR.to_string(),
            start: PAD_START as u32,
            end: PAD_END as u32,
        };

        let id = self
            .tokenizer
            .universe
            .region_to_id
            .get(&region.to_region())
            .unwrap()
            .to_owned();

        Ok(PyTokenizedRegion { region, id })
    }

    #[getter]
    pub fn mask_token(&self) -> Result<PyTokenizedRegion> {
        let region = PyRegion {
            chr: MASK_CHR.to_string(),
            start: MASK_START as u32,
            end: MASK_END as u32,
        };

        let id = self
            .tokenizer
            .universe
            .region_to_id
            .get(&region.to_region())
            .unwrap()
            .to_owned();

        Ok(PyTokenizedRegion { region, id })
    }

    #[getter]
    pub fn universe(&self) -> Result<PyUniverse> {
        let regions = self
            .tokenizer
            .universe
            .regions
            .iter()
            .map(|x| PyRegion::new(x.chr.clone(), x.start, x.end))
            .collect::<Vec<_>>();
        let region_to_id = self
            .tokenizer
            .universe
            .region_to_id
            .iter()
            .map(|(k, v)| (PyRegion::new(k.chr.clone(), k.start, k.end), v.to_owned()))
            .collect::<HashMap<_, _>>();
        let length = self.tokenizer.universe.len();
        Ok(PyUniverse {
            regions,
            region_to_id,
            length,
        })
    }

    pub fn token_to_id(&self, region: &PyRegion) -> usize {
        let region = Region {
            chr: region.chr.to_string(),
            start: region.start,
            end: region.end,
        };
        self.tokenizer
            .universe
            .convert_region_to_id(&region)
            .unwrap() as usize
    }

    pub fn tokens_to_id(&self, regions: &PyList) -> Vec<usize> {
        let ids = regions
            .iter()
            .map(|x| {
                // extract chr, start, end
                // this lets us interface any python object with chr, start, end attributes
                let chr = x.getattr("chr").unwrap().extract::<String>().unwrap();
                let start = x.getattr("start").unwrap().extract::<u32>().unwrap();
                let end = x.getattr("end").unwrap().extract::<u32>().unwrap();

                let region = Region { chr, start, end };

                self.tokenizer
                    .universe
                    .convert_region_to_id(&region)
                    .unwrap() as usize
            })
            .collect::<Vec<_>>();

        ids
    }

    pub fn __len__(&self) -> usize {
        self.tokenizer.universe.len() as usize
    }

    pub fn __repr__(&self) -> String {
        format!(
            "TreeTokenizer({} total regions)",
            self.tokenizer.universe.len()
        )
    }

    ///
    /// Tokenize a list of regions
    ///
    /// # Arguments
    /// - `regions` - a list of regions
    ///
    /// # Returns
    /// A PyTokenizedRegionSet that contains regions, and ids
    pub fn tokenize(&self, regions: &PyList) -> Result<PyTokenizedRegionSet> {
        // attempt to map the list to a vector of regions
        let regions = regions
            .iter()
            .map(|x| {
                // extract chr, start, end
                // this lets us interface any python object with chr, start, end attributes
                let chr = x.getattr("chr").unwrap().extract::<String>().unwrap();
                let start = x.getattr("start").unwrap().extract::<u32>().unwrap();
                let end = x.getattr("end").unwrap().extract::<u32>().unwrap();

                Region { chr, start, end }
            })
            .collect::<Vec<_>>();

        // create RegionSet
        let rs = RegionSet::from(regions);

        // tokenize
        let tokenized_regions = self.tokenizer.tokenize_region_set(&rs);

        // create pytokenizedregionset
        match tokenized_regions {
            Some(tokenized_regions) => {
                let regions = tokenized_regions
                    .into_iter()
                    .map(|x| PyRegion {
                        chr: x.chr,
                        start: x.start,
                        end: x.end,
                    })
                    .collect::<Vec<_>>();

                let ids = tokenized_regions.to_region_ids();

                Ok(PyTokenizedRegionSet::new(regions, ids))
            }
            // return error if tokenized_regions is None
            None => anyhow::bail!("Failed to tokenize regions",),
        }
    }

    pub fn tokenize_bed_file(&self, path: String) -> Result<PyTokenizedRegionSet> {
        let bed_file = Path::new(&path);
        let tokens = self.tokenizer.tokenize_bed_file(bed_file);

        match tokens {
            Some(tokens) => {
                let regions = tokens
                    .into_iter()
                    .map(|x| PyRegion {
                        chr: x.chr,
                        start: x.start,
                        end: x.end,
                    })
                    .collect::<Vec<_>>();

                let ids = tokens.to_region_ids();

                Ok(PyTokenizedRegionSet::new(regions, ids))
            }
            None => anyhow::bail!(format!("Error parsing the bedfile: {}", path)),
        }
    }
}
