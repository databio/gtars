use pyo3::exceptions::{PyIndexError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

use crate::models::PyRegion;
use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::models::ChromosomeStatistics;
use gtars_genomicdist::statistics::GenomicIntervalSetStatistics;
use gtars_genomicdist::IntervalRanges;
use gtars_overlaprs::RegionSetOverlaps;

#[pyclass(name = "ChromosomeStatistics", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyChromosomeStatistics {
    pub chromosome: String,
    pub number_of_regions: u32,
    pub start_nucleotide_position: u32,
    pub end_nucleotide_position: u32,
    pub minimum_region_length: u32,
    pub maximum_region_length: u32,
    pub mean_region_length: f64,
    pub median_region_length: f64,
}

#[pyclass(name = "RegionSet", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyRegionSet {
    pub regionset: RegionSet,
    curr: usize,
    identifier: Option<String>,
    pub strands: Vec<String>,
}

impl PyRegionSet {
    /// Create a PyRegionSet from a RegionSet with default "*" strands
    fn from_regionset(rs: RegionSet) -> Self {
        let n = rs.regions.len();
        Self {
            regionset: rs,
            curr: 0,
            identifier: None,
            strands: vec!["*".to_string(); n],
        }
    }
}

impl From<Vec<PyRegion>> for PyRegionSet {
    fn from(value: Vec<PyRegion>) -> Self {
        let mut rust_regions: Vec<Region> = Vec::new();
        for region in value {
            rust_regions.push(Region {
                chr: region.chr,
                start: region.start,
                end: region.end,
                rest: region.rest,
            })
        }
        let n = rust_regions.len();
        PyRegionSet {
            regionset: RegionSet::from(rust_regions),
            curr: 0,
            identifier: None,
            strands: vec!["*".to_string(); n],
        }
    }
}

#[pymethods]
impl PyRegionSet {
    #[new]
    /// Create a new RegionSet object
    ///
    /// Args:
    ///     path: path to the bed, or bed.gz file
    ///
    /// Returns:
    ///     RegionSet object
    fn py_new(path: &Bound<'_, PyAny>) -> PyResult<Self> {
        let path = path.to_string();
        let regionset =
            RegionSet::try_from(path).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(Self::from_regionset(regionset))
    }

    /// Alternate constructor from a list of PyRegion
    ///
    /// Args:
    ///     regions: a list/vec of PyRegion objects
    ///     strands: optional list of strand strings ("+", "-", "*")
    ///
    /// Returns:
    ///     RegionSet object
    #[staticmethod]
    #[pyo3(signature = (regions, strands = None))]
    fn from_regions(regions: Vec<PyRegion>, strands: Option<Vec<String>>) -> PyResult<Self> {
        let n = regions.len();
        let rust_regions: Vec<Region> = regions.into_iter().map(PyRegion::into_region).collect();
        let strands = strands.unwrap_or_else(|| vec!["*".to_string(); n]);
        if strands.len() != n {
            return Err(PyValueError::new_err(format!(
                "strands length ({}) must match regions length ({})",
                strands.len(),
                n
            )));
        }
        Ok(Self {
            regionset: RegionSet::from(rust_regions),
            curr: 0,
            identifier: None,
            strands,
        })
    }

    /// Create a RegionSet from parallel vectors
    ///
    /// Args:
    ///     chrs: list of chromosome names
    ///     starts: list of start positions (0-based)
    ///     ends: list of end positions (half-open)
    ///     strands: optional list of strand strings ("+", "-", "*")
    ///
    /// Returns:
    ///     RegionSet object
    #[staticmethod]
    #[pyo3(signature = (chrs, starts, ends, strands = None))]
    fn from_vectors(
        chrs: Vec<String>,
        starts: Vec<u32>,
        ends: Vec<u32>,
        strands: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let n = chrs.len();
        if starts.len() != n || ends.len() != n {
            return Err(PyValueError::new_err(
                "chrs, starts, and ends must have the same length",
            ));
        }
        let regions: Vec<Region> = chrs
            .into_iter()
            .zip(starts.into_iter().zip(ends.into_iter()))
            .map(|(chr, (start, end))| Region {
                chr,
                start,
                end,
                rest: None,
            })
            .collect();
        let strands = strands.unwrap_or_else(|| vec!["*".to_string(); n]);
        if strands.len() != n {
            return Err(PyValueError::new_err(format!(
                "strands length ({}) must match regions length ({})",
                strands.len(),
                n
            )));
        }
        Ok(Self {
            regionset: RegionSet::from(regions),
            curr: 0,
            identifier: None,
            strands,
        })
    }

    #[getter]
    fn get_strands(&self) -> Vec<String> {
        self.strands.clone()
    }

    #[getter]
    fn get_identifier(&mut self) -> PyResult<String> {
        if self.identifier.is_none() {
            self.identifier = Some(self.regionset.identifier());
        }
        // Safe: the block above guarantees Some
        self.identifier.clone().ok_or_else(|| {
            pyo3::exceptions::PyRuntimeError::new_err("Failed to compute identifier")
        })
    }

    #[getter]
    fn get_file_digest(&self) -> PyResult<String> {
        Ok(self.regionset.file_digest())
    }

    #[getter]
    fn get_path(&self) -> PyResult<String> {
        let path = self
            .regionset
            .path
            .as_ref()
            .ok_or_else(|| PyValueError::new_err("RegionSet has no associated path"))?;
        let path_str = path
            .to_str()
            .ok_or_else(|| PyValueError::new_err("Path contains invalid UTF-8 characters"))?;
        Ok(path_str.to_string())
    }

    #[getter]
    fn get_header(&self) -> PyResult<Option<String>> {
        Ok(self.regionset.header.clone())
    }

    fn is_empty(&self) -> PyResult<bool> {
        Ok(self.regionset.is_empty())
    }

    fn __repr__(&self) -> String {
        format!("RegionSet with {} regions.", self.regionset.len())
    }

    fn __str__(&self) -> String {
        format!("RegionSet with {} regions.", self.regionset.len())
    }

    fn __len__(&self) -> usize {
        self.regionset.len()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    pub fn __next__(&mut self) -> Option<PyRegion> {
        if self.curr < self.regionset.regions.len() {
            let region = self.regionset.regions[self.curr].clone();
            self.curr += 1;

            Some(PyRegion {
                chr: region.chr,
                start: region.start,
                end: region.end,
                rest: region.rest,
            })
        } else {
            self.curr = 0;
            None
        }
    }

    pub fn __getitem__(&self, indx: isize) -> PyResult<PyRegion> {
        let len = self.regionset.regions.len() as isize;

        let indx = if indx < 0 {
            len + indx // Convert negative index to positive
        } else {
            indx
        };

        if indx < 0 || indx >= len {
            Err(PyIndexError::new_err("Index out of bounds"))
        } else {
            let r = &self.regionset.regions[indx as usize];
            Ok(PyRegion {
                chr: r.chr.clone(),
                start: r.start,
                end: r.end,
                rest: r.rest.clone(),
            })
        }
    }

    fn to_bigbed(
        &self,
        out_path: &Bound<'_, PyAny>,
        chrom_size: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        self.regionset
            .to_bigbed(out_path.to_string(), chrom_size.to_string())
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn to_bed(&self, path: &Bound<'_, PyAny>) -> PyResult<()> {
        self.regionset
            .to_bed(path.to_string())
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn to_bed_gz(&self, path: &Bound<'_, PyAny>) -> PyResult<()> {
        self.regionset
            .to_bed_gz(path.to_string())
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }

    fn sort(&mut self) -> PyResult<()> {
        self.regionset.sort();
        Ok(())
    }

    fn region_widths(&self) -> Vec<u32> {
        self.regionset.region_widths()
    }

    fn widths(&self) -> Vec<u32> {
        self.regionset.calc_widths()
    }

    fn neighbor_distances(&self) -> PyResult<Vec<Option<f64>>> {
        let dists = self
            .regionset
            .calc_neighbor_distances()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(dists
            .into_iter()
            .map(|d| {
                if d == i64::MAX {
                    None
                } else {
                    Some(d as f64)
                }
            })
            .collect())
    }

    fn nearest_neighbors(&self) -> PyResult<Vec<Option<f64>>> {
        let dists = self
            .regionset
            .calc_nearest_neighbors()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(dists
            .into_iter()
            .map(|d| {
                if d == u32::MAX {
                    None
                } else {
                    Some(d as f64)
                }
            })
            .collect())
    }

    #[pyo3(signature = (n_bins = 250))]
    fn distribution<'py>(&self, py: Python<'py>, n_bins: u32) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let bins = self.regionset.region_distribution_with_bins(n_bins);
        let mut sorted_bins: Vec<_> = bins.into_values().collect();
        sorted_bins.sort_by(|a, b| {
            a.chr
                .cmp(&b.chr)
                .then(a.start.cmp(&b.start))
        });
        let mut result = Vec::with_capacity(sorted_bins.len());
        for bin in sorted_bins {
            let dict = PyDict::new(py);
            dict.set_item("chr", &bin.chr)?;
            dict.set_item("start", bin.start)?;
            dict.set_item("end", bin.end)?;
            dict.set_item("n", bin.n)?;
            dict.set_item("rid", bin.rid)?;
            result.push(dict);
        }
        Ok(result)
    }

    fn trim(&self, chrom_sizes: HashMap<String, u32>) -> PyResult<Self> {
        let rs = self.regionset.trim(&chrom_sizes);
        // Trim may drop/reorder regions — strand mapping lost
        Ok(Self::from_regionset(rs))
    }

    fn promoters(&self, upstream: u32, downstream: u32) -> PyResult<Self> {
        let rs = self.regionset.promoters(upstream, downstream);
        // Same regions, same order — preserve strand
        Ok(Self {
            regionset: rs,
            curr: 0,
            identifier: None,
            strands: self.strands.clone(),
        })
    }

    fn reduce(&self) -> PyResult<Self> {
        let rs = self.regionset.reduce();
        // Merge op — strand lost
        Ok(Self::from_regionset(rs))
    }

    fn setdiff(&self, other: &PyRegionSet) -> PyResult<Self> {
        let rs = self.regionset.setdiff(&other.regionset);
        Ok(Self::from_regionset(rs))
    }

    fn pintersect(&self, other: &PyRegionSet) -> PyResult<Self> {
        let rs = self.regionset.pintersect(&other.regionset);
        // Positional — take strand from self
        Ok(Self {
            regionset: rs,
            curr: 0,
            identifier: None,
            strands: self.strands.clone(),
        })
    }

    fn concat(&self, other: &PyRegionSet) -> PyResult<Self> {
        let rs = self.regionset.concat(&other.regionset);
        // Concatenate strand vectors
        let mut strands = self.strands.clone();
        strands.extend(other.strands.iter().cloned());
        Ok(Self {
            regionset: rs,
            curr: 0,
            identifier: None,
            strands,
        })
    }

    fn union(&self, other: &PyRegionSet) -> PyResult<Self> {
        let rs = self.regionset.union(&other.regionset);
        Ok(Self::from_regionset(rs))
    }

    fn jaccard(&self, other: &PyRegionSet) -> f64 {
        self.regionset.jaccard(&other.regionset)
    }

    fn coverage(&self, other: &PyRegionSet) -> f64 {
        self.regionset.coverage(&other.regionset)
    }

    fn overlap_coefficient(&self, other: &PyRegionSet) -> f64 {
        self.regionset.overlap_coefficient(&other.regionset)
    }

    fn mean_region_width(&self) -> f64 {
        self.regionset.mean_region_width()
    }

    fn get_max_end_per_chr(&self) -> HashMap<String, u32> {
        self.regionset.get_max_end_per_chr()
    }

    fn get_nucleotide_length(&self) -> u32 {
        self.regionset.nucleotides_length()
    }

    /// Return a new RegionSet containing only regions that overlap at least
    /// one region in other.
    fn subset_by_overlaps(&self, other: &PyRegionSet) -> PyResult<Self> {
        let rs = self.regionset.subset_by_overlaps(&other.regionset, None);
        Ok(Self::from_regionset(rs))
    }

    /// Return a list of overlap counts, one per region in self.
    fn count_overlaps(&self, other: &PyRegionSet) -> Vec<usize> {
        RegionSetOverlaps::count_overlaps(&self.regionset, &other.regionset, None)
    }

    /// Return a list of booleans indicating whether each region overlaps
    /// any region in other.
    fn any_overlaps(&self, other: &PyRegionSet) -> Vec<bool> {
        self.regionset.any_overlaps(&other.regionset, None)
    }

    /// Return a list of lists of indices into other that overlap each
    /// region in self.
    fn find_overlaps(&self, other: &PyRegionSet) -> Vec<Vec<usize>> {
        RegionSetOverlaps::find_overlaps(&self.regionset, &other.regionset, None)
    }

    fn intersect_all(&self, other: &PyRegionSet) -> PyResult<Self> {
        let rs = self.regionset.intersect_all(&other.regionset);
        Ok(Self::from_regionset(rs))
    }

    fn subtract(&self, other: &PyRegionSet) -> PyResult<Self> {
        let rs = self.regionset.subtract(&other.regionset);
        Ok(Self::from_regionset(rs))
    }

    fn closest(&self, other: &PyRegionSet) -> PyResult<Vec<(usize, usize, i64)>> {
        let result = self.regionset.closest(&other.regionset);
        Ok(result)
    }

    #[pyo3(signature = (max_gap = 0))]
    fn cluster(&self, max_gap: u32) -> PyResult<Vec<u32>> {
        let assignments = self.regionset.cluster(max_gap);
        Ok(assignments)
    }

    fn chromosome_statistics(&self) -> HashMap<String, PyChromosomeStatistics> {
        let mut result: HashMap<String, PyChromosomeStatistics> = HashMap::new();
        let stats: HashMap<String, ChromosomeStatistics> = self.regionset.chromosome_statistics();

        for (key, value) in &stats {
            result.insert(
                key.clone(),
                PyChromosomeStatistics {
                    chromosome: value.chromosome.clone(),
                    number_of_regions: value.number_of_regions,
                    minimum_region_length: value.minimum_region_length,
                    maximum_region_length: value.maximum_region_length,
                    mean_region_length: value.mean_region_length,
                    median_region_length: value.median_region_length,
                    start_nucleotide_position: value.start_nucleotide_position,
                    end_nucleotide_position: value.end_nucleotide_position,
                },
            );
        }

        result
    }
}

#[pymethods]
impl PyChromosomeStatistics {
    #[getter]
    fn get_chromosome(&self) -> &str {
        &self.chromosome
    }

    #[getter]
    fn get_number_of_regions(&self) -> u32 {
        self.number_of_regions
    }

    #[getter]
    fn get_start_nucleotide_position(&self) -> u32 {
        self.start_nucleotide_position
    }

    #[getter]
    fn get_end_nucleotide_position(&self) -> u32 {
        self.end_nucleotide_position
    }

    #[getter]
    fn get_minimum_region_length(&self) -> u32 {
        self.minimum_region_length
    }

    #[getter]
    fn get_maximum_region_length(&self) -> u32 {
        self.maximum_region_length
    }

    #[getter]
    fn get_mean_region_length(&self) -> f64 {
        self.mean_region_length
    }

    #[getter]
    fn get_median_region_length(&self) -> f64 {
        self.median_region_length
    }
}
