use pyo3::exceptions::{PyIndexError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

use crate::models::PyRegion;
use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::models::{
    ChromosomeStatistics, ClusterStats, DensityHomogeneity, DensityVector, SpacingStats,
};
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

/// Python wrapper around `gtars_genomicdist::models::SpacingStats`.
#[pyclass(name = "SpacingStats", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PySpacingStats {
    pub inner: SpacingStats,
}

/// Python wrapper around `gtars_genomicdist::models::ClusterStats`.
#[pyclass(name = "ClusterStats", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyClusterStats {
    pub inner: ClusterStats,
}

/// Python wrapper around `gtars_genomicdist::models::DensityVector`.
#[pyclass(name = "DensityVector", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyDensityVector {
    pub inner: DensityVector,
}

/// Python wrapper around `gtars_genomicdist::models::DensityHomogeneity`.
#[pyclass(name = "DensityHomogeneity", module = "gtars.models")]
#[derive(Clone, Debug)]
pub struct PyDensityHomogeneity {
    pub inner: DensityHomogeneity,
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
    pub fn from_regionset(rs: RegionSet) -> Self {
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

    /// Signed gaps between consecutive regions on each chromosome.
    ///
    /// Output length may be shorter than input region count — chromosomes with
    /// fewer than 2 regions are skipped (no neighbors to measure against). Output
    /// is not aligned 1:1 with input regions.
    fn neighbor_distances(&self) -> PyResult<Vec<i64>> {
        self.regionset
            .calc_neighbor_distances()
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Distance from each region to its nearest neighbor on the same chromosome.
    ///
    /// Output length may be shorter than input region count — chromosomes with
    /// only one region are skipped. Output is not aligned 1:1 with input regions.
    fn nearest_neighbors(&self) -> PyResult<Vec<u32>> {
        self.regionset
            .calc_nearest_neighbors()
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[pyo3(signature = (n_bins = 250, chrom_sizes = None))]
    fn distribution<'py>(
        &self,
        py: Python<'py>,
        n_bins: u32,
        chrom_sizes: Option<HashMap<String, u32>>,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let bins = match chrom_sizes {
            Some(cs) => self.regionset.region_distribution_with_chrom_sizes(n_bins, &cs),
            None => self.regionset.region_distribution_with_bins(n_bins),
        };
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

    /// Break all regions into non-overlapping disjoint pieces.
    ///
    /// Every boundary in the input becomes a boundary in the output.
    /// The result tiles the covered positions exactly, with no overlaps.
    fn disjoin(&self) -> PyResult<Self> {
        let rs = self.regionset.disjoin();
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

    /// Return gaps between regions per chromosome, bounded by chrom sizes.
    ///
    /// Emits leading gaps, inter-region gaps, trailing gaps (to chrom_size),
    /// and full-chromosome gaps for any chromosome in ``chrom_sizes`` that
    /// has no regions at all.
    fn gaps(&self, chrom_sizes: HashMap<String, u32>) -> PyResult<Self> {
        let rs = self.regionset.gaps(&chrom_sizes);
        Ok(Self::from_regionset(rs))
    }

    /// Summary statistics over inter-region spacings. Wraps
    /// ``neighbor_distances()``.
    fn inter_peak_spacing(&self) -> PySpacingStats {
        PySpacingStats {
            inner: self.regionset.calc_inter_peak_spacing(),
        }
    }

    /// Cluster-level summary statistics at a given stitching radius.
    ///
    /// Wraps ``cluster(radius_bp)``.
    ///
    /// ``min_cluster_size`` applies uniformly to every size-dependent
    /// field in the returned ``ClusterStats`` except ``max_cluster_size``
    /// (which is always the biggest cluster in the input). The
    /// arithmetic identity ``n_clusters * mean_cluster_size ==
    /// n_clustered_peaks`` holds at any threshold.
    ///
    /// The **default ``min_cluster_size = 2``** means every field
    /// describes "clusters of at least 2 peaks" — the scientifically
    /// meaningful view for enhancer clustering, super-enhancer
    /// stitching, etc. ``fraction_clustered`` is then the fraction of
    /// peaks with at least one neighbor within ``radius_bp``.
    ///
    /// Pass ``min_cluster_size=1`` to include singletons. Under this
    /// threshold ``mean_cluster_size`` becomes the simple
    /// ``total_peaks / n_clusters`` average, but ``n_clustered_peaks``
    /// degenerates to ``total_peaks`` and ``fraction_clustered`` to
    /// ``1.0`` (every peak is in a cluster of size ≥ 1).
    #[pyo3(signature = (radius_bp, min_cluster_size = 2))]
    fn peak_clusters(&self, radius_bp: u32, min_cluster_size: usize) -> PyClusterStats {
        PyClusterStats {
            inner: self
                .regionset
                .calc_peak_clusters(radius_bp, min_cluster_size),
        }
    }

    /// Dense zero-padded per-window peak count vector, ordered by
    /// karyotypic chromosome order and bin index.
    ///
    /// Args:
    ///     chrom_sizes: mapping of chromosome name → size in bp
    ///     n_bins: target bin count for the **longest chromosome** in
    ///             chrom_sizes (not the total window count). Bin width
    ///             is ``max(chrom_sizes) // n_bins``; shorter chromosomes
    ///             get ``ceil(size / bin_width)`` bins each. Total
    ///             windows returned is typically larger than n_bins.
    ///             See the DensityVector docstring for details on
    ///             per-chromosome bin width and short-contig handling.
    fn density_vector(
        &self,
        chrom_sizes: HashMap<String, u32>,
        n_bins: u32,
    ) -> PyDensityVector {
        PyDensityVector {
            inner: self.regionset.calc_density_vector(&chrom_sizes, n_bins),
        }
    }

    /// Summary statistics over the dense per-window count vector
    /// (variance, CV, Gini). Wraps ``density_vector()``.
    ///
    /// ``n_bins`` is the target bin count for the longest chromosome,
    /// not the total window count — see :py:meth:`density_vector` for
    /// the full semantic. Short contigs in ``chrom_sizes`` each
    /// contribute a narrow single-bin entry which dilutes
    /// ``mean_count``, inflates ``n_windows``, and raises ``gini``.
    fn density_homogeneity(
        &self,
        chrom_sizes: HashMap<String, u32>,
        n_bins: u32,
    ) -> PyDensityHomogeneity {
        PyDensityHomogeneity {
            inner: self
                .regionset
                .calc_density_homogeneity(&chrom_sizes, n_bins),
        }
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

// ── Spatial-arrangement summary getters ─────────────────────────────────

#[pymethods]
impl PySpacingStats {
    #[getter]
    fn n_gaps(&self) -> usize {
        self.inner.n_gaps
    }
    #[getter]
    fn mean(&self) -> f64 {
        self.inner.mean
    }
    #[getter]
    fn median(&self) -> f64 {
        self.inner.median
    }
    #[getter]
    fn std(&self) -> f64 {
        self.inner.std
    }
    #[getter]
    fn iqr(&self) -> f64 {
        self.inner.iqr
    }
    #[getter]
    fn log_mean(&self) -> f64 {
        self.inner.log_mean
    }
    #[getter]
    fn log_std(&self) -> f64 {
        self.inner.log_std
    }
    fn __repr__(&self) -> String {
        format!(
            "SpacingStats(n_gaps={}, mean={:.3}, median={:.3}, std={:.3}, iqr={:.3}, log_mean={:.3}, log_std={:.3})",
            self.inner.n_gaps,
            self.inner.mean,
            self.inner.median,
            self.inner.std,
            self.inner.iqr,
            self.inner.log_mean,
            self.inner.log_std
        )
    }
}

#[pymethods]
impl PyClusterStats {
    #[getter]
    fn radius_bp(&self) -> u32 {
        self.inner.radius_bp
    }
    #[getter]
    fn n_clusters(&self) -> usize {
        self.inner.n_clusters
    }
    #[getter]
    fn n_clustered_peaks(&self) -> usize {
        self.inner.n_clustered_peaks
    }
    #[getter]
    fn mean_cluster_size(&self) -> f64 {
        self.inner.mean_cluster_size
    }
    #[getter]
    fn max_cluster_size(&self) -> usize {
        self.inner.max_cluster_size
    }
    #[getter]
    fn fraction_clustered(&self) -> f64 {
        self.inner.fraction_clustered
    }
    fn __repr__(&self) -> String {
        format!(
            "ClusterStats(radius_bp={}, n_clusters={}, n_clustered_peaks={}, mean_cluster_size={:.3}, max_cluster_size={}, fraction_clustered={:.3})",
            self.inner.radius_bp,
            self.inner.n_clusters,
            self.inner.n_clustered_peaks,
            self.inner.mean_cluster_size,
            self.inner.max_cluster_size,
            self.inner.fraction_clustered
        )
    }
}

#[pymethods]
impl PyDensityVector {
    #[getter]
    fn n_bins(&self) -> u32 {
        self.inner.n_bins
    }
    #[getter]
    fn bin_width(&self) -> u32 {
        self.inner.bin_width
    }
    /// Dense per-window count vector (numpy-compatible via list→array).
    #[getter]
    fn counts(&self) -> Vec<u32> {
        self.inner.counts.clone()
    }
    /// List of (chromosome_name, start_index) per chromosome in karyotypic
    /// order. Slice ``counts[start_index : next_chrom_start_index]`` to
    /// recover per-chromosome subvectors.
    #[getter]
    fn chrom_offsets(&self) -> Vec<(String, usize)> {
        self.inner.chrom_offsets.clone()
    }
    fn __len__(&self) -> usize {
        self.inner.counts.len()
    }
    fn __repr__(&self) -> String {
        format!(
            "DensityVector(n_bins={}, bin_width={}, len={}, chroms={})",
            self.inner.n_bins,
            self.inner.bin_width,
            self.inner.counts.len(),
            self.inner.chrom_offsets.len()
        )
    }
}

#[pymethods]
impl PyDensityHomogeneity {
    #[getter]
    fn bin_width(&self) -> u32 {
        self.inner.bin_width
    }
    #[getter]
    fn n_windows(&self) -> usize {
        self.inner.n_windows
    }
    #[getter]
    fn n_nonzero_windows(&self) -> usize {
        self.inner.n_nonzero_windows
    }
    #[getter]
    fn mean_count(&self) -> f64 {
        self.inner.mean_count
    }
    #[getter]
    fn variance(&self) -> f64 {
        self.inner.variance
    }
    #[getter]
    fn cv(&self) -> f64 {
        self.inner.cv
    }
    #[getter]
    fn gini(&self) -> f64 {
        self.inner.gini
    }
    fn __repr__(&self) -> String {
        format!(
            "DensityHomogeneity(bin_width={}, n_windows={}, n_nonzero={}, mean={:.3}, variance={:.3}, cv={:.3}, gini={:.3})",
            self.inner.bin_width,
            self.inner.n_windows,
            self.inner.n_nonzero_windows,
            self.inner.mean_count,
            self.inner.variance,
            self.inner.cv,
            self.inner.gini
        )
    }
}
