use std::collections::HashMap;

use crate::models::BedEntries;
use gtars_core::models::{Region, RegionSet, RegionSetList};
use gtars_genomicdist::bed_classifier::classify_bed;
use gtars_genomicdist::consensus;
use gtars_genomicdist::interval_ranges::{IntervalRanges, RegionSetListOps};
use gtars_genomicdist::models::RegionBin;
use gtars_genomicdist::statistics::GenomicIntervalSetStatistics;
use wasm_bindgen::prelude::*;

#[wasm_bindgen(js_name = "ChromosomeStatistics")]
#[derive(serde::Serialize)]
pub struct JsChromosomeStatistics {
    chromosome: String,
    number_of_regions: u32,
    minimum_region_length: u32,
    maximum_region_length: u32,
    mean_region_length: f64,
    median_region_length: f64,
    start_nucleotide_position: u32,
    end_nucleotide_position: u32,
}

#[wasm_bindgen(js_class = "ChromosomeStatistics")]
impl JsChromosomeStatistics {
    #[wasm_bindgen(getter)]
    pub fn chromosome(&self) -> String {
        self.chromosome.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn start_nucleotide_position(&self) -> u32 {
        self.start_nucleotide_position
    }
    #[wasm_bindgen(getter)]
    pub fn end_nucleotide_position(&self) -> u32 {
        self.end_nucleotide_position
    }

    #[wasm_bindgen(getter)]
    pub fn number_of_regions(&self) -> u32 {
        self.number_of_regions
    }

    #[wasm_bindgen(getter)]
    pub fn minimum_region_length(&self) -> u32 {
        self.minimum_region_length
    }

    #[wasm_bindgen(getter)]
    pub fn maximum_region_length(&self) -> u32 {
        self.maximum_region_length
    }

    #[wasm_bindgen(getter)]
    pub fn mean_region_length(&self) -> f64 {
        self.mean_region_length
    }

    #[wasm_bindgen(getter)]
    pub fn median_region_length(&self) -> f64 {
        self.median_region_length
    }
}

#[wasm_bindgen(js_name = "RegionDistribution")]
#[derive(serde::Serialize)]
pub struct JsRegionDistribution {
    chr: String,
    start: u32,
    end: u32,
    n: u32,
    rid: u32,
}

#[wasm_bindgen(js_name = "BedClassificationOutput")]
#[derive(serde::Serialize)]
pub struct JsBedClassificationOutput {
    bed_compliance: String,
    data_format: String,
    compliant_columns: usize,
    non_compliant_columns: usize,
}

#[wasm_bindgen(js_class = "BedClassificationOutput")]
impl JsBedClassificationOutput {
    #[wasm_bindgen(getter)]
    pub fn bed_compliance(&self) -> String {
        self.bed_compliance.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn data_format(&self) -> String {
        self.data_format.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn compliant_columns(&self) -> usize {
        self.compliant_columns
    }
    #[wasm_bindgen(getter)]
    pub fn non_compliant_columns(&self) -> usize {
        self.non_compliant_columns
    }
}

#[wasm_bindgen(js_name = "RegionSet")]
pub struct JsRegionSet {
    pub(crate) region_set: RegionSet,
}

#[wasm_bindgen(js_class = "RegionSet")]
impl JsRegionSet {
    #[wasm_bindgen(constructor)]
    pub fn new(regions: &JsValue) -> Result<JsRegionSet, JsValue> {
        let regions: BedEntries = serde_wasm_bindgen::from_value(regions.clone())?;
        let regions: Vec<Region> = regions
            .0
            .into_iter()
            .map(|be| Region {
                chr: be.0,
                start: be.1,
                end: be.2,
                rest: Some(be.3),
            })
            .collect::<Vec<Region>>();
        let mut region_set: RegionSet = RegionSet::from(regions);
        region_set.sort();
        Ok(JsRegionSet { region_set })
    }

    #[wasm_bindgen(getter, js_name = "firstRegion")]
    pub fn get_first(&self) -> String {
        format!("{:#?}", self.region_set.regions.first())
    }

    #[wasm_bindgen(getter, js_name = "numberOfRegions")]
    pub fn get_region_number(&self) -> i32 {
        self.region_set.len() as i32
    }

    #[wasm_bindgen(getter, js_name = "meanRegionWidth")]
    pub fn mean_region_width(&self) -> f64 {
        self.region_set.mean_region_width()
    }

    #[wasm_bindgen(getter, js_name = "nucleotidesLength")]
    pub fn nucleotides_length(&self) -> i32 {
        self.region_set.nucleotides_length() as i32
    }

    #[wasm_bindgen(getter, js_name = "identifier")]
    pub fn identifier(&self) -> String {
        self.region_set.identifier()
    }

    #[wasm_bindgen(js_name = "chromosomeStatistics")]
    pub fn chromosome_statistics(&self) -> Result<JsValue, JsValue> {
        let stats = self.region_set.chromosome_statistics();
        let mut result: HashMap<String, JsChromosomeStatistics> = HashMap::new();

        for (key, value) in stats {
            result.insert(
                key.clone(),
                JsChromosomeStatistics {
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

        serde_wasm_bindgen::to_value(&result).map_err(|e| e.into())
    }

    /// Region distribution across genomic bins.
    ///
    /// When `chrom_sizes` is provided (JS object: `{chr: length, ...}`), per-chromosome
    /// bin sizes are derived from the reference genome (bin_size = chrom_size / n_bins
    /// per chrom). Outputs are comparable across BED files and aligned with reference
    /// genome positions.
    ///
    /// When `chrom_sizes` is null/undefined, bin size is derived from the BED file's
    /// observed max end coordinate — outputs will NOT be comparable across files.
    ///
    /// Returns a JS object `{ bins: [...], outOfRange: {...} }` where `outOfRange`
    /// is populated only when `chrom_sizes` is provided and some regions' midpoints
    /// fell beyond the stated chromosome size (assembly mismatch). `outOfRange` is
    /// always an empty object when `chrom_sizes` is not provided.
    #[wasm_bindgen(js_name = "regionDistribution")]
    pub fn region_distribution(
        &self,
        n_bins: u32,
        chrom_sizes: JsValue,
    ) -> Result<JsValue, JsValue> {
        let (distribution, out_of_range): (HashMap<String, RegionBin>, HashMap<String, u32>) =
            if chrom_sizes.is_null() || chrom_sizes.is_undefined() {
                (
                    self.region_set.region_distribution_with_bins(n_bins),
                    HashMap::new(),
                )
            } else {
                let cs: HashMap<String, u32> = serde_wasm_bindgen::from_value(chrom_sizes)
                    .map_err(|e| JsValue::from_str(&format!("chrom_sizes: {}", e)))?;
                let r = self
                    .region_set
                    .region_distribution_with_chrom_sizes(n_bins, &cs);
                (r.bins, r.out_of_range)
            };

        let mut result_vector: Vec<JsRegionDistribution> = vec![];

        for value in distribution.values() {
            result_vector.push(JsRegionDistribution {
                chr: value.chr.clone(),
                start: value.start,
                end: value.end,
                n: value.n,
                rid: value.rid,
            })
        }

        #[derive(serde::Serialize)]
        struct Output<'a> {
            bins: &'a Vec<JsRegionDistribution>,
            #[serde(rename = "outOfRange")]
            out_of_range: HashMap<String, u32>,
        }
        serde_wasm_bindgen::to_value(&Output {
            bins: &result_vector,
            out_of_range,
        })
        .map_err(|e| e.into())
    }

    #[wasm_bindgen(getter, js_name = "classify")]
    pub fn classify_bed_js(&self) -> Result<JsBedClassificationOutput, JsValue> {
        let output = classify_bed(&self.region_set)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        Ok(JsBedClassificationOutput {
            bed_compliance: output.bed_compliance.clone(),
            data_format: format!("{:#?}", output.data_format),
            compliant_columns: output.compliant_columns,
            non_compliant_columns: output.non_compliant_columns,
        })
    }

    // ── Statistics methods ───────────────────────────────────────────

    #[wasm_bindgen(js_name = "calcWidths")]
    pub fn calc_widths(&self) -> Vec<u32> {
        self.region_set.calc_widths()
    }

    #[wasm_bindgen(js_name = "calcNeighborDistances")]
    pub fn calc_neighbor_distances(&self) -> Result<JsValue, JsValue> {
        let distances = self
            .region_set
            .calc_neighbor_distances()
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        // Convert Vec<i64> to Vec<f64> to avoid BigInt in JS
        let as_f64: Vec<f64> = distances.iter().map(|&d| d as f64).collect();
        serde_wasm_bindgen::to_value(&as_f64).map_err(|e| e.into())
    }

    #[wasm_bindgen(js_name = "calcNearestNeighbors")]
    pub fn calc_nearest_neighbors(&self) -> Result<JsValue, JsValue> {
        let nearest = self
            .region_set
            .calc_nearest_neighbors()
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        serde_wasm_bindgen::to_value(&nearest).map_err(|e| e.into())
    }

    // ── Interval range methods ──────────────────────────────────────

    #[wasm_bindgen(js_name = "trim")]
    pub fn trim(&self, chrom_sizes: &JsValue) -> Result<JsRegionSet, JsValue> {
        let sizes: HashMap<String, u32> = serde_wasm_bindgen::from_value(chrom_sizes.clone())?;
        let trimmed = self.region_set.trim(&sizes);
        Ok(JsRegionSet { region_set: trimmed })
    }

    #[wasm_bindgen(js_name = "promoters")]
    pub fn promoters(&self, upstream: u32, downstream: u32) -> JsRegionSet {
        let result = self.region_set.promoters(upstream, downstream);
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "reduce")]
    pub fn reduce(&self) -> JsRegionSet {
        let result = self.region_set.reduce();
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "setdiff")]
    pub fn setdiff(&self, other: &JsRegionSet) -> JsRegionSet {
        let result = self.region_set.setdiff(&other.region_set);
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "pintersect")]
    pub fn pintersect(&self, other: &JsRegionSet) -> JsRegionSet {
        let result = self.region_set.pintersect(&other.region_set);
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "concat")]
    pub fn concat(&self, other: &JsRegionSet) -> JsRegionSet {
        let result = self.region_set.concat(&other.region_set);
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "union")]
    pub fn union(&self, other: &JsRegionSet) -> JsRegionSet {
        let result = self.region_set.union(&other.region_set);
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "jaccard")]
    pub fn jaccard(&self, other: &JsRegionSet) -> f64 {
        self.region_set.jaccard(&other.region_set)
    }

    #[wasm_bindgen(js_name = "shift")]
    pub fn shift(&self, offset: i64) -> JsRegionSet {
        let result = self.region_set.shift(offset);
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "flank")]
    pub fn flank(&self, width: u32, use_start: bool, both: bool) -> JsRegionSet {
        let result = self.region_set.flank(width, use_start, both);
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "resize")]
    pub fn resize(&self, width: u32, fix: &str) -> JsRegionSet {
        let result = self.region_set.resize(width, fix);
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "narrow")]
    pub fn narrow(
        &self,
        start: Option<u32>,
        end: Option<u32>,
        width: Option<u32>,
    ) -> JsRegionSet {
        let result = self.region_set.narrow(start, end, width);
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "disjoin")]
    pub fn disjoin(&self) -> JsRegionSet {
        let result = self.region_set.disjoin();
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "gaps")]
    pub fn gaps(&self) -> JsRegionSet {
        let result = self.region_set.gaps();
        JsRegionSet { region_set: result }
    }

    #[wasm_bindgen(js_name = "intersect")]
    pub fn intersect(&self, other: &JsRegionSet) -> JsRegionSet {
        let result = self.region_set.intersect(&other.region_set);
        JsRegionSet { region_set: result }
    }

}

// =========================================================================
// JsRegionSetList
// =========================================================================

/// A collection of RegionSets — the gtars equivalent of GRangesList.
#[wasm_bindgen(js_name = "RegionSetList")]
pub struct JsRegionSetList {
    pub(crate) inner: RegionSetList,
}

#[wasm_bindgen(js_class = "RegionSetList")]
impl JsRegionSetList {
    /// Create an empty RegionSetList. Use `add()` to populate it.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        JsRegionSetList {
            inner: RegionSetList::new(Vec::new()),
        }
    }

    /// Add a RegionSet to this list, with an optional name.
    pub fn add(&mut self, set: &JsRegionSet, name: Option<String>) {
        let idx = self.inner.region_sets.len();
        self.inner.region_sets.push(set.region_set.clone());
        let name_val = name.unwrap_or_else(|| format!("set_{}", idx));
        self.inner
            .names
            .get_or_insert_with(Vec::new)
            .push(name_val);
    }

    /// Number of region sets in this list.
    #[wasm_bindgen(getter)]
    pub fn length(&self) -> usize {
        self.inner.len()
    }

    /// Get a region set by 0-based index.
    pub fn get(&self, index: usize) -> Result<JsRegionSet, JsValue> {
        self.inner
            .get(index)
            .map(|rs| JsRegionSet {
                region_set: rs.clone(),
            })
            .ok_or_else(|| {
                JsValue::from_str(&format!(
                    "Index {} out of range (list has {} region sets)",
                    index,
                    self.inner.len()
                ))
            })
    }

    /// Flatten all region sets into a single RegionSet (no merging).
    pub fn concat(&self) -> JsRegionSet {
        JsRegionSet {
            region_set: self.inner.concat(),
        }
    }

    /// Get the names of the region sets, or null if unnamed.
    #[wasm_bindgen(getter)]
    pub fn names(&self) -> JsValue {
        match &self.inner.names {
            Some(names) => serde_wasm_bindgen::to_value(names).unwrap_or(JsValue::NULL),
            None => JsValue::NULL,
        }
    }

    /// Build a RegionSetList directly from arrays of BED entries.
    ///
    /// @param entries - Array of arrays: [[[chr, start, end, rest], ...], ...]
    /// @param names - Optional array of names, one per set
    #[wasm_bindgen(js_name = "fromEntries")]
    pub fn from_entries(entries: &JsValue, names: &JsValue) -> Result<JsRegionSetList, JsValue> {
        let all_entries: Vec<BedEntries> =
            serde_wasm_bindgen::from_value(entries.clone())?;
        let names_vec: Option<Vec<String>> = if names.is_null() || names.is_undefined() {
            None
        } else {
            Some(serde_wasm_bindgen::from_value(names.clone())?)
        };

        let mut sets = Vec::with_capacity(all_entries.len());
        for bed_entries in all_entries {
            let regions: Vec<Region> = bed_entries
                .0
                .into_iter()
                .map(|be| Region {
                    chr: be.0,
                    start: be.1,
                    end: be.2,
                    rest: Some(be.3),
                })
                .collect();
            let mut rs = RegionSet::from(regions);
            rs.sort();
            sets.push(rs);
        }

        Ok(JsRegionSetList {
            inner: RegionSetList {
                region_sets: sets,
                names: names_vec,
                path: None,
            },
        })
    }

    /// Number of regions in the set at the given index.
    #[wasm_bindgen(js_name = "regionCount")]
    pub fn region_count(&self, index: usize) -> Result<u32, JsValue> {
        self.inner.region_count(index)
            .ok_or_else(|| JsValue::from_str("Index out of range"))
    }

    /// Number of overlapping regions between two sets by index.
    #[wasm_bindgen(js_name = "pintersectCount")]
    pub fn pintersect_count(&self, i: usize, j: usize) -> Result<u32, JsValue> {
        self.inner.pintersect_count(i, j)
            .ok_or_else(|| JsValue::from_str("Index out of range"))
    }

    /// Jaccard similarity between two sets by index.
    #[wasm_bindgen(js_name = "jaccardAt")]
    pub fn jaccard_at(&self, i: usize, j: usize) -> Result<f64, JsValue> {
        self.inner.jaccard_at(i, j)
            .ok_or_else(|| JsValue::from_str("Index out of range"))
    }

    /// Union of two sets by index.
    #[wasm_bindgen(js_name = "unionAt")]
    pub fn union_at(&self, i: usize, j: usize) -> Result<JsRegionSet, JsValue> {
        self.inner.union_at(i, j)
            .map(|rs| JsRegionSet { region_set: rs })
            .ok_or_else(|| JsValue::from_str("Index out of range"))
    }

    /// Setdiff of two sets by index (set[i] minus set[j]).
    #[wasm_bindgen(js_name = "setdiffAt")]
    pub fn setdiff_at(&self, i: usize, j: usize) -> Result<JsRegionSet, JsValue> {
        self.inner.setdiff_at(i, j)
            .map(|rs| JsRegionSet { region_set: rs })
            .ok_or_else(|| JsValue::from_str("Index out of range"))
    }

    /// Union of all sets except the one at the given index.
    #[wasm_bindgen(js_name = "unionExcept")]
    pub fn union_except(&self, skip: usize) -> Result<JsRegionSet, JsValue> {
        self.inner.union_except(skip)
            .map(|rs| JsRegionSet { region_set: rs })
            .ok_or_else(|| JsValue::from_str("Index out of range or list too small"))
    }

    /// Compute all N union-except results in O(n) via prefix/suffix.
    /// Returns { union: RegionSet, excepts: RegionSet[] }.
    #[wasm_bindgen(js_name = "bulkUnionExcept")]
    pub fn bulk_union_except(&self) -> Result<JsValue, JsValue> {
        let (full_union, excepts) = self.inner.bulk_union_except()
            .ok_or_else(|| JsValue::from_str("Need at least 2 sets"))?;

        #[derive(serde::Serialize)]
        struct BulkResult {
            union_regions: u32,
            union_nucleotides: u32,
            except_unique: Vec<u32>,
        }

        // For each file, compute setdiff(file_i, union_except_i).len()
        let mut except_unique = Vec::with_capacity(excepts.len());
        for (i, ue) in excepts.iter().enumerate() {
            if let Some(rs) = self.inner.get(i) {
                except_unique.push(rs.setdiff(ue).len() as u32);
            } else {
                except_unique.push(0);
            }
        }

        let result = BulkResult {
            union_regions: full_union.len() as u32,
            union_nucleotides: full_union.nucleotides_length() as u32,
            except_unique,
        };
        serde_wasm_bindgen::to_value(&result).map_err(|e| e.into())
    }

    /// Union of all sets.
    #[wasm_bindgen(js_name = "unionAll")]
    pub fn union_all(&self) -> Result<JsRegionSet, JsValue> {
        self.inner.union_all()
            .map(|rs| JsRegionSet { region_set: rs })
            .ok_or_else(|| JsValue::from_str("Empty list"))
    }

    /// Intersection of all sets.
    #[wasm_bindgen(js_name = "intersectAll")]
    pub fn intersect_all(&self) -> Result<JsRegionSet, JsValue> {
        self.inner.intersect_all()
            .map(|rs| JsRegionSet { region_set: rs })
            .ok_or_else(|| JsValue::from_str("Empty list"))
    }

    /// Compute pairwise Jaccard similarity for all pairs of region sets.
    ///
    /// Returns { matrix: number[][], names: string[] | null }.
    #[wasm_bindgen(js_name = "pairwiseJaccard")]
    pub fn pairwise_jaccard(&self) -> Result<JsValue, JsValue> {
        let sets: Vec<RegionSet> = (0..self.inner.len())
            .filter_map(|i| self.inner.get(i).cloned())
            .collect();
        let n = sets.len();
        let flat = gtars_genomicdist::pairwise_jaccard(&sets);
        let matrix: Vec<Vec<f64>> = (0..n)
            .map(|i| flat[i * n..(i + 1) * n].to_vec())
            .collect();

        #[derive(serde::Serialize)]
        struct PairwiseResult {
            matrix: Vec<Vec<f64>>,
            names: Option<Vec<String>>,
        }

        let result = PairwiseResult {
            matrix,
            names: self.inner.names.clone(),
        };
        serde_wasm_bindgen::to_value(&result).map_err(|e| e.into())
    }
}

/// Builder for computing consensus regions from multiple RegionSet objects.
///
/// Usage from JS:
/// ```js
/// const cb = new ConsensusBuilder();
/// cb.add(rs1);
/// cb.add(rs2);
/// cb.add(rs3);
/// const result = cb.compute(); // [{chr, start, end, count}, ...]
/// ```
#[wasm_bindgen(js_name = "ConsensusBuilder")]
pub struct JsConsensusBuilder {
    sets: Vec<RegionSet>,
}

#[wasm_bindgen(js_class = "ConsensusBuilder")]
impl JsConsensusBuilder {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        JsConsensusBuilder { sets: Vec::new() }
    }

    pub fn add(&mut self, set: &JsRegionSet) {
        self.sets.push(set.region_set.clone());
    }

    pub fn compute(&self) -> Result<JsValue, JsValue> {
        let result = consensus::consensus(&self.sets);

        #[derive(serde::Serialize)]
        struct JsConsensusRegion {
            chr: String,
            start: u32,
            end: u32,
            count: u32,
        }

        let js_result: Vec<JsConsensusRegion> = result
            .into_iter()
            .map(|r| JsConsensusRegion {
                chr: r.chr,
                start: r.start,
                end: r.end,
                count: r.count,
            })
            .collect();

        serde_wasm_bindgen::to_value(&js_result).map_err(|e| e.into())
    }
}
