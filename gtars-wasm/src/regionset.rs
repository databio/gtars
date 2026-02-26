use std::collections::HashMap;

use crate::models::BedEntries;
use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::bed_classifier::classify_bed;
use gtars_genomicdist::consensus;
use gtars_genomicdist::interval_ranges::IntervalRanges;
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

    #[wasm_bindgen(js_name = "regionDistribution")]
    pub fn region_distribution(&self, n_bins: u32) -> Result<JsValue, JsValue> {
        let distribution: HashMap<String, RegionBin> =
            self.region_set.region_distribution_with_bins(n_bins);

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
        serde_wasm_bindgen::to_value(&result_vector).map_err(|e| e.into())
    }

    #[wasm_bindgen(getter, js_name = "classify")]
    pub fn classify_bed_js(&self) -> JsBedClassificationOutput {
        let output = classify_bed(&self.region_set).unwrap();
        JsBedClassificationOutput {
            bed_compliance: output.bed_compliance.clone(),
            data_format: format!("{:#?}", output.data_format),
            compliant_columns: output.compliant_columns,
            non_compliant_columns: output.non_compliant_columns,
        }
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
