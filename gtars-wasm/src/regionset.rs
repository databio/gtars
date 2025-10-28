use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};
use gtars_gd::models::RegionBin;
use gtars_gd::statistics::GenomicIntervalSetStatistics;
use wasm_bindgen::prelude::*;

use crate::models::BedEntries;

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

#[wasm_bindgen(js_name = "RegionSet")]
pub struct JsRegionSet {
    region_set: RegionSet,
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
                rest: None,
            })
            .collect::<Vec<Region>>();
        let mut region_set: RegionSet = RegionSet::from(regions);
        region_set.sort();
        Ok(JsRegionSet { region_set })
    }

    #[wasm_bindgen(getter, js_name = "numberOfRegions")]
    pub fn get_region(&self) -> i32 {
        self.region_set.len() as i32
    }

    #[wasm_bindgen(getter, js_name = "meanRegionWidth")]
    pub fn get_mean_region_width(&self) -> f64 {
        self.region_set.mean_region_width()
    }

    #[wasm_bindgen(getter, js_name = "totalNucleotides")]
    pub fn get_total_nucleotides(&self) -> i32 {
        self.region_set.nucleotides_length() as i32
    }

    #[wasm_bindgen(getter, js_name = "digest")]
    pub fn get_digest(&self) -> String {
        self.region_set.identifier()
    }

    #[wasm_bindgen(js_name = "calculateStatistics")]
    pub fn calculate_statistics(&self) -> Result<JsValue, JsValue> {
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

    #[wasm_bindgen(js_name = "calculateRegionDistribution")]
    pub fn calculate_region_distribution(&self, n_bins: u32) -> Result<JsValue, JsValue> {
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
}
