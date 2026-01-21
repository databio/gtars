use std::collections::HashMap;

use crate::models::BedEntries;
use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::bed_classifier::classify_bed;
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
}
