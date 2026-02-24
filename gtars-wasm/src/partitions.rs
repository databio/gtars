use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::{
    calc_expected_partitions, calc_partitions, genome_partition_list, GeneModel, PartitionList,
    Strand, StrandedRegionSet,
};
use serde::Serialize;
use wasm_bindgen::prelude::*;

use crate::regionset::JsRegionSet;

// ── JsGeneModel ─────────────────────────────────────────────────────

#[wasm_bindgen(js_name = "GeneModel")]
pub struct JsGeneModel {
    model: GeneModel,
}

/// Parse a JS array of region tuples into a StrandedRegionSet.
///
/// Each element is `[chr, start, end]` or `[chr, start, end, strand]`.
/// When 4-element, strand is parsed (`+`/`-`/`*`); when 3-element, unstranded.
fn parse_stranded_regions(js_array: &JsValue) -> Result<StrandedRegionSet, JsValue> {
    let tuples: Vec<Vec<serde_json::Value>> = serde_wasm_bindgen::from_value(js_array.clone())?;
    let mut regions = Vec::with_capacity(tuples.len());
    let mut strands = Vec::with_capacity(tuples.len());

    for tuple in &tuples {
        if tuple.len() < 3 {
            return Err(JsValue::from_str("Each region tuple must have at least 3 elements: [chr, start, end]"));
        }
        let chr = tuple[0]
            .as_str()
            .ok_or_else(|| JsValue::from_str("chr must be a string"))?
            .to_string();
        let start = tuple[1]
            .as_u64()
            .ok_or_else(|| JsValue::from_str("start must be an integer"))? as u32;
        let end = tuple[2]
            .as_u64()
            .ok_or_else(|| JsValue::from_str("end must be an integer"))? as u32;

        let strand = if tuple.len() >= 4 {
            let s = tuple[3]
                .as_str()
                .ok_or_else(|| JsValue::from_str("strand must be a string"))?;
            Strand::from_char(s.chars().next().unwrap_or('.'))
        } else {
            Strand::Unstranded
        };

        regions.push(Region {
            chr,
            start,
            end,
            rest: None,
        });
        strands.push(strand);
    }

    Ok(StrandedRegionSet::new(RegionSet::from(regions), strands))
}

#[wasm_bindgen(js_class = "GeneModel")]
impl JsGeneModel {
    /// Construct a GeneModel from JS arrays of region tuples.
    ///
    /// Each array contains `[chr, start, end]` or `[chr, start, end, strand]` tuples.
    /// `threeUtr` and `fiveUtr` are optional (pass `null`/`undefined`).
    #[wasm_bindgen(constructor)]
    pub fn new(
        genes: &JsValue,
        exons: &JsValue,
        three_utr: &JsValue,
        five_utr: &JsValue,
    ) -> Result<JsGeneModel, JsValue> {
        let genes_srs = parse_stranded_regions(genes)?;
        let exons_srs = parse_stranded_regions(exons)?;

        let three_utr_srs = if three_utr.is_null() || three_utr.is_undefined() {
            None
        } else {
            let srs = parse_stranded_regions(three_utr)?;
            if srs.is_empty() {
                None
            } else {
                Some(srs)
            }
        };

        let five_utr_srs = if five_utr.is_null() || five_utr.is_undefined() {
            None
        } else {
            let srs = parse_stranded_regions(five_utr)?;
            if srs.is_empty() {
                None
            } else {
                Some(srs)
            }
        };

        Ok(JsGeneModel {
            model: GeneModel {
                genes: genes_srs,
                exons: exons_srs,
                three_utr: three_utr_srs,
                five_utr: five_utr_srs,
            },
        })
    }
}

// ── JsPartitionList ─────────────────────────────────────────────────

#[wasm_bindgen(js_name = "PartitionList")]
pub struct JsPartitionList {
    pub(crate) partition_list: PartitionList,
}

#[wasm_bindgen(js_class = "PartitionList")]
impl JsPartitionList {
    /// Build a PartitionList from a GeneModel.
    ///
    /// `chromSizes` is an optional JS object `{"chr1": 249250621, ...}` for
    /// trimming promoters at chromosome boundaries.
    #[wasm_bindgen(js_name = "fromGeneModel")]
    pub fn from_gene_model(
        model: &JsGeneModel,
        core_prom_size: u32,
        prox_prom_size: u32,
        chrom_sizes: &JsValue,
    ) -> Result<JsPartitionList, JsValue> {
        let sizes: Option<HashMap<String, u32>> =
            if chrom_sizes.is_null() || chrom_sizes.is_undefined() {
                None
            } else {
                Some(serde_wasm_bindgen::from_value(chrom_sizes.clone())?)
            };

        let pl = genome_partition_list(
            &model.model,
            core_prom_size,
            prox_prom_size,
            sizes.as_ref(),
        );

        Ok(JsPartitionList { partition_list: pl })
    }
}

// ── Free functions ──────────────────────────────────────────────────

#[derive(Serialize)]
struct PartitionCountEntry {
    name: String,
    count: u32,
}

#[derive(Serialize)]
struct PartitionResultJs {
    partitions: Vec<PartitionCountEntry>,
    total: u32,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ExpectedPartitionRowJs {
    partition: String,
    observed: f64,
    expected: f64,
    log10_oe: f64,
    pvalue: f64,
}

/// Classify query regions into partitions.
#[wasm_bindgen(js_name = "calcPartitions")]
pub fn js_calc_partitions(
    region_set: &JsRegionSet,
    partition_list: &JsPartitionList,
    bp_proportion: bool,
) -> Result<JsValue, JsValue> {
    let result = calc_partitions(&region_set.region_set, &partition_list.partition_list, bp_proportion);

    let js_result = PartitionResultJs {
        partitions: result
            .counts
            .into_iter()
            .map(|(name, count)| PartitionCountEntry { name, count })
            .collect(),
        total: result.total,
    };

    serde_wasm_bindgen::to_value(&js_result).map_err(|e| e.into())
}

/// Compute observed vs expected partition enrichment.
#[wasm_bindgen(js_name = "calcExpectedPartitions")]
pub fn js_calc_expected_partitions(
    region_set: &JsRegionSet,
    partition_list: &JsPartitionList,
    chrom_sizes: &JsValue,
    bp_proportion: bool,
) -> Result<JsValue, JsValue> {
    let sizes: HashMap<String, u32> = serde_wasm_bindgen::from_value(chrom_sizes.clone())?;

    let result = calc_expected_partitions(
        &region_set.region_set,
        &partition_list.partition_list,
        &sizes,
        bp_proportion,
    );

    let rows: Vec<ExpectedPartitionRowJs> = result
        .rows
        .into_iter()
        .map(|r| ExpectedPartitionRowJs {
            partition: r.partition,
            observed: r.observed,
            expected: r.expected,
            log10_oe: r.log10_oe,
            pvalue: r.chi_sq_pval,
        })
        .collect();

    serde_wasm_bindgen::to_value(&rows).map_err(|e| e.into())
}
