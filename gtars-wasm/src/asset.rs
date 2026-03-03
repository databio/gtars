use std::collections::HashMap;

use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::models::{Strand, TssIndex};
use gtars_genomicdist::GenomicDistAnnotation;
use wasm_bindgen::prelude::*;

use crate::partitions::{JsGeneModel, JsPartitionList};
use crate::tss::JsTssIndex;

// ── JsGenomicDistAnnotation ────────────────────────────────────────

#[wasm_bindgen(js_name = "GenomicDistAnnotation")]
pub struct JsGenomicDistAnnotation {
    annotation: GenomicDistAnnotation,
}

#[wasm_bindgen(js_class = "GenomicDistAnnotation")]
impl JsGenomicDistAnnotation {
    /// Load a GenomicDistAnnotation from GDA binary bytes.
    #[wasm_bindgen(js_name = "fromBin")]
    pub fn from_bin(bytes: &[u8]) -> Result<JsGenomicDistAnnotation, JsValue> {
        let annotation = GenomicDistAnnotation::load_bin_from_bytes(bytes)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        Ok(JsGenomicDistAnnotation { annotation })
    }

    /// Get the gene model as a JsGeneModel.
    ///
    /// Note: this clones the gene model data since WASM ownership requires it.
    #[wasm_bindgen(js_name = "geneModel")]
    pub fn gene_model(&self) -> JsGeneModel {
        JsGeneModel::from_inner(self.annotation.gene_model.clone())
    }

    /// Derive a TssIndex from gene starts + strand.
    ///
    /// Plus/Unstranded genes → TSS at start, Minus genes → TSS at end.
    #[wasm_bindgen(js_name = "tssIndex")]
    pub fn tss_index(&self) -> Result<JsTssIndex, JsValue> {
        let model = &self.annotation.gene_model;
        let tss_regions: Vec<Region> = model
            .genes
            .inner
            .regions
            .iter()
            .zip(model.genes.strands.iter())
            .map(|(r, strand)| {
                let tss_pos = match strand {
                    Strand::Minus => r.end.saturating_sub(1),
                    _ => r.start,
                };
                Region {
                    chr: r.chr.clone(),
                    start: tss_pos,
                    end: tss_pos + 1,
                    rest: None,
                }
            })
            .collect();

        let tss_rs = RegionSet {
            regions: tss_regions,
            header: None,
            path: None,
        };

        let index = TssIndex::try_from(tss_rs)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;

        Ok(JsTssIndex::from_inner(index))
    }

    /// Build a PartitionList from the gene model and provided chrom sizes.
    ///
    /// `chrom_sizes` should be a JS object like `{"chr1": 249250621, ...}`.
    #[wasm_bindgen(js_name = "partitionList")]
    pub fn partition_list(
        &self,
        core_prom_size: u32,
        prox_prom_size: u32,
        chrom_sizes: JsValue,
    ) -> Result<JsPartitionList, JsValue> {
        let cs: Option<HashMap<String, u32>> = if chrom_sizes.is_null() || chrom_sizes.is_undefined() {
            None
        } else {
            Some(serde_wasm_bindgen::from_value(chrom_sizes).map_err(|e| {
                JsValue::from_str(&format!("Invalid chromSizes object: {}", e))
            })?)
        };

        let pl = gtars_genomicdist::genome_partition_list(
            &self.annotation.gene_model,
            core_prom_size,
            prox_prom_size,
            cs.as_ref(),
        );
        Ok(JsPartitionList::from_inner(pl))
    }
}
