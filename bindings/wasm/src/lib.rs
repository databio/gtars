mod digests;
mod tokenizers;
mod utils;

use std::str::EncodeUtf16;

use gtars::{
    common::models::Region,
    tokenizers::{
        config::TokenizerType, create_tokenize_core_from_universe, special_tokens::SpecialTokens,
        Tokenizer,
    },
};

use gtars::refget::store;

use wasm_bindgen::prelude::*;
use serde_wasm_bindgen;
use crate::{tokenizers::BedEntry, utils::prepare_universe_from_bed_entries};
use gtars::refget::collection::SequenceCollection;

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
}

#[wasm_bindgen]
pub fn greet(name: &str) {
    alert(&format!("Hello from gtars-js, {}", name));
}

#[wasm_bindgen(js_name="tokenizeBedfileIntoUniverse")]
pub fn tokenize_bed_file_into_universe(
    universe: &JsValue,
    to_tokenize: &JsValue,
    backend: &str,
) -> Result<JsValue, JsValue> {
    // get the backend
    let tokenizer_type = match backend {
        "ailist" => TokenizerType::AiList,
        "bits" => TokenizerType::Bits,
        _ => {
            return Err(JsValue::from_str(&format!(
                "Invalid backend specified: {}",
                backend
            )))
        }
    };

    // destructure the arbitrary js value into something we can work with
    let universe_entries: Vec<BedEntry> = serde_wasm_bindgen::from_value(universe.clone())?;
    let to_tokenize_entries: Vec<BedEntry> = serde_wasm_bindgen::from_value(to_tokenize.clone())?;

    let mut universe = prepare_universe_from_bed_entries(universe_entries);
    let special_tokens = SpecialTokens::default();
    universe.add_special_tokens(&special_tokens);

    let core = create_tokenize_core_from_universe(&universe, tokenizer_type);
    let tokenizer = Tokenizer::new(core, universe, special_tokens);

    let regions_to_tokenize = to_tokenize_entries
        .into_iter()
        .map(|e| Region {
            chr: e.chr,
            start: e.start,
            end: e.end,
            rest: None,
        })
        .collect::<Vec<Region>>();
    let tokens = tokenizer.tokenize(&regions_to_tokenize);
    if tokens.is_err() {
        return Err(JsValue::from_str("Error tokenizing regions."));
    }
    let tokens = tokens.unwrap();

    serde_wasm_bindgen::to_value(&tokens)
        .map_err(|e| JsValue::from_str(&format!("Error serializing tokens: {}", e)))
}


//////// GLOBALREFGETSTORE ADDITIONS

#[wasm_bindgen]
pub struct WasmRefgetStore {
    inner: store::GlobalRefgetStore,
}

#[wasm_bindgen]
impl WasmRefgetStore {
    #[wasm_bindgen(constructor)]
    pub fn new() -> WasmRefgetStore {
        WasmRefgetStore {
            inner: store::GlobalRefgetStore::new(store::StorageMode::Raw),
        }
    }
    
    #[wasm_bindgen]
    pub fn add_fasta(&mut self, fasta_bytes: &[u8]) -> Result<(), JsValue> {
        self.inner
            .import_fasta_bytes(fasta_bytes)
            .map_err(|e| JsValue::from_str(&format!("Failed to add FASTA: {}", e)))
    }

}

#[wasm_bindgen]
pub struct SequenceCollectionWasm {
    inner: SequenceCollection,
}

#[wasm_bindgen]
impl SequenceCollectionWasm {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            inner: SequenceCollection::from_records(vec![]),
        }
    }

    #[wasm_bindgen]
    pub fn from_fasta_bytes(fasta_bytes: &[u8]) -> Result<SequenceCollectionWasm, JsValue> {
        let collection = SequenceCollection::from_fasta_bytes(fasta_bytes)
            .map_err(|e| JsValue::from_str(&format!("Failed to parse FASTA: {}", e)))?;

        Ok(SequenceCollectionWasm { inner: collection })
    }

    #[wasm_bindgen(getter)]
    pub fn sequences(&self) -> JsValue {
        // Convert sequences to JS-compatible format
        let sequences: Vec<_> = self.inner.sequences.iter().map(|seq| {
            serde_json::json!({
            "metadata": {
                "name": seq.metadata.name,
                "length": seq.metadata.length,
                "sha512t24u": seq.metadata.sha512t24u,
                "md5": seq.metadata.md5,
                "alphabet": seq.metadata.alphabet
            }
        })
        }).collect();

        serde_wasm_bindgen::to_value(&sequences).unwrap_or(JsValue::NULL)
    }

}