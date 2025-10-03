use wasm_bindgen::prelude::*;

use gtars_tokenizers::tokenizer::Tokenizer;
use gtars_tokenizers::special_tokens::SpecialTokens;
use gtars_tokenizers::utils::create_tokenize_core_from_universe;
use gtars_tokenizers::config::TokenizerType;

use crate::models::BedEntries;
use crate::utils::prepare_universe_from_bed_entries;

#[wasm_bindgen(js_name = "Tokenizer")]
pub struct JsTokenizer {
    internal: Tokenizer
}

#[wasm_bindgen(js_class = "Tokenizer")]

impl JsTokenizer {

    #[wasm_bindgen(constructor)]
    pub fn new(universe: &JsValue, backend: &str) -> Result<JsTokenizer, JsValue> {
        let universe_regions: BedEntries = serde_wasm_bindgen::from_value(universe.to_owned())?;
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


        let mut universe = prepare_universe_from_bed_entries(universe_regions);
        let special_tokens = SpecialTokens::default();
        
        universe.add_special_tokens(&special_tokens);

        let core = create_tokenize_core_from_universe(&universe, tokenizer_type);
        let tokenizer = Tokenizer::new(core, universe, special_tokens);

        Ok(JsTokenizer { internal:tokenizer })
    }
}
