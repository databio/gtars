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

    pub fn tokenize(&self, regions: &JsValue) -> Result<JsValue, JsValue> {
        let regions: BedEntries = serde_wasm_bindgen::from_value(regions.clone())?;
        let tokens = self.internal.tokenize(&regions);
        serde_wasm_bindgen::to_value(&tokens).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    pub fn encode(&self, tokens: &JsValue) -> Result<JsValue, JsValue> {
        let encoded = self.internal.encode(tokens);
        serde_wasm_bindgen::to_value(&encoded).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    pub fn decode(&self, ids: &JsValue) -> Result<JsValue, JsValue> {
        let decoded = self.internal.decode(ids);
        serde_wasm_bindgen::to_value(&decoded).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "convertIdsToTokens")]
    pub fn convert_ids_to_tokens(&self, ids: &JsValue) -> Result<JsValue, JsValue> {
        let tokens = self.internal.convert_ids_to_tokens(ids);
        serde_wasm_bindgen::to_value(&tokens).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "convertTokensToIds")]
    pub fn convert_tokens_to_ids(&self, tokens: &JsValue) -> Result<JsValue, JsValue> {
        let ids = self.internal.convert_tokens_to_ids(tokens);
        serde_wasm_bindgen::to_value(&ids).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    #[wasm_bindgen(getter, js_name = "unkToken")]
    pub fn unk_token(&self) -> String {
        self.internal.get_unk_token().to_string()
    }

    #[wasm_bindgen(getter, js_name = "padToken")]
    pub fn pad_token(&self) -> String {
        self.internal.get_pad_token().to_string()
    }

    #[wasm_bindgen(getter, js_name = "maskToken")]
    pub fn mask_token(&self) -> String {
        self.internal.get_mask_token().to_string()
    }

    #[wasm_bindgen(getter, js_name = "clsToken")]
    pub fn cls_token(&self) -> String {
        self.internal.get_cls_token().to_string()
    }

    #[wasm_bindgen(getter, js_name = "bosToken")]
    pub fn bos_token(&self) -> String {
        self.internal.get_bos_token().to_string()
    }

    #[wasm_bindgen(getter, js_name = "eosToken")]
    pub fn eos_token(&self) -> String {
        self.internal.get_eos_token().to_string()
    }

    #[wasm_bindgen(getter, js_name = "sepToken")]
    pub fn sep_token(&self) -> String {
        self.internal.get_sep_token().to_string()
    }

    #[wasm_bindgen(getter, js_name = "padTokenId")]
    pub fn pad_token_id(&self) -> u32 {
        self.internal.get_pad_token_id()
    }

    #[wasm_bindgen(getter, js_name = "maskTokenId")]
    pub fn mask_token_id(&self) -> u32 {
        self.internal.get_mask_token_id()
    }

    #[wasm_bindgen(getter, js_name = "clsTokenId")]
    pub fn cls_token_id(&self) -> u32 {
        self.internal.get_cls_token_id()
    }

    #[wasm_bindgen(getter, js_name = "bosTokenId")]
    pub fn bos_token_id(&self) -> u32 {
        self.internal.get_bos_token_id()
    }

    #[wasm_bindgen(getter, js_name = "eosTokenId")]
    pub fn eos_token_id(&self) -> u32 {
        self.internal.get_eos_token_id()
    }

    #[wasm_bindgen(getter, js_name = "sepTokenId")]
    pub fn sep_token_id(&self) -> u32 {
        self.internal.get_sep_token_id()
    }

    #[wasm_bindgen(getter, js_name = "unkTokenId")]
    pub fn unk_token_id(&self) -> u32 {
        self.internal.get_unk_token_id()
    }

    #[wasm_bindgen(getter, js_name = "vocabSize")]
    pub fn vocab_size(&self) -> usize {
        self.internal.get_vocab_size()
    }

    #[wasm_bindgen(js_name = "specialTokensMap")]
    pub fn special_tokens_map(&self) -> Result<JsValue, JsValue> {
        let map = self.internal.get_special_tokens();
        serde_wasm_bindgen::to_value(&map).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    #[wasm_bindgen(js_name = "getVocab")]
    pub fn get_vocab(&self) -> Result<JsValue, JsValue> {
        let vocab = self.internal.get_vocab();
        serde_wasm_bindgen::to_value(&vocab).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    #[wasm_bindgen(getter)]
    pub fn length(&self) -> usize {
        self.internal.get_universe().len()
    }
}
