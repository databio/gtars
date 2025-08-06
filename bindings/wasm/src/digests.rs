use wasm_bindgen::prelude::*;

use gtars::refget::fasta::digest_fasta;

#[wasm_bindgen(js_name="digestFasta")]
pub fn digest_fasta_from_string(fasta: &JsValue) -> Result<JsValue, JsValue> {

    let fasta_string: String = serde_wasm_bindgen::from_value(fasta.clone())?;

    match digest_fasta(&fasta_string) {
        Ok(seq_col) => {
            Ok(JsValue::from_str(&seq_col.digest))
        },
        Err(err) => {
            Err(JsValue::from_str(&format!("There was an error digesting the fasta file: {}", err)))
        }
    }
}