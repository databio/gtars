use wasm_bindgen::prelude::*;

use refget::digest::{md5, sha512t24u};

#[wasm_bindgen(js_name="digestFasta")]
pub fn digest_fasta_from_string(fasta: &JsValue, algorithm: &JsValue) -> Result<JsValue, JsValue> {

    let fasta_string: String = serde_wasm_bindgen::from_value(fasta.clone())?;
    let alg: String = serde_wasm_bindgen::from_value(algorithm.clone())?;

    let digest = match alg.as_str() {
        "md5" => {
            md5(fasta_string)
        },
        "sha512t24u" => {
            sha512t24u(fasta_string)
        },
        _ => {
            return Err(JsValue::from_str(&format!(
                "Invalid algorithm specified: {}",
                alg
            )))
        }
    };
   
    Ok(JsValue::from_str(&digest))
}