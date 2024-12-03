use extendr_api::prelude::*;

use gtars::io::{read_tokens_from_gtok, write_tokens_to_gtok};

/// Write tokens to a gtok file
/// @export
/// @param filename A string representing the path to the gtok file.
#[extendr(r_name = "write_tokens_to_gtok")]
pub fn r_write_tokens_to_gtok(filename: String, tokens: Vec<i32>) {
    let tokens: Vec<u32> = tokens.into_iter().map(|t| t as u32).collect();
    let _ = write_tokens_to_gtok(&filename, &tokens);
}

/// Write tokens to a gtok file
/// @export
/// @param filename A string representing the path to the gtok file.
#[extendr(r_name = "read_tokens_from_gtok")]
pub fn r_read_tokens_from_gtok(filename: String) -> Vec<i32> {
    read_tokens_from_gtok(&filename)
        .unwrap()
        .into_iter()
        .map(|gtok| gtok as i32)
        .collect()
}

extendr_module! {
    mod io;
    fn r_read_tokens_from_gtok;
    fn r_write_tokens_to_gtok;
}
