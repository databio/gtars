use gtars::common::utils::{generate_id_to_region_string_map, generate_region_string_to_id_map};
use gtars::tokenizers::universe::Universe;

use crate::tokenizers::BedEntry;

pub fn set_panic_hook() {
    // When the `console_error_panic_hook` feature is enabled, we can call the
    // `set_panic_hook` function at least once during initialization, and then
    // we will get better error messages if our code ever panics.
    //
    // For more details see
    // https://github.com/rustwasm/console_error_panic_hook#readme
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}

// very basic wrapper that is specific to wasm
// since we can do I/O in a browser and must pass regions
// to the function in memory
pub fn prepare_universe_from_bed_entries(region_entries: Vec<BedEntry>) -> Universe {
    let regions = region_entries
        .into_iter()
        .map(|e| format!("{}:{}-{}", e.chr, e.start, e.end))
        .collect::<Vec<String>>();

    let region_to_id = generate_region_string_to_id_map(&regions);
    let id_to_region = generate_id_to_region_string_map(&regions);

    Universe {
        regions,
        region_to_id,
        id_to_region,
        names: None,
        scores: None,
        special_tokens: None,
    }
}
