use std::io::Write;
use std::fs::File;

use pyo3::prelude::*;

use genimtools::vocab::create_count_map;

#[pyfunction]
pub fn prune_universe(universe: &str, data: &str, min_count: Option<u32>, output: Option<&str>) -> PyResult<()> {
    
    let min_count = min_count.unwrap_or(1);
    let output = output.unwrap_or("output.bed");

    let cnt_map = create_count_map(data, universe).unwrap();

    // create output file
    let mut file = File::create(output).unwrap();
    for (region, cnt) in cnt_map {
        if cnt < min_count {
            // skip this region
            continue;
        }
        let line = format!("{}\t{}\t{}\n", region.chr, region.start, region.end);

        // write to file
        file.write_all(line.as_bytes())?;
    }

    Ok(())
}