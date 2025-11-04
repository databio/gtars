Attempting to replicate IGD in Rust from C:
https://github.com/databio/IGD
https://academic.oup.com/bioinformatics/article/37/1/118/6050710

Current manual test:

Input: /home/drc/IGD_TEST/bedfiles/
Output: /home/drc/IGD_TEST/output/

Full command:

Create
```
cargo run igd create --output /home/drc/IGD_TEST/output/ --filelist /home/drc/IGD_TEST/bedfiles/
```

temp comparison
```
cargo run igd create --output /home/drc/IGD_TEST_2/igd_rust_output/ --filelist /home/drc/IGD_TEST_2/source_bedfiles/
```

Search
```
cargo run igd search --database /home/drc/IGD_TEST_2/igd_rust_output/igd_database.igd --query /home/drc/IGD_TEST_2/query_bed_file/igd_query_test.bed
```


## USING BLOOM FILTERS

You must enable the `bloom` feature for igd.

### Create Bloom Filters from a Directory of BED Files

```rust
use gtars_igd::igdbloom::{process_bed_directory, load_bloom_directory};
use gtars_tokenizers::tokenizer::Tokenizer;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Directory containing your .bed files
    let input_dir = "/path/to/bed/files";
    
    // Directory where .bloom files will be saved
    let output_dir = "/path/to/bloom/output";
    
    // Bloom filter parameters
    let num_of_items = 10000;           // Expected number of tokens per file
    let false_positive_rate = 0.01;      // 1% false positive rate
    
    // Create tokenizer from a sample BED file
    let sample_bed = "/path/to/sample.bed"; // Preferably this should be a UNIVERSE BED file
    let tokenizer = Tokenizer::from_auto(sample_bed)?;
    
    // Process all BED files and create bloom filters
    let processed_count = process_bed_directory(
        &tokenizer,
        input_dir,
        output_dir,
        num_of_items,
        false_positive_rate,
    )?;
    
    println!("Created {} bloom filters", processed_count);
    
    Ok(())
}
```

### Load and Query Bloom Filters

```rust
use gtars_igd::igdbloom::load_bloom_directory;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Load all bloom filters from directory into memory
    let bloom_dir = "/path/to/bloom/output";
    let filters = load_bloom_directory(bloom_dir)?;
    
    println!("Loaded {} bloom filters", filters.len());
    
    // Query a specific filter
    if let Some(filter) = filters.get("my_file.bed") { // filters have the same name as their original BED file
        let token = "chr1:100-200".to_string();
        if filter.check(&token) {
            println!("Token found in my_file.bed");
        }
    }
    
    // Query all filters to find which files contain a token
    let query_token = "chr1:1000-2000".to_string();
    for (filename, filter) in filters.iter() {
        if filter.check(&query_token) {
            println!("Token found in: {}", filename);
        }
    }
    
    Ok(())
}
```

### Add to Cargo.toml

```toml
[dependencies]
gtars-igd = { version = "*", features = ["bloom"] }
gtars-tokenizers = "*"
```



