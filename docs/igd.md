# gtars-igd

Implementation of IGD (Integrated Genome Database) - a high-performance genomic interval search tool.

## Features

- Fast interval queries (logarithmic time complexity)
- Memory-efficient indexing
- Binary format for persistence
- Support for multiple genomes

## Usage

Ensure `gtars` is compiled with igd:
`cargo build --release --all-features` or `cargo build --release --features igd`


### Building an Index
```shell
gtars igd create --output /home/igd_output/ --filelist /home/my_bedfiles/
```

### Querying with a single bed file
```shell
gtars igd search --database my_igd_database.igd --query my_query.bed
```

## Performance

- Index creation: O(n log n)
- Query time: O(log n + k) where k is number of results
- Memory usage: ~30% of input file size

## File Format

IGD uses a binary format optimized for:

- Fast loading
- Compact storage
- Memory-mapped access