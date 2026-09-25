# gtars-io

I/O focused utilities for gtars, providing efficient file format parsers and stream processing.

## Features

- BED file parsing and writing
- Fragment file I/O
- Compression support (gzip, bgzip)
- Stream processing for large files
- Memory-efficient iterators

## File Formats

### BED Files
```rust
use gtars_io::BedReader;

let reader = BedReader::from_path("regions.bed")?;
for region in reader {
    println!("{:?}", region?);
}
```

### Fragment Files
```rust
use gtars_io::FragmentReader;

let fragments = FragmentReader::from_path("fragments.tsv.gz")?;
for fragment in fragments {
    let frag = fragment?;
    println!("{}\t{}\t{}", frag.chr, frag.start, frag.end);
}
```

## Compression

Automatic detection and handling of compressed files:

- `.gz` - gzip compression
- `.bgz` - bgzip compression
- Uncompressed files

## Performance

- Buffered I/O for optimal throughput
- Zero-copy parsing where possible
- Parallel decompression support