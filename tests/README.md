# Testing gtars

## Basic Commands

```bash
# Run all tests
cargo test

# Run tests for a specific package/crate
cargo test -p gtars-refget

# Run tests matching a name pattern
cargo test refget
cargo test test_import_fasta

# Combine filters
cargo test -p gtars-refget test_import
```

## Useful Options

```bash
# Show output from tests
cargo test -- --nocapture

# Run ignored tests
cargo test -- --ignored

# List available tests without running
cargo test -- --list
```
