# gtars-refget

Rust implementation of GA4GH refget sequence collection functions.

## Features

- GA4GH refget protocol support
- Sequence digest computation (SHA512t24u and others)
- Sequence collection management
- FASTA file reading and writing
- Efficient sequence storage and retrieval

## Modules

### digest
Compute and verify sequence digests:
```rust
use gtars_refget::digest;

// Compute SHA512t24u digest for a sequence
let seq = b"ACGTACGTACGT";
let digest = digest::sha512t24u(seq);
```

### collection
Manage sequence collections:
```rust
use gtars_refget::collection::SequenceCollection;

// Create and work with sequence collections
let mut collection = SequenceCollection::new();
// Add sequences, compute digests, etc.
```

### store
Efficient sequence storage indexed by digest:
```rust
use gtars_refget::store::SequenceStore;

// Store and retrieve sequences by their digest
let store = SequenceStore::new();
```

### fasta
Read and write FASTA files:
```rust
use gtars_refget::fasta;

// Parse FASTA files
// Write sequences to FASTA format
```

### alphabet
Define sequence alphabets (DNA, protein, ASCII):
```rust
use gtars_refget::alphabet;

// Work with different sequence alphabets
```

## Standards Compliance

Implements GA4GH refget specification for sequence collections and digests.