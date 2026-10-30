# gtars Modules

gtars is organized as a workspace of independent Rust crates, each providing specific functionality for genomic interval analysis. All modules share common dependencies and infrastructure through `gtars-core`.

📖 **See [Module Vision](module-vision.md) for the conceptual organization and design principles.**

## Core Infrastructure

- **[gtars-core](core.md)** - Fundamental data structures (`Region`, `RegionSet`, `RegionSetList`, `Interval`, `Fragment`, `CoordinateMode`), interval set operations, and utilities
- **[gtars-io](io.md)** - I/O operations and file format parsers

## Genomic Analysis

- **[gtars-overlaprs](overlaprs.md)** - High-performance interval overlap operations
- **[gtars-uniwig](uniwig.md)** - Coverage computation from BED/BAM files
- **[gtars-igd](igd.md)** - Fast genomic interval search
- **[gtars-scoring](scoring.md)** - Fragment overlap scoring
- **[gtars-fragsplit](fragsplit.md)** - Fragment file splitting for pseudobulk

## Distribution Analysis & Enrichment

- **[gtars-genomicdist](genomicdist.md)** - Rust port of R GenomicDistributions: summary stats, partitions, signal matrices, consensus
- **[gtars-lola](lola.md)** - LOLA (Locus Overlap Analysis) enrichment testing built on IGD + genomicdist

## Machine Learning

- **[gtars-tokenizers](tokenizers.md)** - Genomic region tokenizers for ML

## Data Access

- **[gtars-refget](refget.md)** - GA4GH refget protocol implementation
- **[gtars-bbcache](bbcache.md)** - BED file caching for bedbase.org
- **gtars-vrs** - GA4GH VRS (Variation Representation Specification) allele digest computation

## Language Bindings

- **[gtars-python](python-overview.md)** - Python API bindings
- **[gtars-wasm](wasm.md)** - WebAssembly bindings
- **[gtars-r](r.md)** - R bindings
- **gtars-node** - Native Node.js bindings (refget)
- **[gtars-cli](cli.md)** - Command-line interface