# Gtars command line interface (CLI)

We've packaged some gtars functionality into a command line interface (CLI) for easy use in shell scripts and terminal environments. The CLI provides access to various genomic tools and utilities built on top of the core gtars Rust library.

## Installation

To facilitate fast compilation and broad compatibility, the CLI is feature-gated. You can enable only the tools you need to minimize dependencies, binary size, and build time.

There are two main ways to install the gtars CLI:

### Via Cargo (Rust's package manager)
```bash
cargo install gtars-cli --features "uniwig overlaprs igd bbcache scoring fragsplit blast"
```

### From source
```bash
git clone https://github.com/databio/gtars.git
cd gtars
cargo install --path gtars-cli --features "uniwig overlaprs igd bbcache scoring fragsplit blast"
```

You can also enable all features with:
```bash
cargo install gtars-cli --all-features
```

## Available Commands

The CLI provides the following subcommands (availability depends on features enabled during compilation):

### igd
Build and query IGD (Integrated Genome Database) indexes:
```bash
gtars igd create --input regions.bed --output index.igd
gtars igd query --database index.igd --query chr1:1000-2000
```

### overlaprs
Compute overlaps between genomic intervals:
```bash
gtars overlaprs --input1 regions1.bed --input2 regions2.bed
```

### uniwig
Generate coverage tracks from BED/BAM files:
```bash
gtars uniwig --input reads.bam --output coverage.bw
```

### bbcache
Cache and manage BED files from bedbase.org:
```bash
gtars bbcache get --id GSM123456
gtars bbcache list
```

### scoring
Score fragment overlaps against a reference:
```bash
gtars scoring --fragments frags.tsv.gz --universe peaks.bed --output scores.txt
```

### fragsplit
Split fragment files by cell barcodes or clusters:
```bash
gtars fragsplit --fragments frags.tsv.gz --barcodes clusters.csv --output-dir splits/
```

### blast
BLAST over a RefgetStore: build databases straight from a store (much faster than `makeblastdb`), search a store, and report where hits come from. The verbs are those of the `rgmakeblastdb` (`build`, `verify`, `view`, `nstats`) and `rgblast` (`search`, `report`, `origin`, `reps`) binaries, with the same flags:
```bash
gtars blast build --store ./store --collection ucsc:hg38 --out hg38
gtars blast verify --db hg38 --store ./store
gtars blast view --db hg38 --per-collection
gtars blast nstats --store ./store --all
gtars blast search --store ./store --all -- -query q.fa -outfmt 6 -out hits.tsv
gtars blast search --engine blastn --db hg38 -- -query q.fa -outfmt 6
gtars blast report hits.tsv --query q.fa
gtars blast origin --store ./store --all --query q.fa
gtars blast reps --store ./store --all
```
See [gtars-rgblast](blast.md) for every verb and option. `search` and `origin` run the in-process engine `rgblastn` on the store, which needs no BLAST+, or with `--engine blastn --db BASE` run BLAST+ `blastn` (from `--blastn PATH` or the `PATH`) on a database you built. Every command that takes `--store` also takes the fetch flags `--cache-dir`, `--resume` and `--fetch-jobs` for remote stores. Exit codes: 0 success, 1 runtime error, 2 usage error or missing engine, 130 cancelled.

## Global Options

```bash
gtars --help              # Show help
gtars --version           # Show version
gtars <command> --help    # Show command-specific help
```

## Building with Specific Features

To build the CLI with specific tools:
```bash
cargo build --release --features "uniwig,overlaprs,igd"
```