# gtars-core

Core library providing the fundamental data structures and utilities that every other gtars crate builds on. If you're working directly with genomic regions in Rust, everything starts here.

!!! info "What's new"
    Recent additions to `gtars-core` that pre-existing doc snippets don't cover:

    - **`RegionSetList`** — a `GRangesList`-style collection of named `RegionSet`s with `concat`, `iter`, and set-level identifier computation. Consumed by `gtars-genomicdist` and `gtars-lola`.
    - **`CoordinateMode`** enum — selects BED (0-based half-open) vs. GRanges (1-based closed) conventions for midpoint calculations.
    - **`Fragment`** — fragment-file record type for single-cell ATAC-seq workflows.
    - **`Interval<I, T>`** — generic interval with payload, used internally by overlap indexes.
    - **Typed error enum** `RegionSetError` — replaces the previous panicking parse paths.
    - **`to_polars`** + `dataframe` feature flag — zero-copy conversion to a Polars DataFrame.
    - **`bigbed`** and **`http`** feature flags — optional bigBed writing and URL-backed `RegionSet::try_from`.
    - **Interval set operations**: `reduce`, `union`, `trim`, `promoters`, `gaps`, and friends are now inherent methods on `RegionSet`, with `setdiff`/`intersect`/`jaccard` in the `IntervalSetOps` trait. `SortedRegionSet` also lives here now (re-exported by `gtars-genomicdist`).

## Core data types

The `gtars_core::models` module re-exports the core types at the top level, so you typically import them directly:

```rust
use gtars_core::models::{
    Region, RegionSet, RegionSetList, SortedRegionSet, Interval, Fragment, CoordinateMode,
    IntervalSetOps,
};
```

### `Region`

A single genomic interval. `start` and `end` are 0-based half-open by BED convention; `rest` carries any trailing tab-separated fields from the original line verbatim.

```rust
use gtars_core::models::Region;

let region = Region {
    chr: "chr1".to_string(),
    start: 1000,
    end: 2000,
    rest: None,
};

assert_eq!(region.width(), 1000);
println!("{}", region.as_string()); // chr1\t1000\t2000
```

Methods:

| method | returns | notes |
|---|---|---|
| `width()` | `u32` | `end - start` |
| `as_string()` | `String` | tab-separated BED line |
| `digest()` | `String` | MD5 digest of `"chr,start,end"` |
| `mid_point()` | `u32` | `start + width() / 2` (BED/floor) |
| `mid_point_with_mode(mode)` | `u32` | BED or GRanges convention — see `CoordinateMode` below |
| `distance_to(&other)` | `i64` | gap in bp between the two regions; 0 if they overlap |

`Region` implements `Display` (as tab-separated text), `Clone`, `Debug`, `Eq`, `Hash`, and — under the `serde` feature — `Serialize`/`Deserialize`.

### `RegionSet`

An ordered collection of `Region`s. `RegionSet::try_from` accepts `&Path`, `&str`, `String`, or `PathBuf`, auto-detects gzip by magic bytes (not file extension), and with the `http` feature will fetch from URLs. Loading from a file always sorts by `(chr, start)`; `RegionSet::from(Vec<Region>)` keeps the order you give it.

```rust
use gtars_core::models::RegionSet;
use std::path::Path;

// From a local BED (or BED.gz) file
let rs = RegionSet::try_from(Path::new("peaks.bed"))?;

assert!(!rs.is_empty());
println!("{} regions, {} bp total", rs.len(), rs.nucleotides_length());

for region in &rs {
    println!("{}: {}-{}", region.chr, region.start, region.end);
}
# Ok::<(), gtars_core::errors::RegionSetError>(())
```

Key methods:

**Construction**

- `RegionSet::try_from(path)` — accepts `&Path`, `&str`, `String`, `PathBuf`. Sorts on load.
- `RegionSet::from(regions: Vec<Region>)` — in-memory construction.
- `RegionSet::from(bytes: &[u8])` — parse from an in-memory byte slice (no gzip handling).

**Iteration**

- `for region in &rs { ... }` — `IntoIterator` is implemented for `&RegionSet`.
- `iter_chroms()` → unique chromosomes in insertion order (post-sort).
- `iter_chr_regions(chr)` → all regions on a specific chromosome.

**Summaries**

- `len()`, `is_empty()`, `nucleotides_length()` — count and total bp.
- `region_widths()` → `Vec<u32>`.
- `mean_region_width()` → `f64`, rounded to 2 decimals.
- `get_max_end_per_chr()` → `HashMap<String, u32>`.
- `calc_mid_points()` → `HashMap<String, Vec<u32>>` (BED convention).
- `calc_mid_points_with_mode(CoordinateMode)` → same, with mode control.

**Identifiers**

- `identifier()` — MD5 digest over the first-three-column layout; the canonical BEDbase identifier.
- `file_digest()` — MD5 digest over the full serialized file content.

**I/O**

- `to_bed(path)` / `to_bed_gz(path)` — write plain or gzipped BED.
- `to_bigbed(path, chrom_sizes)` — under the `bigbed` feature, write a bigBed file.
- `to_polars()` — under the `dataframe` feature, return a `PolarsResult<DataFrame>`.

**Mutation**

- `sort()` — in-place sort by `(chr, start)`.

**Interval set operations**

GRanges-style operations are inherent methods on `RegionSet`; each returns a new `RegionSet` (with `rest: None`):

- `reduce()`, `concat(&other)`, `union(&other)` (plus consuming `concat_into` / `union_into`).
- `trim(&chrom_sizes)`, `gaps(&chrom_sizes)`, `promoters(upstream, downstream)`.
- `shift(offset)`, `flank(width, use_start, both)`, `resize(width, fix)`, `narrow(start, end, width)`.
- `pintersect(&other)`, `disjoin()`, `cluster(max_gap)` → `Vec<u32>`, `closest(&other)` → `Vec<(usize, usize, i64)>`.

Two-set operations live in the `IntervalSetOps` trait (also implemented by the overlap indexes in `gtars-overlaprs`): `setdiff`, `intersect`, `jaccard`, `coverage`, `overlap_coefficient`, and `closest`. Bring it into scope with `use gtars_core::models::IntervalSetOps;`.

`SortedRegionSet` is a newtype (`SortedRegionSet(pub RegionSet)`) that guarantees `(chr, start)` order. `SortedRegionSet::new(rs)` sorts in place without cloning.

### `RegionSetList`

A collection of `RegionSet`s — the gtars equivalent of Bioconductor's `GRangesList`. This is the type that downstream crates (genomicdist, lola) use to pass multiple region sets across FFI boundaries without paying N × clone costs.

```rust
use gtars_core::models::{RegionSet, RegionSetList};

let peaks1 = RegionSet::try_from("peaks_rep1.bed")?;
let peaks2 = RegionSet::try_from("peaks_rep2.bed")?;
let peaks3 = RegionSet::try_from("peaks_rep3.bed")?;

// Unnamed — names are optional
let rsl = RegionSetList::new(vec![peaks1, peaks2, peaks3]);

// Or with explicit names
let rsl = RegionSetList::with_names(
    rsl.region_sets,
    vec!["rep1".into(), "rep2".into(), "rep3".into()],
);

// Iterate
for rs in &rsl {
    println!("{} regions", rs.len());
}

// Flatten all regions into a single RegionSet (no merge/dedup)
let combined: RegionSet = rsl.concat();

// Stable identifier over the full set (order-independent)
let id = rsl.identifier();
# Ok::<(), gtars_core::errors::RegionSetError>(())
```

`RegionSetList::try_from` also accepts:

- A path to a **bedset manifest file** — a text file listing one BED path per line (`read_bedset_file` under the hood).
- A `Vec<&Path>`, `Vec<&str>`, `Vec<String>`, or `Vec<PathBuf>` — each is loaded as its own `RegionSet`.

`concat()` flattens without merging; if you need a reduced union, call `.reduce()` on the result.

Key methods: `new`, `with_names`, `add`, `get(i)`, `iter`, `len`, `is_empty`, `concat`, `identifier`.

### `CoordinateMode`

Switches between BED (0-based half-open, floor division) and GRanges (1-based closed, banker's rounding) conventions for midpoint calculation. BED is the default and matches Python/numpy conventions; GRanges exists for exact bit-compatibility with R GenomicDistributions output.

```rust
use gtars_core::models::{Region, CoordinateMode};

let r = Region { chr: "chr1".into(), start: 100, end: 206, rest: None };

assert_eq!(r.mid_point_with_mode(CoordinateMode::Bed), 153);      // 100 + 53
assert_eq!(r.mid_point_with_mode(CoordinateMode::GRanges), 152);  // banker's rounding
```

For widths `w` where `w % 4 == 2`, BED and GRanges midpoints differ by 1 bp — this affects approximately 2.6% of typical peak distance calculations. Use `GRanges` mode only when you need byte-for-byte parity with R output.

### `Fragment`

A fragment-file record (chromatin accessibility / scATAC-seq). Implements `FromStr` for parsing lines from 10x-style `fragments.tsv` files and `From<Fragment> for Region` for dropping the barcode/support metadata when you only care about coordinates.

```rust
use gtars_core::models::{Fragment, Region};
use std::str::FromStr;

let f = Fragment::from_str("chr1\t1000\t1200\tAAACCTGAGAAACCAT-1\t3")?;
assert_eq!(f.barcode, "AAACCTGAGAAACCAT-1");
assert_eq!(f.read_support, 3);

let as_region: Region = f.into();
# Ok::<(), anyhow::Error>(())
```

### `Interval<I, T>`

A generic `[start, end)` range with a payload `T`, parameterized over any unsigned integer type. This is primarily a building block for overlap indexes (consumed by `gtars-overlaprs` and `gtars-igd`); most user code should prefer `Region`/`RegionSet`.

```rust
use gtars_core::models::Interval;

let a: Interval<u32, usize> = Interval { start: 10, end: 50, val: 0 };
let b: Interval<u32, usize> = Interval { start: 40, end: 80, val: 1 };

assert!(a.overlap(b.start, b.end));
assert_eq!(a.intersect(&b), 10); // overlap width in bp
```

## Error handling

Parse and I/O errors from `RegionSet::try_from` come back as the typed `RegionSetError` enum (no panics on malformed input):

```rust
use gtars_core::errors::RegionSetError;
use gtars_core::models::RegionSet;

match RegionSet::try_from("not_a_real_file.bed") {
    Ok(rs) => println!("loaded {} regions", rs.len()),
    Err(RegionSetError::FileReadError(msg)) => eprintln!("read failed: {msg}"),
    Err(RegionSetError::RegionParseError(msg)) => eprintln!("parse failed: {msg}"),
    Err(RegionSetError::EmptyRegionSet(path)) => eprintln!("no regions in {path}"),
    Err(e) => eprintln!("other: {e}"),
}
```

Variants: `FileReadError`, `InvalidPathOrUrl`, `InvalidBedbaseIdentifier`, `BedbaseFetchError`, `RegionParseError`, `EmptyRegionSet`, `HttpFeatureDisabled`, `BigBedError`, and a transparent `Io` wrapper around `std::io::Error`.

## Feature flags

| flag | effect |
|---|---|
| *(default)* | Pure in-memory BED reading/writing, no optional dependencies. |
| `serde` | Derives `Serialize`/`Deserialize` on `Region`, `RegionSet`. |
| `http` | Enables `RegionSet::try_from(&Path)` to fetch from HTTP(S) URLs via `ureq`. Without this feature, non-file paths return `HttpFeatureDisabled`. |
| `dataframe` | Enables `RegionSet::to_polars()` (pulls in `polars`). Required transitively by `gtars-genomicdist`'s `bedclassifier` feature. |
| `bigbed` | Enables `RegionSet::to_bigbed()` via `bigtools`. |

Enable them in `Cargo.toml` like so:

```toml
[dependencies]
gtars-core = { version = "0.5", features = ["serde", "dataframe"] }
```

## Available modules

- **`models`** — all core data types (`Region`, `RegionSet`, `RegionSetList`, `SortedRegionSet`, `Interval`, `Fragment`, `CoordinateMode`) and the `IntervalSetOps` trait. Import them from `gtars_core::models`.
- **`errors`** — `RegionSetError` enum.
- **`utils`** — readers, file-type detection, chromosome-sizes parsing, and `Region` ↔ id hash-map helpers:
    - `get_dynamic_reader(&Path)` / `get_dynamic_reader_w_stdin(&str)` — transparent gzip/stdin handling (gzip detected by magic bytes).
    - `get_dynamic_reader_from_url(&Path)` — under the `http` feature.
    - `get_file_info(&Path) -> FileInfo` — detect type (BED, BAM, NARROWPEAK, UNKNOWN) and gzip.
    - `parse_bedlike_file(line)` → `Option<(chr, start, end)>` from a single line.
    - `get_chrom_sizes(path)` → `HashMap<String, u32>`.
    - `read_bedset_file(path)` → `Vec<String>` of BED paths from a bedset manifest.
    - `generate_region_to_id_map` / `generate_id_to_region_map` and string variants — stable id assignment for tokenizer vocabularies.
    - `remove_all_extensions(&Path)` → stem with *all* extensions stripped (handles `.bed.gz`).
    - `chrom_karyotype_key(chr)` → sort key for karyotype order (`chr1`…`chr22`, `chrX`, `chrY`, `chrM`, then others).
- **`consts`** — column-name constants (`CHR_COL_NAME`, `START_COL_NAME`, `END_COL_NAME`, `DELIMITER`) and file-extension constants (`BED_FILE_EXTENSION`, `BAM_FILE_EXTENSION`, `GZ_FILE_EXTENSION`, `IGD_FILE_EXTENSION`, `GTOK_EXT`).

## Where to go next

- **[Core models tour](regionSet.md)** — a cross-language (Python + Rust) walkthrough of `Region`, `RegionSet`, and friends.
- **[gtars-overlaprs](overlaprs.md)** — high-performance overlap queries that operate on `RegionSet`.
- **[gtars-genomicdist](genomicdist.md)** — the `GenomicIntervalSetStatistics` trait extends `RegionSet` with R GenomicDistributions–style summary stats, plus partitions, signal, and consensus.
- **[gtars-lola](lola.md)** — LOLA enrichment built on top of the IGD index and `RegionSetList`.
