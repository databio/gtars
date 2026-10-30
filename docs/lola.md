# gtars-lola

Rust port of the R [LOLA](http://lola.databio.org/) (Locus Overlap Analysis) package: region-set enrichment testing against a database of reference region sets, using Fisher's exact test on 2×2 contingency tables.

Given a **user set** of regions (your peaks), a **universe** (the background you consider sampleable — typically a union of publicly available peak sets), and a **region database** (LOLA core or custom), `gtars-lola` computes for every database set the probability that your user set's overlap with that database set is larger (or smaller) than expected under the null.

The crate is built on top of:

- **[gtars-igd](igd.md)** — the interval overlap index. A LOLA database is exposed as a single `Igd` index over all its files, so one call can score a user set against thousands of reference sets simultaneously.
- **[gtars-core](core.md)** — all region types come from here; `RegionSet::concat` and `RegionSet::disjoin` are used to build restricted universes, and `RegionSetList` is used as the FFI-friendly return type for extracting region sets from a database.

## Installation

```toml
[dependencies]
gtars-lola = "0.2"
```

No feature flags — the crate has a single default configuration.

## Module layout

`gtars-lola` does not re-export symbols at the crate root; import from the submodules directly:

- **`gtars_lola::database`** — `RegionDB`, `CollectionAnno`, `RegionSetAnno`
- **`gtars_lola::enrichment`** — `run_lola`, `ContingencyTable` impls
- **`gtars_lola::universe`** — `check_universe_appropriateness`, `redefine_user_sets`, `build_restricted_universe`, `UniverseReport`, `UserSetReport`
- **`gtars_lola::output`** — `annotate_results`, `apply_fdr_correction`, `results_to_columns`, `write_results_tsv`, `LolaColumnar`
- **`gtars_lola::models`** — `Direction`, `LolaConfig`, `ContingencyTable`, `LolaResult`
- **`gtars_lola::errors`** — `LolaError`

## End-to-end workflow

```rust
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use gtars_core::models::RegionSet;
use gtars_lola::database::RegionDB;
use gtars_lola::enrichment::run_lola;
use gtars_lola::models::LolaConfig;
use gtars_lola::output::{annotate_results, apply_fdr_correction, write_results_tsv};

// 1. Load the region database from a LOLA-format folder
let db = RegionDB::from_lola_folder(
    Path::new("LOLACore/hg38"),
    None,   // collections filter (None = all)
    None,   // per-collection file limit (None = all)
)?;

// 2. Load the user set(s) and universe
let user_sets = vec![
    RegionSet::try_from("peaks_condition_a.bed")?,
    RegionSet::try_from("peaks_condition_b.bed")?,
];
let universe = RegionSet::try_from("universe.bed")?;

// 3. Run Fisher's exact test for every (user_set, db_set) pair
let mut results = run_lola(
    &db.igd,
    &user_sets,
    &universe,
    &LolaConfig::default(),
)?;

// 4. Attach database metadata (collection, cell type, tissue, …)
annotate_results(&mut results, &db);

// 5. Benjamini-Hochberg FDR correction, per user set
apply_fdr_correction(&mut results);

// 6. Write results to TSV (format matches R LOLA's writeCombinedEnrichment)
let mut out = BufWriter::new(File::create("lola_results.tsv")?);
write_results_tsv(&mut out, &results)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

The steps correspond 1:1 to the R LOLA workflow: `loadRegionDB` → `runLOLA` → `getLolaResults` → `writeCombinedEnrichment`.

## Loading a region database

`RegionDB::from_lola_folder` expects a directory laid out in the standard R LOLA format:

```text
db_path/
├── collection1/
│   ├── collection.txt      # collector, date, source, description (TSV)
│   ├── index.txt           # per-file annotations (TSV or CSV)
│   └── regions/
│       ├── file1.bed
│       └── file2.bed
├── collection2/
│   ├── collection.txt
│   ├── index.txt
│   └── regions/
│       └── …
```

`collection.txt` and `index.txt` are both auto-detected TSV or CSV (R's `fread` convention). Recognized columns in `index.txt`: `filename`, `cellType`, `description`, `tissue`, `dataSource`, `antibody`, `treatment`. Unknown columns are ignored; missing optional columns are preserved as `None`.

### Filters and limits

```rust
// Load only specific collections
let db = RegionDB::from_lola_folder(
    path,
    Some(&["encode_segmentation", "roadmap_epigenomics"]),
    None,
)?;

// Limit files per collection (useful for smoke tests)
let db = RegionDB::from_lola_folder(path, None, Some(10))?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

### In-memory construction (for API integrations)

If you're serving region sets from an API rather than a folder, build the `RegionDB` directly from a pre-constructed `Igd`:

```rust
use gtars_lola::database::{RegionDB, RegionSetAnno};

let db = RegionDB::from_igd_with_regions(igd, region_sets, region_annos);
```

This is the path BEDbase uses to expose LOLA over the web — the IGD is built server-side from the bedset collection and shipped to the client.

### Convenience accessors

- **`db.num_region_sets()`** — total region sets in the database.
- **`db.list_region_sets(collections)`** → `Vec<String>` — filenames, optionally filtered by collection.
- **`db.get_region_set(filenames, collections)`** → `Vec<&RegionSet>` — look up by filename.
- **`db.get_region_set_list(&[usize])`** → `RegionSetList` — extract indexed region sets as a named `RegionSetList` (for FFI, plotting, or downstream analysis). Out-of-bounds indices are silently skipped.
- **`RegionDB::merge(a, b)`** — combine two databases into one; rebuilds the IGD.

## Checking and shaping the universe

LOLA results are only meaningful if your universe **contains** your user set: every user region should overlap at least one universe region, with no many-to-many mappings. The `universe` module provides three diagnostics/fixers:

### `check_universe_appropriateness`

```rust
use gtars_lola::universe::check_universe_appropriateness;
use gtars_igd::igd::Igd;

// Pre-build the universe IGD once — you can reuse it across calls
let universe_igd = Igd::from_single_region_set(&universe);

let report = check_universe_appropriateness(&user_sets, &universe_igd);

for r in &report.user_set_reports {
    println!(
        "user set {}: {:.1}% coverage ({} of {} regions), {} many-to-many",
        r.user_set_index,
        r.coverage * 100.0,
        r.regions_in_universe,
        r.total_regions,
        r.many_to_many_count,
    );
    for w in &r.warnings {
        eprintln!("  ⚠ {w}");
    }
}
```

Warnings fire when coverage is under 50% (severe) or under 90% (moderate), or when any user region overlaps more than one universe region (many-to-many).

### `redefine_user_sets`

Rewrite each user set in terms of universe regions: for every user region, find the universe regions it overlaps and emit those as the new user set. This eliminates many-to-many mapping artifacts and is the Rust equivalent of R LOLA's `redefineUserSets()`.

```rust
use gtars_lola::universe::redefine_user_sets;

let universe_igd = Igd::from_single_region_set(&universe);
let redefined: Vec<RegionSet> =
    redefine_user_sets(&user_sets, &universe, &universe_igd);

// Now run enrichment on `redefined` instead of the raw user_sets
let results = run_lola(&db.igd, &redefined, &universe, &LolaConfig::default())?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

### `build_restricted_universe`

For differential enrichment analysis — build a universe that is exactly the union of all user sets, disjoined into non-overlapping pieces. This is R LOLA's `disjoin(unlist(userSets))`.

```rust
use gtars_lola::universe::build_restricted_universe;

let restricted = build_restricted_universe(&user_sets);
// Use `restricted` as the universe in a differential run
```

## Running enrichment

`run_lola(igd, user_sets, universe, config)` is the core engine. For every `(user_set, db_set)` pair it builds a 2×2 contingency table and runs a one-sided Fisher exact test:

```text
                   In DB set    Not in DB set
In user set           a              c
Not in user set       b              d
```

- **a** = overlap count between user set and DB set (the *support* — shows up as `support` in results).
- **b** = universe-DB overlap − user-DB overlap.
- **c** = user set size − user-DB overlap.
- **d** = universe size − a − b − c.

### `LolaConfig`

```rust
pub struct LolaConfig {
    pub min_overlap: i32,         // default 1
    pub direction: Direction,     // default Enrichment
}

pub enum Direction {
    Enrichment, // P(X ≥ a), alternative = "greater"
    Depletion,  // P(X ≤ a), alternative = "less"
}
```

```rust
use gtars_lola::models::{LolaConfig, Direction};

let config = LolaConfig {
    min_overlap: 1,
    direction: Direction::Enrichment,
};
```

### p-values and odds ratios

The p-value is computed from the hypergeometric distribution using the survival function (`sf`) rather than `1 - cdf` to avoid catastrophic cancellation when the tail probability is very small. It is reported as `-log10(p)` (capped at ~322 via a `1e-322` floor, matching R LOLA).

The odds ratio is the **conditional maximum likelihood estimate** — the noncentrality parameter of Fisher's noncentral hypergeometric distribution — solved by Brent's method. This matches `fisher.test()$estimate` in R exactly, not the simple `(a·d)/(b·c)` point estimate.

### Negative contingency values

If your user set contains regions *outside* the universe, the contingency table can produce negative `b`, `c`, or `d`. `run_lola` matches R LOLA's behavior: it prints a warning to stderr, stores the signed values on the `LolaResult` row, and returns `p_value_log = 0.0`, `odds_ratio = NaN` for that row. To avoid this entirely, pre-process your inputs with `redefine_user_sets`.

## Results

Each `(user_set, db_set)` pair produces a `LolaResult`:

```rust
pub struct LolaResult {
    pub user_set: usize,              // 0-based user set index
    pub db_set: usize,                // 0-based db set index
    pub p_value_log: f64,             // -log10(p), ≥ 0
    pub odds_ratio: f64,              // CMLE, matches R fisher.test()$estimate
    pub support: u64,                 // contingency table cell 'a'
    pub rnk_pv: usize,                // rank by p-value (1-based, ascending)
    pub rnk_or: usize,                // rank by odds ratio (1-based, descending)
    pub rnk_sup: usize,                // rank by support (1-based, descending)
    pub max_rnk: usize,                // max(rnk_pv, rnk_or, rnk_sup)
    pub mean_rnk: f64,                 // mean of the three ranks
    pub b: i64, pub c: i64, pub d: i64, // signed contingency values
    pub q_value: Option<f64>,          // BH-adjusted p (None until apply_fdr_correction)
    pub filename: String,
    pub collection: Option<String>,
    pub description: Option<String>,
    pub cell_type: Option<String>,
    pub tissue: Option<String>,
    pub antibody: Option<String>,
    pub treatment: Option<String>,
    pub data_source: Option<String>,
    pub db_set_size: u64,
}
```

`run_lola` returns results sorted by descending `p_value_log`, then ascending `mean_rnk` — identical to R LOLA's output order. Ranks are assigned per user set so the rank columns within a user set are independent.

## Annotation and FDR

Two post-processing steps, both in `gtars_lola::output`:

- **`annotate_results(&mut results, &db)`** — fills in `collection`, `description`, `cell_type`, `tissue`, `antibody`, `treatment`, `data_source`, and `db_set_size` from the database. `description` is truncated to 80 chars to match R LOLA.
- **`apply_fdr_correction(&mut results)`** — Benjamini-Hochberg q-value computation, applied **independently per user set**. Writes to the `q_value: Option<f64>` field.

## Output

### TSV (R LOLA–compatible)

```rust
use gtars_lola::output::write_results_tsv;
use std::fs::File;
use std::io::BufWriter;

let mut out = BufWriter::new(File::create("lola_results.tsv")?);
write_results_tsv(&mut out, &results)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

The header exactly matches R LOLA's `writeCombinedEnrichment`:

```text
userSet  dbSet  collection  pValueLog  oddsRatio  support
rnkPV  rnkOR  rnkSup  maxRnk  meanRnk  b  c  d
description  cellType  tissue  antibody  treatment  dataSource
filename  qValue  size
```

`userSet` and `dbSet` are 1-based in the TSV (R convention), even though they are 0-based in memory.

### Columnar (for FFI bindings)

`results_to_columns(&results) -> LolaColumnar` pivots the row-oriented `Vec<LolaResult>` into parallel column vectors. This is what the Python, WASM, and R bindings consume to build native data frames without reimplementing the row→column transpose:

```rust
use gtars_lola::output::{results_to_columns, LolaColumnar};

let cols: LolaColumnar = results_to_columns(&results);
println!("{} rows", cols.p_value_log.len());
// cols.p_value_log[i], cols.odds_ratio[i], cols.cell_type[i], … all align on row i
```

Every field of `LolaResult` has a corresponding `Vec<T>` on `LolaColumnar`; the ordering is preserved from the input slice.

## Direct contingency-table use

If you want to compute a single Fisher test outside of the full `run_lola` pipeline:

```rust
use gtars_lola::models::{ContingencyTable, Direction};

let ct = ContingencyTable { a: 45, b: 155, c: 155, d: 9645 };

let p     = ct.fisher_pvalue(Direction::Enrichment);
let log_p = ct.p_value_log(Direction::Enrichment);
let or    = ct.odds_ratio();

println!("p = {p:.2e}  (-log10 = {log_p:.2})  odds = {or:.2}");
```

`fisher_pvalue`, `p_value_log`, and `odds_ratio` are `impl ContingencyTable` methods — the same functions `run_lola` calls internally.

## Errors

```rust
pub enum LolaError {
    EmptyUniverse,                // universe has 0 regions
    EmptyDatabase,                // IGD has 0 files
    Other(anyhow::Error),         // wrapped I/O etc.
}
```

`run_lola` returns early with `EmptyUniverse` / `EmptyDatabase` so you can distinguish configuration problems from per-row negative-contingency warnings (which print to stderr and continue).

## Where to go next

- **[gtars-igd](igd.md)** — the overlap index that backs `RegionDB.igd`. Worth reading if you're tuning universe construction or wrapping your own database.
- **[gtars-core](core.md)** — `RegionSet` and `RegionSetList`, the input/output types throughout this page (plus the `concat`/`disjoin` operations used by `build_restricted_universe`).
- **R LOLA** ([docs](http://lola.databio.org/)) — the reference implementation this port targets.
