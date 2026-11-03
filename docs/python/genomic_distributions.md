# `gtars.genomic_distributions`

Python bindings for [`gtars-genomicdist`](../genomicdist.md) — the Rust port of R GenomicDistributions plus extras. This submodule contains the **free functions** that operate on regions-plus-something-else (a reference genome, a partition list, a signal matrix). Methods that operate only on a `RegionSet` (summary stats, interval algebra, peak clustering) live directly on `RegionSet` in [`gtars.models`](models.md).

!!! info "Where to look for what"
    - Per-region summaries, nearest-neighbor distances, region distribution bins, set algebra → `RegionSet` methods in [`gtars.models`](models.md).
    - GC content / dinucleotide frequencies → **this module** (need a reference genome).
    - Partition classification → **this module** (need a `PartitionList`).
    - Signal matrix overlap → **this module** (need a `SignalMatrix`).
    - Consensus regions across replicates → **this module**.
    - LOLA enrichment → [`gtars.lola`](lola.md).

The Rust reference page, [gtars-genomicdist](../genomicdist.md), has the algorithmic detail, caveats about bin-width semantics, divergences from R GenomicDistributions, and so on. This page is a Python-focused reference.

## Sequence statistics

These both take a `RegionSet` and a reference genome (`GenomeAssembly` for in-memory or `BinaryGenomeAssembly` for mmap-backed `.fab` files — see [`gtars.models`](models.md#genomeassembly)).

### `calc_gc_content`

```python
from gtars.models import RegionSet, BinaryGenomeAssembly
from gtars.genomic_distributions import calc_gc_content

peaks = RegionSet("peaks.bed")
genome = BinaryGenomeAssembly("hg38.fab")

gc = calc_gc_content(peaks, genome, ignore_unk_chroms=True)
# gc: list[float] — one value per region, 0.0–1.0
```

With `ignore_unk_chroms=False` (default), the call raises if any region lives on a chromosome missing from the assembly. Set it to `True` to skip those regions silently.

### `calc_dinucl_freq`

```python
from gtars.genomic_distributions import calc_dinucl_freq

result = calc_dinucl_freq(
    peaks, genome,
    raw_counts=False,       # False (default) → percentages 0–100; True → integer counts
    ignore_unk_chroms=True,
)

# result is a dict:
result["region_labels"]   # list[str] — "chr_start_end" per region
result["dinucleotides"]   # list[str] — 16 dinucleotide names in canonical order
result["frequencies"]     # list[list[float]] — one row per region, 16 values per row
```

The 16 columns in `frequencies` correspond one-for-one with the 16 names in `dinucleotides`. For pooled global counts across all regions, run with `raw_counts=True` and sum each column of the `frequencies` matrix.

This matches R GenomicDistributions `calcDinuclFreq` column-for-column, including the canonical order.

## Partition classification

Builds on the `PartitionList` type from [`gtars.models`](models.md#partitionlist).

### `calc_partitions`

```python
from gtars.models import RegionSet, PartitionList
from gtars.genomic_distributions import calc_partitions

peaks = RegionSet("peaks.bed")
pl = PartitionList.from_gtf("gencode.v47.gtf.gz", core_prom=100, prox_prom=2000)

result = calc_partitions(peaks, pl, bp_proportion=False)

# result is a dict with:
result["partition"]   # list[str] — partition names + "intergenic"
result["count"]       # list[int] — region (or bp) count per partition
result["total"]       # int — sum of all counts
```

- `bp_proportion=False` (priority mode): each query region is assigned to the first partition it overlaps. Mutually exclusive. The `intergenic` bucket collects everything not matched.
- `bp_proportion=True`: for each partition, sum the overlapping base pairs from the query set. Not mutually exclusive — a region that straddles two partitions contributes bp to each (matching R behavior).

### `calc_expected_partitions`

Observed vs expected enrichment with chi-square test, given chromosome sizes:

```python
from gtars.genomic_distributions import calc_expected_partitions

chrom_sizes = {"chr1": 248_956_422, "chr2": 242_193_529, ...}

result = calc_expected_partitions(
    peaks, pl, chrom_sizes, bp_proportion=False,
)

# result is a dict with parallel lists:
result["partition"]   # list[str]
result["observed"]    # list[float]
result["expected"]    # list[float]
result["log10OE"]     # list[float] — log10(observed / expected)
result["pvalue"]      # list[float] — chi-square goodness-of-fit p-value
```

!!! warning "p-values will differ slightly from R"
    `calc_expected_partitions` uses a goodness-of-fit `(O-E)²/E` formula. R's `chisq.test()` computes a 2×2 test of independence with optional Yates correction, so p-values will not match R GenomicDistributions output byte-for-byte. The `log10OE` column is directly comparable.

## Signal matrix overlap

### `calc_summary_signal`

Overlap a query region set against a pre-loaded `SignalMatrix` and compute Tukey boxplot statistics per condition:

```python
from gtars.models import RegionSet, SignalMatrix
from gtars.genomic_distributions import calc_summary_signal

peaks = RegionSet("peaks.bed")
sm = SignalMatrix.from_tsv("signal_matrix.tsv")
# or: sm = SignalMatrix.load_bin("signal_matrix.sigm")

result = calc_summary_signal(peaks, sm)

result["condition_names"]   # list[str] — columns of the signal matrix
result["region_labels"]     # list[str] — "chr_start_end" per matched query
result["signal_matrix"]     # list[list[float]] — row per matched query, col per condition (MAX)
result["matrix_stats"]      # list of dicts — Tukey boxplot stats per condition
```

Each entry in `matrix_stats` has: `condition`, `lower_whisker`, `lower_hinge`, `median`, `upper_hinge`, `upper_whisker`.

For each query region, the function finds all overlapping rows in the signal matrix and takes the **max** per condition — answering "what's the peak signal at this region in each cell type?".

## Consensus regions

### `consensus`

Given N region sets (e.g. replicates), produce the reduced union with per-region replicate support count:

```python
from gtars.models import RegionSet
from gtars.genomic_distributions import consensus

reps = [
    RegionSet("rep1.bed"),
    RegionSet("rep2.bed"),
    RegionSet("rep3.bed"),
]

cons = consensus(reps)
# cons: list of dicts with keys: chr, start, end, count

# Keep regions present in ≥ 2/3 replicates
robust = [c for c in cons if c["count"] >= 2]
print(f"{len(robust)} robust regions")
```

The `count` field is how many input sets overlap each union region — not the number of replicate *regions*, but the number of input *sets*.

## Utilities

### `median_abs_distance`

Median of absolute values, filtering sentinels. Pairs well with `TssIndex.feature_distances`:

```python
from gtars.models import TssIndex, RegionSet
from gtars.genomic_distributions import median_abs_distance

tss = TssIndex("gencode.v47.tss.bed")
peaks = RegionSet("peaks.bed")

signed = tss.feature_distances(peaks)        # list[float | None]
# Drop None values (regions on chromosomes not in the TSS index)
clean = [d for d in signed if d is not None]

median = median_abs_distance(clean)          # Optional[float]
```

Returns `None` if the input is empty or contains only missing/sentinel values.

## End-to-end example

```python
from gtars.models import (
    RegionSet, BinaryGenomeAssembly, PartitionList, SignalMatrix, TssIndex,
)
from gtars.genomic_distributions import (
    calc_gc_content, calc_partitions, calc_expected_partitions,
    calc_summary_signal, consensus,
)

# 1. Load inputs once
peaks = RegionSet("peaks.bed")
genome = BinaryGenomeAssembly("hg38.fab")
pl = PartitionList.from_gtf("gencode.v47.gtf.gz", core_prom=100, prox_prom=2000)
chrom_sizes = {"chr1": 248_956_422, "chr2": 242_193_529}

# 2. Peak-level statistics — live on RegionSet itself
print(f"{len(peaks)} regions, mean width {peaks.mean_region_width():.1f}")
print(f"Per-chromosome: {peaks.chromosome_statistics()}")

# 3. GC content per region
gc = calc_gc_content(peaks, genome, ignore_unk_chroms=True)
print(f"Mean GC: {sum(gc) / len(gc):.3f}")

# 4. Partition enrichment (observed vs expected)
exp = calc_expected_partitions(peaks, pl, chrom_sizes)
for p, o, e, lo, pv in zip(
    exp["partition"], exp["observed"], exp["expected"], exp["log10OE"], exp["pvalue"]
):
    print(f"  {p:<15} obs={o:>6.0f}  exp={e:>10.1f}  log10(O/E)={lo:+.2f}  p={pv:.2e}")

# 5. Nearest-TSS distances
tss = TssIndex("gencode.v47.tss.bed")
tss_dists = tss.calc_tss_distances(peaks)
print(f"Median distance to nearest TSS: {sorted(tss_dists)[len(tss_dists)//2]:,}")

# 6. Consensus across replicates
reps = [RegionSet(f"rep{i}.bed") for i in (1, 2, 3)]
cons = consensus(reps)
robust = [c for c in cons if c["count"] >= 2]
print(f"{len(robust)} of {len(cons)} union regions present in ≥2 replicates")
```

## See also

- **[`gtars.models`](models.md)** — `RegionSet`, `GenomeAssembly`, `PartitionList`, `SignalMatrix`, and the peak-stats methods attached to `RegionSet`.
- **[`gtars.lola`](lola.md)** — LOLA enrichment testing.
- **[gtars-genomicdist](../genomicdist.md)** — full Rust API reference, algorithmic notes, and caveats that apply identically in Python.
