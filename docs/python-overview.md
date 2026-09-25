# gtars-python

Python bindings for gtars — native Python API for genomic interval analysis, built on the Rust core.

## Installation

```bash
pip install gtars
```

For development installs from source, `uv pip install -e ./gtars-python` from the repo root.

## Submodules

The `gtars` Python package exposes several submodules, each wrapping a gtars Rust crate:

| submodule | Rust crate | purpose | dedicated doc |
|---|---|---|---|
| `gtars.models` | `gtars-core` + `gtars-genomicdist::models` | `Region`, `RegionSet`, `RegionSetList`, reference genomes, gene models, partition lists, signal matrices — plus all the per-`RegionSet` statistics and interval algebra methods | [models](python/models.md) |
| `gtars.genomic_distributions` | `gtars-genomicdist` | GC content, dinucleotide freq, partition classification, signal summaries, consensus regions (free functions that need a reference genome or a partition list) | [genomic_distributions](python/genomic_distributions.md) |
| `gtars.lola` | `gtars-lola` | LOLA enrichment testing: `RegionDB`, `run_lola`, universe helpers | [lola](python/lola.md) |
| `gtars.tokenizers` | `gtars-tokenizers` | Genomic region tokenizers for ML | *(see in-package docs)* |
| `gtars.refget` | `gtars-refget` | GA4GH refget protocol client + local store | [digests](python/digests.md), [refget API](python/refget-api.md), [RefgetStore](python/refget-store.md) |
| `gtars.reftx` | `gtars-refget` (transcripts) | Transcript store, MANE index, and transcript-to-genome coordinate mapping | *(see in-package docs)* |
| `gtars.vrs` | `gtars-vrs` | GA4GH VRS allele digests, ids, and normalization; HGVS parsing in `gtars.vrs.hgvs` | *(see in-package docs)* |
| `gtars.utils` | `gtars-io` | Read and write `.gtok` token files | *(see in-package docs)* |

## Quick start

```python
from gtars.models import RegionSet

# Load a BED file (local path, or URL if the http feature is enabled)
rs = RegionSet("peaks.bed")

print(len(rs))                       # region count
print(rs.identifier)                 # canonical bedbase MD5 id
print(rs.mean_region_width())        # average width

# Iterate
for region in rs:
    print(region.chr, region.start, region.end)

# Set algebra
merged = rs.reduce()                 # merge overlapping
promoters = rs.promoters(2000, 0)    # 2kb upstream
```

See [`gtars.models`](python/models.md) for the full `RegionSet` surface — including `reduce`, `setdiff`, `union`, `intersect_all`, `jaccard`, `coverage`, `neighbor_distances`, `nearest_neighbors`, and more.

## Extended example — genomic distributions

```python
from gtars.models import RegionSet, BinaryGenomeAssembly, PartitionList
from gtars.genomic_distributions import calc_gc_content, calc_expected_partitions

peaks = RegionSet("peaks.bed")
genome = BinaryGenomeAssembly("hg38.fab")
pl = PartitionList.from_gtf("gencode.v47.gtf.gz", core_prom=100, prox_prom=2000)
chrom_sizes = {"chr1": 248_956_422, "chr2": 242_193_529}

# GC content per region
gc = calc_gc_content(peaks, genome, ignore_unk_chroms=True)
print(f"Mean GC: {sum(gc) / len(gc):.3f}")

# Observed vs expected partition enrichment
result = calc_expected_partitions(peaks, pl, chrom_sizes)
for p, o, e, lo, pv in zip(
    result["partition"], result["observed"],
    result["expected"], result["log10OE"], result["pvalue"],
):
    print(f"  {p:<15} obs={o:>6.0f}  exp={e:>10.1f}  log10(O/E)={lo:+.2f}  p={pv:.2e}")
```

## Extended example — LOLA enrichment

```python
import pandas as pd
from gtars.models import RegionSet
from gtars.lola import RegionDB, run_lola

def to_tuples(rs):
    return [(r.chr, r.start, r.end) for r in rs]

db = RegionDB.from_folder("LOLACore/hg38")

user_sets = [to_tuples(RegionSet("peaks.bed"))]
universe = to_tuples(RegionSet("universe.bed"))

results = run_lola(user_sets, universe, db)
df = pd.DataFrame(results).sort_values("pValueLog", ascending=False)
print(df.head(20))
```

See [`gtars.lola`](python/lola.md) for universe preparation helpers, `RegionDB` accessors, and the full result schema.

## Performance notes

- Loading a `RegionSet` from a file sorts it; repeated operations don't re-sort. `RegionSet.from_regions` / `from_vectors` keep the input order.
- Overlap queries (`count_overlaps`, `find_overlaps`, `subset_by_overlaps`) go through an AIList index for O(log n) per-query lookups.
- `BinaryGenomeAssembly` uses `mmap` for zero-copy sequence access — use it over `GenomeAssembly` for large genomes unless you specifically need the in-memory version.
- Most numeric-result methods return Python lists; for very large outputs, consider iterating or using numpy wrappers at the call site.
