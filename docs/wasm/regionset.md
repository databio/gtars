# RegionSet & models (Wasm)

The `@databio/gtars` Wasm bindings expose gtars' core interval-analysis types directly to JavaScript and TypeScript, with no server round-trip. This page covers everything in the `regionset` Wasm module — `RegionSet`, `RegionSetList`, `ConsensusBuilder`, and the result types for statistics and BED classification.

See the Rust reference pages for the underlying semantics:

- [gtars-core](../core.md) — `Region`, `RegionSet`, `RegionSetList`
- [gtars-genomicdist](../genomicdist.md) — the stats, interval algebra, and consensus/classifier methods that `RegionSet` gains in Wasm

## Import

```ts
import init, {
  RegionSet,
  RegionSetList,
  ConsensusBuilder,
} from '@databio/gtars';

await init();  // initializes the wasm module; top-level await in ESM environments
```

## `RegionSet`

A sorted collection of genomic regions constructed from an array of `[chr, start, end, rest]` tuples. Construction sorts by `(chr, start)`.

### Construction

```ts
// BED-like entries: [chr, start, end, rest]
// `rest` is any trailing BED metadata — pass "" if you have BED3 data
const entries: [string, number, number, string][] = [
  ['chr1', 100, 200, 'peak1'],
  ['chr1', 300, 400, 'peak2'],
  ['chr2', 500, 600, 'peak3'],
];

const rs = new RegionSet(entries);
```

### Properties

```ts
rs.numberOfRegions    // number
rs.meanRegionWidth    // number
rs.nucleotidesLength  // number (total bp)
rs.identifier         // string — canonical MD5 id over first 3 columns
rs.firstRegion        // string — debug repr of the first region
rs.classify           // BedClassificationOutput (see below)
```

### Statistics

```ts
// Per-chromosome summary stats — returns { [chr]: ChromosomeStatistics }
const stats = rs.chromosomeStatistics();

// Region widths (end - start)
const widths: number[] = rs.calcWidths();

// Gaps between consecutive regions (positive gaps only; overlapping pairs are dropped)
const gaps: number[] = rs.calcNeighborDistances();

// Per-region min-neighbor distance
const nn: number[] = rs.calcNearestNeighbors();

// Single-linkage clusters at a stitching radius; returns a cluster id per region
const clusterIds: number[] = rs.cluster(5_000);
```

!!! warning "Output length for `calcNeighborDistances` / `calcNearestNeighbors`"
    Both methods **skip chromosomes with only one region** (matching R GenomicDistributions). The returned array is therefore **not aligned 1:1 with input regions**. `calcNearestNeighbors` is shorter than `rs.numberOfRegions` whenever any chromosome has a single peak; `calcNeighborDistances` returns one value per positive gap, so it is always shorter. No sentinel values are emitted.

### Region distribution

```ts
// Bin counts for plotting; chrom_sizes is optional
const chromSizes = { chr1: 248956422, chr2: 242193529 };

// With chrom_sizes → reference-aligned bins (comparable across files)
const dist = rs.regionDistribution(250, chromSizes);

// Without chrom_sizes → bins derived from observed max end (NOT comparable)
const distLocal = rs.regionDistribution(250, null);
// Array of { chr, start, end, n, rid }
```

!!! warning "`n_bins` is a target, not a total"
    When `chrom_sizes` is provided, `n_bins` is the bin count for the **longest** chromosome in `chrom_sizes`. Bin width is derived from that, and every chromosome is tiled at the same bp width, so the total bin count across chromosomes can substantially exceed `n_bins`. No chromosome gets more than `n_bins` bins: a midpoint in the leftover tail is clamped into the last bin, which stretches to the chromosome end. Regions on chromosomes missing from `chrom_sizes`, or whose midpoint falls past the stated chromosome size, are skipped.

### Interval set algebra

All operations return a new `RegionSet`. Metadata (`rest` fields) is dropped by operations that merge or synthesize new intervals.

```ts
const chromSizes = { chr1: 248956422 };

const trimmed   = rs.trim(chromSizes);
const merged    = rs.reduce();
const disjoint  = rs.disjoin();
const proms     = rs.promoters(2000, 0);
const shifted   = rs.shift(100);
const resized   = rs.resize(500, 'center');
const flanked   = rs.flank(1000, true, false);
const narrowed  = rs.narrow(100, 200, null);
const gapped    = rs.gaps(chromSizes);

const diff      = rs.setdiff(other);
const pInter    = rs.pintersect(other);
const inter     = rs.intersect(other);
const combined  = rs.concat(other);
const unioned   = rs.union(other);

const jac: number = rs.jaccard(other);
```

## `RegionSetList`

A collection of named `RegionSet`s — the gtars equivalent of Bioconductor's `GRangesList`. Provides indexed pairwise operations without copying whole `RegionSet`s on every call.

### Construction

```ts
// Empty builder + add()
const rsl = new RegionSetList();
rsl.add(rs1, 'rep1');
rsl.add(rs2, 'rep2');
rsl.add(rs3, 'rep3');

// Or build directly from BED entries for multiple sets at once
const rsl2 = RegionSetList.fromEntries(
  [
    [['chr1', 100, 200, ''], ['chr1', 300, 400, '']],  // set 1
    [['chr1', 150, 250, '']],                          // set 2
    [['chr2', 500, 600, '']],                          // set 3
  ],
  ['rep1', 'rep2', 'rep3'],   // names; pass null to leave unnamed
);
```

### Accessors

```ts
rsl.length                    // number of sets
const rs = rsl.get(0);        // RegionSet at index (throws on out-of-range)
rsl.names                     // string[] or null
const flat = rsl.concat();    // flatten into a single RegionSet (no merge)
```

### Indexed pair operations

These let you compute operations on pairs by index without shuttling full `RegionSet`s across the JS/Wasm boundary — useful for building N×N analysis grids:

```ts
const n = rsl.regionCount(0);               // regions in set 0
const pCount = rsl.pintersectCount(0, 1);   // pairwise intersect count
const jac    = rsl.jaccardAt(0, 1);         // Jaccard between sets 0 and 1
const un     = rsl.unionAt(0, 1);           // union as a new RegionSet
const diff   = rsl.setdiffAt(0, 1);         // set 0 minus set 1

const except = rsl.unionExcept(2);          // union of everything but index 2
```

### Bulk operations

O(n) N-way operations via prefix/suffix arrays.

```ts
// Union of all sets
const unionAll = rsl.unionAll();

// Intersection of all sets
const interAll = rsl.intersectAll();

// Prefix/suffix-based leave-one-out: all N "union except i" results in O(n) unions
const bulk = rsl.bulkUnionExcept();
// { union_regions, union_nucleotides, except_unique: number[] }
// except_unique[i] = number of regions unique to set i vs. the union of all others
```

### Pairwise Jaccard matrix

```ts
const result = rsl.pairwiseJaccard();
// { matrix: number[][], names: string[] | null }
//
// matrix is symmetric with 1.0 on the diagonal — one row per set, one col per set.
```

## `ConsensusBuilder`

Builder pattern for consensus region analysis — given N input region sets, compute the reduced union annotated with per-region replicate support.

```ts
import { ConsensusBuilder } from '@databio/gtars';

const cb = new ConsensusBuilder();
cb.add(rep1);
cb.add(rep2);
cb.add(rep3);

const consensus = cb.compute();
// Array of { chr, start, end, count }

// Keep regions present in ≥ 2/3 replicates
const robust = consensus.filter(r => r.count >= 2);
console.log(`${robust.length} robust regions`);
```

`count` is the number of input **sets** (not regions) that overlap each union region.

## Result types

### `ChromosomeStatistics`

Returned by `rs.chromosomeStatistics()` → `{ [chr: string]: ChromosomeStatistics }`.

| field | type |
|---|---|
| `chromosome` | `string` |
| `number_of_regions` | `number` |
| `start_nucleotide_position` | `number` — leftmost start |
| `end_nucleotide_position` | `number` — rightmost end |
| `minimum_region_length`, `maximum_region_length` | `number` |
| `mean_region_length`, `median_region_length` | `number` |

The values in the returned object are plain JS objects with these fields (serialized from Rust), not class instances.

### `RegionDistribution`

Entries in the array returned by `rs.regionDistribution(n_bins, chrom_sizes)`.

| field | type |
|---|---|
| `chr` | `string` |
| `start`, `end`, `n`, `rid` | `number` |

`n` is the count of regions whose midpoint falls in the bin; `rid` is the bin's row index within its chromosome.

### `BedClassificationOutput`

Returned by the `rs.classify` getter. Identifies the BED/ENCODE subtype based on column analysis (uses the `bedclassifier` feature, enabled by default in the Wasm build).

| field | type |
|---|---|
| `bed_compliance` | `string` |
| `data_format` | `string` — e.g. `"UcscBed"`, `"EncodeNarrowPeak"`, etc. |
| `compliant_columns`, `non_compliant_columns` | `number` |

## See also

- **[wasm/partitions](partitions.md)** — `GeneModel`, `PartitionList`, `calcPartitions` / `calcExpectedPartitions`.
- **[wasm/signal](signal.md)** — `SignalMatrix` and `calcSummarySignal`.
- **[wasm/lola](lola.md)** — `LolaRegionDB` and `runLOLA`.
- **[wasm/overlappers](overlappers.md)** — low-level overlap engine.
- **[gtars-core](../core.md)** and **[gtars-genomicdist](../genomicdist.md)** — Rust reference with algorithmic details and caveats that apply identically in Wasm.
