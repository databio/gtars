# Partitions (Wasm)

Wasm bindings for the gtars-genomicdist partition system — classify query regions into promoter / UTR / exon / intron / intergenic buckets given a gene model. Everything is in-memory; there is no GTF loader on the JS side (because GTF loading requires filesystem access), so you either build a `GeneModel` from pre-parsed region arrays or lift one out of a `GenomicDistAnnotation` served from the API.

See [gtars-genomicdist → Partitions](../genomicdist.md#partitions) for the full semantics — priority order, mutually-exclusive vs. bp-proportion modes, R GenomicDistributions divergences, and the chi-square p-value caveat.

## Import

```ts
import init, {
  GeneModel,
  GenomicDistAnnotation,
  PartitionList,
  calcPartitions,
  calcExpectedPartitions,
} from '@databio/gtars';

await init();
```

## `GeneModel`

Constructed from four JS arrays of region tuples: `genes`, `exons`, and optional `three_utr` / `five_utr`. Each region tuple is `[chr, start, end]` or `[chr, start, end, strand]`. When strand is provided, `+`/`-` become stranded; any other character becomes `Unstranded`.

```ts
const genes = [
  ['chr1', 1000, 5000, '+'],
  ['chr1', 8000, 12000, '-'],
];
const exons = [
  ['chr1', 1000, 2000, '+'],
  ['chr1', 4500, 5000, '+'],
  ['chr1', 11500, 12000, '-'],
];
const threeUtr = [['chr1', 4800, 5000, '+']];
const fiveUtr  = [['chr1', 1000, 1100, '+']];

// Any of three_utr / five_utr can be null if you don't have UTR annotations
const model = new GeneModel(genes, exons, threeUtr, fiveUtr);
```

!!! note "No GTF loader in Wasm"
    The Rust `GeneModel::from_gtf` and `GeneModel::from_bed_files` methods are **not** available in Wasm — they require filesystem access. In the browser, typically:

    - Fetch a pre-built `.gda` (Genomic Dist Annotation) binary from the BEDbase API and decode it with `GenomicDistAnnotation.fromBin` (see below), **or**
    - Parse regions yourself from JSON / TSV over the network and pass them into the `GeneModel` constructor.

### From a `.gda` binary

`GenomicDistAnnotation.fromBin(bytes)` loads a packed annotation. From it you can get a `GeneModel` (`geneModel()`), a `PartitionList` (`partitionList(coreSize, proxSize, chromSizes)`, with `chromSizes` optional), or a `TssIndex` (`tssIndex()`).

```ts
const bytes = new Uint8Array(await (await fetch('/api/annotations/hg38.gda')).arrayBuffer());
const gda = GenomicDistAnnotation.fromBin(bytes);
const pl = gda.partitionList(100, 2000, null);  // or pass a chromSizes object
```

## `PartitionList`

Build an ordered, priority-based partition list from a `GeneModel`. Partition priority order is hard-coded: `promoterCore` → `promoterProx` → `threeUTR` → `fiveUTR` → `exon` → `intron` (with `intergenic` added implicitly by `calcPartitions`).

```ts
const chromSizes = { chr1: 248956422 };

const pl = PartitionList.fromGeneModel(
  model,
  100,        // core promoter size (bp upstream of gene start)
  2000,       // proximal promoter size
  chromSizes, // optional — pass null to skip boundary trimming
);
```

Pass `chromSizes: null` if you don't need promoter trimming at chromosome boundaries.

## `calcPartitions`

Classify query regions into partition categories. Two modes:

- **Priority mode** (`bpProportion = false`): each query region is assigned to the first partition it overlaps; everything else goes to `intergenic`. Mutually exclusive, one region → one partition.
- **BP proportion mode** (`bpProportion = true`): for each partition, computes the total overlapping base pairs of the query set. **Not** mutually exclusive — a region that straddles two partitions contributes bp to each.

```ts
const peaks: RegionSet = /* ... */;

// Priority mode (region counts)
const result = calcPartitions(peaks, pl, false);
// {
//   partitions: [{ name: 'promoterCore', count: 42 }, ...],
//   total: 523
// }

for (const { name, count } of result.partitions) {
  const pct = (count / result.total) * 100;
  console.log(`${name.padEnd(15)} ${String(count).padStart(6)}  (${pct.toFixed(1)}%)`);
}
```

## `calcExpectedPartitions`

Observed vs. expected partition enrichment with a chi-square test. Requires chromosome sizes to compute the expected fraction based on each partition's share of the genome.

```ts
const chromSizes = { chr1: 248956422, chr2: 242193529 };

const expected = calcExpectedPartitions(peaks, pl, chromSizes, false);
// Array of { partition, observed, expected, log10Oe, pvalue }

for (const row of expected) {
  console.log(
    `${row.partition.padEnd(15)} ` +
    `obs=${row.observed.toFixed(0).padStart(6)} ` +
    `exp=${row.expected.toFixed(1).padStart(10)} ` +
    `log10(O/E)=${row.log10Oe >= 0 ? '+' : ''}${row.log10Oe.toFixed(2)} ` +
    `p=${row.pvalue.toExponential(2)}`
  );
}
```

Result fields use camelCase: `partition`, `observed`, `expected`, `log10Oe`, `pvalue`.

!!! warning "p-values will differ slightly from R GenomicDistributions"
    `calcExpectedPartitions` uses a goodness-of-fit `(O-E)²/E` formula. R's `chisq.test()` computes a 2×2 test of independence with optional Yates correction, so p-values will not match R output byte-for-byte. The `log10Oe` column is directly comparable.

## End-to-end example

```ts
import init, {
  RegionSet,
  GeneModel,
  PartitionList,
  calcExpectedPartitions,
} from '@databio/gtars';

await init();

// 1. Build query RegionSet from BED-like entries
const peaks = new RegionSet([
  ['chr1', 1500, 1800, ''],
  ['chr1', 4600, 4900, ''],
  ['chr1', 11700, 12000, ''],
]);

// 2. Build GeneModel — in a real app you'd fetch annotations from an API
const model = new GeneModel(
  [['chr1', 1000, 5000, '+'], ['chr1', 8000, 12000, '-']],
  [['chr1', 1000, 2000, '+'], ['chr1', 4500, 5000, '+'], ['chr1', 11500, 12000, '-']],
  null,
  null,
);

const chromSizes = { chr1: 248956422 };
const pl = PartitionList.fromGeneModel(model, 100, 2000, chromSizes);

// 3. Compute observed vs expected partition enrichment
const enriched = calcExpectedPartitions(peaks, pl, chromSizes, false);
console.table(enriched);
```

## See also

- **[wasm/regionset](regionset.md)** — the `RegionSet` type used as the query input.
- **[wasm/signal](signal.md)** — signal matrix overlap.
- **[gtars-genomicdist → Partitions](../genomicdist.md#partitions)** — Rust reference with the priority-order and mutually-exclusive-vs-bp-proportion semantics.
