# Signal matrix (Wasm)

Wasm bindings for the gtars-genomicdist signal matrix overlap — overlay query regions on a region × condition matrix of signal values (e.g. a peak × cell-type ChIP intensity matrix), aggregate by MAX per query region, and compute Tukey boxplot statistics per condition.

See [gtars-genomicdist → Signal matrix overlap](../genomicdist.md#signal-matrix-overlap) for the full algorithmic detail.

## Import

```ts
import init, {
  SignalMatrix,
  calcSummarySignal,
} from '@databio/gtars';

await init();
```

## `SignalMatrix`

A region × condition matrix of signal values, loaded from either the packed binary format or built directly from JS data.

### `SignalMatrix.fromBin`

Load from a packed binary buffer (the `.sigm` format produced by the Rust `SignalMatrix::save_bin` or the CLI). This is the fast path — pass the bytes of a `.sigm` file served from your API to `SignalMatrix.fromBin`.

```ts
const response = await fetch('/api/signal-matrix/my_matrix.sigm');
const bytes = new Uint8Array(await response.arrayBuffer());

const sm = SignalMatrix.fromBin(bytes);
```

### From JS arrays

For the rare case where you have the signal data in memory already (e.g. parsed from a TSV fetched separately), you can build a `SignalMatrix` directly via the `new SignalMatrix(...)` constructor:

```ts
const regionIds = [
  'chr1_100_200',
  'chr1_300_400',
  'chr2_500_600',
];
const conditionNames = ['K562', 'HeLa', 'GM12878'];

// Flat row-major array: row i, condition j = values[i * nConditions + j]
const values = new Float64Array([
  1.2, 0.8, 2.3,   // region 0
  0.4, 1.5, 0.9,   // region 1
  3.1, 2.7, 1.8,   // region 2
]);

const sm = new SignalMatrix(
  regionIds,
  conditionNames,
  values,
  3,  // nRegions
  3,  // nConditions
);
```

!!! note "Region IDs must be `chr_start_end`"
    The constructor parses each `regionIds` entry by splitting on `_` and expects exactly three parts: chromosome, start, end. IDs with more or fewer underscores error out. This matches the R GenomicDistributions convention used by `signal_matrix.tsv` files.

## `calcSummarySignal`

Overlap a query region set against a `SignalMatrix`, take the MAX signal per query region per condition, and compute Tukey boxplot statistics per condition.

```ts
const peaks = new RegionSet([
  ['chr1', 150, 250, ''],
  ['chr2', 550, 580, ''],
]);

const result = calcSummarySignal(peaks, sm);

console.log(result.conditionNames);
// ['K562', 'HeLa', 'GM12878']

console.log(result.signalMatrix);
// Array of { region, values } — per-query-region max signal vector
// [
//   { region: "chr1_150_250", values: [1.2, 0.8, 2.3] },
//   { region: "chr2_550_580", values: [3.1, 2.7, 1.8] },
// ]

console.log(result.matrixStats);
// Per-condition Tukey stats
// [
//   { condition: 'K562',    lowerWhisker, lowerHinge, median, upperHinge, upperWhisker },
//   { condition: 'HeLa',    ... },
//   { condition: 'GM12878', ... },
// ]
```

### Result schema

**Top level:**
- `signalMatrix: { region: string, values: number[] }[]` — one entry per query region that matched at least one signal row. `region` is the query region label in `chr_start_end` form; `values` are the per-condition max signals.
- `matrixStats: ConditionStats[]` — one entry per condition, in the order of `conditionNames`.
- `conditionNames: string[]` — column labels, same as the condition names stored in the input `SignalMatrix`.

**`ConditionStats`** (per condition, standard Tukey 5-number summary):

| field | type |
|---|---|
| `condition` | `string` |
| `lowerWhisker` | `number` |
| `lowerHinge` | `number` — Q1 |
| `median` | `number` |
| `upperHinge` | `number` — Q3 |
| `upperWhisker` | `number` |

## End-to-end example

```ts
import init, { RegionSet, SignalMatrix, calcSummarySignal } from '@databio/gtars';

await init();

// 1. Fetch a packed binary signal matrix from the API
const response = await fetch('/api/signal/encode_k562_hela.sigm');
const sm = SignalMatrix.fromBin(new Uint8Array(await response.arrayBuffer()));

// 2. Build a query RegionSet from user-provided peaks
const peaks = new RegionSet(userPeakEntries);

// 3. Compute summary signal
const result = calcSummarySignal(peaks, sm);

// 4. Render per-condition boxplot
for (const stats of result.matrixStats) {
  renderBoxplot(stats.condition, stats);
}
```

## See also

- **[wasm/regionset](regionset.md)** — the `RegionSet` type used as the query input.
- **[wasm/partitions](partitions.md)** — partition classification.
- **[gtars-genomicdist → Signal matrix overlap](../genomicdist.md#signal-matrix-overlap)** — Rust reference for the underlying algorithm and the packed `.sigm` binary format.
