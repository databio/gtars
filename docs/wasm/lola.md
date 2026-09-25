# LOLA enrichment (Wasm)

Wasm bindings for [gtars-lola](../lola.md) — LOLA (Locus Overlap Analysis) enrichment testing, runnable entirely in the browser. Given a user set of regions, a universe, and a region database (built in memory from API-served BED data), computes Fisher's exact test enrichment of the user set against every database region set.

See the [Rust gtars-lola reference](../lola.md) for the full statistical details — Fisher's exact test via hypergeometric survival function, CMLE odds ratio, contingency table semantics, and compatibility notes with R LOLA. This page is a Wasm-focused binding reference.

!!! note "In-memory only — no folder loading"
    The Wasm binding's `LolaRegionDB` is built entirely from in-memory data (`{regions, name}` entries passed in at construction). The Rust `RegionDB::from_lola_folder` loader is not available in Wasm — it requires a filesystem. In practice you'd fetch the region database from an API and pass the regions directly into the constructor, or build it from a pre-indexed IGD on the server side.

## Import

```ts
import init, {
  LolaRegionDB,
  runLOLA,
  checkUniverseAppropriateness,
} from '@databio/gtars';

await init();
```

## `LolaRegionDB`

### Construction

Construct from an array of `{ regions, name }` objects, where `regions` is an array of `[chr, start, end]` tuples and `name` is the filename-like identifier:

```ts
const dbEntries = [
  {
    name: 'encode_k562_ctcf.bed',
    regions: [
      ['chr1', 100, 200],
      ['chr1', 400, 500],
      ['chr2', 1000, 1100],
    ],
  },
  {
    name: 'encode_hela_ctcf.bed',
    regions: [
      ['chr1', 150, 250],
      ['chr2', 1050, 1150],
    ],
  },
];

const db = new LolaRegionDB(dbEntries);
```

The constructor builds an IGD overlap index internally and wraps it with minimal annotation (filename only). To attach richer per-file metadata — cell type, tissue, antibody, etc. — build the annotated `RegionDB` on the Rust side and serve it through a different transport (e.g. a custom binary format), or expose a method that accepts `RegionSetAnno` objects alongside the regions.

### Accessors

```ts
db.numRegionSets              // number
db.listRegionSets()           // string[] — filenames
db.collectionAnno             // Array of { collectionname, collector, date, source, description } (empty for an in-memory DB)

// Extract region sets by 0-based index as a RegionSetList
// (pass null for "all sets"; names are populated from filenames)
const rsl = db.getRegionSets(null);
const rsl2 = db.getRegionSets([0, 5, 12]);
```

`getRegionSets()` returns a [`RegionSetList`](regionset.md#regionsetlist) with `.names` populated from the database filenames.

## `runLOLA`

Runs Fisher's exact test for every `(user_set, db_set)` pair and returns a column-oriented object suitable for dropping directly into a DataFrame or table.

```ts
runLOLA(
  userSets: Array<Array<[string, number, number]>>,
  universe: Array<[string, number, number]>,
  regionDb: LolaRegionDB,
  minOverlap?: number,       // default 1
  direction?: string,         // "enrichment" (default) | "depletion"
): LolaResultColumnar
```

**Parameters:**
- `userSets` — **array of user sets**, where each user set is an array of `[chr, start, end]` tuples. Pass a single-element outer array `[peaks]` for a one-user-set analysis.
- `universe` — a single array of `[chr, start, end]` tuples representing the background.
- `regionDb` — a `LolaRegionDB`.
- `minOverlap` — minimum bp overlap to count as overlapping (default 1).
- `direction` — `"enrichment"` (default, P(X ≥ a), alternative "greater") or `"depletion"` (P(X ≤ a), alternative "less"). `"less"` is accepted as an alias for depletion; any other value (including `"greater"`) runs enrichment.

**Returns** an object with parallel arrays (one entry per `(user_set, db_set)` pair) using camelCase keys:

| key | type | meaning |
|---|---|---|
| `userSet` | `number[]` | 0-based user-set index |
| `dbSet` | `number[]` | 0-based db-set index |
| `collection` | `(string \| null)[]` | per-file annotation, usually null in Wasm |
| `pValueLog` | `number[]` | `-log10(p)` from Fisher's exact test, capped at ~322 |
| `oddsRatio` | `number[]` | CMLE odds ratio (matches R `fisher.test()$estimate`) |
| `support` | `number[]` | overlap count between user set and db set |
| `rnkPv`, `rnkOr`, `rnkSup` | `number[]` | 1-based per-metric ranks within the user set |
| `maxRnk` | `number[]` | max of the three ranks |
| `meanRnk` | `number[]` | mean of the three ranks |
| `b`, `c`, `d` | `number[]` | signed contingency values (can be negative) |
| `qValue` | `(number \| null)[]` | BH-adjusted p-value (applied automatically inside `runLOLA`) |
| `description`, `cellType`, `tissue`, `antibody`, `treatment`, `dataSource` | `(string \| null)[]` | per-file metadata |
| `filename` | `string[]` | db set file name |
| `size` | `number[]` | number of regions in the db set |

Like the Python binding, Wasm `runLOLA` **auto-applies** `annotate_results` + `apply_fdr_correction` internally — you don't need to call them separately. Rows come back sorted by descending `pValueLog`, then ascending `meanRnk` (matching R LOLA output order).

### Usage example

```ts
import init, {
  LolaRegionDB,
  runLOLA,
  RegionSet,
} from '@databio/gtars';

await init();

// Helper to convert a RegionSet to the tuple form LOLA expects
function toTuples(rs: RegionSet): Array<[string, number, number]> {
  // Iterate rs and pull chr/start/end — or keep the original BedEntries
  // if you still have them in scope from constructing the RegionSet.
  // ...
}

// 1. Build the region database from API-fetched data
const dbEntries = await fetchDatabaseRegions();  // [{ name, regions }, ...]
const db = new LolaRegionDB(dbEntries);

// 2. Prepare user sets and universe as tuple arrays
const userSets = [userPeakTuples];
const universe = universeTuples;

// 3. Run enrichment
const results = runLOLA(userSets, universe, db, 1, 'enrichment');

// 4. Pivot to row-oriented view for rendering
const n = results.pValueLog.length;
for (let i = 0; i < n; i++) {
  console.log(
    `${results.filename[i]}  ` +
    `p=${results.pValueLog[i].toFixed(2)}  ` +
    `OR=${results.oddsRatio[i].toFixed(2)}  ` +
    `support=${results.support[i]}`
  );
}
```

!!! note "Negative contingency values"
    If your user set contains regions *outside* the universe, the contingency table produces negative `b`/`c`/`d`. The binding matches R LOLA's behavior — it logs a warning, stores the signed values, and returns `pValueLog = 0.0`, `oddsRatio = NaN` for that row. Pre-process with `checkUniverseAppropriateness` (below) to detect this before running enrichment.

## `checkUniverseAppropriateness`

Diagnostic check — for each user set, reports what fraction of user-set regions overlap the universe and whether there are many-to-many mappings.

```ts
const report = checkUniverseAppropriateness(userSets, universe);
// Returns a column-oriented object:
// {
//   userSet: number[],
//   totalRegions: number[],
//   regionsInUniverse: number[],
//   coverage: number[],        // 0.0 to 1.0 per user set
//   manyToMany: number[],      // count of user regions overlapping >1 universe region
//   warnings: string[]         // free-form warning messages, pooled across user sets
// }

const n = report.userSet.length;
for (let i = 0; i < n; i++) {
  console.log(
    `user set ${report.userSet[i]}: ` +
    `${(report.coverage[i] * 100).toFixed(1)}% coverage, ` +
    `${report.manyToMany[i]} many-to-many`
  );
}
for (const w of report.warnings) {
  console.warn('⚠', w);
}
```

Warnings fire when coverage drops below 50% (severe) or 90% (moderate), or when any user region overlaps more than one universe region.

!!! note "Rust-only helpers not exposed in Wasm"
    The Rust universe module has two additional helpers — `redefine_user_sets` (rewrite user sets in terms of universe regions) and `build_restricted_universe` (disjoint union of all user sets, for differential analysis). Neither is currently exposed in the Wasm binding. If you need them, you can replicate the logic client-side using `RegionSetList` operations from [wasm/regionset](regionset.md), or request a binding.

## See also

- **[wasm/regionset](regionset.md)** — `RegionSet`, `RegionSetList`, and the `ConsensusBuilder` used upstream of LOLA analyses.
- **[gtars-lola](../lola.md)** — full Rust API reference, including the contingency table math, p-value computation, and CMLE odds ratio details.
