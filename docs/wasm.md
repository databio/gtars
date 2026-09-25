# gtars Wasm/JS

WebAssembly bindings for gtars, enabling genomic interval analysis directly in browsers and Node — no server round-trip required. The package is `@databio/gtars`; the underlying Rust crate lives in `gtars-wasm/` (Cargo package name `gtars-js`).

## Features

- Zero-installation genomic tools running client-side.
- Full `RegionSet` interval algebra and summary statistics.
- Genomic distribution analysis: partitions, signal matrices, consensus regions.
- LOLA enrichment testing with in-memory region databases.
- TypeScript-friendly — the published package ships `.d.ts` declarations generated from the Rust source.

## Installation

```bash
npm install @databio/gtars
```

## Quick start

Every Wasm entry point needs the module to be initialized once before use. In ESM environments with top-level await, that's a single line:

```ts
import init, { RegionSet, Overlapper } from '@databio/gtars';

await init();

const peaks = new RegionSet([
  ['chr1', 100, 200, 'peak1'],
  ['chr2', 150, 250, 'peak2'],
  ['chr3', 300, 400, 'peak3'],
]);

console.log(`${peaks.numberOfRegions} regions, mean width ${peaks.meanRegionWidth}`);
```

## Submodule pages

| page | covers |
|---|---|
| [**Overlappers**](wasm/overlappers.md) | `Overlapper` — low-level interval overlap engine with AIList / BITS backends. |
| [**RegionSet & models**](wasm/regionset.md) | `RegionSet`, `RegionSetList`, `ConsensusBuilder`, `ChromosomeStatistics`, `RegionDistribution`, `BedClassificationOutput`, and all the per-`RegionSet` statistics and interval algebra methods. |
| [**Partitions**](wasm/partitions.md) | `GeneModel`, `PartitionList`, `calcPartitions`, `calcExpectedPartitions` — classify regions into promoter / UTR / exon / intron / intergenic. |
| [**Signal matrix**](wasm/signal.md) | `SignalMatrix` and `calcSummarySignal` — overlay region sets on region × condition signal matrices. |
| [**LOLA enrichment**](wasm/lola.md) | `LolaRegionDB`, `runLOLA`, `checkUniverseAppropriateness` — Fisher's exact test enrichment against a database of reference region sets. |

The package also exports a few things without their own page yet: `Tokenizer`, `TssIndex` (`calcTssDistances`, `calcFeatureDistances`), `GenomicDistAnnotation` (loads a `.gda` buffer), a streaming BED parser (`bedParserNew` / `bedParserUpdate` / `bedParserFinish`), refget helpers, and HGVS/VCF to VRS functions.

## Limitations vs. the Rust/Python APIs

Wasm is a constrained runtime — a few crate features are not exposed:

- **No filesystem loaders.** `GeneModel::from_gtf`, `GeneModel::from_bed_files`, `SignalMatrix::from_tsv`, and `RegionDB::from_lola_folder` are not available. Everything is constructed from in-memory JS data or from packed binary buffers (e.g. `.sigm`, `.gda`) fetched over the network.
- **No `redefine_user_sets` / `build_restricted_universe`.** These LOLA universe helpers are Rust-only at the moment; use `checkUniverseAppropriateness` for diagnostics and replicate the set-algebra logic client-side with `RegionSetList` operations if you need the rewriting behavior.
- **Refget ships in the same package.** Functions like `digestSeqcol`, `sequenceDigest`, and the `RefgetStore` class are exported from `@databio/gtars`; `RemoteRefgetStore` is under the `@databio/gtars/remote` subpath. The refget binding docs are still being written.

## Integration example — React

```tsx
import { useEffect, useState } from 'react';
import init, { RegionSet, ConsensusBuilder } from '@databio/gtars';

function ConsensusComponent({ replicates }: { replicates: Array<[string, number, number, string][]> }) {
  const [ready, setReady] = useState(false);
  const [robust, setRobust] = useState<number | null>(null);

  useEffect(() => {
    init().then(() => setReady(true));
  }, []);

  useEffect(() => {
    if (!ready) return;
    const cb = new ConsensusBuilder();
    for (const rep of replicates) {
      cb.add(new RegionSet(rep));
    }
    const cons = cb.compute();
    setRobust(cons.filter((r) => r.count >= 2).length);
  }, [ready, replicates]);

  if (!ready) return <p>Loading gtars…</p>;
  return <p>{robust} regions present in ≥ 2 replicates</p>;
}
```

## Performance notes

- `RegionSet` construction sorts on load. Repeated operations don't re-sort.
- Set operations (`jaccard`, `setdiff`, `intersect`, etc.) merge each input and then do a linear sweep per chromosome. `Overlapper` and the partition/signal functions use an AIList index for overlap lookups.
- Statistics methods are implemented in Rust, so they're substantially faster than pure-JS equivalents for typical peak-set sizes.
- For very large numeric outputs, the Wasm glue converts between Rust `Vec<u32>` and JS arrays — if you're hitting GC pressure, consider batching calls or keeping `RegionSet` handles around instead of returning raw arrays.
