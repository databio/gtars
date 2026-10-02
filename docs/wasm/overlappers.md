# Overlappers (Wasm)

Wasm binding for [gtars-overlaprs](../overlaprs.md): a low-level interval overlap engine. Build an `Overlapper` once from a set of reference regions (the "universe"), then ask which of those regions each query region overlaps. For whole-region-set operations (intersections, unions, statistics), use [`RegionSet`](regionset.md) instead.

## Import

```ts
import init, { Overlapper } from '@databio/gtars';

await init();
```

## Construction

```ts
const universe = [
  ['chr1', 100, 200, 'a'],
  ['chr1', 150, 250, 'b'],
  ['chr2', 300, 400, 'c'],
];

const overlapper = new Overlapper(universe, 'ailist');
```

- `universe`: an array of `[chr, start, end, name]` entries. The binding reads each entry as a four-element tuple, so include a `name` (any string; it is not used).
- `backend`: the index to build for each chromosome, either `'ailist'` (augmented interval list) or `'bits'` (binary interval search). Any other value throws `Invalid backend specified`.

## Methods

### `get_backend()`

Returns the backend name the overlapper was built with (`'ailist'` or `'bits'`).

### `find(regions)`

Takes an array of query regions in the same `[chr, start, end, name]` form and returns one flat array of the universe regions they overlap, each as `[chr, start, end]`. Overlaps from all queries are concatenated in query order, so a universe region that overlaps two queries appears twice. Queries on a chromosome that is not in the universe contribute nothing.

```ts
const hits = overlapper.find([['chr1', 180, 220, 'q1']]);
// [['chr1', 100, 200], ['chr1', 150, 250]]
```
