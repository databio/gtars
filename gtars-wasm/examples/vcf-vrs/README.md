# VCF → VRS browser harness (OPFS + Web Worker scaffold)

A **starting-point scaffold** for an in-browser tool that converts every variant
in a VCF into a GA4GH VRS allele identifier (`ga4gh:VA.<digest>`), entirely
client-side, using the WASM build of `gtars` (`gtars-js`).

This example is deliberately a **harness**: it implements the full, real
JavaScript / Web Worker / OPFS plumbing — the hard, browser-specific part — with
one clearly-marked seam where the not-yet-built WASM batch API will plug in. The
plan for that WASM work is `gtars_browser_vcf_vrs_plan_v1.md`.

## The core idea

> **WASM has no I/O.** JavaScript reads bytes from a browser store and copies
> them into WASM linear memory.

`wasm32` has a hard **4 GB linear-memory ceiling**, and a decoded human genome is
~3 GB. So the reference genome is held **ENCODED** (e.g. a ~750 MB `.2bit`) and
decoded **on the fly** per variant — we never materialize the decoded genome.

## Build the WASM bundle

From the `gtars-wasm/` crate directory:

```bash
# Install once: https://rustwasm.github.io/wasm-pack/installer/
wasm-pack build --target web
```

This emits a `pkg/` directory (`gtars_js.js`, `gtars_js_bg.wasm`, …) next to
`Cargo.toml`. Both `main.js` and `worker.js` import `../../pkg/gtars_js.js`.

## Run the demo

ES modules + Workers must be served over HTTP (not `file://`). From
`gtars-wasm/`:

```bash
python3 -m http.server 8080
# then open http://localhost:8080/examples/vcf-vrs/
```

- No build step beyond `wasm-pack`, no `npm install`, no frameworks.
- The Worker is loaded as `type: "module"` so it can `import` the wasm bundle —
  this requires a reasonably modern browser (see caveats).
- **Cross-origin isolation is NOT required.** This harness does not use
  `SharedArrayBuffer` or wasm threads, so no `COOP`/`COEP` headers are needed.
  (If a future version adds wasm threading you'd then need
  `Cross-Origin-Opener-Policy: same-origin` +
  `Cross-Origin-Embedder-Policy: require-corp`.)

### Steps in the UI

1. **Download / verify genome** — streams the genome URL into OPFS once, then
   reuses it. UCSC's `hg38.2bit` (~750 MB) is the default; any reachable
   encoded genome URL works (CORS must allow the fetch).
2. **Drop a VCF** (`.vcf` or `.vcf.gz`) onto the drop zone.
3. **Compute** runs in the Worker; watch the progress bar, then **Download TSV**.

## What genuinely works now vs. what is stubbed

### ✅ Real and working (this is the value of the scaffold)

- **OPFS download-once layer** (`worker.js → prepareGenome`):
  `navigator.storage.persist()`, `estimate()` headroom check, and a streamed
  `fetch` whose response body is written into an OPFS file **chunk-by-chunk**
  via a sync access handle — peak JS memory stays at ~one network chunk, not
  750 MB. Re-runs detect the cached file and skip the download.
- **Web Worker** doing all OPFS + compute work, because OPFS *synchronous access
  handles* (`createSyncAccessHandle()`) are Worker-only.
- **Both genome feed modes**, as real `handle.read(buffer, { at })` call shapes:
  - **Mode A — load-whole:** one `read(whole, { at: 0 })` of the entire encoded
    genome to build a resident store once.
  - **Mode B — byte-range:** per-variant positioned `read(slice, { at: byteStart })`,
    the local analog of an HTTP Range request.
- **Drag-drop / file-picker VCF** on the main thread, read to an `ArrayBuffer`,
  with **gzip** handled by `DecompressionStream("gzip")`.
- **Zero-copy handoff:** the VCF `ArrayBuffer` is **transferred** to the Worker.
- **VCF parsing** (header skip, first 5 columns, multi-allelic ALT split).
- **Progress + results streaming:** the Worker posts compute progress every N
  variants and results in batches; the main thread renders progress bars, a
  preview table, and assembles a downloadable **TSV**
  (`chrom  pos  ref  alt  vrs_id`). The Worker yields to the event loop so a
  **1M-variant** file does not freeze the UI.
- **Quota / eviction messaging** via `persist()` + `estimate()`.

### ⛔ Stubbed (pending `gtars_browser_vcf_vrs_plan_v1.md`)

Every stub is marked with a `TODO:` pointing at the plan.

- **`computeVrsIdsForVcf(...)` — the seam.** Currently returns a clearly-labeled
  placeholder id (`STUB:VA.<fnv-hash>`) per variant. It will be replaced by the
  WASM batch entry `vcf_to_vrs_ids(store, vcfText, onResult)`.
- **Resident `RefgetStore`** built once from the encoded genome bytes
  (`feedWholeGenome` returns the raw bytes today; the real path constructs a
  WASM `RefgetStore` and reuses it across all variants).
- **Byte-window math** (`byteRangeForBases`) — returns a placeholder window. The
  real path uses the WASM-safe encoder helpers `byte_range_for_bases` /
  `decode_substring_from_bytes_at_offset` (which may not be exported to JS yet).
  The `handle.read({ at })` *call shape* around it is real.

> The single-variant `hgvs_to_vrs_id` / `parse_hgvs` entries DO exist today (see
> `../hgvs/`). `parse_hgvs` is loaded on the main thread here as a smoke test.

## File map

| File         | Thread | Role |
|--------------|--------|------|
| `index.html` | main   | Drop zone, genome-URL input, progress bars, results table, TSV button. Loads the wasm bundle (main thread) and spawns the Worker. |
| `main.js`    | main   | Drag-drop, Worker orchestration, progress UI, TSV download. |
| `worker.js`  | worker | OPFS download-once, sync access handle, both feed modes, the `computeVrsIdsForVcf` seam + stub, progress/result messaging. |

## Browser-support caveats

- **OPFS** (`navigator.storage.getDirectory`) and **sync access handles** are
  available in current Chromium and Safari; Firefox shipped OPFS sync access
  handles relatively recently — use an up-to-date browser. `main.js` warns if
  OPFS is missing.
- **Module Workers** (`new Worker(url, { type: "module" })`) are required to
  `import` the wasm bundle in the Worker; supported in current browsers.
- **`DecompressionStream("gzip")`** is broadly available in current browsers but
  not ancient ones; plain `.vcf` always works.
- **Persistent storage** (`navigator.storage.persist()`) may be granted or
  denied silently based on the browser's site-engagement heuristics; without it
  the cached genome can be evicted under storage pressure (the UI says so).
- **CORS:** the genome URL must permit a cross-origin `fetch` from your serving
  origin. UCSC's `hgdownload.soe.ucsc.edu` generally does; mirror it yourself if
  not.
