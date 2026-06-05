# VCF → VRS in the browser (OPFS + Web Worker)

An in-browser tool that converts every variant in a VCF into a GA4GH VRS allele
identifier (`ga4gh:VA.<digest>`), **entirely client-side** with zero server
compute, using the WASM build of `gtars` (`gtars-js`). All parsing,
normalization, and digesting happen in WASM; the only network traffic is a
one-time fetch of the reference genome, which is then cached locally in OPFS.

## The core idea

> **WASM has no I/O.** JavaScript fetches bytes (network or OPFS) and copies them
> into WASM linear memory.

`wasm32` has a hard **4 GB linear-memory ceiling**, and a decoded human genome is
~3 GB. So the reference genome is held **ENCODED** (~750 MB, 2-bit packed) in a
resident `RefgetStore` and decoded **on the fly** per variant — the decoded
genome is never materialized.

## Genome source: a served gtars refgetstore

The genome comes from a directory that serves a **gtars refgetstore** (encoded
mode) over HTTP:

```
{base}/rgstore.json                      manifest: path templates, index names, mode
{base}/sequences.rgsi                    TSV index: name, length, alphabet, sha512t24u, md5
{base}/sequences/<ab>/<digest>.seq       per-sequence ENCODED bytes (no header)
```

Build one from FASTA with `gtars refget build <fasta> -o <dir>` and serve `<dir>`.
The `.seq` bytes are exactly the encoder output, so the Worker feeds them straight
into the WASM store via `add_encoded_sequence` — no decode needed in JS, and the
chromosome digests (the VRS sequence accessions) come from the index, not from
re-hashing.

## Build the WASM bundle

From the `gtars-wasm/` crate directory:

```bash
# Install once: https://rustwasm.github.io/wasm-pack/installer/
wasm-pack build --target web
```

This emits a `pkg/` directory (`gtars_js.js`, `gtars_js_bg.wasm`, …) next to
`Cargo.toml`. Both `main.js` and `worker.js` import `../../pkg/gtars_js.js`.

## Run

ES modules + Workers must be served over HTTP (not `file://`). From `gtars-wasm/`:

```bash
python3 -m http.server 8080
# then open http://localhost:8080/examples/vcf-vrs/
```

- No build step beyond `wasm-pack`, no `npm install`, no frameworks.
- The Worker is loaded as `type: "module"` so it can `import` the wasm bundle.
- **Cross-origin isolation is NOT required** — no `SharedArrayBuffer`/wasm threads,
  so no `COOP`/`COEP` headers needed.

### Steps in the UI

1. **Download / verify genome** — paste the refgetstore base URL. The Worker reads
   only the small `rgstore.json` + `sequences.rgsi` (no sequence bytes yet).
   Chromosomes are fetched **lazily**: when you drop a VCF, only the chromosomes it
   references are downloaded (and cached in OPFS for reuse). A whole-genome store has
   hundreds of contigs; a VCF touches a few, so this avoids downloading the rest.
2. **Drop a VCF** (`.vcf` or `.vcf.gz`) onto the drop zone.
3. **Compute** runs `vcf_to_vrs_ids` in the Worker; watch the progress, then
   **Download TSV** (`chrom  pos  ref  alt  vrs_id`).

## How it fits together

| File         | Thread | Role |
|--------------|--------|------|
| `index.html` | main   | Drop zone, refgetstore-URL input, progress bars, results table, TSV button. |
| `main.js`    | main   | Drag-drop, Worker orchestration, progress UI, batched-results rendering, TSV download. |
| `worker.js`  | worker | Genome provisioning (manifest/index → OPFS-cached `.seq` → resident `RefgetStore`), then `vcf_to_vrs_ids` over the dropped VCF, streaming results/progress back. |

The WASM surface used here (all real): the `RefgetStore` class
(`add_sequence` / `add_encoded_sequence` / `add_alias`) and the batch entry
`vcf_to_vrs_ids(store, vcfText, (chrom, pos, ref, alt, vrsId) => …)`. The
single-variant `hgvs_to_vrs_id` / `parse_hgvs` entries also exist (see `../hgvs/`);
`parse_hgvs` is loaded on the main thread here as a smoke test.

## Notes & current limits

- **Lazy per-chromosome loading.** Only the chromosomes a VCF references are
  downloaded and held resident (each is tens of MB encoded; cached in OPFS after the
  first touch). A whole-genome VCF still pulls the full primary set; a memory-frugal
  **byte-range** mode (read only the window each variant needs via
  `byte_range_for_bases` + a positioned OPFS `read`) would shrink even that, and is a
  future optimization.
- **Chromosome-name matching.** A VCF `CHROM` resolves to a store sequence by exact
  name or a `chr`/no-`chr` toggle. Names the store doesn't have are skipped and
  logged. Cross-naming (e.g. `chr20` ↔ `NC_000020.11`) needs the store's alias
  namespaces, not yet wired here.
- **Whole-file VCF in memory.** The dropped VCF is decoded to one string and run
  in a single `vcf_to_vrs_ids` call; progress still streams out per result batch.
  Streaming-parse of truly huge files is a future refinement.
- **A normalization error aborts the run.** Per-row error collection (skip-and-
  report) is noted in the plan but not yet implemented.

## Browser-support caveats

- **OPFS** (`navigator.storage.getDirectory`) and **sync access handles**
  (`createSyncAccessHandle()`, Worker-only) require a current Chromium/Safari or a
  recent Firefox. `main.js` warns if OPFS is missing.
- **Module Workers** (`new Worker(url, { type: "module" })`) are required.
- **`DecompressionStream("gzip")`** handles `.vcf.gz`; plain `.vcf` always works.
- **Persistent storage** (`navigator.storage.persist()`) may be granted/denied
  silently; without it the cached genome can be evicted under pressure.
- **CORS:** the refgetstore origin must permit a cross-origin `fetch` from your
  serving origin.
