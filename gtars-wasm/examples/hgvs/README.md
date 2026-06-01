# HGVS → VRS browser demo

A minimal browser tool that converts a genomic (`g.`) HGVS string into a GA4GH
VRS allele identifier (`ga4gh:VA.<digest>`) entirely client-side, using the
WASM build of `gtars-vrs`.

Scope (v1): genomic (`g.`) variants only. Transcript (`c.`/`n.`) variants need
the mmap-backed `gtars-reftx` store, which is not yet WASM-safe.

## Build the WASM bundle

From the `gtars-wasm/` crate directory:

```bash
# Install once: https://rustwasm.github.io/wasm-pack/installer/
wasm-pack build --target web
```

This compiles `gtars-js` to `wasm32-unknown-unknown` and emits a `pkg/`
directory (`gtars_js.js`, `gtars_js_bg.wasm`, …) next to `Cargo.toml`. The demo
page imports `../../pkg/gtars_js.js`.

> The HGVS entry point pulls only the in-memory refget store — no filesystem or
> HTTP features — so the bundle stays browser-safe.

## Run the demo

WASM modules must be served over HTTP (not `file://`). From `gtars-wasm/`:

```bash
python3 -m http.server 8080
# then open http://localhost:8080/examples/hgvs/
```

Enter an HGVS expression (e.g. `chrF:g.6C>T`), the sequence name/accession it
references (must match the HGVS reference, e.g. `chrF`), and the reference bases
— either pasted directly or fetched from a URL. Click **Compute VRS id**.

## How it works

```js
import init, { hgvs_to_vrs_id, parse_hgvs } from "../../pkg/gtars_js.js";
await init();

const bases = await (await fetch(seqUrl)).text(); // JS owns the async fetch
const id = hgvs_to_vrs_id("chrF:g.6C>T", "chrF", bases.replace(/\s+/g, ""));
// → "ga4gh:VA.<digest>"
```

- `parse_hgvs(hgvs)` — pure parse-only summary `{ accession, gene,
  reference_type }`; throws on a parse error. Mirrors the Python
  `gtars.vrs.hgvs.parse_hgvs` surface.
- `hgvs_to_vrs_id(hgvs, sequenceName, sequenceBases)` — builds an in-memory
  refget store from the supplied bases, then runs the readonly HGVS→VRS bridge
  (`hgvs_str_to_vrs_id_readonly`) with the g.-only `NoTranscriptProvider`.

### Reference bases and canonical digests

To get a **canonical** VRS id the supplied bases must be the *full* bases of the
referenced sequence, so the computed refget digest (`SQ.<sha512t24u>`) matches
the true reference. Supplying only a window changes the digest and the
coordinates. A future revision can fetch just the covering bytes via HTTP
`Range` and decode them with the WASM-safe `byte_range_for_bases` /
`decode_substring_from_bytes_at_offset` encoder helpers (route (b) in the
implementation plan).
