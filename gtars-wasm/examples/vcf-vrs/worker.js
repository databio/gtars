// worker.js — the heavy-lifting Worker for the VCF → VRS harness.
//
// Why a Worker at all:
//   1. OPFS *synchronous access handles* (createSyncAccessHandle) are only
//      available in Workers. They give us fast, seekable read()/write() over a
//      file in the Origin-Private File System — the local analog of HTTP Range
//      requests against a remote genome.
//   2. Processing a 1M-variant VCF must not block the main thread; doing it here
//      keeps the UI responsive.
//
// The core mental model: WASM HAS NO I/O. This Worker reads bytes from OPFS and
// copies the relevant window into WASM linear memory. wasm32 has a hard 4 GB
// linear-memory ceiling, so the reference genome (e.g. ~750 MB 2bit-encoded)
// must be held ENCODED and decoded on the fly — we never materialize the ~3 GB
// decoded genome.
//
// Build: wasm-pack build --target web   (emits ../../pkg/gtars_js.js)

// The Worker loads its OWN wasm instance (workers don't share the main thread's
// linear memory). The resident RefgetStore / batch entry points do not exist
// yet, so we only need the bundle's default init() here (see the STUB seam
// below); once the batch API lands, import RefgetStore / vcf_to_vrs_ids too.
import init from "../../pkg/gtars_js.js";

const GENOME_FILENAME = "reference.2bit"; // name of the cached file in OPFS
const PROGRESS_EVERY = 5000; // emit compute progress every N variants
const RESULT_BATCH = 2000; // flush results to main thread every N rows

let wasmReady = false;

function post(msg, transfer) {
  self.postMessage(msg, transfer || []);
}
function workerLog(text, error = false) {
  post({ type: "log", text, error });
}

self.addEventListener("message", async (ev) => {
  const msg = ev.data;
  try {
    if (!wasmReady) {
      await init();
      wasmReady = true;
      workerLog("wasm initialized in worker.");
    }
    if (msg.type === "prepare-genome") {
      await prepareGenome(msg.url);
    } else if (msg.type === "process-vcf") {
      await processVcf(msg.buffer, msg.isGz, msg.name);
    }
  } catch (err) {
    post({ type: "error", text: String(err && err.stack ? err.stack : err) });
  }
});

// ===========================================================================
// (1) OPFS download-once layer
// ===========================================================================
// On first run: request persistence, check headroom, stream the genome into an
// OPFS file chunk-by-chunk (never buffering 750 MB in JS memory). On later runs,
// detect the existing file and skip the download.
async function prepareGenome(url) {
  // Ask the browser to make this origin's storage persistent so the cached
  // genome is not silently evicted under storage pressure. This may prompt the
  // user or be granted/denied silently depending on engagement heuristics.
  let persisted = false;
  if (navigator.storage && navigator.storage.persist) {
    persisted = await navigator.storage.persist();
  }

  // Report headroom so the UI can warn before a 750 MB download.
  if (navigator.storage && navigator.storage.estimate) {
    const est = await navigator.storage.estimate();
    post({
      type: "storage-estimate",
      usage: est.usage,
      quota: est.quota,
      persisted,
    });
    const headroom = (est.quota || 0) - (est.usage || 0);
    if (headroom > 0 && headroom < 1.2 * 1024 * 1024 * 1024) {
      workerLog(
        `Only ${(headroom / 1e9).toFixed(2)} GB free — a full genome may not fit.`,
        true
      );
    }
  }

  const root = await navigator.storage.getDirectory();

  // Detect an already-cached genome and skip the download.
  const existing = await tryGetFile(root, GENOME_FILENAME);
  if (existing && existing.size > 0) {
    post({ type: "genome-ready", cached: true, size: existing.size });
    return;
  }

  // Stream the response body into OPFS chunk-by-chunk.
  workerLog(`Fetching ${url} …`);
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`fetch ${url} → HTTP ${resp.status}`);
  if (!resp.body) throw new Error("Response has no readable body to stream.");

  const total = Number(resp.headers.get("Content-Length")) || null;

  // Create the destination file handle. We write via a sync access handle so we
  // can append at a running offset without holding the whole file in memory.
  const fileHandle = await root.getFileHandle(GENOME_FILENAME, { create: true });
  const access = await fileHandle.createSyncAccessHandle();
  let offset = 0;
  try {
    access.truncate(0);
    const reader = resp.body.getReader();
    // Each `value` is a Uint8Array chunk straight off the network. We write it
    // immediately and drop it — peak JS memory stays at ~one chunk, not 750 MB.
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      access.write(value, { at: offset });
      offset += value.byteLength;
      post({ type: "download-progress", received: offset, total });
    }
    access.flush();
  } finally {
    access.close();
  }

  post({ type: "genome-ready", cached: false, size: offset });
  workerLog(`Wrote ${offset} bytes to OPFS:/${GENOME_FILENAME}.`);
}

async function tryGetFile(dir, name) {
  try {
    const fh = await dir.getFileHandle(name, { create: false });
    return await fh.getFile();
  } catch {
    return null; // NotFoundError → not cached yet
  }
}

// ===========================================================================
// (2) Two genome feed modes (both shown as commented, correct API shapes)
// ===========================================================================
// We open the cached OPFS genome with a synchronous access handle (Worker-only).
// From there, two strategies feed bytes into the resident wasm store.
async function openGenomeAccessHandle() {
  const root = await navigator.storage.getDirectory();
  const fileHandle = await root.getFileHandle(GENOME_FILENAME, { create: false });
  // Sync access handle: fast, seekable read(buffer, {at}) / write — only valid
  // inside a Worker. One exclusive handle per file at a time.
  return await fileHandle.createSyncAccessHandle();
}

// ----- Mode A: load-whole -------------------------------------------------
// Read the ENTIRE encoded genome into one buffer and hand it to wasm to build a
// resident store once. Fine for the ENCODED genome (~750 MB) — that fits under
// the 4 GB wasm32 ceiling. NEVER decode the whole thing (the ~3 GB decoded form
// would not be safe to hold).
function feedWholeGenome(access) {
  const size = access.getSize();
  const whole = new Uint8Array(size);
  // Single seekable read of the entire encoded file at offset 0.
  access.read(whole, { at: 0 });

  // TODO(gtars_browser_vcf_vrs_plan_v1.md): hand `whole` to the not-yet-built
  // wasm `RefgetStore` constructor to build a resident, encoded genome store:
  //
  //   import { RefgetStore } from "../../pkg/gtars_js.js";
  //   const store = RefgetStore.from_2bit_bytes(whole); // builds resident store
  //   return store;
  //
  // For now we just return the raw bytes so the shape is exercised.
  return whole;
}

// ----- Mode B: byte-range -------------------------------------------------
// Per-variant, compute a byte window covering the needed bases and read ONLY
// that slice — the local analog of an HTTP Range request. This avoids ever
// holding the whole genome and is how you'd scale to many contigs.
function readGenomeByteRange(access, byteStart, byteLen) {
  const slice = new Uint8Array(byteLen);
  // The real, load-bearing call shape: a seekable positioned read.
  access.read(slice, { at: byteStart });
  return slice;
}

// Window math for Mode B. The real implementation needs the wasm-safe encoder
// helpers (byte_range_for_bases / decode_substring_from_bytes_at_offset) which
// may not be wired for JS yet, so this is STUBBED.
function byteRangeForBases(/* chrom, start, end, contigByteOffsets */) {
  // TODO(gtars_browser_vcf_vrs_plan_v1.md): replace with the wasm export
  //   byte_range_for_bases(contig, start, end) -> { byteStart, byteLen }
  // which accounts for the 2bit header, per-contig offsets, and the 4-bases/byte
  // packing. Returning a tiny placeholder window so the read() shape is real.
  return { byteStart: 0, byteLen: 64 };
}

// ===========================================================================
// (3) VCF parsing + compute
// ===========================================================================
async function processVcf(buffer, isGz, name) {
  let text;
  if (isGz) {
    // DecompressionStream('gzip') is a fine, dependency-free gunzip in modern
    // browsers. Wrapping the ArrayBuffer in a Response gives us .text() cheaply.
    // (For a true 1M-line file you'd stream-decompress + stream-parse instead of
    // decompressing to one big string; that optimization is noted, not done.)
    const ds = new DecompressionStream("gzip");
    const stream = new Blob([buffer]).stream().pipeThrough(ds);
    text = await new Response(stream).text();
  } else {
    text = new TextDecoder().decode(buffer);
  }

  const variants = parseVcf(text);
  workerLog(`Parsed ${variants.length} variants from ${name}.`);

  // Open the genome handle so the two feed modes are genuinely exercised. If the
  // genome isn't cached yet we proceed anyway (the compute is stubbed).
  let access = null;
  try {
    access = await openGenomeAccessHandle();
    // Mode A demo: build the resident store once (returns raw bytes for now).
    const wholeOrStore = feedWholeGenome(access);
    workerLog(
      `Mode A: read whole encoded genome (${wholeOrStore.byteLength} bytes) for resident store.`
    );
  } catch (e) {
    workerLog(
      `Genome not available in OPFS (${e}); running compute stub without bytes. ` +
        `Run "Download / verify genome" first for the real path.`,
      true
    );
  }

  try {
    await computeVrsIdsForVcf(variants, access);
  } finally {
    if (access) access.close(); // release the exclusive sync handle
  }
}

// Minimal VCF parser: skips headers (#), splits the first five columns. Handles
// multi-allelic ALT by emitting one record per comma-separated allele.
function parseVcf(text) {
  const out = [];
  let lineStart = 0;
  const len = text.length;
  // Manual line scan avoids building a giant array of all lines at once.
  while (lineStart < len) {
    let nl = text.indexOf("\n", lineStart);
    if (nl === -1) nl = len;
    const line = text.slice(lineStart, nl);
    lineStart = nl + 1;
    if (line.length === 0 || line.charCodeAt(0) === 35 /* '#' */) continue;
    const cols = line.split("\t");
    if (cols.length < 5) continue;
    const [chrom, pos, , ref, altField] = cols;
    for (const alt of altField.split(",")) {
      out.push([chrom, pos, ref, alt]);
    }
  }
  return out;
}

// ===========================================================================
// (4) THE SEAM — where the real wasm batch API will plug in.
// ===========================================================================
// TODO(gtars_browser_vcf_vrs_plan_v1.md): replace the stub body below with the
// not-yet-built wasm batch entry:
//
//   import { RefgetStore, vcf_to_vrs_ids } from "../../pkg/gtars_js.js";
//   // store built once in Mode A (feedWholeGenome) and reused across variants:
//   vcf_to_vrs_ids(store, vcfText, (chrom, pos, ref, alt, vrsId) => { ... });
//
// The Rust side builds a RESIDENT encoded genome once and decodes windows on the
// fly per variant (Mode B byte ranges), staying under the 4 GB wasm32 ceiling.
//
// Until that exists, this STUB returns a clearly-derived placeholder id per
// variant so the OPFS/Worker/progress/TSV wiring runs end-to-end.
async function computeVrsIdsForVcf(variants, access) {
  const total = variants.length;
  let batch = [];

  for (let i = 0; i < total; i++) {
    const [chrom, pos, ref, alt] = variants[i];

    // ---- STUB compute (NOT real VRS) -------------------------------------
    // In the real path this is where we'd, per variant:
    //   (Mode B) const win = byteRangeForBases(chrom, start, end, offsets);
    //            const enc = readGenomeByteRange(access, win.byteStart, win.byteLen);
    //   then call the wasm batch/single entry to get the canonical
    //   ga4gh:VA.<digest>. We exercise the byte-range read shape on a few
    //   variants so the OPFS read path is genuinely visited.
    if (access && i < 3) {
      const win = byteRangeForBases(chrom, pos, pos);
      const enc = readGenomeByteRange(access, win.byteStart, win.byteLen);
      void enc; // bytes are read but not yet decoded (stub)
    }
    const vrsId = stubVrsId(chrom, pos, ref, alt);
    // ----------------------------------------------------------------------

    batch.push([chrom, pos, ref, alt, vrsId]);

    if (batch.length >= RESULT_BATCH) {
      post({ type: "results", rows: batch });
      batch = [];
    }
    if ((i + 1) % PROGRESS_EVERY === 0) {
      post({ type: "compute-progress", processed: i + 1, total });
      // Yield to the event loop so the Worker can flush messages and stays
      // responsive on huge inputs.
      await Promise.resolve();
    }
  }

  if (batch.length) post({ type: "results", rows: batch });
  post({ type: "compute-progress", processed: total, total });
  post({ type: "compute-done" });
}

// STUB id generator. Clearly NOT a real VRS digest — a deterministic placeholder
// derived from the variant so the table/TSV are populated and reproducible.
// TODO(gtars_browser_vcf_vrs_plan_v1.md): delete once the wasm batch entry lands.
function stubVrsId(chrom, pos, ref, alt) {
  const key = `${chrom}:${pos}:${ref}:${alt}`;
  let h = 0x811c9dc5; // FNV-1a 32-bit
  for (let i = 0; i < key.length; i++) {
    h ^= key.charCodeAt(i);
    h = Math.imul(h, 0x01000193);
  }
  const hex = (h >>> 0).toString(16).padStart(8, "0");
  return `STUB:VA.${hex}`;
}
