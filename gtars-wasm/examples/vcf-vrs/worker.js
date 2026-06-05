// worker.js — the heavy-lifting Worker for the VCF → VRS harness.
//
// Why a Worker at all:
//   1. OPFS *synchronous access handles* (createSyncAccessHandle) are only
//      available in Workers. They give us fast, seekable read()/write() over a
//      file in the Origin-Private File System — the local cache for the genome.
//   2. Processing a 1M-variant VCF must not block the main thread; doing it here
//      keeps the UI responsive (progress still streams out via postMessage even
//      while the wasm compute loop runs synchronously).
//
// The core mental model: WASM HAS NO I/O. This Worker fetches bytes (network or
// OPFS) and copies them into WASM linear memory. wasm32 has a hard 4 GB linear-
// memory ceiling, so the reference genome (~750 MB 2bit-encoded for a human
// genome) is held ENCODED and decoded on the fly per variant — we never
// materialize the ~3 GB decoded genome.
//
// Genome source: a served **gtars refgetstore** directory, i.e. a base URL that
// exposes:
//   {base}/rgstore.json                     — manifest (templates, index names, mode)
//   {base}/sequences.rgsi                    — TSV: name, length, alphabet, sha512t24u, md5
//   {base}/sequences/<first2>/<digest>.seq   — per-sequence ENCODED bytes (no header)
// The .seq bytes are exactly `encode_sequence` output, so they go straight into
// the wasm store via `add_encoded_sequence`.
//
// Build: wasm-pack build --target web   (emits ../../pkg/gtars_js.js)

import init, { RefgetStore, vcf_to_vrs_ids } from "../../pkg/gtars_js.js";

const OPFS_DIR = "refget"; // OPFS subdir holding cached .seq files
const PROGRESS_EVERY = 5000; // emit compute progress every N results
const RESULT_BATCH = 2000; // flush results to the main thread every N rows

let wasmReady = false;
// The resident genome store, built once during prepare-genome and reused across
// every process-vcf call.
let residentStore = null;

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
// (1) Provision the genome: manifest → per-sequence encoded bytes → resident
//     wasm store. .seq files are cached in OPFS so the big download happens once.
// ===========================================================================
async function prepareGenome(baseUrl) {
  const base = baseUrl.replace(/\/+$/, ""); // strip trailing slashes

  // Ask for persistent storage so the cached genome is not evicted under
  // pressure, and report headroom before a potentially large download.
  let persisted = false;
  if (navigator.storage && navigator.storage.persist) {
    persisted = await navigator.storage.persist();
  }
  if (navigator.storage && navigator.storage.estimate) {
    const est = await navigator.storage.estimate();
    post({ type: "storage-estimate", usage: est.usage, quota: est.quota, persisted });
  }

  // --- manifest ---
  workerLog(`Fetching ${base}/rgstore.json …`);
  const manifest = await fetchJson(`${base}/rgstore.json`);
  const template = manifest.seqdata_path_template || "sequences/%s2/%s.seq";
  const indexName = manifest.sequence_index || "sequences.rgsi";
  workerLog(`Store mode=${manifest.mode}, index=${indexName}.`);

  // --- sequence index ---
  const rgsiText = await (await fetchOk(`${base}/${indexName}`)).text();
  const seqs = parseRgsi(rgsiText);
  if (seqs.length === 0) throw new Error("sequences.rgsi listed no sequences.");
  workerLog(`Index lists ${seqs.length} sequences.`);

  // --- download-once each .seq into OPFS, then build the resident store ---
  const root = await navigator.storage.getDirectory();
  const dir = await root.getDirectoryHandle(OPFS_DIR, { create: true });

  const store = new RefgetStore();
  let totalBytes = 0;
  let cachedCount = 0;

  for (let i = 0; i < seqs.length; i++) {
    const s = seqs[i];
    const seqUrl = `${base}/${expandTemplate(template, s.sha512t24u)}`;
    const cacheName = `${s.sha512t24u}.seq`;

    const { bytes, cached } = await ensureCached(dir, cacheName, seqUrl);
    if (cached) cachedCount++;
    totalBytes += bytes.byteLength;
    post({ type: "download-progress", received: totalBytes, total: null });

    // Encoded alphabets go in as packed bytes; a Raw/ASCII store would carry
    // plain bases, which we ingest via the digesting path instead.
    if (isEncodedAlphabet(s.alphabet)) {
      store.add_encoded_sequence(s.name, s.sha512t24u, s.md5, s.length, bytes, s.alphabet);
    } else {
      store.add_sequence(s.name, bytes);
    }
    workerLog(`Loaded ${s.name} (${s.length} bp, ${s.alphabet}).`);
  }

  residentStore = store;
  post({
    type: "genome-ready",
    cached: cachedCount === seqs.length,
    size: totalBytes,
  });
  workerLog(
    `Resident store built: ${seqs.length} sequences, ${totalBytes} encoded bytes ` +
      `(${cachedCount} already cached in OPFS).`
  );
}

// Fetch a .seq from OPFS if present, else download it and write it to OPFS.
// Returns { bytes: Uint8Array, cached: boolean }.
async function ensureCached(dir, name, url) {
  const existing = await tryReadOpfs(dir, name);
  if (existing && existing.byteLength > 0) {
    return { bytes: existing, cached: true };
  }
  const resp = await fetchOk(url);
  const bytes = new Uint8Array(await resp.arrayBuffer());
  await writeOpfs(dir, name, bytes);
  return { bytes, cached: false };
}

async function tryReadOpfs(dir, name) {
  try {
    const fh = await dir.getFileHandle(name, { create: false });
    const access = await fh.createSyncAccessHandle();
    try {
      const size = access.getSize();
      const buf = new Uint8Array(size);
      access.read(buf, { at: 0 });
      return buf;
    } finally {
      access.close();
    }
  } catch {
    return null; // NotFoundError → not cached yet
  }
}

async function writeOpfs(dir, name, bytes) {
  const fh = await dir.getFileHandle(name, { create: true });
  const access = await fh.createSyncAccessHandle();
  try {
    access.truncate(0);
    access.write(bytes, { at: 0 });
    access.flush();
  } finally {
    access.close();
  }
}

// ===========================================================================
// (2) VCF → VRS compute, against the resident store.
// ===========================================================================
async function processVcf(buffer, isGz, name) {
  if (!residentStore) {
    throw new Error('Genome not loaded. Click "Download / verify genome" first.');
  }

  let text;
  if (isGz) {
    // DecompressionStream('gzip') is a dependency-free gunzip in modern browsers.
    const ds = new DecompressionStream("gzip");
    const stream = new Blob([buffer]).stream().pipeThrough(ds);
    text = await new Response(stream).text();
  } else {
    text = new TextDecoder().decode(buffer);
  }
  workerLog(`Decoded ${name} (${text.length} chars). Running VRS compute…`);

  // The wasm batch entry runs the whole file in one synchronous call; the
  // callback fires once per result. We accumulate rows and post them in batches,
  // and emit progress every PROGRESS_EVERY results. postMessage from inside the
  // loop still reaches the main thread promptly, so the UI updates live even
  // though this Worker thread is busy.
  let batch = [];
  let processed = 0;

  const onResult = (chrom, pos, ref, alt, vrsId) => {
    batch.push([chrom, pos, ref, alt, vrsId]);
    processed++;
    if (batch.length >= RESULT_BATCH) {
      post({ type: "results", rows: batch });
      batch = [];
    }
    if (processed % PROGRESS_EVERY === 0) {
      post({ type: "compute-progress", processed, total: null });
    }
  };

  const count = vcf_to_vrs_ids(residentStore, text, onResult);

  if (batch.length) post({ type: "results", rows: batch });
  post({ type: "compute-progress", processed: count, total: count });
  post({ type: "compute-done" });
  workerLog(`Computed ${count} VRS ids for ${name}.`);
}

// ===========================================================================
// (3) refgetstore format helpers
// ===========================================================================

// Expand a seqdata path template (e.g. "sequences/%s2/%s.seq") for a digest.
// %s2 → first 2 chars, %s4 → first 4 chars, %s → full digest. Order matters:
// replace the longer keys before the %s catch-all.
function expandTemplate(template, digest) {
  return template
    .replaceAll("%s2", digest.slice(0, 2))
    .replaceAll("%s4", digest.slice(0, 4))
    .replaceAll("%s", digest);
}

// Parse sequences.rgsi: tab-separated, '#' header/comment lines skipped.
// Columns: name, length, alphabet, sha512t24u, md5, [description].
function parseRgsi(text) {
  const out = [];
  for (const line of text.split("\n")) {
    if (line.length === 0 || line.charCodeAt(0) === 35 /* '#' */) continue;
    const c = line.split("\t");
    if (c.length < 5) continue;
    out.push({
      name: c[0],
      length: parseInt(c[1], 10),
      alphabet: c[2],
      sha512t24u: c[3],
      md5: c[4],
    });
  }
  return out;
}

// The wasm add_encoded_sequence accepts the bit-packed alphabets; ASCII/raw and
// anything else falls back to the digesting add_sequence path.
function isEncodedAlphabet(alphabet) {
  const a = alphabet.toLowerCase();
  return a === "dna2bit" || a === "dna3bit" || a === "dnaio" || a === "dnaiupac";
}

// ===========================================================================
// (4) fetch helpers
// ===========================================================================
async function fetchOk(url) {
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`fetch ${url} → HTTP ${resp.status}`);
  return resp;
}
async function fetchJson(url) {
  return (await fetchOk(url)).json();
}
