// worker.js — the heavy-lifting Worker for the VCF → VRS harness.
//
// Why a Worker at all:
//   1. OPFS *synchronous access handles* (createSyncAccessHandle) are only
//      available in Workers. They give fast, seekable read()/write() over a
//      file in the Origin-Private File System — our local genome cache.
//   2. Processing a 1M-variant VCF must not block the main thread; doing it here
//      keeps the UI responsive (progress streams out via postMessage even while
//      the wasm compute loop runs synchronously).
//
// The core mental model: WASM HAS NO I/O. This Worker fetches bytes (network or
// OPFS) and copies them into WASM linear memory. wasm32 has a hard 4 GB linear-
// memory ceiling, so the reference is held ENCODED (2-bit packed) and decoded on
// the fly per variant; the decoded genome is never materialized.
//
// LAZY loading: a whole-genome refgetstore has hundreds of contigs (alts, decoys,
// HLA, patches) and full chromosomes. A given VCF references only a few. So we
// fetch the small manifest+index up front, but download a chromosome's encoded
// .seq ONLY when a dropped VCF actually references it — caching it in OPFS so
// later runs (and later files) reuse it without re-downloading.
//
// Genome source: a served **gtars refgetstore** directory (base URL) exposing:
//   {base}/rgstore.json                      manifest (templates, index names, mode)
//   {base}/sequences.rgsi                     TSV: name, length, alphabet, sha512t24u, md5
//   {base}/sequences/<ab>/<digest>.seq        per-sequence ENCODED bytes (no header)
//
// Build: wasm-pack build --target web   (emits ../../pkg/gtars_js.js)

import init, { RefgetStore, vcf_to_vrs_ids, parse_hgvs } from "../../pkg/gtars_js.js";

const OPFS_DIR = "refget"; // OPFS subdir holding cached .seq files
const PROGRESS_EVERY = 5000; // emit compute progress every N results
const RESULT_BATCH = 2000; // flush results to the main thread every N rows

let wasmReady = false;
// Genome handle, populated by prepare-genome and reused across process-vcf calls.
let genome = null; // { base, template, byName: Map(name -> seqMeta) }
let residentStore = null; // wasm RefgetStore; grows as chromosomes are loaded
let loaded = null; // Set of sha512t24u digests already ingested into the store

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
    } else if (msg.type === "process-hgvs") {
      await processHgvs(msg.lines);
    } else if (msg.type === "list-cache") {
      await listCache();
    } else if (msg.type === "delete-genome") {
      await deleteGenome(msg.base);
    } else if (msg.type === "delete-unattributed") {
      await deleteUnattributed();
    } else if (msg.type === "purge-cache") {
      await purgeCache();
    }
  } catch (err) {
    post({ type: "error", text: String(err && err.stack ? err.stack : err) });
  }
});

// ===========================================================================
// (1) Prepare: fetch ONLY the small manifest + index. No sequence bytes yet —
//     chromosomes are downloaded lazily when a VCF references them.
// ===========================================================================
async function prepareGenome(baseUrl) {
  const base = baseUrl.replace(/\/+$/, ""); // strip trailing slashes

  let persisted = false;
  if (navigator.storage && navigator.storage.persist) {
    persisted = await navigator.storage.persist();
  }
  if (navigator.storage && navigator.storage.estimate) {
    const est = await navigator.storage.estimate();
    post({ type: "storage-estimate", usage: est.usage, quota: est.quota, persisted });
  }

  workerLog(`Fetching ${base}/rgstore.json …`);
  const manifest = await fetchJson(`${base}/rgstore.json`);
  const template = manifest.seqdata_path_template || "sequences/%s2/%s.seq";
  const indexName = manifest.sequence_index || "sequences.rgsi";

  const rgsiText = await (await fetchOk(`${base}/${indexName}`)).text();
  const byName = new Map();
  for (const s of parseRgsi(rgsiText)) byName.set(s.name, s);
  if (byName.size === 0) throw new Error("sequences.rgsi listed no sequences.");

  genome = { base, template, byName };
  residentStore = new RefgetStore();
  loaded = new Set();

  // Record this genome's identity + its full sequence-digest set, so the local
  // store panel can group cached chromosomes by genome (one row per genome).
  await registerGenome(base, byName);

  // size:null — nothing downloaded yet; sequences load on demand.
  post({ type: "genome-ready", cached: true, size: null, sequences: byName.size });
  workerLog(
    `Index ready: ${byName.size} sequences available. Chromosomes download on ` +
      `demand when a VCF references them (saved locally for reuse).`
  );
}

// Resolve a VCF chromosome name to a store sequence: exact match, else the
// chr-prefix toggle (chr20 <-> 20). Returns the seqMeta or null.
function resolveChrom(name) {
  if (genome.byName.has(name)) return genome.byName.get(name);
  const alt = name.startsWith("chr") ? name.slice(3) : `chr${name}`;
  return genome.byName.get(alt) || null;
}

// Ensure the sequence for `seqMeta` is resident in the wasm store, downloading
// its encoded .seq (OPFS-cached) on first use. Returns bytes downloaded (0 if
// already loaded). Registers `vcfName` as an alias so the engine resolves the
// exact name the VCF used.
async function ensureChromLoaded(dir, seqMeta, vcfName) {
  if (loaded.has(seqMeta.sha512t24u)) {
    if (vcfName !== seqMeta.name) residentStore.add_alias(vcfName, seqMeta.sha512t24u);
    return 0;
  }
  const cacheName = `${seqMeta.sha512t24u}.seq`;
  const seqUrl = `${genome.base}/${expandTemplate(genome.template, seqMeta.sha512t24u)}`;
  const { bytes, cached } = await ensureCached(dir, cacheName, seqUrl);

  if (isEncodedAlphabet(seqMeta.alphabet)) {
    residentStore.add_encoded_sequence(
      seqMeta.name,
      seqMeta.sha512t24u,
      seqMeta.md5,
      seqMeta.length,
      bytes,
      seqMeta.alphabet
    );
  } else {
    residentStore.add_sequence(seqMeta.name, bytes);
  }
  if (vcfName !== seqMeta.name) residentStore.add_alias(vcfName, seqMeta.sha512t24u);
  loaded.add(seqMeta.sha512t24u);
  workerLog(
    `Loaded ${seqMeta.name} (${seqMeta.length} bp, ${seqMeta.alphabet})` +
      `${cached ? " from local storage" : " — downloaded"}.`
  );
  return cached ? 0 : bytes.byteLength;
}

// ===========================================================================
// (2) Process a VCF: scan its chromosomes, lazily load just those, then compute.
// ===========================================================================
async function processVcf(buffer, isGz, name) {
  if (!genome) {
    throw new Error('Genome not prepared. Click "Download / verify genome" first.');
  }

  let text;
  if (isGz) {
    const ds = new DecompressionStream("gzip");
    const stream = new Blob([buffer]).stream().pipeThrough(ds);
    text = await new Response(stream).text();
  } else {
    text = new TextDecoder().decode(buffer);
  }

  // One pass to find the distinct chromosomes this VCF references, so we only
  // download those (not the whole genome).
  const chroms = scanVcfChroms(text);
  workerLog(`VCF references ${chroms.size} distinct chromosome(s).`);

  const root = await navigator.storage.getDirectory();
  const dir = await root.getDirectoryHandle(OPFS_DIR, { create: true });

  let downloaded = 0;
  const missing = [];
  for (const c of chroms) {
    const meta = resolveChrom(c);
    if (!meta) {
      missing.push(c);
      continue;
    }
    downloaded += await ensureChromLoaded(dir, meta, c);
    post({ type: "download-progress", received: downloaded, total: null });
  }
  if (missing.length) {
    workerLog(
      `Skipping ${missing.length} chromosome(s) not in the store: ` +
        `${missing.slice(0, 10).join(", ")}${missing.length > 10 ? " …" : ""}`,
      true
    );
  }

  workerLog(`Running VRS compute over ${name} …`);
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

// Convert a batch of HGVS expressions (one per line) against the resident
// genome, lazily loading whichever chromosomes they reference. Emits results as
// [hgvs, vrs_id] rows (bad/unconvertible lines get an "ERROR: …" id).
async function processHgvs(lines) {
  if (!genome) {
    throw new Error('Genome not prepared. Click "Download / verify genome" first.');
  }
  const root = await navigator.storage.getDirectory();
  const dir = await root.getDirectoryHandle(OPFS_DIR, { create: true });

  // Resolve + load referenced sequences once (the accession in each HGVS).
  const accessions = new Set();
  for (const line of lines) {
    try {
      const p = parse_hgvs(line);
      if (p && p.accession) accessions.add(p.accession);
    } catch {
      /* unparseable; the per-line conversion below reports it */
    }
  }
  let downloaded = 0;
  const missing = [];
  for (const acc of accessions) {
    const meta = resolveChrom(acc);
    if (!meta) {
      missing.push(acc);
      continue;
    }
    downloaded += await ensureChromLoaded(dir, meta, acc);
    post({ type: "download-progress", received: downloaded, total: null });
  }
  if (missing.length) {
    workerLog(
      `Skipping HGVS referencing sequences not in the store: ${missing.slice(0, 10).join(", ")}`,
      true
    );
  }

  let batch = [];
  let processed = 0;
  for (const line of lines) {
    let row;
    try {
      row = [line, residentStore.hgvs_to_vrs_id(line)];
    } catch (e) {
      row = [line, `ERROR: ${String(e && e.message ? e.message : e)}`];
    }
    batch.push(row);
    processed++;
    if (batch.length >= RESULT_BATCH) {
      post({ type: "results", rows: batch });
      batch = [];
    }
    if (processed % PROGRESS_EVERY === 0) {
      post({ type: "compute-progress", processed, total: lines.length });
    }
  }
  if (batch.length) post({ type: "results", rows: batch });
  post({ type: "compute-progress", processed: lines.length, total: lines.length });
  post({ type: "compute-done" });
  workerLog(`Converted ${lines.length} HGVS expression(s).`);
}

// Distinct CHROM values (column 0) across all non-header VCF rows.
function scanVcfChroms(text) {
  const set = new Set();
  let lineStart = 0;
  const len = text.length;
  while (lineStart < len) {
    let nl = text.indexOf("\n", lineStart);
    if (nl === -1) nl = len;
    if (text.charCodeAt(lineStart) !== 35 /* '#' */ && nl > lineStart) {
      const tab = text.indexOf("\t", lineStart);
      if (tab !== -1 && tab < nl) set.add(text.slice(lineStart, tab));
    }
    lineStart = nl + 1;
  }
  return set;
}

// ===========================================================================
// (2b) Local genome store: inventory grouped BY GENOME, not by file.
//
// The on-disk cache is just per-sequence `<digest>.seq` files. To present one
// row per genome, we keep a tiny registry (`genomes.json`) mapping each genome's
// base URL to its label and full set of sequence digests, written when a genome
// is verified. The panel then reports, per genome, how many of its sequences are
// cached and their total size.
// ===========================================================================
const REGISTRY_FILE = "genomes.json";

async function storageEstimate() {
  if (navigator.storage && navigator.storage.estimate) {
    const e = await navigator.storage.estimate();
    return { usage: e.usage, quota: e.quota };
  }
  return { usage: null, quota: null };
}

function deriveLabel(base) {
  const parts = base.split("/").filter((p) => p && !p.endsWith(":"));
  return parts[parts.length - 1] || base;
}

async function readRegistry(dir) {
  try {
    const fh = await dir.getFileHandle(REGISTRY_FILE, { create: false });
    const f = await fh.getFile();
    return JSON.parse(await f.text());
  } catch {
    return {};
  }
}

async function writeRegistry(dir, reg) {
  const fh = await dir.getFileHandle(REGISTRY_FILE, { create: true });
  const access = await fh.createSyncAccessHandle();
  try {
    const bytes = new TextEncoder().encode(JSON.stringify(reg));
    access.truncate(0);
    access.write(bytes, { at: 0 });
    access.flush();
  } finally {
    access.close();
  }
}

async function registerGenome(base, byName) {
  const root = await navigator.storage.getDirectory();
  const dir = await root.getDirectoryHandle(OPFS_DIR, { create: true });
  const reg = await readRegistry(dir);
  const digests = [...new Set([...byName.values()].map((m) => m.sha512t24u))];
  reg[base] = { label: deriveLabel(base), digests };
  await writeRegistry(dir, reg);
}

// List cached chromosomes grouped by genome. Emits one entry per registered
// genome (cached/total counts + bytes on disk) plus an "other" bucket for any
// cached files not attributed to a known genome (e.g. cached before the genome
// was registered — verifying that genome will attribute them).
async function listCache() {
  const est = await storageEstimate();
  const root = await navigator.storage.getDirectory();
  let dir;
  try {
    dir = await root.getDirectoryHandle(OPFS_DIR, { create: false });
  } catch {
    post({ type: "genome-list", genomes: [], other: { count: 0, bytes: 0 }, totalBytes: 0, ...est });
    return;
  }

  // Sizes of every cached .seq on disk.
  const sizes = new Map();
  let totalBytes = 0;
  for await (const [name, handle] of dir.entries()) {
    if (handle.kind !== "file" || !name.endsWith(".seq")) continue;
    const f = await handle.getFile();
    sizes.set(name, f.size);
    totalBytes += f.size;
  }

  const reg = await readRegistry(dir);
  const attributed = new Set();
  const genomes = [];
  for (const [base, info] of Object.entries(reg)) {
    let cachedCount = 0;
    let cachedBytes = 0;
    for (const d of info.digests) {
      const fn = `${d}.seq`;
      if (sizes.has(fn)) {
        cachedCount++;
        cachedBytes += sizes.get(fn);
        attributed.add(fn);
      }
    }
    genomes.push({
      base,
      label: info.label || base,
      cachedCount,
      totalCount: info.digests.length,
      cachedBytes,
      active: genome && genome.base === base,
    });
  }

  let otherCount = 0;
  let otherBytes = 0;
  for (const [fn, sz] of sizes) {
    if (!attributed.has(fn)) {
      otherCount++;
      otherBytes += sz;
    }
  }

  genomes.sort((a, b) => b.cachedBytes - a.cachedBytes);
  post({ type: "genome-list", genomes, other: { count: otherCount, bytes: otherBytes }, totalBytes, ...est });
}

// Delete one genome's cached chromosomes. Sequences shared with another
// registered genome are kept; the genome's registry entry is removed.
async function deleteGenome(base) {
  const root = await navigator.storage.getDirectory();
  const dir = await root.getDirectoryHandle(OPFS_DIR, { create: false });
  const reg = await readRegistry(dir);
  const info = reg[base];
  if (!info) {
    await listCache();
    return;
  }
  const keep = new Set();
  for (const [b, i] of Object.entries(reg)) {
    if (b !== base) for (const d of i.digests) keep.add(d);
  }
  let removed = 0;
  for (const d of info.digests) {
    if (keep.has(d)) continue;
    try {
      await dir.removeEntry(`${d}.seq`);
      removed++;
    } catch {
      /* not cached */
    }
    if (loaded) loaded.delete(d);
  }
  delete reg[base];
  await writeRegistry(dir, reg);
  workerLog(`Deleted genome "${info.label || base}" (${removed} chromosome file(s)) from local storage.`);
  await listCache();
}

// Delete every cached .seq not claimed by a registered genome (e.g. files left
// by an earlier download before the genome was registered). Verifying that
// genome instead would attribute them; this just frees the disk.
async function deleteUnattributed() {
  const root = await navigator.storage.getDirectory();
  const dir = await root.getDirectoryHandle(OPFS_DIR, { create: false });
  const reg = await readRegistry(dir);
  const known = new Set();
  for (const info of Object.values(reg)) for (const d of info.digests) known.add(`${d}.seq`);

  const toDelete = [];
  for await (const [name, handle] of dir.entries()) {
    if (handle.kind === "file" && name.endsWith(".seq") && !known.has(name)) toDelete.push(name);
  }
  let removed = 0;
  for (const name of toDelete) {
    try {
      await dir.removeEntry(name);
      removed++;
      if (loaded) loaded.delete(name.replace(/\.seq$/, ""));
    } catch {
      /* ignore */
    }
  }
  workerLog(`Deleted ${removed} unattributed chromosome file(s) from local storage.`);
  await listCache();
}

async function purgeCache() {
  const root = await navigator.storage.getDirectory();
  try {
    await root.removeEntry(OPFS_DIR, { recursive: true });
  } catch {
    /* dir may not exist yet */
  }
  if (loaded) loaded.clear();
  workerLog("Deleted all genomes from local storage.");
  await listCache();
}

// ===========================================================================
// (3) OPFS cache: read a .seq if present, else download and write it.
// ===========================================================================
async function ensureCached(dir, name, url) {
  const existing = await tryReadOpfs(dir, name);
  if (existing && existing.byteLength > 0) return { bytes: existing, cached: true };
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
// (4) refgetstore format + fetch helpers
// ===========================================================================
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
    out.push({ name: c[0], length: parseInt(c[1], 10), alphabet: c[2], sha512t24u: c[3], md5: c[4] });
  }
  return out;
}

function isEncodedAlphabet(alphabet) {
  const a = alphabet.toLowerCase();
  return a === "dna2bit" || a === "dna3bit" || a === "dnaio" || a === "dnaiupac";
}

async function fetchOk(url) {
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`fetch ${url} → HTTP ${resp.status}`);
  return resp;
}
async function fetchJson(url) {
  return (await fetchOk(url)).json();
}
