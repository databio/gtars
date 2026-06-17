// remote-refget-store.js — a lean, byte-range + OPFS-caching refget client for
// the browser, shipped as part of the @databio/gtars package.
//
// ---------------------------------------------------------------------------
// What this is
// ---------------------------------------------------------------------------
// A `RemoteRefgetStore` reads reference sequence out of a *served* gtars
// refgetstore directory WITHOUT downloading whole chromosomes. It issues HTTP
// `Range:` requests for just the bytes a region needs, decodes them with the
// wasm `encodedByteRange` / `decodeEncodedRange` primitives, and caches fetched
// byte-windows in OPFS (Origin Private File System) so repeated / overlapping
// reads are served from local disk with no network.
//
// This is the browser twin of the native remote-range read path in gtars
// (`open_remote_range` / `get_substring_from_remote`): wasm has no HTTP client,
// so JavaScript does the `fetch`; the wasm primitives do the bit-packing byte
// math and the on-the-fly decode. WASM never has to hold (or even download) a
// whole encoded chromosome to read a small region.
//
// ---------------------------------------------------------------------------
// Three-flow model (mirrors the native bindings, adapted for the browser)
// ---------------------------------------------------------------------------
//   Flow 1 — whole-sequence (resident): `prefetchSequence(name)` downloads the
//            entire encoded `.seq` once (OPFS-cached) for hot/dense access. The
//            classic "load the chromosome" path; best when an input touches most
//            of a contig.
//   Flow 2 — byte-range (DEFAULT): `getSubstring(name, start, end)` fetches only
//            the covering bytes via an HTTP `Range:` request and decodes them.
//            Lean: a few KB instead of 16–60 MB per chromosome for sparse reads.
//   Flow 3 — cached windows: every byte-range fetch is persisted in OPFS keyed
//            by (digest, byteStart, byteEnd); overlapping later reads are served
//            locally. `getSubstring` transparently prefers a covering cached
//            window (whole-sequence cache included) before going to the network.
//
// ---------------------------------------------------------------------------
// Store layout (served gtars refgetstore)
// ---------------------------------------------------------------------------
//   {base}/rgstore.json                       manifest (path templates, indexes)
//   {base}/collections.rgci                   collection index
//   {base}/collections/<coll-digest>.rgsi     a collection's sequence index
//   {base}/sequences.rgsi                     (flat-store) the sequence index
//   {base}/sequences/<ab>/<seq-digest>.seq    per-sequence ENCODED bytes (no header)
//
// ---------------------------------------------------------------------------
// Usage (inside a Web Worker — OPFS sync handles are Worker-only)
// ---------------------------------------------------------------------------
//   import init, { encodedByteRange, decodeEncodedRange }
//       from "@databio/gtars";
//   import { RemoteRefgetStore } from "@databio/gtars/remote";
//   await init();
//   const store = new RemoteRefgetStore({
//     baseUrl: "https://example.org/mystore",
//     wasm: { encodedByteRange, decodeEncodedRange },
//   });
//   await store.openCollection(collectionDigest);   // loads that genome's index
//   const bases = await store.getSubstring("chr1", 1000, 1100); // 100 bp, ~26 B fetched

// Default OPFS subdirectory holding cached sequence bytes/windows.
const DEFAULT_OPFS_DIR = "refget";
// Encoded alphabets whose `.seq` bytes are bit-packed and byte-addressable.
const ENCODED_ALPHABETS = new Set(["dna2bit", "dna3bit", "dnaio", "dnaiupac"]);
// When fetching a byte-range, pad the request by this many BASES on each side and
// snap to the byte grid, so a run of nearby small reads (a sorted VCF) reuses one
// fetched/cached window instead of issuing a request per variant. Tunable.
const DEFAULT_WINDOW_PAD_BASES = 4096;

/**
 * A byte-range + OPFS-caching reader for a served gtars refgetstore.
 *
 * The genome is never downloaded whole (unless you ask via `prefetchSequence`);
 * `getSubstring` pulls only the covering bytes and decodes them with the wasm
 * primitives. Construct one per store; call `openCollection` (or
 * `loadSequenceIndex` for a flat store) to bind a genome's sequence index, then
 * read regions with `getSubstring`.
 */
export class RemoteRefgetStore {
  /**
   * @param {object} opts
   * @param {string} opts.baseUrl                 Base URL of the served store.
   * @param {object} opts.wasm                    The wasm module (or the two fns).
   *        Must expose `encodedByteRange(start, end, alphabet) -> [bs, be]` and
   *        `decodeEncodedRange(bytes, byteOffset, start, end, alphabet) -> string`.
   * @param {string} [opts.opfsDir]               OPFS subdir for the cache.
   * @param {number} [opts.windowPadBases]        Base padding per range fetch.
   * @param {string} [opts.seqTemplate]           Override seq path template.
   * @param {string} [opts.collTemplate]          Override collection path template.
   * @param {(text:string, isError?:boolean)=>void} [opts.log]  Optional logger.
   */
  constructor(opts) {
    if (!opts || !opts.baseUrl) throw new Error("RemoteRefgetStore: baseUrl is required");
    const wasm = opts.wasm;
    if (!wasm || typeof wasm.encodedByteRange !== "function" || typeof wasm.decodeEncodedRange !== "function") {
      throw new Error(
        "RemoteRefgetStore: opts.wasm must expose encodedByteRange and decodeEncodedRange " +
          "(import them from @databio/gtars and call init() first)"
      );
    }
    this.base = opts.baseUrl.replace(/\/+$/, "");
    this.wasm = wasm;
    this.opfsDir = opts.opfsDir || DEFAULT_OPFS_DIR;
    this.windowPadBases = opts.windowPadBases ?? DEFAULT_WINDOW_PAD_BASES;
    this.log = opts.log || (() => {});

    // Filled by prepare()/openCollection(); both have sane defaults.
    this.seqTemplate = opts.seqTemplate || "sequences/%s2/%s.seq";
    this.collTemplate = opts.collTemplate || "collections/%s.rgsi";

    // Sequence index of the currently-open genome: name -> seqMeta and the same
    // keyed by digest. seqMeta = { name, length, alphabet, sha512t24u, md5 }.
    this.byName = new Map();
    this.byDigest = new Map();

    // In-memory cache of OPFS-resident byte windows per digest, so repeated reads
    // in one session skip even the OPFS open. Map(digest -> Array<{bs, be, bytes}>).
    // The whole-sequence cache is just a window with bs=0, be=fileLength.
    this._windows = new Map();
    this._opfsDirHandle = null;
  }

  // -------------------------------------------------------------------------
  // Store / genome setup
  // -------------------------------------------------------------------------

  /**
   * Fetch the store manifest (rgstore.json) and adopt its path templates. Safe
   * to skip if you passed templates to the constructor or use a known layout.
   */
  async prepare() {
    const manifest = await this._fetchJson(`${this.base}/rgstore.json`);
    this.seqTemplate = manifest.seqdata_path_template || this.seqTemplate;
    this.collTemplate = manifest.collections_path_template || this.collTemplate;
    return manifest;
  }

  /**
   * Bind a specific collection (genome/assembly) by its digest: loads that
   * collection's sequence index so names resolve. After this, `getSubstring`
   * works for any sequence in the collection.
   * @param {string} collectionDigest
   */
  async openCollection(collectionDigest) {
    const url = `${this.base}/${expandTemplate(this.collTemplate, collectionDigest)}`;
    const text = await (await this._fetchOk(url)).text();
    this._adoptSequences(parseRgsi(text));
    this.collectionDigest = collectionDigest;
    return this.byName.size;
  }

  /**
   * Load a flat store's top-level `sequences.rgsi` (no collection layer). Use
   * this for stores that serve a single, un-collected sequence index.
   * @param {string} [indexName]
   */
  async loadSequenceIndex(indexName = "sequences.rgsi") {
    const text = await (await this._fetchOk(`${this.base}/${indexName}`)).text();
    this._adoptSequences(parseRgsi(text));
    return this.byName.size;
  }

  /** Register a pre-parsed sequence list (e.g. an index fetched by the caller). */
  setSequences(seqMetas) {
    this._adoptSequences(seqMetas);
    return this.byName.size;
  }

  _adoptSequences(seqMetas) {
    this.byName = new Map();
    this.byDigest = new Map();
    for (const s of seqMetas) {
      this.byName.set(s.name, s);
      this.byDigest.set(s.sha512t24u, s);
    }
    if (this.byName.size === 0) throw new Error("RemoteRefgetStore: index listed no sequences.");
  }

  /** seqMeta for a sequence name, or null. EXACT match — no chr-prefix guessing. */
  resolve(name) {
    return this.byName.get(name) || null;
  }

  // -------------------------------------------------------------------------
  // Reads
  // -------------------------------------------------------------------------

  /**
   * Read bases `[start, end)` (0-based, half-open) of sequence `name`, fetching
   * only the covering bytes via an HTTP `Range:` request and decoding them. The
   * fetched window is padded + OPFS-cached so overlapping later reads are local.
   *
   * @param {string} name   Sequence name (must be in the open index).
   * @param {number} start  0-based inclusive start base.
   * @param {number} end    0-based exclusive end base.
   * @returns {Promise<string>} The decoded ASCII bases.
   */
  async getSubstring(name, start, end) {
    const meta = this.resolve(name);
    if (!meta) throw new Error(`RemoteRefgetStore: sequence "${name}" not in the open genome`);
    if (start < 0 || end > meta.length || start > end) {
      throw new Error(`RemoteRefgetStore: range [${start}, ${end}) out of bounds for ${name} (len ${meta.length})`);
    }
    if (start === end) return "";

    const alphabet = meta.alphabet;
    if (!ENCODED_ALPHABETS.has(alphabet.toLowerCase())) {
      // Plain (unencoded) `.seq`: byte offset == base offset, no bit packing.
      const bytes = await this._readBytes(meta, start, end);
      return new TextDecoder().decode(bytes);
    }

    // Exact covering byte range for the requested bases.
    const [bs, be] = this.wasm.encodedByteRange(start, end, alphabet);

    // Pad the byte window (snapped to the encoded grid) so nearby reads coalesce.
    const [padBs, padBe] = this._paddedByteRange(start, end, meta);

    // Serve from an OPFS-cached window if one already covers [bs, be); else fetch
    // the padded window, cache it, and decode the exact requested bases from it.
    const window = await this._ensureWindow(meta, bs, be, padBs, padBe);
    return this.wasm.decodeEncodedRange(window.bytes, window.bs, start, end, alphabet);
  }

  /**
   * Flow 1: download the ENTIRE encoded `.seq` for `name` once (OPFS-cached) and
   * keep it resident, so subsequent reads of any region are served locally with
   * zero further network. Use for dense access to a contig. Returns the bytes
   * downloaded (0 if already cached). The bytes themselves remain wasm's to
   * decode per-read via `getSubstring`.
   */
  async prefetchSequence(name) {
    const meta = this.resolve(name);
    if (!meta) throw new Error(`RemoteRefgetStore: sequence "${name}" not in the open genome`);
    const fileLen = encodedByteLength(meta);
    // A whole-file window (bs=0) covers every later read.
    const existing = this._coveringWindow(meta.sha512t24u, 0, fileLen);
    if (existing) return 0;

    const cacheName = `${meta.sha512t24u}.seq`;
    const dir = await this._dir();
    const onDisk = await tryReadOpfs(dir, cacheName);
    let bytes;
    let downloaded = 0;
    if (onDisk && onDisk.byteLength > 0) {
      bytes = onDisk;
    } else {
      const url = this._seqUrl(meta);
      bytes = new Uint8Array(await (await this._fetchOk(url)).arrayBuffer());
      await writeOpfs(dir, cacheName, bytes);
      downloaded = bytes.byteLength;
    }
    this._addWindow(meta.sha512t24u, 0, bytes.byteLength, bytes);
    this.log(`Prefetched ${meta.name} (${meta.length} bp)${downloaded ? " — downloaded" : " from local storage"}.`);
    return downloaded;
  }

  /**
   * Return the WHOLE encoded `.seq` bytes for `name`, downloading once if needed
   * (OPFS-cached). This is the bridge to the synchronous wasm VRS compute path:
   * `vcf_to_vrs_ids` / `hgvs_to_vrs_id` cannot `await fetch()` mid-loop, so they
   * need the sequence resident in the wasm `RefgetStore`; hand these bytes to
   * `store.add_encoded_sequence(...)`. Equivalent to `prefetchSequence` but
   * returns the bytes (and an `alreadyCached` flag) for ingestion.
   *
   * NOTE: this downloads the whole chromosome. True byte-range VRS compute is
   * blocked on a missing wasm windowed-ingestion API (`add_encoded_window` +
   * roll-edge signaling); see the byte-range VRS plan. Until that lands, the
   * compute path is whole-sequence; `getSubstring` is the lean byte-range read
   * for direct sequence retrieval.
   *
   * @returns {Promise<{bytes: Uint8Array, cached: boolean, meta: SeqMeta}>}
   */
  async getSequenceBytes(name) {
    const meta = this.resolve(name);
    if (!meta) throw new Error(`RemoteRefgetStore: sequence "${name}" not in the open genome`);
    const cacheName = `${meta.sha512t24u}.seq`;
    const dir = await this._dir();
    const onDisk = await tryReadOpfs(dir, cacheName);
    if (onDisk && onDisk.byteLength > 0) {
      this._addWindow(meta.sha512t24u, 0, onDisk.byteLength, onDisk);
      return { bytes: onDisk, cached: true, meta };
    }
    const url = this._seqUrl(meta);
    const bytes = new Uint8Array(await (await this._fetchOk(url)).arrayBuffer());
    await writeOpfs(dir, cacheName, bytes);
    this._addWindow(meta.sha512t24u, 0, bytes.byteLength, bytes);
    return { bytes, cached: false, meta };
  }

  // -------------------------------------------------------------------------
  // Window cache (in-memory + OPFS), keyed by (digest, byteStart, byteEnd)
  // -------------------------------------------------------------------------

  // Return an in-memory window covering [bs, be), or null.
  _coveringWindow(digest, bs, be) {
    const list = this._windows.get(digest);
    if (!list) return null;
    for (const w of list) {
      if (w.bs <= bs && w.be >= be) return w;
    }
    return null;
  }

  _addWindow(digest, bs, be, bytes) {
    let list = this._windows.get(digest);
    if (!list) {
      list = [];
      this._windows.set(digest, list);
    }
    // Drop any existing window now subsumed by this one, to bound the list.
    const kept = list.filter((w) => !(bs <= w.bs && be >= w.be));
    kept.push({ bs, be, bytes });
    this._windows.set(digest, kept);
  }

  // Ensure a window covering encoded bytes [coverBs, coverBe) is resident,
  // fetching the padded range [padBs, padBe) from OPFS or network if needed.
  // Returns the window { bs, be, bytes } (bs is the absolute byte offset of
  // bytes[0], as decodeEncodedRange's byteOffset expects).
  async _ensureWindow(meta, coverBs, coverBe, padBs, padBe) {
    const digest = meta.sha512t24u;

    // 1) In-memory window already covers the needed bytes?
    let w = this._coveringWindow(digest, coverBs, coverBe);
    if (w) return w;

    const dir = await this._dir();

    // 2) Whole-sequence file cached in OPFS? (Flow 1 result, or a prior prefetch.)
    const wholeName = `${digest}.seq`;
    const whole = await tryReadOpfs(dir, wholeName);
    if (whole && whole.byteLength >= coverBe) {
      this._addWindow(digest, 0, whole.byteLength, whole);
      return this._coveringWindow(digest, coverBs, coverBe);
    }

    // 3) A previously-cached byte window in OPFS that covers [coverBs, coverBe)?
    const cachedWin = await this._readCachedWindow(dir, digest, coverBs, coverBe);
    if (cachedWin) {
      this._addWindow(digest, cachedWin.bs, cachedWin.be, cachedWin.bytes);
      return this._coveringWindow(digest, coverBs, coverBe);
    }

    // 4) Fetch the padded window over HTTP Range, cache it, return it.
    const bytes = await this._fetchRange(meta, padBs, padBe);
    this._addWindow(digest, padBs, padBs + bytes.byteLength, bytes);
    await this._writeCachedWindow(dir, digest, padBs, padBs + bytes.byteLength, bytes);
    return this._coveringWindow(digest, coverBs, coverBe);
  }

  // OPFS filename for a cached byte window. Windows are stored as separate small
  // files so reads never load the whole chromosome. The whole-file cache uses the
  // bare `<digest>.seq` name (shared with prefetchSequence / the legacy path).
  _windowFileName(digest, bs, be) {
    return `${digest}.${bs}-${be}.win`;
  }

  // Find and read a cached window file for `digest` that covers [coverBs, coverBe).
  // Window files are named `<digest>.<bs>-<be>.win`; we scan the directory for a
  // covering one. (Directory scans are cheap relative to a network round-trip.)
  async _readCachedWindow(dir, digest, coverBs, coverBe) {
    const prefix = `${digest}.`;
    try {
      for await (const [fname, handle] of dir.entries()) {
        if (handle.kind !== "file" || !fname.startsWith(prefix) || !fname.endsWith(".win")) continue;
        const range = parseWindowName(fname, digest);
        if (!range) continue;
        if (range.bs <= coverBs && range.be >= coverBe) {
          const bytes = await tryReadOpfs(dir, fname);
          if (bytes && bytes.byteLength > 0) return { bs: range.bs, be: range.bs + bytes.byteLength, bytes };
        }
      }
    } catch {
      /* directory iteration unsupported / empty — fall through to network */
    }
    return null;
  }

  async _writeCachedWindow(dir, digest, bs, be, bytes) {
    try {
      await writeOpfs(dir, this._windowFileName(digest, bs, be), bytes);
    } catch (e) {
      // Caching is an optimization; a write failure must not break the read.
      this.log(`window cache write failed for ${digest} [${bs},${be}): ${e}`, true);
    }
  }

  // -------------------------------------------------------------------------
  // Byte fetching
  // -------------------------------------------------------------------------

  // Pad the requested bases by windowPadBases on each side, clamp to the
  // sequence, and convert to a byte range snapped to the encoding grid.
  _paddedByteRange(start, end, meta) {
    const pStart = Math.max(0, start - this.windowPadBases);
    const pEnd = Math.min(meta.length, end + this.windowPadBases);
    const [bs, be] = this.wasm.encodedByteRange(pStart, pEnd, meta.alphabet);
    return [bs, be];
  }

  // HTTP Range fetch of encoded bytes [bs, be). Returns the fetched Uint8Array.
  // Falls back to a full GET if the server ignores Range (200 instead of 206) and
  // slices locally, so correctness never depends on Range support.
  async _fetchRange(meta, bs, be) {
    const url = this._seqUrl(meta);
    const resp = await fetch(url, { headers: { Range: `bytes=${bs}-${be - 1}` } });
    if (!resp.ok) throw new Error(`fetch ${url} (Range ${bs}-${be - 1}) → HTTP ${resp.status}`);
    const buf = new Uint8Array(await resp.arrayBuffer());
    if (resp.status === 206) return buf;
    // Server ignored Range and returned the whole file: slice the window out.
    return buf.subarray(bs, be);
  }

  // Plain (unencoded) byte read: base offset == byte offset.
  async _readBytes(meta, start, end) {
    return this._fetchRange(meta, start, end);
  }

  _seqUrl(meta) {
    return `${this.base}/${expandTemplate(this.seqTemplate, meta.sha512t24u)}`;
  }

  // -------------------------------------------------------------------------
  // OPFS + fetch helpers
  // -------------------------------------------------------------------------

  async _dir() {
    if (this._opfsDirHandle) return this._opfsDirHandle;
    const root = await navigator.storage.getDirectory();
    this._opfsDirHandle = await root.getDirectoryHandle(this.opfsDir, { create: true });
    return this._opfsDirHandle;
  }

  async _fetchOk(url) {
    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`fetch ${url} → HTTP ${resp.status}`);
    return resp;
  }

  async _fetchJson(url) {
    return (await this._fetchOk(url)).json();
  }
}

// ===========================================================================
// Module-level helpers (shared by the methods above)
// ===========================================================================

/** Expand a refgetstore path template against a digest. */
function expandTemplate(template, digest) {
  return template
    .replaceAll("%s2", digest.slice(0, 2))
    .replaceAll("%s4", digest.slice(0, 4))
    .replaceAll("%s", digest);
}

/** Bytes of an encoded sequence: ceil(length * bits/symbol / 8). */
function encodedByteLength(meta) {
  const a = meta.alphabet.toLowerCase();
  const bits = a === "dna2bit" ? 2 : a === "dna3bit" ? 3 : 4; // dnaio/dnaiupac = 4
  return Math.ceil((meta.length * bits) / 8);
}

/** Parse `<digest>.<bs>-<be>.win` into { bs, be } (the recorded byte range). */
function parseWindowName(fname, digest) {
  const tail = fname.slice(digest.length + 1, -4); // strip "<digest>." and ".win"
  const dash = tail.indexOf("-");
  if (dash < 0) return null;
  const bs = parseInt(tail.slice(0, dash), 10);
  const be = parseInt(tail.slice(dash + 1), 10);
  if (Number.isNaN(bs) || Number.isNaN(be)) return null;
  return { bs, be };
}

/**
 * Parse a sequences.rgsi index: tab-separated, '#' comment/header lines skipped.
 * Columns: name, length, alphabet, sha512t24u, md5, [description].
 */
export function parseRgsi(text) {
  const out = [];
  for (const line of text.split("\n")) {
    if (line.length === 0 || line.charCodeAt(0) === 35 /* '#' */) continue;
    const c = line.split("\t");
    if (c.length < 5) continue;
    out.push({ name: c[0], length: parseInt(c[1], 10), alphabet: c[2], sha512t24u: c[3], md5: c[4] });
  }
  return out;
}

/** Read an OPFS file fully into a Uint8Array, or null if it doesn't exist. */
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

/** Write bytes to an OPFS file (truncating any existing content). */
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
