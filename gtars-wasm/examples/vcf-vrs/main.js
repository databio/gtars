// main.js — main-thread orchestration for the VCF → VRS browser harness.
//
// Responsibilities (all REAL, not stubbed):
//   - Load the wasm bundle on the main thread for pure helpers (parse_hgvs).
//   - Drag-drop / file-picker a VCF, read it to an ArrayBuffer.
//   - Spawn the Worker and drive it via postMessage; the Worker owns OPFS +
//     (stubbed) wasm compute, because OPFS sync access handles are Worker-only.
//   - Render download + compute progress and accumulate results into a
//     downloadable TSV.
//
// Build: wasm-pack build --target web   (emits ../../pkg/gtars_js.js)

// The main-thread wasm import is used only for the pure parse helper; the
// heavy lifting happens in the Worker, which loads its own wasm instance.
import init, { parse_hgvs } from "../../pkg/gtars_js.js";

// ---------------------------------------------------------------------------
// DOM handles
// ---------------------------------------------------------------------------
const $ = (id) => document.getElementById(id);
const genomeUrlInput = $("genomeUrl");
const prepBtn = $("prepBtn");
const dlBar = $("dlBar");
const storageInfo = $("storageInfo");
const dropZone = $("drop");
const filePicker = $("filePicker");
const computeBar = $("computeBar");
const computeInfo = $("computeInfo");
const downloadBtn = $("downloadBtn");
const resultsBody = $("results").querySelector("tbody");
const logEl = $("log");
const cacheRefreshBtn = $("cacheRefreshBtn");
const cachePurgeBtn = $("cachePurgeBtn");
const cacheInfo = $("cacheInfo");
const cacheBody = $("cacheTable").querySelector("tbody");

const MAX_PREVIEW_ROWS = 50;

function log(msg, isError = false) {
  const line = `${new Date().toLocaleTimeString()}  ${msg}\n`;
  logEl.textContent += line;
  logEl.scrollTop = logEl.scrollHeight;
  if (isError) console.error(msg);
  else console.log(msg);
}

function setBar(bar, fraction, done = false) {
  bar.style.width = `${Math.max(0, Math.min(1, fraction)) * 100}%`;
  bar.classList.toggle("done", done);
}

// ---------------------------------------------------------------------------
// Results accumulation + TSV download
// ---------------------------------------------------------------------------
// Rows are kept as compact arrays [chrom, pos, ref, alt, vrs_id]. For a 1M
// variant file this is the only thing the main thread holds (the Worker streams
// results in batches), so we keep it lean and only render a small preview.
const TSV_HEADER = ["chrom", "pos", "ref", "alt", "vrs_id"];
let resultRows = [];
let previewCount = 0;

function addResults(rows) {
  for (const r of rows) {
    resultRows.push(r);
    if (previewCount < MAX_PREVIEW_ROWS) {
      const tr = document.createElement("tr");
      for (const cell of r) {
        const td = document.createElement("td");
        td.textContent = cell;
        tr.appendChild(td);
      }
      resultsBody.appendChild(tr);
      previewCount++;
    }
  }
}

function buildTsv() {
  // Build the TSV lazily at download time. Joining ~1M short rows is fine in
  // one pass; if it ever isn't, switch to a streamed Blob via a generator.
  const lines = [TSV_HEADER.join("\t")];
  for (const r of resultRows) lines.push(r.join("\t"));
  return lines.join("\n") + "\n";
}

downloadBtn.addEventListener("click", () => {
  const blob = new Blob([buildTsv()], { type: "text/tab-separated-values" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "vcf_vrs.tsv";
  a.click();
  URL.revokeObjectURL(url);
});

// ---------------------------------------------------------------------------
// Worker setup
// ---------------------------------------------------------------------------
// type: "module" so the Worker can `import` the wasm ES module bundle.
const worker = new Worker(new URL("./worker.js", import.meta.url), {
  type: "module",
});

let computeStartTs = 0;

worker.addEventListener("message", (ev) => {
  const msg = ev.data;
  switch (msg.type) {
    case "log":
      log(`[worker] ${msg.text}`, !!msg.error);
      break;

    case "download-progress": {
      // msg.received / msg.total (bytes). total may be null if the server
      // didn't send Content-Length, in which case we show received only.
      if (msg.total) {
        setBar(dlBar, msg.received / msg.total);
        storageInfo.textContent =
          `Downloading… ${fmtBytes(msg.received)} / ${fmtBytes(msg.total)}`;
      } else {
        storageInfo.textContent = `Downloading… ${fmtBytes(msg.received)}`;
      }
      break;
    }

    case "genome-ready":
      setBar(dlBar, 1, true);
      storageInfo.textContent =
        `Index ready: ${msg.sequences ?? "?"} sequences available. ` +
        `Chromosomes download on demand when you drop a VCF.`;
      prepBtn.disabled = false;
      log(`Genome index ready (${msg.sequences ?? "?"} sequences).`);
      requestCacheList();
      break;

    case "storage-estimate":
      storageInfo.textContent =
        `Storage: ${fmtBytes(msg.usage)} used / ${fmtBytes(msg.quota)} quota` +
        (msg.persisted ? " — persisted ✓" : " — NOT persisted (may be evicted)");
      break;

    case "compute-progress": {
      const frac = msg.total ? msg.processed / msg.total : 0;
      setBar(computeBar, frac);
      const elapsed = (performance.now() - computeStartTs) / 1000;
      const rate = msg.processed / Math.max(elapsed, 0.001);
      computeInfo.textContent =
        `${msg.processed.toLocaleString()}` +
        (msg.total ? ` / ${msg.total.toLocaleString()}` : "") +
        ` variants  (${rate.toFixed(0)}/s)`;
      break;
    }

    case "results":
      // Batched results: array of [chrom, pos, ref, alt, vrs_id].
      addResults(msg.rows);
      break;

    case "compute-done":
      setBar(computeBar, 1, true);
      downloadBtn.disabled = resultRows.length === 0;
      computeInfo.textContent =
        `Done: ${resultRows.length.toLocaleString()} variants in ` +
        `${((performance.now() - computeStartTs) / 1000).toFixed(1)}s.`;
      log(`Compute done: ${resultRows.length} results.`);
      requestCacheList(); // newly-downloaded chromosomes show up
      break;

    case "cache-list":
      renderCache(msg);
      break;

    case "error":
      log(`ERROR: ${msg.text}`, true);
      break;
  }
});

worker.addEventListener("error", (ev) => {
  log(`Worker error: ${ev.message} (${ev.filename}:${ev.lineno})`, true);
});

// ---------------------------------------------------------------------------
// OPFS cache panel: inventory cached chromosomes, delete one, or purge all.
// ---------------------------------------------------------------------------
function requestCacheList() {
  worker.postMessage({ type: "list-cache" });
}

function renderCache(msg) {
  cacheBody.innerHTML = "";
  const used = msg.usage != null ? fmtBytes(msg.usage) : "?";
  const quota = msg.quota != null ? fmtBytes(msg.quota) : "?";
  cacheInfo.textContent =
    `${msg.entries.length} chromosome(s) saved locally, ${fmtBytes(msg.totalBytes)}` +
    `  ·  browser storage used: ${used} / ${quota}`;

  if (msg.entries.length === 0) {
    const tr = document.createElement("tr");
    const td = document.createElement("td");
    td.colSpan = 4;
    td.className = "hint";
    td.textContent = "No chromosomes downloaded yet.";
    tr.appendChild(td);
    cacheBody.appendChild(tr);
    return;
  }

  for (const e of msg.entries) {
    const tr = document.createElement("tr");
    const label = e.label || "(unknown — prepare a genome to map names)";
    const shortFile = e.name.length > 20 ? e.name.slice(0, 12) + "…" + e.name.slice(-7) : e.name;
    const cells = [label + (e.resident ? " · resident" : ""), shortFile, fmtBytes(e.size)];
    for (const c of cells) {
      const td = document.createElement("td");
      td.textContent = c;
      tr.appendChild(td);
    }
    const tdBtn = document.createElement("td");
    const del = document.createElement("button");
    del.textContent = "Delete";
    del.addEventListener("click", () => {
      if (confirm(`Delete cached ${label} (${fmtBytes(e.size)})?`)) {
        worker.postMessage({ type: "delete-cache", name: e.name });
      }
    });
    tdBtn.appendChild(del);
    tr.appendChild(tdBtn);
    cacheBody.appendChild(tr);
  }
}

cacheRefreshBtn.addEventListener("click", requestCacheList);
cachePurgeBtn.addEventListener("click", () => {
  if (confirm("Delete ALL downloaded chromosomes from local storage? They will re-download on next use.")) {
    worker.postMessage({ type: "purge-cache" });
  }
});

// ---------------------------------------------------------------------------
// Genome prep button → ask Worker to download-once into OPFS
// ---------------------------------------------------------------------------
prepBtn.addEventListener("click", () => {
  const url = genomeUrlInput.value.trim();
  if (!url) {
    log("Please provide a genome URL.", true);
    return;
  }
  prepBtn.disabled = true;
  setBar(dlBar, 0);
  storageInfo.textContent = "Requesting persistent storage…";
  log(`Preparing genome from ${url}`);
  worker.postMessage({ type: "prepare-genome", url });
});

// ---------------------------------------------------------------------------
// Drag-drop / file-picker VCF on the MAIN thread
// ---------------------------------------------------------------------------
dropZone.addEventListener("click", () => filePicker.click());
filePicker.addEventListener("change", () => {
  if (filePicker.files.length) handleFile(filePicker.files[0]);
});

["dragenter", "dragover"].forEach((evt) =>
  dropZone.addEventListener(evt, (e) => {
    e.preventDefault();
    dropZone.classList.add("drag");
  })
);
["dragleave", "drop"].forEach((evt) =>
  dropZone.addEventListener(evt, (e) => {
    e.preventDefault();
    dropZone.classList.remove("drag");
  })
);
dropZone.addEventListener("drop", (e) => {
  const file = e.dataTransfer?.files?.[0];
  if (file) handleFile(file);
});

async function handleFile(file) {
  // Reset result state for a fresh run.
  resultRows = [];
  previewCount = 0;
  resultsBody.innerHTML = "";
  downloadBtn.disabled = true;
  setBar(computeBar, 0);
  computeInfo.textContent = "";

  const isGz = /\.gz$/i.test(file.name);
  log(`Loaded "${file.name}" (${fmtBytes(file.size)})${isGz ? " — gzip" : ""}.`);

  // Read the whole File into an ArrayBuffer. For very large VCFs a streamed
  // approach (file.stream()) would be better, but reading-to-buffer keeps the
  // scaffold simple; gzip VCFs are compact, plain VCFs are still typically
  // manageable. The Worker is what keeps the UI responsive during compute.
  const buf = await file.arrayBuffer();

  computeStartTs = performance.now();
  // Transfer the ArrayBuffer to the Worker (zero-copy ownership handoff): after
  // this `buf` is detached/unusable on the main thread, which is exactly what
  // we want — the Worker now owns those bytes.
  worker.postMessage(
    { type: "process-vcf", buffer: buf, isGz, name: file.name },
    [buf]
  );
  log("VCF bytes transferred to Worker for processing.");
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
function fmtBytes(n) {
  if (n == null) return "?";
  const u = ["B", "KB", "MB", "GB", "TB"];
  let i = 0;
  while (n >= 1024 && i < u.length - 1) {
    n /= 1024;
    i++;
  }
  return `${n.toFixed(i === 0 ? 0 : 1)} ${u[i]}`;
}

// ---------------------------------------------------------------------------
// Main-thread wasm init (pure helpers only)
// ---------------------------------------------------------------------------
await init();
log("Main-thread wasm initialized (parse_hgvs available).");

// Tiny sanity demo of the existing pure entry point — proves the bundle loads.
try {
  const parsed = parse_hgvs("chrF:g.6C>T");
  log(`parse_hgvs("chrF:g.6C>T") → ${JSON.stringify(parsed)}`);
} catch (e) {
  log(`parse_hgvs demo failed: ${e}`, true);
}

if (!("storage" in navigator) || !navigator.storage.getDirectory) {
  log(
    "WARNING: This browser can't store genomes locally (no OPFS support). " +
      "The local genome store and compute will not function.",
    true
  );
} else {
  // Show whatever is already cached from a previous session (names resolve once
  // a genome index is prepared).
  requestCacheList();
}
