import { MSG, FILTER_DEFAULTS } from "./protocol.js?v=20261023-1";
import { SEQ_RE, describeGroups, pairSequenceFiles } from "./pairing.js?v=20261023-1";

const ASSET_VERSION = "20261023-1";

const MAX_INDEX_BYTES = 2 * 1024 * 1024 * 1024;

const INDEX_URL =
  "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panhuman-1.k31w21.pidx";
const INDEX_FILENAME = "panhuman-1.k31w21.pidx";
const INDEX_DISPLAY_NAME = "Index: panhuman-1";

const $ = (id) => document.getElementById(id);
const indexStatus = $("index-status");
const indexInput = $("index-input");
const indexName = $("index-name");
const indexInfo = $("index-info");
const indexProgress = $("index-progress");
const indexRetry = $("index-retry");
const indexDownload = $("index-download");
const dirZone = $("dir-zone");
const dirInput = $("dir-input");
const dirSummary = $("dir-summary");
const modeEl = $("mode");
const absThresholdEl = $("abs-threshold");
const relThresholdEl = $("rel-threshold");
const filterBtn = $("filter-btn");
const resetBtn = $("reset-btn");
const progressWrap = $("progress-wrap");
const dehostProgress = $("dehost-progress");
const dehostLabel = $("dehost-label");
const fileListEl = $("file-list");
const statusEl = $("status");
const purgeIndexEl = $("purge-index");

modeEl.value = FILTER_DEFAULTS.deplete ? "deplete" : "search";
absThresholdEl.value = FILTER_DEFAULTS.absThreshold;
relThresholdEl.value = FILTER_DEFAULTS.relThreshold;

let selectedFiles = [];
let fileRows = [];
let totalBytes = 0;
let isFiltering = false;
let indexLoaded = false;
let indexFilename = "";
let workerReady = false;
let pendingIndexBuffer = null;
let indexFailed = false;
let indexNeedsDownload = false;

function setStatus(msg, kind) {
  statusEl.textContent = msg;
  statusEl.className = kind || "";
}

function humanBytes(n) {
  const units = ["B", "KB", "MB", "GB", "TB"];
  let i = 0;
  while (n >= 1024 && i < units.length - 1) {
    n /= 1024;
    i++;
  }
  return `${n.toFixed(i === 0 ? 0 : 1)}${units[i]}`;
}

function humanBases(bp) {
  const units = ["bp", "kbp", "Mbp", "Gbp", "Tbp"];
  let i = 0;
  while (bp >= 1000 && i < units.length - 1) {
    bp /= 1000;
    i++;
  }
  return `${bp.toFixed(i === 0 ? 0 : 1)}${units[i]}`;
}

function setButtonDisabled(btn, disabled) {
  btn.disabled = disabled;
  btn.classList.toggle("pulse-glow", !disabled);
}

function readyStatus() {
  if (!indexLoaded) return "";
  if (!selectedFiles.length) return "Index loaded. Select sequences to filter";
  return "Press 'Filter' to begin";
}

function updateFilterEnabled() {
  setButtonDisabled(filterBtn, !indexLoaded || isFiltering || !selectedFiles.length);
}

function filterParams() {
  const abs = parseInt(absThresholdEl.value, 10);
  const rel = parseFloat(relThresholdEl.value);
  return {
    deplete: modeEl.value === "deplete",
    absThreshold: Number.isNaN(abs) ? FILTER_DEFAULTS.absThreshold : abs,
    relThreshold: Number.isNaN(rel) ? FILTER_DEFAULTS.relThreshold : rel,
  };
}

const GZ_RE = /\.gz$/i;
const gzName = (name) => (GZ_RE.test(name) ? name : name + ".gz");

function setSelection(files) {
  const seqs = files.filter((f) => SEQ_RE.test(f.name));
  try {
    selectedFiles = pairSequenceFiles(seqs);
  } catch (err) {
    selectedFiles = [];
    totalBytes = 0;
    dirZone.classList.remove("loaded");
    dirSummary.textContent = "Pairing error";
    setStatus(err?.message || String(err), "error");
    renderFileList();
    updateFilterEnabled();
    return;
  }
  totalBytes = selectedFiles.reduce((sum, group) => sum + group.size, 0);

  if (seqs.length === 0) {
    dirZone.classList.remove("loaded");
    dirSummary.textContent = files.length
      ? "No fastq/fasta files found in selection"
      : "";
  } else {
    dirZone.classList.add("loaded");
    const pairs = selectedFiles.filter((group) => group.kind === "paired").length;
    const singles = selectedFiles.length - pairs;
    const parts = [];
    if (pairs) parts.push(`${pairs} pair${pairs === 1 ? "" : "s"}`);
    if (singles) parts.push(`${singles} unpaired`);
    dirSummary.textContent = `${seqs.length} sequence file${
      seqs.length === 1 ? "" : "s"
    } selected (${parts.join(", ")}), ${humanBytes(totalBytes)} in total`;
  }
  setStatus(readyStatus());
  progressWrap.style.display = "none";
  renderFileList();
  updateFilterEnabled();
}

function renderFileList() {
  fileListEl.innerHTML = "";
  fileRows = selectedFiles.map((group) => {
    const li = document.createElement("li");
    const name = document.createElement("span");
    name.className = "name";
    name.textContent = group.label;
    const st = document.createElement("span");
    st.className = "st";
    li.append(name, st);
    fileListEl.appendChild(li);
    return { li, st };
  });
  fileListEl.style.display = fileRows.length ? "" : "none";
}

const worker = new Worker(`./worker.js?v=${ASSET_VERSION}`, { type: "module" });

let onOutputBatch = null;
let onFilterProgress = null;
let onWorkerError = null;

worker.onmessage = (e) => {
  const m = e.data;
  switch (m.type) {
    case MSG.READY:
      workerReady = true;
      if (pendingIndexBuffer) {
        const buf = pendingIndexBuffer;
        pendingIndexBuffer = null;
        worker.postMessage({ type: MSG.LOAD_INDEX, data: buf }, [buf]);
      }
      break;
    case MSG.INDEX_LOADED:
      indexLoaded = true;
      indexFailed = false;
      indexNeedsDownload = false;
      indexStatus.classList.remove("failed", "needs-download");
      indexStatus.classList.add("loaded");
      indexProgress.style.display = "none";
      indexRetry.hidden = true;
      indexDownload.hidden = true;
      clearIndexLoading();
      indexInfo.textContent = m.info;
      setStatus(readyStatus());
      updateFilterEnabled();
      break;
    case MSG.PROGRESS:
      if (onFilterProgress) onFilterProgress(m);
      break;
    case MSG.OUTPUT_CHUNK_BATCH:
      if (onOutputBatch) onOutputBatch(m);
      break;
    case MSG.ERROR:
      if (onWorkerError) onWorkerError(m.message);
      else {
        clearIndexLoading();
        indexInfo.textContent = "";
        setStatus("Error: " + m.message, "error");
      }
      break;
    default:
      break;
  }
};

worker.postMessage({ type: MSG.INIT });

function filterSingleStream(file, params, hooks) {
  return new ReadableStream({
    start() {
      onFilterProgress = hooks.onProgress;
      worker.postMessage({ type: MSG.FILTER, data: { file, ...params } });
    },
    pull(controller) {
      return new Promise((resolve, reject) => {
        onWorkerError = (message) => {
          onWorkerError = null;
          onOutputBatch = null;
          reject(new Error(message));
        };
        onOutputBatch = (m) => {
          onOutputBatch = null;
          for (const buf of m.chunks) controller.enqueue(new Uint8Array(buf));
          if (m.done) {
            onWorkerError = null;
            onFilterProgress = null;
            if (hooks.onStats) hooks.onStats(m.stats, m.elapsed);
            controller.close();
          }
          resolve();
        };
        worker.postMessage({ type: MSG.PULL });
      });
    },
    cancel() {
      onOutputBatch = null;
      onFilterProgress = null;
      onWorkerError = null;
    },
  });
}

function filterGroup(group, params, hooks) {
  if (group.kind === "paired") {
    return filterPairedStreams(group, params, hooks);
  }
  return { stream: filterSingleStream(group.file, params, hooks) };
}

function filterPairedStreams(group, params, hooks) {
  const state = {
    queues: { r1: [], r2: [] },
    done: false,
    pulling: null,
    started: false,
  };

  const clearHandlers = () => {
    onWorkerError = null;
    onOutputBatch = null;
  };

  const pump = () => {
    if (state.done || state.pulling) return state.pulling || Promise.resolve();
    state.pulling = new Promise((resolve, reject) => {
      onWorkerError = (message) => {
        clearHandlers();
        reject(new Error(message));
      };
      onOutputBatch = (m) => {
        clearHandlers();
        for (const buf of m.chunksR1 || []) state.queues.r1.push(new Uint8Array(buf));
        for (const buf of m.chunksR2 || []) state.queues.r2.push(new Uint8Array(buf));
        if (m.done) {
          state.done = true;
          onFilterProgress = null;
          if (hooks.onStats) hooks.onStats(m.stats, m.elapsed);
        }
        resolve();
      };
      worker.postMessage({ type: MSG.PULL });
    }).finally(() => {
      state.pulling = null;
    });
    return state.pulling;
  };

  const makeStream = (mate) =>
    new ReadableStream({
      start() {
        if (!state.started) {
          state.started = true;
          onFilterProgress = hooks.onProgress;
          worker.postMessage({
            type: MSG.FILTER,
            data: { file1: group.file1, file2: group.file2, ...params },
          });
        }
      },
      async pull(controller) {
        while (!state.queues[mate].length && !state.done) {
          await pump();
        }
        const chunk = state.queues[mate].shift();
        if (chunk) {
          controller.enqueue(chunk);
        } else if (state.done) {
          controller.close();
        }
      },
      cancel() {
        clearHandlers();
        onFilterProgress = null;
      },
    });

  return { streamR1: makeStream("r1"), streamR2: makeStream("r2") };
}

indexStatus.addEventListener("click", () => {
  if (isFiltering) return;
  if (indexNeedsDownload) downloadAndLoadIndex();
  else if (indexFailed) indexInput.click();
});
indexStatus.addEventListener("dragover", (e) => {
  if (!(indexNeedsDownload || indexFailed)) return;
  e.preventDefault();
  if (!isFiltering) indexStatus.classList.add("dragover");
});
indexStatus.addEventListener("dragleave", () => indexStatus.classList.remove("dragover"));
indexStatus.addEventListener("drop", (e) => {
  if (!(indexNeedsDownload || indexFailed)) return;
  e.preventDefault();
  indexStatus.classList.remove("dragover");
  if (isFiltering) return;
  if (e.dataTransfer.files.length > 0) loadIndexFile(e.dataTransfer.files[0]);
});
indexInput.addEventListener("change", () => {
  if (indexInput.files.length > 0) loadIndexFile(indexInput.files[0]);
});

async function loadIndexFile(file) {
  if (file.size > MAX_INDEX_BYTES) {
    setStatus(
      `Index file too large for browser use (${(file.size / 1024 / 1024).toFixed(0)}MB). ` +
        `Maximum supported size is ${MAX_INDEX_BYTES / 1024 / 1024 / 1024}GB`,
      "error"
    );
    return;
  }
  indexLoaded = false;
  indexFailed = false;
  indexNeedsDownload = false;
  indexFilename = file.name;
  indexStatus.classList.remove("loaded", "failed", "needs-download", "dragover");
  indexRetry.hidden = true;
  indexDownload.hidden = true;
  updateFilterEnabled();
  indexName.textContent = file.name;
  showIndexLoading("Loading");
  setStatus(`Loading index (${(file.size / 1024 / 1024).toFixed(0)}MB)…`);
  const buf = await file.arrayBuffer();
  sendIndexToWorker(buf);
}

const IDB_NAME = "deacon-wasm";
const IDB_STORE = "index";

function openIdxDB() {
  return new Promise((resolve, reject) => {
    const req = indexedDB.open(IDB_NAME, 1);
    req.onupgradeneeded = () => req.result.createObjectStore(IDB_STORE);
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
}

function idbRequest(makeReq) {
  return new Promise((resolve, reject) => {
    const req = makeReq();
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
}

async function getCachedIndex(url) {
  const db = await openIdxDB();
  try {
    return await idbRequest(() =>
      db.transaction(IDB_STORE, "readonly").objectStore(IDB_STORE).get(url)
    );
  } finally {
    db.close();
  }
}

async function putCachedIndex(url, blob) {
  if (navigator.storage?.persist) {
    try {
      await navigator.storage.persist();
    } catch (_) {}
  }
  const db = await openIdxDB();
  try {
    await idbRequest(() =>
      db.transaction(IDB_STORE, "readwrite").objectStore(IDB_STORE).put(blob, url)
    );
  } finally {
    db.close();
  }
}

async function clearCachedIndex(url) {
  const db = await openIdxDB();
  try {
    await idbRequest(() =>
      db.transaction(IDB_STORE, "readwrite").objectStore(IDB_STORE).delete(url)
    );
  } finally {
    db.close();
  }
}

async function pruneStaleIndexes(keepUrl) {
  const db = await openIdxDB();
  try {
    const store = db.transaction(IDB_STORE, "readwrite").objectStore(IDB_STORE);
    const keys = await idbRequest(() => store.getAllKeys());
    await Promise.all(
      keys
        .filter((key) => key !== keepUrl)
        .map((key) => idbRequest(() => store.delete(key)))
    );
  } finally {
    db.close();
  }
}

function showIndexLoading(label) {
  indexInfo.textContent = label;
  indexInfo.classList.add("loading-dots");
}

function clearIndexLoading() {
  indexInfo.classList.remove("loading-dots");
}

function sendIndexToWorker(arrayBuffer) {
  if (workerReady) {
    worker.postMessage({ type: MSG.LOAD_INDEX, data: arrayBuffer }, [arrayBuffer]);
  } else {
    pendingIndexBuffer = arrayBuffer;
  }
}

async function downloadIndex(url, onProgress) {
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`HTTP ${resp.status} ${resp.statusText}`);
  const total = Number(resp.headers.get("content-length")) || 0;
  const reader = resp.body.getReader();
  const chunks = [];
  let received = 0;
  for (;;) {
    const { value, done } = await reader.read();
    if (done) break;
    chunks.push(value);
    received += value.length;
    onProgress(received, total);
  }
  return new Blob(chunks);
}

function showIndexProgress(received, total) {
  clearIndexLoading();
  indexProgress.style.display = "";
  if (total) {
    const pct = (received / total) * 100;
    indexProgress.removeAttribute("indeterminate");
    indexProgress.value = pct;
    indexInfo.textContent = `${humanBytes(received)} / ${humanBytes(total)} (${pct.toFixed(0)}%)`;
  } else {
    indexProgress.removeAttribute("value");
    indexInfo.textContent = `${humanBytes(received)} downloaded`;
  }
}

function showIndexNeedsDownload() {
  indexNeedsDownload = true;
  indexFailed = false;
  clearIndexLoading();
  indexProgress.style.display = "none";
  indexRetry.hidden = true;
  indexDownload.hidden = false;
  indexDownload.classList.remove("pulse-glow");
  void indexDownload.offsetWidth;
  indexDownload.classList.add("pulse-glow");
  indexName.textContent = "";
  indexInfo.textContent =
    "No Deacon index cached. Click to download or drop an .idx/.pidx file.";
  indexStatus.classList.remove("loaded", "failed");
  indexStatus.classList.add("needs-download");
  setStatus("");
  updateFilterEnabled();
}

function showIndexFailure(err) {
  console.error(err);
  indexFailed = true;
  indexNeedsDownload = false;
  clearIndexLoading();
  indexProgress.style.display = "none";
  indexName.textContent = INDEX_DISPLAY_NAME;
  indexInfo.textContent =
    `Download failed (${err?.message || err}). ` +
    "Click to select a local .idx/.pidx file, drag one here, or retry.";
  indexRetry.hidden = false;
  indexDownload.hidden = true;
  indexStatus.classList.remove("loaded", "needs-download");
  indexStatus.classList.add("failed");
  setStatus("");
  updateFilterEnabled();
}

async function initIndex() {
  indexFilename = INDEX_FILENAME;
  indexFailed = false;
  indexNeedsDownload = false;
  indexName.textContent = INDEX_DISPLAY_NAME;
  indexRetry.hidden = true;
  indexDownload.hidden = true;
  indexStatus.classList.remove("failed", "needs-download", "loaded", "dragover");
  pruneStaleIndexes(INDEX_URL).catch((e) =>
    console.warn("deacon-wasm: failed to prune stale cached indexes", e)
  );
  try {
    const cached = await getCachedIndex(INDEX_URL);
    if (cached) {
      showIndexLoading("Loading cached index");
      setStatus("Loading cached index…");
      sendIndexToWorker(await cached.arrayBuffer());
      return;
    }
    showIndexNeedsDownload();
  } catch (err) {
    showIndexFailure(err);
  }
}

async function downloadAndLoadIndex() {
  indexNeedsDownload = false;
  indexFailed = false;
  indexFilename = INDEX_FILENAME;
  indexName.textContent = INDEX_DISPLAY_NAME;
  indexRetry.hidden = true;
  indexDownload.hidden = true;
  indexStatus.classList.remove("needs-download", "failed", "dragover");
  updateFilterEnabled();
  try {
    indexInfo.textContent = "Downloading index…";
    setStatus("Downloading index…");
    showIndexProgress(0, 0);
    const blob = await downloadIndex(INDEX_URL, showIndexProgress);
    await putCachedIndex(INDEX_URL, blob).catch((e) =>
      console.warn("deacon-wasm: failed to cache index", e)
    );
    indexProgress.style.display = "none";
    showIndexLoading("Loading index");
    setStatus("Loading index…");
    sendIndexToWorker(await blob.arrayBuffer());
  } catch (err) {
    showIndexFailure(err);
  }
}

indexDownload.addEventListener("click", (e) => {
  e.stopPropagation();
  downloadAndLoadIndex();
});

indexRetry.addEventListener("click", (e) => {
  e.stopPropagation();
  downloadAndLoadIndex();
});

purgeIndexEl.addEventListener("click", async (e) => {
  e.preventDefault();
  if (isFiltering) return;
  indexLoaded = false;
  worker.postMessage({ type: MSG.RESET });
  updateFilterEnabled();
  showIndexNeedsDownload();
  setStatus("Clearing index cache…");
  try {
    await clearCachedIndex(INDEX_URL);
    setStatus("Index cache cleared");
  } catch (err) {
    console.warn("deacon-wasm: failed to clear cached index", err);
    setStatus("Couldn't clear the index cache: " + (err?.message || err), "error");
  }
});

dirZone.addEventListener("click", () => {
  if (!isFiltering) dirInput.click();
});

dirInput.addEventListener("change", () => {
  setSelection(Array.from(dirInput.files));
});

dirZone.addEventListener("dragover", (e) => {
  e.preventDefault();
  if (!isFiltering) dirZone.classList.add("dragover");
});
dirZone.addEventListener("dragleave", () => dirZone.classList.remove("dragover"));
dirZone.addEventListener("drop", async (e) => {
  e.preventDefault();
  dirZone.classList.remove("dragover");
  if (isFiltering) return;
  const items = Array.from(e.dataTransfer.items)
    .map((it) => (it.webkitGetAsEntry ? it.webkitGetAsEntry() : null))
    .filter(Boolean);
  if (items.length) {
    const files = await collectEntries(items);
    setSelection(files);
  } else {
    setSelection(Array.from(e.dataTransfer.files || []));
  }
});

async function collectEntries(entries, prefix = "") {
  const out = [];
  for (const entry of entries) {
    if (entry.isFile) {
      const file = await new Promise((res, rej) => entry.file(res, rej));
      try {
        Object.defineProperty(file, "webkitRelativePath", {
          value: prefix + entry.name,
          configurable: true,
        });
      } catch (_) {
      }
      out.push(file);
    } else if (entry.isDirectory) {
      const reader = entry.createReader();
      const children = await readAllEntries(reader);
      out.push(...(await collectEntries(children, prefix + entry.name + "/")));
    }
  }
  return out;
}

function readAllEntries(reader) {
  return new Promise((resolve, reject) => {
    const all = [];
    const next = () =>
      reader.readEntries((batch) => {
        if (!batch.length) resolve(all);
        else {
          all.push(...batch);
          next();
        }
      }, reject);
    next();
  });
}

resetBtn.addEventListener("click", () => {
  if (isFiltering) return;
  dirInput.value = "";
  setSelection([]);
  dirSummary.textContent = "";
  dirZone.classList.remove("loaded");
});

filterBtn.addEventListener("click", filterAll);

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

async function filterAll() {
  const groups = selectedFiles;
  if (!groups.length) return;
  const params = filterParams();

  console.log("deacon-wasm: filter params", {
    index: indexFilename,
    deplete: params.deplete,
    abs_threshold: params.absThreshold,
    rel_threshold: params.relThreshold,
    prefix_length: FILTER_DEFAULTS.prefixLength,
  });

  isFiltering = true;
  updateFilterEnabled();
  resetBtn.disabled = true;
  setStatus("Filtering…");

  progressWrap.style.display = "block";
  dehostProgress.value = 0;
  dehostLabel.textContent = "";
  fileRows.forEach(({ li, st }) => {
    li.className = "";
    st.textContent = "queued";
  });
  const rows = fileRows;

  let filteredBefore = 0;
  let totalBasesIn = 0;
  let completed = 0;
  try {
    for (let i = 0; i < groups.length; i++) {
      const group = groups[i];
      rows[i].li.className = "active";
      rows[i].st.textContent = "filtering…";

      let curProcessed = 0;
      let fileStats = null;

      const renderRow = () => {
        const dpct = group.size ? Math.min(100, (curProcessed / group.size) * 100) : 100;
        rows[i].st.textContent = `filtered ${dpct.toFixed(0)}%`;
      };
      const onProgress = (m) => {
        curProcessed = m.bytesProcessed || 0;
        const pct = totalBytes ? ((filteredBefore + curProcessed) / totalBytes) * 100 : 0;
        dehostProgress.value = pct;
        dehostLabel.textContent =
          `Filtering ${i + 1} of ${groups.length}: ` +
          `${humanBytes(filteredBefore + curProcessed)} / ${humanBytes(totalBytes)} ` +
          `(${pct.toFixed(1)}%)`;
        renderRow();
      };
      const onStats = (stats) => {
        fileStats = stats;
      };

      const outputs = filterGroup(group, params, { onProgress, onStats });
      if (group.kind === "paired") {
        const [blob1, blob2] = await Promise.all([
          new Response(outputs.streamR1).blob(),
          new Response(outputs.streamR2).blob(),
        ]);
        downloadBlob(blob1, gzName(group.file1.name));
        downloadBlob(blob2, gzName(group.file2.name));
      } else {
        const blob = await new Response(outputs.stream).blob();
        downloadBlob(blob, gzName(group.file.name));
      }

      filteredBefore += group.size;
      totalBasesIn += Number(fileStats?.basesIn || 0);
      rows[i].li.className = "done";
      const readsIn = Number(fileStats?.readsIn || 0);
      const readsOut = Number(fileStats?.readsOut || 0);
      rows[i].st.textContent = readsIn
        ? `done. ${readsOut.toLocaleString()}/${readsIn.toLocaleString()} records kept`
        : "done";
      completed++;
    }

    dehostProgress.value = 100;
    const done = describeGroups(groups.slice(0, completed));
    dehostLabel.textContent = `Processed ${humanBases(totalBasesIn)} of input across ${done}.`;
    setStatus(`Filtering complete: ${done} downloaded`, "success");
  } catch (err) {
    const idx = completed;
    if (rows[idx]) {
      rows[idx].li.className = "failed";
      rows[idx].st.textContent = "failed";
    }
    console.error(err);
    setStatus("Filtering failed: " + (err?.message || err), "error");
  } finally {
    isFiltering = false;
    resetBtn.disabled = false;
    updateFilterEnabled();
  }
}

initIndex();
