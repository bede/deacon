import { MSG, FILTER_DEFAULTS } from "./protocol.js?v=20261023-1";

// Off-main-thread WASM dehosting; pull-driven (one batch per MSG.PULL) so memory stays bounded
let wasm = null;
let index = null;
const ASSET_VERSION = "20261023-1";
const OUTPUT_BATCH_BYTES = 4 * 1024 * 1024; // 4 MiB

// wasm-bindgen throws Result errors as plain strings with no .message; normalise so the reason survives
function errText(err) {
  if (err == null) return "Unknown error";
  if (typeof err === "string") return err;
  return err.message || String(err);
}

// Per-file state, set up by MSG.FILTER and drained by MSG.PULL
let active = null;

function isGzipFilename(name) {
  return /\.(fastq|fq|fasta|fa)\.gz$/i.test(name || "");
}

function isSeqFilename(name) {
  return /\.(fastq|fq|fasta|fa)(\.gz)?$/i.test(name || "");
}

// Per-session threshold overrides from the FILTER payload; fall back to FILTER_DEFAULTS
function resolveParams(data) {
  return {
    deplete: data.deplete ?? FILTER_DEFAULTS.deplete,
    absThreshold: data.absThreshold ?? FILTER_DEFAULTS.absThreshold,
    relThreshold: data.relThreshold ?? FILTER_DEFAULTS.relThreshold,
  };
}

function startFilter(data) {
  const params = resolveParams(data);
  if (data.file1 && data.file2) {
    return startPairedFilter(data.file1, data.file2, params);
  }
  return startSingleFilter(data.file, params);
}

function startSingleFilter(file, params) {
  if (!index) throw new Error("No index loaded");
  if (!file || typeof file.stream !== "function") throw new Error("Missing sequence file");
  if (!isSeqFilename(file.name)) {
    throw new Error("Only FASTA/FASTQ files are supported (.fasta/.fa/.fastq/.fq, optionally .gz)");
  }

  const isGz = isGzipFilename(file.name);
  const session = new wasm.FilterSession(
    index,
    params.deplete,
    params.absThreshold,
    params.relThreshold,
    isGz, // decompress_input
    true, // compress_output (always gzip; uncompressed inputs get a .gz name in app.js)
    FILTER_DEFAULTS.rename,
    FILTER_DEFAULTS.discardQuality
  );

  active = {
    mode: "single",
    session,
    reader: file.stream().getReader(),
    isGz,
    totalBytes: file.size || 0,
    processedBytes: 0,
    pendingBuffers: [],
    pendingBytes: 0,
    finished: false,
    lastProgressTs: 0,
    t0: performance.now(),
  };
}

function startPairedFilter(file1, file2, params) {
  if (!index) throw new Error("No index loaded");
  for (const file of [file1, file2]) {
    if (!file || typeof file.stream !== "function") throw new Error("Missing paired sequence file");
    if (!isSeqFilename(file.name)) {
      throw new Error("Only FASTA/FASTQ files are supported (.fasta/.fa/.fastq/.fq, optionally .gz)");
    }
  }

  const r1Gz = isGzipFilename(file1.name);
  const r2Gz = isGzipFilename(file2.name);
  const session = new wasm.PairedFilterSession(
    index,
    params.deplete,
    params.absThreshold,
    params.relThreshold,
    r1Gz, // decompress_r1
    r2Gz, // decompress_r2
    true, // compress_r1 (always gzip; uncompressed inputs get a .gz name in app.js)
    true, // compress_r2 (always gzip; uncompressed inputs get a .gz name in app.js)
    FILTER_DEFAULTS.rename,
    FILTER_DEFAULTS.discardQuality
  );

  active = {
    mode: "paired",
    session,
    readerR1: file1.stream().getReader(),
    readerR2: file2.stream().getReader(),
    doneR1: false,
    doneR2: false,
    isGz: r1Gz || r2Gz,
    totalBytes: (file1.size || 0) + (file2.size || 0),
    processedBytes: 0,
    pendingBuffersR1: [],
    pendingBuffersR2: [],
    pendingBytes: 0,
    finished: false,
    lastProgressTs: 0,
    t0: performance.now(),
  };
}

// Append a chunk of filtered output, transferring its buffer
function appendOutput(chunk) {
  if (!chunk || chunk.length === 0) return;
  const buffer =
    chunk.byteOffset === 0 && chunk.byteLength === chunk.buffer.byteLength
      ? chunk.buffer
      : chunk.slice().buffer;
  active.pendingBuffers.push(buffer);
  active.pendingBytes += buffer.byteLength;
}

function appendPairedOutput(out) {
  for (const [chunk, key] of [[out.r1, "pendingBuffersR1"], [out.r2, "pendingBuffersR2"]]) {
    if (!chunk || chunk.length === 0) continue;
    const buffer =
      chunk.byteOffset === 0 && chunk.byteLength === chunk.buffer.byteLength
        ? chunk.buffer
        : chunk.slice().buffer;
    active[key].push(buffer);
    active.pendingBytes += buffer.byteLength;
  }
}

function postProgress() {
  // Running read counts for the live host tier; mid-stream stats() may be unsupported
  let readsIn = 0;
  let readsOut = 0;
  try {
    const s = active.session.stats();
    readsIn = Number(s?.readsIn || 0);
    readsOut = Number(s?.readsOut || 0);
  } catch (_) {}
  self.postMessage({
    type: MSG.PROGRESS,
    bytesProcessed: active.processedBytes,
    bytesTotal: active.totalBytes,
    progressCompressed: active.isGz,
    readsIn,
    readsOut,
  });
}

// Pump input until one batch is ready or input ends, then post OUTPUT_CHUNK_BATCH (done + stats at EOF)
async function pump() {
  const a = active;
  if (a.mode === "paired") return pumpPaired();

  while (a.pendingBytes < OUTPUT_BATCH_BYTES && !a.finished) {
    const { value, done } = await a.reader.read();
    if (done) {
      appendOutput(a.session.finish());
      a.finished = true;
      break;
    }
    if (!value || value.length === 0) continue;
    a.processedBytes += value.length;
    appendOutput(a.session.push_chunk(value));

    const now = performance.now();
    if (now - a.lastProgressTs >= 150) {
      postProgress();
      a.lastProgressTs = now;
    }
  }

  const chunks = a.pendingBuffers;
  a.pendingBuffers = [];
  a.pendingBytes = 0;

  if (a.finished) {
    const stats = a.session.stats();
    const elapsed = ((performance.now() - a.t0) / 1000).toFixed(2);
    postProgress();
    self.postMessage(
      {
        type: MSG.OUTPUT_CHUNK_BATCH,
        chunks,
        done: true,
        stats,
        elapsed,
        bytesProcessed: a.processedBytes,
        bytesTotal: a.totalBytes,
        progressCompressed: a.isGz,
      },
      chunks
    );
    a.session.free();
    active = null;
  } else {
    self.postMessage({ type: MSG.OUTPUT_CHUNK_BATCH, chunks, done: false }, chunks);
  }
}

async function pumpPaired() {
  const a = active;
  while (a.pendingBytes < OUTPUT_BATCH_BYTES && !a.finished) {
    if (!a.doneR1) {
      const { value, done } = await a.readerR1.read();
      if (done) {
        a.doneR1 = true;
        appendPairedOutput(a.session.finish_r1());
      } else if (value && value.length > 0) {
        a.processedBytes += value.length;
        appendPairedOutput(a.session.push_r1(value));
      }
    }

    if (!a.doneR2) {
      const { value, done } = await a.readerR2.read();
      if (done) {
        a.doneR2 = true;
        appendPairedOutput(a.session.finish_r2());
      } else if (value && value.length > 0) {
        a.processedBytes += value.length;
        appendPairedOutput(a.session.push_r2(value));
      }
    }

    if (a.doneR1 && a.doneR2) a.finished = true;

    const now = performance.now();
    if (now - a.lastProgressTs >= 150) {
      postProgress();
      a.lastProgressTs = now;
    }
  }

  const chunksR1 = a.pendingBuffersR1;
  const chunksR2 = a.pendingBuffersR2;
  a.pendingBuffersR1 = [];
  a.pendingBuffersR2 = [];
  a.pendingBytes = 0;

  if (a.finished) {
    const stats = a.session.stats();
    const elapsed = ((performance.now() - a.t0) / 1000).toFixed(2);
    postProgress();
    self.postMessage(
      {
        type: MSG.OUTPUT_CHUNK_BATCH,
        chunksR1,
        chunksR2,
        done: true,
        stats,
        elapsed,
        bytesProcessed: a.processedBytes,
        bytesTotal: a.totalBytes,
        progressCompressed: a.isGz,
      },
      [...chunksR1, ...chunksR2]
    );
    a.session.free();
    active = null;
  } else {
    self.postMessage({ type: MSG.OUTPUT_CHUNK_BATCH, chunksR1, chunksR2, done: false }, [
      ...chunksR1,
      ...chunksR2,
    ]);
  }
}

function disposeActive() {
  if (!active) return;
  try {
    if (active.mode === "paired") {
      active.readerR1.cancel();
      active.readerR2.cancel();
    } else {
      active.reader.cancel();
    }
  } catch (_) {}
  try {
    active.session.free();
  } catch (_) {}
  active = null;
}

self.onmessage = async function (e) {
  const { type, data } = e.data;

  if (type === MSG.INIT) {
    try {
      const mod = await import(`./pkg/deacon_wasm.js?v=${ASSET_VERSION}`);
      const wasmUrl = new URL(`./pkg/deacon_wasm_bg.wasm?v=${ASSET_VERSION}`, import.meta.url);
      await mod.default({ module_or_path: wasmUrl });
      wasm = mod;
      self.postMessage({ type: MSG.READY });
    } catch (err) {
      self.postMessage({ type: MSG.ERROR, message: "Failed to initialize WASM: " + errText(err) });
    }
    return;
  }

  if (type === MSG.LOAD_INDEX) {
    try {
      index = new wasm.WasmIndex(new Uint8Array(data));
      self.postMessage({ type: MSG.INDEX_LOADED, info: index.info() });
    } catch (err) {
      // wasm error already carries a "Failed to load index: …" prefix; surface as-is
      self.postMessage({ type: MSG.ERROR, message: errText(err) });
    }
    return;
  }

  if (type === MSG.RESET) {
    disposeActive();
    index = null;
    self.postMessage({ type: MSG.RESET_DONE });
    return;
  }

  if (type === MSG.FILTER) {
    try {
      disposeActive();
      startFilter(data);
    } catch (err) {
      self.postMessage({ type: MSG.ERROR, message: "Filtering failed: " + errText(err) });
    }
    return;
  }

  if (type === MSG.PULL) {
    if (!active) {
      self.postMessage({ type: MSG.ERROR, message: "No active filter session" });
      return;
    }
    try {
      await pump();
    } catch (err) {
      disposeActive();
      self.postMessage({ type: MSG.ERROR, message: "Filtering failed: " + errText(err) });
    }
    return;
  }
};
