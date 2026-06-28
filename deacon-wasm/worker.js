import { MSG, STAGE } from "./protocol.js?v=20260302-3";

// Web Worker for off-main-thread WASM filtering
let wasm = null;
let index = null;
const ASSET_VERSION = "20260302-3";
const OUTPUT_BATCH_BYTES = 4 * 1024 * 1024; // 4 MiB

function isGzipFilename(name) {
  return /\.(fastq|fq|fasta|fa)\.gz$/i.test(name || "");
}

function isSeqFilename(name) {
  return /\.(fastq|fq|fasta|fa)(\.gz)?$/i.test(name || "");
}

async function streamFilterFile(file, opts) {
  const { deplete, absThreshold, relThreshold } = opts;

  if (!index) {
    throw new Error("No index loaded");
  }
  if (!file || typeof file.stream !== "function") {
    throw new Error("Missing sequence file");
  }
  if (!isSeqFilename(file.name)) {
    throw new Error("Streaming mode supports FASTA/FASTQ files only (.fasta/.fa/.fastq/.fq, optionally .gz)");
  }

  const isGz = isGzipFilename(file.name);
  const totalBytes = file.size || 0;

  const session = new wasm.FilterSession(
    index, deplete, absThreshold, relThreshold,
    isGz,  // decompress_input
    isGz   // compress_output (match input format)
  );

  const reader = file.stream().getReader();
  let processedBytes = 0;
  let lastProgressTs = 0;
  let pendingOutputBuffers = [];
  let pendingOutputBytes = 0;

  const flushPendingOutput = (force = false) => {
    if (!force && pendingOutputBytes < OUTPUT_BATCH_BYTES) return;
    if (pendingOutputBuffers.length === 0) return;
    const transferList = pendingOutputBuffers;
    self.postMessage(
      {
        type: MSG.OUTPUT_CHUNK_BATCH,
        chunks: transferList,
      },
      transferList
    );
    pendingOutputBuffers = [];
    pendingOutputBytes = 0;
  };

  const postOutputChunk = (chunk) => {
    if (!chunk || chunk.length === 0) return;
    const transfer =
      chunk.byteOffset === 0 && chunk.byteLength === chunk.buffer.byteLength
        ? chunk.buffer
        : chunk.slice().buffer;
    pendingOutputBuffers.push(transfer);
    pendingOutputBytes += transfer.byteLength;
    flushPendingOutput(false);
  };

  try {
    while (true) {
      const { value, done } = await reader.read();
      if (done) break;
      if (!value || value.length === 0) continue;

      processedBytes += value.length;
      const outChunk = session.push_chunk(value);
      postOutputChunk(outChunk);

      const now = performance.now();
      if (now - lastProgressTs >= 150) {
        self.postMessage({
          type: MSG.PROGRESS,
          bytesProcessed: processedBytes,
          bytesTotal: totalBytes,
          progressCompressed: isGz,
        });
        lastProgressTs = now;
      }
    }

    self.postMessage({ type: MSG.STAGE, stage: STAGE.FINALIZING });
    const tail = session.finish();
    postOutputChunk(tail);
    flushPendingOutput(true);

    const stats = session.stats();

    return {
      stats,
      bytesProcessed: processedBytes,
      bytesTotal: totalBytes,
      progressCompressed: isGz,
    };
  } finally {
    reader.releaseLock();
  }
}

async function streamFilterPairedFiles(file1, file2, opts) {
  const { deplete, absThreshold, relThreshold } = opts;

  if (!index) {
    throw new Error("No index loaded");
  }
  for (const file of [file1, file2]) {
    if (!file || typeof file.stream !== "function") {
      throw new Error("Missing paired sequence file");
    }
    if (!isSeqFilename(file.name)) {
      throw new Error("Streaming mode supports FASTA/FASTQ files only (.fasta/.fa/.fastq/.fq, optionally .gz)");
    }
  }

  const r1Gz = isGzipFilename(file1.name);
  const r2Gz = isGzipFilename(file2.name);
  const totalBytes = (file1.size || 0) + (file2.size || 0);
  const session = new wasm.PairedFilterSession(
    index, deplete, absThreshold, relThreshold,
    r1Gz, r2Gz,
    r1Gz, r2Gz
  );

  const readerR1 = file1.stream().getReader();
  const readerR2 = file2.stream().getReader();
  let doneR1 = false;
  let doneR2 = false;
  let processedBytes = 0;
  let lastProgressTs = 0;
  let pendingR1 = [];
  let pendingR2 = [];
  let pendingBytes = 0;

  const appendOutput = (out) => {
    for (const [chunk, pending] of [[out.r1, pendingR1], [out.r2, pendingR2]]) {
      if (!chunk || chunk.length === 0) continue;
      const transfer =
        chunk.byteOffset === 0 && chunk.byteLength === chunk.buffer.byteLength
          ? chunk.buffer
          : chunk.slice().buffer;
      pending.push(transfer);
      pendingBytes += transfer.byteLength;
    }
  };

  const flushPendingOutput = (force = false) => {
    if (!force && pendingBytes < OUTPUT_BATCH_BYTES) return;
    if (pendingR1.length === 0 && pendingR2.length === 0) return;
    const transferList = [...pendingR1, ...pendingR2];
    self.postMessage(
      {
        type: MSG.OUTPUT_CHUNK_BATCH,
        chunksR1: pendingR1,
        chunksR2: pendingR2,
      },
      transferList
    );
    pendingR1 = [];
    pendingR2 = [];
    pendingBytes = 0;
  };

  const postProgress = () => {
    self.postMessage({
      type: MSG.PROGRESS,
      bytesProcessed: processedBytes,
      bytesTotal: totalBytes,
      progressCompressed: r1Gz || r2Gz,
    });
  };

  try {
    while (!doneR1 || !doneR2) {
      if (!doneR1) {
        const { value, done } = await readerR1.read();
        if (done) {
          doneR1 = true;
          appendOutput(session.finish_r1());
        } else if (value && value.length > 0) {
          processedBytes += value.length;
          appendOutput(session.push_r1(value));
        }
      }

      if (!doneR2) {
        const { value, done } = await readerR2.read();
        if (done) {
          doneR2 = true;
          appendOutput(session.finish_r2());
        } else if (value && value.length > 0) {
          processedBytes += value.length;
          appendOutput(session.push_r2(value));
        }
      }

      const now = performance.now();
      if (now - lastProgressTs >= 150) {
        postProgress();
        lastProgressTs = now;
      }
      flushPendingOutput(false);
    }

    self.postMessage({ type: MSG.STAGE, stage: STAGE.FINALIZING });
    flushPendingOutput(true);
    return {
      stats: session.stats(),
      bytesProcessed: processedBytes,
      bytesTotal: totalBytes,
      progressCompressed: r1Gz || r2Gz,
    };
  } finally {
    readerR1.releaseLock();
    readerR2.releaseLock();
    session.free();
  }
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
      self.postMessage({ type: MSG.ERROR, message: "Failed to initialize WASM: " + err.message });
    }
    return;
  }

  if (type === MSG.LOAD_INDEX) {
    try {
      const bytes = new Uint8Array(data);
      index = new wasm.WasmIndex(bytes);
      const info = index.info();
      self.postMessage({ type: MSG.INDEX_LOADED, info });
    } catch (err) {
      self.postMessage({ type: MSG.ERROR, message: "Failed to load index: " + err.message });
    }
    return;
  }

  if (type === MSG.RESET) {
    index = null;
    self.postMessage({ type: MSG.RESET_DONE });
    return;
  }

  if (type === MSG.FILTER) {
    try {
      const { file, file1, file2, deplete, absThreshold, relThreshold } = data;
      const t0 = performance.now();
      const result = file1 && file2
        ? await streamFilterPairedFiles(file1, file2, { deplete, absThreshold, relThreshold })
        : await streamFilterFile(file, { deplete, absThreshold, relThreshold });

      const {
        stats,
        bytesProcessed,
        bytesTotal,
        progressCompressed,
      } = result;
      const elapsed = ((performance.now() - t0) / 1000).toFixed(2);
      self.postMessage({
        type: MSG.FILTERED_DONE,
        elapsed,
        stats,
        bytesProcessed,
        bytesTotal,
        progressCompressed,
      });
    } catch (err) {
      self.postMessage({ type: MSG.ERROR, message: "Filtering failed: " + err.message });
    }
    return;
  }
};
