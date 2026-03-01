// Web Worker for off-main-thread WASM filtering
let wasm = null;
let index = null;
const ASSET_VERSION = "20260301-5";
const STREAM_CHUNK_BYTES = 2 * 1024 * 1024; // 2 MiB (browser may choose different chunk sizes)

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

  const session = new wasm.WasmFilterSession(
    index, deplete, absThreshold, relThreshold,
    isGz,  // decompress_input
    isGz   // compress_output (match input format)
  );

  const reader = file.stream().getReader();
  const outputChunks = [];
  let processedBytes = 0;
  let lastProgressTs = 0;

  try {
    while (true) {
      const { value, done } = await reader.read();
      if (done) break;
      if (!value || value.length === 0) continue;

      processedBytes += value.length;
      const outChunk = session.push_chunk(value);
      if (outChunk && outChunk.length > 0) {
        outputChunks.push(outChunk);
      }

      const now = performance.now();
      if (now - lastProgressTs >= 150) {
        self.postMessage({
          type: "progress",
          bytesProcessed: processedBytes,
          bytesTotal: totalBytes,
          progressCompressed: isGz,
          chunkHint: STREAM_CHUNK_BYTES,
        });
        lastProgressTs = now;
      }
    }

    const outputCompressed = session.output_compressed();
    const tail = session.finish();
    if (tail && tail.length > 0) outputChunks.push(tail);

    const stats = session.stats();

    const blob = new Blob(outputChunks, { type: "application/octet-stream" });
    const output = new Uint8Array(await blob.arrayBuffer());

    return {
      output,
      stats,
      bytesProcessed: processedBytes,
      bytesTotal: totalBytes,
      progressCompressed: isGz,
      outputCompressed,
    };
  } finally {
    reader.releaseLock();
  }
}

self.onmessage = async function (e) {
  const { type, data } = e.data;

  if (type === "init") {
    try {
      const mod = await import(`./pkg/deacon_wasm.js?v=${ASSET_VERSION}`);
      const wasmUrl = new URL(`./pkg/deacon_wasm_bg.wasm?v=${ASSET_VERSION}`, import.meta.url);
      await mod.default({ module_or_path: wasmUrl });
      wasm = mod;
      self.postMessage({ type: "ready" });
    } catch (err) {
      self.postMessage({ type: "error", message: "Failed to initialize WASM: " + err.message });
    }
    return;
  }

  if (type === "load_index") {
    try {
      const bytes = new Uint8Array(data);
      index = new wasm.WasmIndex(bytes);
      const info = index.info();
      self.postMessage({ type: "index_loaded", info });
    } catch (err) {
      self.postMessage({ type: "error", message: "Failed to load index: " + err.message });
    }
    return;
  }

  if (type === "reset") {
    index = null;
    self.postMessage({ type: "reset_done" });
    return;
  }

  if (type === "filter") {
    try {
      const { file, deplete, absThreshold, relThreshold } = data;
      const t0 = performance.now();
      const result = await streamFilterFile(file, { deplete, absThreshold, relThreshold });

      const {
        output,
        stats,
        bytesProcessed,
        bytesTotal,
        progressCompressed,
        outputCompressed,
      } = result;
      const elapsed = ((performance.now() - t0) / 1000).toFixed(2);
      self.postMessage(
        {
          type: "filtered",
          result: output.buffer,
          elapsed,
          stats,
          bytesProcessed,
          bytesTotal,
          progressCompressed,
          outputCompressed,
        },
        [output.buffer]
      );
    } catch (err) {
      self.postMessage({ type: "error", message: "Filtering failed: " + err.message });
    }
    return;
  }
};
