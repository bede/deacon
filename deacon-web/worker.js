// Web Worker for off-main-thread WASM filtering
let wasm = null;
let index = null;
const ASSET_VERSION = "20260301-3";
const STREAM_CHUNK_BYTES = 2 * 1024 * 1024; // 2 MiB (browser may choose different chunk sizes)
const MAX_COMPRESSED_OUTPUT_STREAM_BYTES = 512 * 1024 * 1024; // Avoid unstable browser compression behavior on very large jobs.

function isGzipFilename(name) {
  return /\.(fastq|fq)\.gz$/i.test(name || "");
}

function isFastqFilename(name) {
  return /\.(fastq|fq)(\.gz)?$/i.test(name || "");
}

async function runLegacyFilter(file, opts) {
  const { deplete, absThreshold, relThreshold } = opts;
  const buf = await file.arrayBuffer();
  const bytes = new Uint8Array(buf);
  const { output, stats } = wasm.filter_with_stats(index, bytes, deplete, absThreshold, relThreshold);
  return {
    output,
    stats,
    bytesProcessed: bytes.length,
    bytesTotal: bytes.length,
    progressCompressed: true,
    outputCompressed: true,
    mode: "legacy",
  };
}

async function streamFilterFile(file, opts) {
  const { deplete, absThreshold, relThreshold } = opts;

  if (!index) {
    throw new Error("No index loaded");
  }
  if (!file || typeof file.stream !== "function") {
    throw new Error("Missing sequence file");
  }
  if (!isFastqFilename(file.name)) {
    return runLegacyFilter(file, opts);
  }

  let inputStream = file.stream();
  let totalBytes = file.size || 0;
  let progressCompressed = true;

  if (isGzipFilename(file.name)) {
    if (typeof DecompressionStream === "undefined") {
      return runLegacyFilter(file, opts);
    }
    inputStream = inputStream.pipeThrough(new DecompressionStream("gzip"));
    // Decompressed length is not known ahead of time.
    totalBytes = 0;
    progressCompressed = false;
  }

  let outputWriter = null;
  let outputBytesPromise = null;
  const plainChunks = [];
  let outputCompressed = false;
  let outputWriterClosed = false;

  const shouldCompressOutput =
    typeof CompressionStream !== "undefined" &&
    totalBytes > 0 &&
    totalBytes <= MAX_COMPRESSED_OUTPUT_STREAM_BYTES;

  if (shouldCompressOutput) {
    const compression = new CompressionStream("gzip");
    outputWriter = compression.writable.getWriter();
    outputBytesPromise = new Response(compression.readable).arrayBuffer();
    outputCompressed = true;
  }

  const session = new wasm.WasmFilterSession(index, deplete, absThreshold, relThreshold);
  const reader = inputStream.getReader();

  let processedBytes = 0;
  let lastProgressTs = 0;

  try {
    while (true) {
      const { value, done } = await reader.read();
      if (done) {
        break;
      }
      if (!value || value.length === 0) {
        continue;
      }

      processedBytes += value.length;
      const outChunk = session.push_chunk(value);
      if (outChunk && outChunk.length > 0) {
        if (outputWriter) {
          await outputWriter.write(outChunk);
        } else {
          plainChunks.push(outChunk);
        }
      }

      const now = performance.now();
      if (now - lastProgressTs >= 150) {
        self.postMessage({
          type: "progress",
          bytesProcessed: processedBytes,
          bytesTotal: totalBytes,
          progressCompressed,
          chunkHint: STREAM_CHUNK_BYTES,
        });
        lastProgressTs = now;
      }
    }

    const tail = session.finish();
    if (tail && tail.length > 0) {
      if (outputWriter) {
        await outputWriter.write(tail);
      } else {
        plainChunks.push(tail);
      }
    }

    const stats = session.stats();

    let output;
    if (outputWriter) {
      await outputWriter.close();
      outputWriterClosed = true;
      const compressed = await outputBytesPromise;
      output = new Uint8Array(compressed);
    } else {
      const blob = new Blob(plainChunks, { type: "application/octet-stream" });
      output = new Uint8Array(await blob.arrayBuffer());
    }

    return {
      output,
      stats,
      bytesProcessed: processedBytes,
      bytesTotal: totalBytes,
      progressCompressed,
      outputCompressed,
      mode: "stream",
    };
  } finally {
    reader.releaseLock();
    if (outputWriter) {
      try {
        if (!outputWriterClosed) {
          await outputWriter.close();
          outputWriterClosed = true;
        }
      } catch (_) {
        // Ignore close failures after errors/cancel.
      }
    }
  }
}

self.onmessage = async function (e) {
  const { type, data } = e.data;

  if (type === "init") {
    try {
      const mod = await import(`./pkg/deacon_web.js?v=${ASSET_VERSION}`);
      const wasmUrl = new URL(`./pkg/deacon_web_bg.wasm?v=${ASSET_VERSION}`, import.meta.url);
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
      const {
        output,
        stats,
        bytesProcessed,
        bytesTotal,
        progressCompressed,
        outputCompressed,
      } = await streamFilterFile(
        file,
        { deplete, absThreshold, relThreshold }
      );
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
