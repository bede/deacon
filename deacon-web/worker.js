// Web Worker for off-main-thread WASM filtering
let wasm = null;

self.onmessage = async function (e) {
  const { type, data } = e.data;

  if (type === "init") {
    try {
      const { default: init, load_index, filter_fastq, get_index_info } = await import("./pkg/deacon_web.js");
      await init();
      wasm = { load_index, filter_fastq, get_index_info };
      self.postMessage({ type: "ready" });
    } catch (err) {
      self.postMessage({ type: "error", message: "Failed to initialize WASM: " + err.message });
    }
    return;
  }

  if (type === "load_index") {
    try {
      const bytes = new Uint8Array(data);
      wasm.load_index(bytes);
      const info = wasm.get_index_info();
      self.postMessage({ type: "index_loaded", info });
    } catch (err) {
      self.postMessage({ type: "error", message: "Failed to load index: " + err.message });
    }
    return;
  }

  if (type === "filter") {
    try {
      const { input, filename, deplete, absThreshold, relThreshold } = data;
      const bytes = new Uint8Array(input);
      const isGzipped = filename.endsWith(".gz");
      const t0 = performance.now();
      const result = wasm.filter_fastq(bytes, isGzipped, deplete, absThreshold, relThreshold);
      const elapsed = ((performance.now() - t0) / 1000).toFixed(2);
      self.postMessage(
        { type: "filtered", result: result.buffer, elapsed },
        [result.buffer]
      );
    } catch (err) {
      self.postMessage({ type: "error", message: "Filtering failed: " + err.message });
    }
    return;
  }
};
