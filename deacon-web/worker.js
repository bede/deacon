// Web Worker for off-main-thread WASM filtering
let wasm = null;
let index = null;
const ASSET_VERSION = "20260227-1";

self.onmessage = async function (e) {
  const { type, data } = e.data;

  if (type === "init") {
    try {
      const mod = await import(`./pkg/deacon_web.js?v=${ASSET_VERSION}`);
      const wasmUrl = new URL(`./pkg/deacon_web_bg.wasm?v=${ASSET_VERSION}`, import.meta.url);
      await mod.default(wasmUrl);
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

  if (type === "filter") {
    try {
      const { input, deplete, absThreshold, relThreshold } = data;
      const bytes = new Uint8Array(input);
      const t0 = performance.now();
      const { output, stats } = wasm.filter_with_stats(index, bytes, deplete, absThreshold, relThreshold);
      const elapsed = ((performance.now() - t0) / 1000).toFixed(2);
      self.postMessage(
        { type: "filtered", result: output.buffer, elapsed, stats },
        [output.buffer]
      );
    } catch (err) {
      self.postMessage({ type: "error", message: "Filtering failed: " + err.message });
    }
    return;
  }
};
