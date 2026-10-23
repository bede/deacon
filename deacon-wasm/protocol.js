export const MSG = Object.freeze({
  INIT: "init",
  READY: "ready",
  ERROR: "error",
  LOAD_INDEX: "load_index",
  INDEX_LOADED: "index_loaded",
  RESET: "reset",
  RESET_DONE: "reset_done",
  FILTER: "filter",
  PULL: "pull",
  PROGRESS: "progress",
  STAGE: "stage",
  OUTPUT_CHUNK_BATCH: "output_chunk_batch",
});

export const STAGE = Object.freeze({
  FINALIZING: "finalizing",
});

// Deacon CLI filter defaults (-a 2 -r 0.01, deplete); shared by worker and main thread
export const FILTER_DEFAULTS = Object.freeze({
  deplete: true,
  absThreshold: 2,
  relThreshold: 0.01,
  prefixLength: 0,
  rename: false,
  discardQuality: false,
});
