import test from "node:test";
import assert from "node:assert/strict";

import { MSG, STAGE } from "./protocol.js";

test("protocol message types are stable and unique", () => {
  assert.deepEqual(MSG, {
    INIT: "init",
    READY: "ready",
    ERROR: "error",
    LOAD_INDEX: "load_index",
    INDEX_LOADED: "index_loaded",
    RESET: "reset",
    RESET_DONE: "reset_done",
    FILTER: "filter",
    PROGRESS: "progress",
    STAGE: "stage",
    OUTPUT_CHUNK_BATCH: "output_chunk_batch",
    FILTERED_DONE: "filtered_done",
  });

  const values = Object.values(MSG);
  assert.equal(new Set(values).size, values.length);
});

test("stage constants are stable", () => {
  assert.deepEqual(STAGE, {
    FINALIZING: "finalizing",
  });
});

