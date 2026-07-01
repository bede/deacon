// Benchmarks the actual compiled WASM bundle under Node. See AGENTS.md.

import { readFile } from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";

function parseArgs(argv) {
  const args = new Map();
  for (let i = 0; i < argv.length; i += 2) {
    const key = argv[i];
    const value = argv[i + 1];
    if (!key?.startsWith("--") || value === undefined) {
      throw new Error(`Invalid arguments near: ${key ?? "<end>"}`);
    }
    args.set(key.slice(2), value);
  }
  return args;
}

function requireArg(args, name) {
  const value = args.get(name);
  if (!value) {
    throw new Error(`Missing required argument: --${name}`);
  }
  return value;
}

const args = parseArgs(process.argv.slice(2));
const pkgDir = path.resolve(requireArg(args, "pkg"));
const indexPath = path.resolve(requireArg(args, "index"));
const readsPath = path.resolve(requireArg(args, "reads"));
const chunkBytes = Number(args.get("chunk-kb") ?? 256) * 1024;
const iters = Number(args.get("iters") ?? 3);

const wasmModule = await import(
  pathToFileURL(path.join(pkgDir, "deacon_wasm.js")).href
);
if (typeof wasmModule.default === "function") {
  const wasmBytes = await readFile(path.join(pkgDir, "deacon_wasm_bg.wasm"));
  await wasmModule.default({ module_or_path: wasmBytes });
}

const indexBytes = await readFile(indexPath);
const index = new wasmModule.WasmIndex(new Uint8Array(indexBytes));

// Read the reads file into memory once so file I/O is excluded from timing.
const reads = new Uint8Array(await readFile(readsPath));
const decompressInput = readsPath.endsWith(".gz");
const totalMB = reads.length / (1024 * 1024);

console.log(`index=${index.info()}`);
console.log(
  `reads=${path.basename(readsPath)} bytes=${reads.length} gz=${decompressInput} ` +
    `chunk_kb=${chunkBytes / 1024} iters=${iters}`,
);

const secsPerIter = [];
for (let iter = 0; iter < iters; iter++) {
  const session = new wasmModule.FilterSession(
    index,
    false, // deplete
    2, // abs_threshold
    0.01, // rel_threshold
    decompressInput,
    false, // compress_output
    false, // rename
    false, // output_fasta
  );

  const start = process.hrtime.bigint();
  let bytesOut = 0;
  for (let off = 0; off < reads.length; off += chunkBytes) {
    const end = Math.min(off + chunkBytes, reads.length);
    bytesOut += session.push_chunk(reads.subarray(off, end)).length;
  }
  bytesOut += session.finish().length;
  const secs = Number(process.hrtime.bigint() - start) / 1e9;
  secsPerIter.push(secs);

  const stats = session.stats();
  session.free();
  console.log(
    `iter=${iter} secs=${secs.toFixed(3)} MB/s=${(totalMB / secs).toFixed(1)} ` +
      `reads_in=${stats.readsIn} reads_out=${stats.readsOut} bytes_out=${bytesOut}`,
  );
}

const meanSecs = secsPerIter.reduce((a, b) => a + b, 0) / secsPerIter.length;
const bestSecs = Math.min(...secsPerIter);
console.log(
  `RESULT mean_secs=${meanSecs.toFixed(3)} best_MBps=${(totalMB / bestSecs).toFixed(1)} ` +
    `chunk_kb=${chunkBytes / 1024}`,
);
