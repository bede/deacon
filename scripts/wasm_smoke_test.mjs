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
  return path.resolve(value);
}

function assertEqual(actual, expected, message) {
  if (actual !== expected) {
    throw new Error(`${message}: expected ${expected}, got ${actual}`);
  }
}

const args = parseArgs(process.argv.slice(2));
const pkgDir = requireArg(args, "pkg");
const indexPath = requireArg(args, "index");
const readsPath = requireArg(args, "reads");

const wasmModule = await import(pathToFileURL(path.join(pkgDir, "deacon_wasm.js")).href);

if (typeof wasmModule.default === "function") {
  const wasmBytes = await readFile(path.join(pkgDir, "deacon_wasm_bg.wasm"));
  await wasmModule.default({ module_or_path: wasmBytes });
}

const [indexBytes, readsBytes] = await Promise.all([
  readFile(indexPath),
  readFile(readsPath),
]);

const index = new wasmModule.WasmIndex(new Uint8Array(indexBytes));
const session = new wasmModule.FilterSession(
  index,
  false, // deplete
  2,     // abs_threshold
  0.01,  // rel_threshold
  false, // decompress_input
  false, // compress_output
  false, // rename
  false, // output_fasta
);

const outputChunks = [];
for (let offset = 0; offset < readsBytes.length; offset += 137) {
  const chunk = new Uint8Array(readsBytes.subarray(offset, offset + 137));
  const out = session.push_chunk(chunk);
  if (out.length > 0) {
    outputChunks.push(out);
  }
}

const tail = session.finish();
if (tail.length > 0) {
  outputChunks.push(tail);
}

const stats = session.stats();
const outputBytes = outputChunks.reduce((sum, chunk) => sum + chunk.length, 0);

assertEqual(stats.readsIn, 151, "readsIn");
assertEqual(stats.readsOut, 145, "readsOut");
assertEqual(stats.basesIn, 299126, "basesIn");
assertEqual(stats.basesOut, 297446, "basesOut");
assertEqual(outputBytes, 598267, "output byte count");

const decoder = new TextDecoder();
const firstChunk = outputChunks[0] ?? new Uint8Array();
const firstLine = decoder.decode(firstChunk).split("\n", 1)[0];
assertEqual(firstLine, "@NC_045512.2__1_1", "first output header");

console.log(
  [
    "WASM runtime smoke test passed",
    `index=${index.info()}`,
    `reads_in=${stats.readsIn}`,
    `reads_out=${stats.readsOut}`,
    `bytes_out=${outputBytes}`,
  ].join(" "),
);
