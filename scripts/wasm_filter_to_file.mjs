import { createReadStream, createWriteStream } from "node:fs";
import { readFile } from "node:fs/promises";
import { once } from "node:events";
import path from "node:path";
import { pathToFileURL } from "node:url";

function parseArgs(argv) {
  const args = new Map();
  const flags = new Set();
  for (let i = 0; i < argv.length; i++) {
    const key = argv[i];
    if (!key?.startsWith("--")) {
      throw new Error(`Invalid arguments near: ${key ?? "<end>"}`);
    }
    const name = key.slice(2);
    if (name === "rename" || name === "fasta") {
      flags.add(name);
      continue;
    }
    const value = argv[++i];
    if (value === undefined) {
      throw new Error(`Invalid arguments near: ${key}`);
    }
    args.set(name, value);
  }
  args.flags = flags;
  return args;
}

function requireArg(args, name) {
  const value = args.get(name);
  if (!value) {
    throw new Error(`Missing required argument: --${name}`);
  }
  return path.resolve(value);
}

async function writeChunk(stream, chunk) {
  if (!chunk || chunk.length === 0) {
    return;
  }
  if (!stream.write(chunk)) {
    await once(stream, "drain");
  }
}

const args = parseArgs(process.argv.slice(2));
const pkgDir = requireArg(args, "pkg");
const indexPath = requireArg(args, "index");
const readsPath = requireArg(args, "reads");
const outputPath = requireArg(args, "output");
const rename = args.flags.has("rename");
const outputFasta = args.flags.has("fasta");

const wasmModule = await import(pathToFileURL(path.join(pkgDir, "deacon_wasm.js")).href);

if (typeof wasmModule.default === "function") {
  const wasmBytes = await readFile(path.join(pkgDir, "deacon_wasm_bg.wasm"));
  await wasmModule.default({ module_or_path: wasmBytes });
}

const indexBytes = await readFile(indexPath);
const index = new wasmModule.WasmIndex(new Uint8Array(indexBytes));
const session = new wasmModule.FilterSession(
  index,
  false, // deplete
  2,     // abs_threshold
  0.01,  // rel_threshold
  true,  // decompress_input
  false, // compress_output
  rename,
  outputFasta,
);

const input = createReadStream(readsPath, { highWaterMark: 256 * 1024 });
const output = createWriteStream(outputPath);

let bytesIn = 0;
let bytesOut = 0;

for await (const chunk of input) {
  bytesIn += chunk.length;
  const out = session.push_chunk(
    new Uint8Array(chunk.buffer, chunk.byteOffset, chunk.byteLength),
  );
  bytesOut += out.length;
  await writeChunk(output, out);
}

const tail = session.finish();
bytesOut += tail.length;
await writeChunk(output, tail);

output.end();
await once(output, "finish");

const stats = session.stats();

console.log(
  [
    "WASM file filter complete",
    `index=${index.info()}`,
    `bytes_in=${bytesIn}`,
    `bytes_out=${bytesOut}`,
    `reads_in=${stats.readsIn}`,
    `reads_out=${stats.readsOut}`,
  ].join(" "),
);
