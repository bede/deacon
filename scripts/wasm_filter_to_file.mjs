import { createReadStream, createWriteStream } from "node:fs";
import { readFile } from "node:fs/promises";
import { once } from "node:events";
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

const wasmModule = await import(pathToFileURL(path.join(pkgDir, "deacon_wasm.js")).href);

if (typeof wasmModule.default === "function") {
  const wasmBytes = await readFile(path.join(pkgDir, "deacon_wasm_bg.wasm"));
  await wasmModule.default({ module_or_path: wasmBytes });
}

const indexBytes = await readFile(indexPath);
const index = new wasmModule.WasmIndex(new Uint8Array(indexBytes));
const session = new wasmModule.FilterSession(
  index,
  false,
  2,
  0.01,
  true,
  false,
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
