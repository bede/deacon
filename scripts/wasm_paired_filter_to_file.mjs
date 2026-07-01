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
    if (name === "deplete" || name === "rename" || name === "fasta") {
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
  if (!chunk || chunk.length === 0) return;
  if (!stream.write(chunk)) await once(stream, "drain");
}

const args = parseArgs(process.argv.slice(2));
const pkgDir = requireArg(args, "pkg");
const indexPath = requireArg(args, "index");
const reads1Path = requireArg(args, "reads1");
const reads2Path = requireArg(args, "reads2");
const output1Path = requireArg(args, "output1");
const output2Path = requireArg(args, "output2");
const deplete = args.flags.has("deplete");
const rename = args.flags.has("rename");
const outputFasta = args.flags.has("fasta");
const absThreshold = Number(args.get("abs-threshold") ?? 2);
const relThreshold = Number(args.get("rel-threshold") ?? 0.01);

const wasmModule = await import(pathToFileURL(path.join(pkgDir, "deacon_wasm.js")).href);

if (typeof wasmModule.default === "function") {
  const wasmBytes = await readFile(path.join(pkgDir, "deacon_wasm_bg.wasm"));
  await wasmModule.default({ module_or_path: wasmBytes });
}

const indexBytes = await readFile(indexPath);
const index = new wasmModule.WasmIndex(new Uint8Array(indexBytes));
const session = new wasmModule.PairedFilterSession(
  index,
  deplete,
  absThreshold,
  relThreshold,
  reads1Path.endsWith(".gz"), // decompress_r1
  reads2Path.endsWith(".gz"), // decompress_r2
  false, // compress_r1
  false, // compress_r2
  rename,
  outputFasta,
);

const input1 = createReadStream(reads1Path, { highWaterMark: 256 * 1024 });
const input2 = createReadStream(reads2Path, { highWaterMark: 256 * 1024 });
const output1 = createWriteStream(output1Path);
const output2 = createWriteStream(output2Path);

async function writePaired(out) {
  await Promise.all([writeChunk(output1, out.r1), writeChunk(output2, out.r2)]);
}

let bytesIn = 0;
const iter1 = input1[Symbol.asyncIterator]();
const iter2 = input2[Symbol.asyncIterator]();
let done1 = false;
let done2 = false;

while (!done1 || !done2) {
  if (!done1) {
    const { value, done } = await iter1.next();
    if (done) {
      done1 = true;
      await writePaired(session.finish_r1());
    } else {
      bytesIn += value.length;
      await writePaired(session.push_r1(new Uint8Array(value.buffer, value.byteOffset, value.byteLength)));
    }
  }

  if (!done2) {
    const { value, done } = await iter2.next();
    if (done) {
      done2 = true;
      await writePaired(session.finish_r2());
    } else {
      bytesIn += value.length;
      await writePaired(session.push_r2(new Uint8Array(value.buffer, value.byteOffset, value.byteLength)));
    }
  }
}

output1.end();
output2.end();
await Promise.all([once(output1, "finish"), once(output2, "finish")]);

const stats = session.stats();
session.free();

console.log(
  [
    "WASM paired file filter complete",
    `index=${index.info()}`,
    `bytes_in=${bytesIn}`,
    `reads_in=${stats.readsIn}`,
    `reads_out=${stats.readsOut}`,
  ].join(" "),
);
