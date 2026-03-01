import { readFileSync } from 'fs';
import pkg from './pkg-node/deacon_web.js';
const { initSync, WasmIndex, filter } = pkg;

// Init WASM
const wasmBytes = readFileSync('./pkg-node/deacon_web_bg.wasm');
initSync({ module: wasmBytes });

// Load index
const indexData = readFileSync('../data/chm13v2.k31w61.idx');
const index = new WasmIndex(new Uint8Array(indexData));
console.log('Index:', index.info());

// Load reads (compressed)
const readsData = readFileSync('../data/HG02334.100MB.fastq.gz');
console.log(`Reads: ${(readsData.length / 1024 / 1024).toFixed(1)} MB compressed`);

// Test with deplete=true, abs=1, rel=0.01 (WASM defaults)
console.log('\n--- deplete=true, abs=1, rel=0.01 (WASM UI defaults) ---');
let t0 = performance.now();
let result = filter(index, new Uint8Array(readsData), true, 1, 0.01);
console.log(`Output: ${(result.length / 1024 / 1024).toFixed(1)} MB, took ${((performance.now() - t0) / 1000).toFixed(2)}s`);

// Test with deplete=true, abs=1, rel=0.05 (native test params)
console.log('\n--- deplete=true, abs=1, rel=0.05 ---');
t0 = performance.now();
result = filter(index, new Uint8Array(readsData), true, 1, 0.05);
console.log(`Output: ${(result.length / 1024 / 1024).toFixed(1)} MB, took ${((performance.now() - t0) / 1000).toFixed(2)}s`);

// Test with deplete=false, abs=1, rel=0.01
console.log('\n--- deplete=false, abs=1, rel=0.01 ---');
t0 = performance.now();
result = filter(index, new Uint8Array(readsData), false, 1, 0.01);
console.log(`Output: ${(result.length / 1024 / 1024).toFixed(1)} MB, took ${((performance.now() - t0) / 1000).toFixed(2)}s`);
