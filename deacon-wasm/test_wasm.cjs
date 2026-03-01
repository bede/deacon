const { readFileSync, writeFileSync } = require('fs');
const { WasmIndex, filter } = require('./pkg-node/deacon_web.js');

// Load index
const indexData = readFileSync('../data/chm13v2.k31w61.idx');
const index = new WasmIndex(new Uint8Array(indexData));
console.log('Index:', index.info());

// Load reads (compressed)
const readsData = readFileSync('../data/HG02334.100MB.fastq.gz');
console.log(`Reads: ${(readsData.length / 1024 / 1024).toFixed(1)} MB compressed`);

// Test with deplete=true, abs=1, rel=0.01 (WASM UI defaults)
console.log('\n--- deplete=true, abs=1, rel=0.01 (WASM UI defaults) ---');
let t0 = performance.now();
let result1 = filter(index, new Uint8Array(readsData), true, 1, 0.01);
console.log(`Output: ${(result1.length / 1024 / 1024).toFixed(1)} MB, took ${((performance.now() - t0) / 1000).toFixed(2)}s`);
writeFileSync('out_wasm_deplete_r001.fq.gz', result1);

// Test with deplete=false, abs=1, rel=0.01
console.log('\n--- deplete=false, abs=1, rel=0.01 ---');
t0 = performance.now();
let result2 = filter(index, new Uint8Array(readsData), false, 1, 0.01);
console.log(`Output: ${(result2.length / 1024 / 1024).toFixed(1)} MB, took ${((performance.now() - t0) / 1000).toFixed(2)}s`);
writeFileSync('out_wasm_search_r001.fq.gz', result2);

// Test with deplete=true, abs=1, rel=0.05
console.log('\n--- deplete=true, abs=1, rel=0.05 ---');
t0 = performance.now();
let result3 = filter(index, new Uint8Array(readsData), true, 1, 0.05);
console.log(`Output: ${(result3.length / 1024 / 1024).toFixed(1)} MB, took ${((performance.now() - t0) / 1000).toFixed(2)}s`);
writeFileSync('out_wasm_deplete_r005.fq.gz', result3);
