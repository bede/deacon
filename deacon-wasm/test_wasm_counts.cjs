const { readFileSync } = require('fs');
const { WasmIndex, filter } = require('./pkg-node/deacon_web.js');

// Load index
const indexData = readFileSync('../data/chm13v2.k31w61.idx');
const index = new WasmIndex(new Uint8Array(indexData));

// Use a small uncompressed FASTQ to remove any decompression variable
// Create a tiny test: just the first few records
const zlib = require('zlib');
const compressed = readFileSync('../data/HG02334.100MB.fastq.gz');
const decompressed = zlib.gunzipSync(compressed);
console.log(`Decompressed: ${(decompressed.length / 1024 / 1024).toFixed(1)} MB`);

// Filter on decompressed input
console.log('\n--- Uncompressed input, deplete=true, a=1, r=0.01 ---');
let result = filter(index, new Uint8Array(decompressed), true, 1, 0.01);
console.log('Done');

// Also try on compressed
console.log('\n--- Compressed input, deplete=true, a=1, r=0.01 ---');
result = filter(index, new Uint8Array(compressed), true, 1, 0.01);
console.log('Done');
