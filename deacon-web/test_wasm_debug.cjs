const { readFileSync } = require('fs');
const { WasmIndex, filter } = require('./pkg-node/deacon_web.js');

const indexData = readFileSync('../data/chm13v2.k31w61.idx');
const index = new WasmIndex(new Uint8Array(indexData));
const readsData = readFileSync('../data/HG02334.100MB.fastq.gz');

// Capture console.log from WASM
const origLog = console.log;

console.log('=== deplete=true ===');
filter(index, new Uint8Array(readsData), true, 1, 0.01);

console.log('\n=== deplete=false ===');
filter(index, new Uint8Array(readsData), false, 1, 0.01);
