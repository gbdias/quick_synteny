#!/usr/bin/env node
// Pipeline entry point for bin/chain.js: chains one comparison (two genomes'
// hit tables, or one genome against itself with --self) and writes
// links.tsv (bin/chain.js's OUTPUTS header). Any parameter left out is
// auto-tuned by SYNCHAIN.autoParams, the same function the interactive page
// uses for its initial control values.
import fs from 'node:fs';
import v8 from 'node:v8';
import zlib from 'node:zlib';
import {createRequire} from 'node:module';
import {parseArgs} from 'node:util';

// Maglev, V8's mid-tier optimizing compiler, miscompiles chain.js's hot
// loops when node runs under x86-64 emulation on an arm64 host (an amd64
// container on Apple Silicon): SYNCHAIN.meanBestPositive came back 0 instead
// of 0.95 in about 1 run in 20, and almost always once the function got hot,
// which silently re-tunes minPositive to its 0.3 floor. TurboFan alone and
// the interpreter compute it correctly (160/160 under the same load, flag
// set here). Off before chain.js is loaded, so none of it is ever
// Maglev-compiled; the cost natively is nil at these run times.
v8.setFlagsFromString('--no-maglev');

const require = createRequire(import.meta.url);
const SYNCHAIN = require('./chain.js');

const {values: args} = parseArgs({
    options: {
        query_hits: {type: 'string'},
        subject_hits: {type: 'string'},
        self: {type: 'boolean', default: false},
        min_identity: {type: 'string'},
        max_gap: {type: 'string'},
        min_block: {type: 'string'},
        max_hit_rank: {type: 'string'},
        max_lookback: {type: 'string'},
        gap_penalty: {type: 'string'},
        out: {type: 'string'},
        help: {type: 'boolean', short: 'h', default: false},
    },
});

const USAGE = `usage: chain_blocks.mjs --query_hits Q.hits.tsv.gz (--subject_hits S.hits.tsv.gz | --self)
                       [--min_identity F] [--max_gap N] [--min_block N] [--max_hit_rank N]
                       [--max_lookback N] [--gap_penalty F] --out links.tsv`;

if (args.help) { console.log(USAGE); process.exit(0); }
if (!args.query_hits || !args.out || (!args.self && !args.subject_hits)) {
    console.error(USAGE);
    process.exit(2);
}

const readRows = (path) => SYNCHAIN.parseHitsTsv(zlib.gunzipSync(fs.readFileSync(path)).toString('utf8'));
const genomes = [{name: 'target', rows: readRows(args.query_hits)}];
if (!args.self) { genomes.push({name: 'reference', rows: readRows(args.subject_hits)}); }
const {tables} = SYNCHAIN.makeTables(genomes);
const A = tables[0], B = args.self ? tables[0] : tables[1];

const auto = SYNCHAIN.autoParams(A, args.self ? undefined : B);
// Every hit carries a miniprot Positive= above 0, so with hits on both sides
// the weaker mean best-hit positive can only be in (0, 1]. Anything else is
// a miscomputation (see --no-maglev above), not data -- fail loudly rather
// than chain with auto-tuned parameters derived from it. Exit status 3 is
// what the CHAIN_* processes retry on.
const w = auto.weakerMeanPositive;
if (A.n > 0 && B.n > 0 && !(w > 0 && w <= 1)) {
    console.error(`[chain_blocks] ERROR: weaker mean best-hit positive computed as ${w} from `
        + `${A.n}${args.self ? '' : ` + ${B.n}`} hits -- impossible for real hits, so the node `
        + `runtime miscomputed it (e.g. an amd64 container emulated on an arm64 host)`);
    process.exit(3);
}
const num = (v, parse) => (v === undefined ? undefined : parse(v));
const given = {
    minPositive: num(args.min_identity, Number),
    maxGap: num(args.max_gap, (v) => parseInt(v, 10)),
    minBlock: num(args.min_block, (v) => parseInt(v, 10)),
    maxHitRank: num(args.max_hit_rank, (v) => parseInt(v, 10)),
    maxLookback: num(args.max_lookback, (v) => parseInt(v, 10)),
    gapPenalty: num(args.gap_penalty, Number),
};
const params = {...auto};
for (const [k, v] of Object.entries(given)) { if (v !== undefined) { params[k] = v; } }
params.selfMode = args.self;

const result = SYNCHAIN.chain(A, B, params);
fs.writeFileSync(args.out, SYNCHAIN.blocksToTsv(result));

const shown = ['minPositive', 'maxGap', 'minBlock', 'maxHitRank', 'maxLookback', 'gapPenalty']
    .map((k) => `${k}=${k === 'minPositive' ? params[k].toFixed(4) : params[k]}${given[k] === undefined ? ' (auto)' : ''}`)
    .join(' ');
const s = result.stats;
console.error(`[chain_blocks] ${args.self ? 'self' : 'cross'}: ${shown}; weaker mean best-hit positive `
    + `${auto.weakerMeanPositive.toFixed(4)}`);
console.error(`[chain_blocks] ${A.n}${args.self ? '' : ` + ${B.n}`} hits -> ${s.nAnchors} anchors -> `
    + `${s.nChains} chains -> ${s.nBlocks} block(s) in ${s.ms.total.toFixed(0)} ms -> ${args.out}`);
