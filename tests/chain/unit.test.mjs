// Behavioural unit tests for bin/chain.js, including regressions for the
// defects found in the old bin/build_synteny_blocks.py chainer (removed; see
// git history): inversions reported as '+', gene-sparse
// genomes yielding no blocks, and a tandem family aborting its whole
// chromosome pair.
import test from 'node:test';
import assert from 'node:assert/strict';
import zlib from 'node:zlib';
import {SYNCHAIN, tablesFrom, P, syntheticPair} from './helpers.mjs';

const collinear = (n, spacing, prefix = 'p', start = 1_000_000) =>
    Array.from({length: n}, (_, i) => ({chrom: 'A1', pos: start + i * spacing, prot: `${prefix}${i}`}));
const onB = (hits, map = (h) => h) => hits.map((h) => ({...map(h), chrom: 'B1'}));

test('perfect inversion is one block with orientation "-"', () => {
    const a = collinear(30, 50_000);
    const b = a.map((h, i) => ({chrom: 'B1', pos: 1_000_000 + (29 - i) * 50_000, prot: h.prot}));
    const [A, B] = tablesFrom(a, b);
    const r = SYNCHAIN.chain(A, B, P());
    assert.equal(r.blocks.length, 1);
    assert.equal(r.blocks[0].orientation, '-');
    assert.equal(r.blocks[0].nAnchors, 30);
});

test('collinear genes 1 Mb apart still chain (gene-rank distance, not bp)', () => {
    const a = collinear(200, 1_000_000);
    const [A, B] = tablesFrom(a, onB(a));
    const r = SYNCHAIN.chain(A, B, P());
    assert.equal(r.blocks.length, 1);
    assert.equal(r.blocks[0].nAnchors, 200);
    assert.equal(r.blocks[0].orientation, '+');
});

test('a tandem family next to a real run neither forms a block nor hides the run', () => {
    const a = [], b = [];
    for (let fam = 0; fam < 2; fam++) {
        for (let c = 0; c < 6; c++) { a.push({chrom: 'A1', pos: 5_000_000 + fam * 70_000 + c * 11_000, prot: `fam${fam}`}); }
        b.push({chrom: 'B1', pos: 5_000_000 + fam * 70_000, prot: `fam${fam}`});
    }
    for (let i = 0; i < 8; i++) {
        a.push({chrom: 'A1', pos: 50_000_000 + i * 40_000, prot: `g${i}`});
        b.push({chrom: 'B1', pos: 50_000_000 + i * 40_000, prot: `g${i}`});
    }
    // unmatched genes between the family and the run: distance is counted in
    // genes, so without these the two would be adjacent however many bp apart
    for (let i = 0; i < 10; i++) {
        a.push({chrom: 'A1', pos: 20_000_000 + i * 40_000, prot: `ua${i}`});
        b.push({chrom: 'B1', pos: 20_000_000 + i * 40_000, prot: `ub${i}`});
    }
    const [A, B] = tablesFrom(a, b);
    const r = SYNCHAIN.chain(A, B, P({maxGap: 5}));
    assert.deepEqual(r.blocks.map((x) => x.nAnchors), [8]);
    assert.ok(r.blocks[0].qStart >= 50_000_000 - 500);
});

test('isoforms at one locus count once', () => {
    const a = [], b = [];
    for (let g = 0; g < 6; g++) {
        for (const iso of ['a', 'b']) {
            a.push({chrom: 'A1', pos: 1_000_000 + g * 30_000, prot: `g${g}${iso}`});
            b.push({chrom: 'B1', pos: 2_000_000 + g * 30_000, prot: `g${g}${iso}`});
        }
    }
    const [A, B] = tablesFrom(a, b);
    const r = SYNCHAIN.chain(A, B, P());
    assert.deepEqual(r.blocks.map((x) => x.nAnchors), [6]);
});

test('self mode: duplicated segment found, no diagonal, no tandem-neighbour blocks', () => {
    const hits = [];
    for (let i = 0; i < 10; i++) {
        hits.push({chrom: 'A1', pos: 1_000_000 + i * 20_000, prot: `d${i}`});
        hits.push({chrom: 'A2', pos: 3_000_000 + i * 20_000, prot: `d${i}`, rank: 2, positive: 0.9});
    }
    // tandem pairs: each protein also hits the very next locus on A1
    for (let i = 0; i < 10; i++) {
        hits.push({chrom: 'A1', pos: 8_000_000 + i * 40_000, prot: `t${i}`});
        hits.push({chrom: 'A1', pos: 8_000_000 + i * 40_000 + 20_000, prot: `t${i}`, rank: 2});
    }
    const [A] = tablesFrom(hits, null, ['A1', 'A2']);
    const r = SYNCHAIN.chain(A, A, P({selfMode: true}));
    assert.equal(r.blocks.length, 1);
    const blk = r.blocks[0];
    assert.equal(A.chromNames[blk.qChrom], 'A1');
    assert.equal(A.chromNames[blk.sChrom], 'A2');
    assert.equal(blk.nAnchors, 10);
});

test('maxHitRank=1 ignores secondary hits', () => {
    const a = collinear(20, 30_000);
    const b = onB(a).map((h) => ({...h, rank: 2}));
    const [A, B] = tablesFrom(a, b);
    assert.equal(SYNCHAIN.chain(A, B, P()).blocks.length, 1);
    assert.equal(SYNCHAIN.chain(A, B, P({maxHitRank: 1})).blocks.length, 0);
});

test('minPositive filters hits (and so loci) below the floor', () => {
    const a = collinear(12, 30_000).map((h, i) => ({...h, positive: i < 6 ? 0.95 : 0.4}));
    const [A, B] = tablesFrom(a, onB(a));
    assert.deepEqual(SYNCHAIN.chain(A, B, P({minPositive: 0.5})).blocks.map((x) => x.nAnchors), [6]);
    assert.deepEqual(SYNCHAIN.chain(A, B, P({minPositive: 0.3})).blocks.map((x) => x.nAnchors), [12]);
});

const {A: SA, B: SB} = syntheticPair({nGenes: 4000, nChrom: 6, seed: 7});

test('minBlock is a pure filter: blocks at T == blocks at a lower T filtered to >= T', () => {
    const lo = SYNCHAIN.chain(SA, SB, P({minPositive: 0.3, minBlock: 3}));
    for (const T of [5, 8, 15]) {
        const hi = SYNCHAIN.chain(SA, SB, P({minPositive: 0.3, minBlock: T}));
        const filtered = {...lo, blocks: lo.blocks.filter((b) => b.nAnchors >= T)};
        assert.equal(SYNCHAIN.blocksToTsv(hi), SYNCHAIN.blocksToTsv(filtered), `T=${T}`);
    }
});

test('stage cache never leaks results across parameter changes', () => {
    const ch = SYNCHAIN.createChainer(SA, SB);
    const sets = [P({minPositive: 0.6}), P({minPositive: 0.6, maxGap: 10}), P({minPositive: 0.3, maxHitRank: 1}),
                  P({minPositive: 0.6, minBlock: 12}), P({minPositive: 0.6, gapPenalty: 0.05}), P({minPositive: 0.6})];
    for (const s of sets) {
        assert.equal(SYNCHAIN.blocksToTsv(ch.run(s)), SYNCHAIN.blocksToTsv(SYNCHAIN.chain(SA, SB, s)), JSON.stringify(s));
    }
});

test('deterministic: identical input gives identical output', () => {
    const again = syntheticPair({nGenes: 4000, nChrom: 6, seed: 7});
    const p = P({minPositive: 0.3, gapPenalty: 0.02});
    assert.equal(SYNCHAIN.blocksToTsv(SYNCHAIN.chain(SA, SB, p)), SYNCHAIN.blocksToTsv(SYNCHAIN.chain(again.A, again.B, p)));
});

test('gap penalty splits a chain at a large gap', () => {
    const a = [...collinear(10, 30_000), ...collinear(10, 30_000, 'q', 20_000_000)];
    // 15 unmatched loci between the two runs on both genomes
    const filler = (chrom) => Array.from({length: 15}, (_, i) => ({chrom, pos: 10_000_000 + i * 30_000, prot: `x${chrom}${i}`}));
    const [A, B] = tablesFrom([...a, ...filler('A1')], [...onB(a), ...filler('B1')]);
    assert.deepEqual(SYNCHAIN.chain(A, B, P()).blocks.map((x) => x.nAnchors), [20]);
    // joining across the 16-rank step scores f(prev) + 1 - 15*pen: it still
    // beats starting over (1) until pen exceeds (10 + 1 - 1) / 15 = 0.667
    assert.deepEqual(SYNCHAIN.chain(A, B, P({gapPenalty: 0.6})).blocks.map((x) => x.nAnchors), [20]);
    assert.deepEqual(SYNCHAIN.chain(A, B, P({gapPenalty: 0.7})).blocks.map((x) => x.nAnchors), [10, 10]);
});

test('every anchor lands in exactly one chain', () => {
    const ch = SYNCHAIN.createChainer(SA, SB);
    const r = ch.run(P({minPositive: 0.3, minBlock: 2}));
    assert.equal(r.stats.nChains >= r.blocks.length, true);
    const seen = new Uint8Array(r.anchors.n);
    for (const b of r.blocks) { for (const m of b.members) { assert.equal(seen[m], 0); seen[m] = 1; } }
});

test('autoParams follows the PARAMETERS rules', () => {
    const a = collinear(10, 30_000).map((h) => ({...h, positive: 0.7}));
    const [A, B] = tablesFrom(a, onB(a).map((h) => ({...h, positive: 0.95})));
    const p = SYNCHAIN.autoParams(A, B);
    assert.equal(p.minPositive, 0.7);
    assert.equal(p.minBlock, 5);
    assert.equal(SYNCHAIN.autoParams(B).minPositive, 0.9);
    assert.equal(SYNCHAIN.autoParams(B).minBlock, 15);
});

// ---- embedded payload (bin/chain.js's EMBEDDED PAYLOAD header) round trip, with a reference encoder

function encodeColumn(arr, dtype) {
    const bytes = Buffer.from(arr.buffer, arr.byteOffset, arr.byteLength);
    return {dtype, data: zlib.gzipSync(bytes).toString('base64')};
}

function encodeGenome(t) {
    const startDelta = new Uint32Array(t.n), length = new Uint32Array(t.n);
    let prevChrom = -1, prev = 0;
    for (let i = 0; i < t.n; i++) {
        if (t.chrom[i] !== prevChrom) { prevChrom = t.chrom[i]; prev = 0; }
        startDelta[i] = t.start[i] - prev; prev = t.start[i];
        length[i] = t.end[i] - t.start[i];
    }
    return {n: t.n, nLoci: t.nLoci, cols: {
        chrom: encodeColumn(t.chrom, 'u16'), startDelta: encodeColumn(startDelta, 'u32'),
        length: encodeColumn(length, 'u32'), strand: encodeColumn(t.strand, 'u8'),
        positive: encodeColumn(t.positive, 'u16'), rank: encodeColumn(t.rank, 'u8'),
        prot: encodeColumn(t.prot, 'u32'), locus: encodeColumn(t.locus, 'u32'),
    }};
}

test('decodePayload reproduces the tables exactly', async () => {
    const {proteins, tables: [A, B]} = SYNCHAIN.makeTables([
        {name: 'target', rows: syntheticRows(SA)}, {name: 'reference', rows: syntheticRows(SB)}]);
    const payload = {version: 1, proteins: zlib.gzipSync(proteins.join('\n')).toString('base64'),
                     genomes: {target: encodeGenome(A), reference: encodeGenome(B)}};
    const d = await SYNCHAIN.decodePayload(payload, {
        target: {names: A.chromNames, sizes: A.chromNames.map(() => 0)},
        reference: {names: B.chromNames, sizes: B.chromNames.map(() => 0)}});
    assert.deepEqual(d.proteins, proteins);
    for (const [orig, dec] of [[A, d.target], [B, d.reference]]) {
        for (const k of ['chrom', 'start', 'end', 'strand', 'positive', 'rank', 'prot', 'locus']) {
            assert.deepEqual(Array.from(dec[k]), Array.from(orig[k]), k);
        }
        assert.equal(dec.nLoci, orig.nLoci);
    }
    const p = P({minPositive: 0.6});
    assert.equal(SYNCHAIN.blocksToTsv(SYNCHAIN.chain(d.target, d.reference, p)), SYNCHAIN.blocksToTsv(SYNCHAIN.chain(A, B, p)));
});

// rebuild TSV-style rows from a table, to feed makeTables in the round trip
function syntheticRows(t) {
    const rows = [];
    for (let i = 0; i < t.n; i++) {
        rows.push([t.chromNames[t.chrom[i]], String(t.start[i]), String(t.end[i]), t.strand[i] ? '+' : '-',
                   (t.positive[i] / 10000).toFixed(4), (t.positive[i] / 10000).toFixed(4), '1000', String(t.rank[i]),
                   `p${String(t.prot[i]).padStart(6, '0')}`, String(t.locus[i])]);
    }
    return rows;
}
