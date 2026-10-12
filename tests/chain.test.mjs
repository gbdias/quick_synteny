// bin/chain.js on small synthetic genomes: the cases its CHAINING RULES are
// about. No data or Nextflow needed: `node --test tests/chain.test.mjs`.
import assert from 'node:assert/strict';
import fs from 'node:fs';
import test from 'node:test';
import {createRequire} from 'node:module';

const require = createRequire(import.meta.url);
const SYNCHAIN = require('../bin/chain.js');

// one hit row per [protein, locus, positive, rank], on chr1, 2 kb apart by locus
const rows = (hits) => hits.map(([p, l, pos = 0.95, rank = 1]) =>
    ['chr1', String(l * 2000), String(l * 2000 + 1500), '+', String(pos), String(pos), '500', String(rank), p, String(l)]);

// target locus i <-> reference locus refOf(i), one protein each, for i < n
function genomes(n, refOf, extra = []) {
    const name = (i) => `p${String(i).padStart(4, '0')}`;
    const t = [], r = [];
    for (let i = 0; i < n; i++) { t.push([name(i), i]); r.push([name(i), refOf(i)]); }
    for (const [i, l, pos, rank] of extra) { r.push([name(i), l, pos, rank]); }
    r.sort((a, b) => a[1] - b[1]);
    const {tables} = SYNCHAIN.makeTables([{name: 'target', rows: rows(t)}, {name: 'reference', rows: rows(r)}]);
    return tables;
}

const blocks = (tables, params) => SYNCHAIN.chain(tables[0], tables[1], {...SYNCHAIN.DEFAULTS, ...params}).blocks
    .map((b) => `${b.orientation}${b.nAnchors}`).sort();

test('a nested inversion is its own block, inside the bridging block', () => {
    // 100 genes, 36-65 inverted (loci 35..64)
    const tables = genomes(100, (i) => (i >= 35 && i <= 64 ? 99 - i : i));
    // the + chain bridges through the inversion's middle gene at either gap,
    // and the other 29 inverted genes stay one - block (VERSION 1 gave
    // -5, -5 at 25 and lost the inversion at 31)
    for (const maxGap of [25, 31]) {
        assert.deepEqual(blocks(tables, {maxGap, minBlock: 5}), ['+71', '-29'], `maxGap ${maxGap}`);
    }
});

test('tandem copies collapse to one anchor instead of a parallel chain', () => {
    // each protein also hits the next reference locus, weaker: without the
    // collapse that's a second, parallel 29-anchor block
    const extra = Array.from({length: 29}, (_, i) => [i, i + 1, 0.8, 2]);
    assert.deepEqual(blocks(genomes(30, (i) => i, extra), {minBlock: 5}), ['+30']);
});

test('min block is a pure filter, and blocks keep to the rules', () => {
    // a collinear run, an inversion, a translocated piece, and scattered noise
    let seed = 7;
    const rand = () => { seed = (seed * 48271) % 2147483647; return seed / 2147483647; };
    const refOf = (i) => (i < 60 ? i : i < 90 ? 149 - i : i < 120 ? i + 40 : i - 30);
    const extra = Array.from({length: 40}, () => [Math.floor(rand() * 160), Math.floor(rand() * 160), 0.75, 2]);
    const tables = genomes(160, refOf, extra);
    const all = SYNCHAIN.chain(tables[0], tables[1], {...SYNCHAIN.DEFAULTS, minBlock: 2});
    for (const t of [3, 5, 10, 20]) {
        const at = SYNCHAIN.chain(tables[0], tables[1], {...SYNCHAIN.DEFAULTS, minBlock: t});
        assert.deepEqual(at.blocks.map((b) => Array.from(b.members)), all.blocks.filter((b) => b.nAnchors >= t).map((b) => Array.from(b.members)));
    }
    const an = all.anchors, seen = new Set();
    for (const b of all.blocks) {
        const m = Array.from(b.members);
        for (const x of m) { assert.ok(!seen.has(x), 'an anchor in two blocks'); seen.add(x); }
        for (let k = 1; k < m.length; k++) {
            const dA = an.rA[m[k]] - an.rA[m[k - 1]];
            const dB = b.orientation === '+' ? an.rB[m[k]] - an.rB[m[k - 1]] : an.rB[m[k - 1]] - an.rB[m[k]];
            assert.ok(dA >= 1 && dA <= SYNCHAIN.DEFAULTS.maxGap && dB >= 1 && dB <= SYNCHAIN.DEFAULTS.maxGap, 'a step breaks rule 4');
        }
    }
});

test("the page's payload version is chain.js's", () => {
    const py = fs.readFileSync(new URL('../bin/plot_synteny_interactive.py', import.meta.url), 'utf8');
    const m = py.match(/VERSION: decodePayload refuses any other\s*\n\s*'version': (\d+)/);
    assert.ok(m, "plot_synteny_interactive.py's payload version not found");
    assert.equal(Number(m[1]), SYNCHAIN.VERSION);
});
