// End-to-end on the real GFF pairs in test_data/ (gitignored; skipped when
// absent): GFF -> bin/extract_hits.py -> bin/chain.js, plus parity with the
// bin/chain_blocks.mjs CLI. Coverage/depth floors come from the prototype
// run of 2026-09-23 (docs/plans/C1_chain_core.md).
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import zlib from 'node:zlib';
import {execFileSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {SYNCHAIN, coverage} from './helpers.mjs';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const TMP = fs.mkdtempSync(path.join(os.tmpdir(), 'chain-real-'));

const PAIRS = [
    // the same assembly against itself: every region should be one diagonal block
    {name: 'axolotl', minCoverage: 0.95, depth: [1.0, 1.1]},
    // A. thaliana vs the allotetraploid A. suecica: ~2 subgenome copies plus
    // Arabidopsis' own ancient duplications
    {name: 'thaliana_vs_suecica', minCoverage: 0.99, depth: [2.0, 3.0]},
];

function extract(pair, role) {
    const out = path.join(TMP, `${pair}.${role}.hits.tsv.gz`);
    execFileSync('python3', [path.join(ROOT, 'bin/extract_hits.py'),
        '--gff', path.join(ROOT, 'test_data', pair, `${role}.gff`),
        '--chrom_sizes', path.join(ROOT, 'test_data', pair, `${role}.chrom.sizes`), '--out', out],
        {stdio: ['ignore', 'ignore', 'pipe']});
    return out;
}

const readSizes = (p) => fs.readFileSync(p, 'utf8').trim().split('\n').map((l) => l.split('\t'));

for (const pair of PAIRS) {
    const dir = path.join(ROOT, 'test_data', pair.name);
    const present = fs.existsSync(path.join(dir, 'target.gff')) && fs.existsSync(path.join(dir, 'comparison.gff'));
    test(`${pair.name}: coverage, depth, speed, CLI parity`, {skip: present ? false : 'test_data absent'}, () => {
        const files = {target: extract(pair.name, 'target'), comparison: extract(pair.name, 'comparison')};
        const sizesT = readSizes(path.join(dir, 'target.chrom.sizes'));
        const sizesC = readSizes(path.join(dir, 'comparison.chrom.sizes'));
        const rows = (f) => SYNCHAIN.parseHitsTsv(zlib.gunzipSync(fs.readFileSync(f)).toString('utf8'));
        const {tables: [A, B]} = SYNCHAIN.makeTables([
            {name: 'target', rows: rows(files.target), chromNames: sizesT.map((r) => r[0])},
            {name: 'reference', rows: rows(files.comparison), chromNames: sizesC.map((r) => r[0])}]);
        const genome = sizesT.reduce((s, r) => s + +r[1], 0);

        const auto = SYNCHAIN.autoParams(A, B);
        const ch = SYNCHAIN.createChainer(A, B);
        const times = [];
        let r;
        for (let i = 0; i < 3; i++) {
            // fresh chainer each time so every stage is actually recomputed
            const t = performance.now();
            r = SYNCHAIN.createChainer(A, B).run({...auto, minBlock: 5});
            times.push(performance.now() - t);
        }
        const cov = coverage(r.blocks, A.chromNames.length);
        const inv = r.blocks.filter((b) => b.orientation === '-').length;
        ch.run({...auto, minBlock: 5});
        const tFilter = performance.now();
        ch.run({...auto, minBlock: 12});
        const filterMs = performance.now() - tFilter;
        console.log(`  ${pair.name}: minPositive ${auto.minPositive.toFixed(3)}, ${r.stats.nAnchors} anchors -> `
            + `${r.blocks.length} blocks (${inv} '-'), coverage ${(100 * cov.union / genome).toFixed(2)} %, `
            + `depth ${cov.depth.toFixed(2)}, full recompute ${Math.min(...times).toFixed(0)} ms `
            + `(anchors ${r.stats.ms.anchors.toFixed(0)} / chains ${r.stats.ms.chains.toFixed(0)} / `
            + `blocks ${r.stats.ms.blocks.toFixed(0)}), minBlock-only change on a warm chainer ${filterMs.toFixed(1)} ms`);
        assert.ok(cov.union / genome >= pair.minCoverage, 'coverage');
        assert.ok(cov.depth >= pair.depth[0] && cov.depth <= pair.depth[1], 'depth');
        assert.ok(Math.min(...times) <= 200, 'speed');

        // the pipeline CLI must write exactly what the library produces
        const out = path.join(TMP, `${pair.name}.links.tsv`);
        execFileSync('node', [path.join(ROOT, 'bin/chain_blocks.mjs'), '--query_hits', files.target,
            '--subject_hits', files.comparison, '--min_block', '5', '--out', out], {stdio: ['ignore', 'ignore', 'pipe']});
        assert.equal(fs.readFileSync(out, 'utf8'), SYNCHAIN.blocksToTsv(r));
    });
}
