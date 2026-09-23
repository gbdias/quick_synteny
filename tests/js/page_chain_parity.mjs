#!/usr/bin/env node
// Runs a generated synteny page's own client-side chaining path (SYN.init ->
// SYN.startChainer -> SYN.requestChain -> SYN.applyChainResult) in node:vm
// and checks that its "download blocks TSV" export is byte-identical to
// `node bin/chain_blocks.mjs` run with the same parameters against the same
// hit tables -- the two are supposed to be the exact same chainer
// (bin/chain.js), just inlined into the page instead of run as a CLI, so
// any difference here is a real bug in how the page drives it (payload
// decoding, param wiring, or the min-block-size post-filter), not a
// legitimate discrepancy.
//
// Usage: node tests/js/page_chain_parity.mjs <page.html> <target.hits.tsv.gz> <comparison.hits.tsv.gz>
//
// The page's own scripts run in a node:vm context with a `window` stub (see
// makeSandbox), a synchronous stand-in for Worker (see makeFakeWorkerEnv --
// SYN.startChainer's `new Worker(URL.createObjectURL(new Blob([...])))` path
// still runs for real, just against a Worker that executes its source
// in-process instead of on a real thread), and Bokeh sources/widgets as
// plain {data, change: {emit(){}}, selected: {indices: []}} objects (or
// {value}/{active}/{title:{text}} for widgets/figures) -- no real Bokeh
// runtime needed, same idiom tests/js/page_geometry.mjs already uses.

import { readFileSync, unlinkSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import path from 'node:path';
import vm from 'node:vm';

function fail(msg) {
    console.error(`FAIL: ${msg}`);
    process.exitCode = 1;
}

function assert(cond, msg) {
    if (!cond) fail(msg);
}

// ---- extraction: pull the raw chain.js copy + the main executable script
// (chain.js + SHARED_JS + SYN.data/SYN.hitsPayload) out of a generated page --

function allScriptBodies(html) {
    return [...html.matchAll(/<script([^>]*)>([\s\S]*?)<\/script>/g)]
        .map((m) => ({ attrs: m[1], body: m[2] }));
}

function extractPage(pagePath) {
    const html = readFileSync(pagePath, 'utf8');
    const scripts = allScriptBodies(html);
    const synchainSrcTag = scripts.find((s) => /id="synchain-src"/.test(s.attrs));
    if (!synchainSrcTag) throw new Error(`${pagePath}: no <script type="text/plain" id="synchain-src"> found`);
    const mainScript = scripts[scripts.length - 1];
    if (mainScript === synchainSrcTag) throw new Error(`${pagePath}: synchain-src is the LAST script -- expected the executable chain.js+SHARED_JS block last`);
    return { chainJsRaw: synchainSrcTag.body, mainScript: mainScript.body };
}

// ---- sandbox: node:vm context with window/document/Worker/Blob/URL stubs --

// Blob subclass that remembers its constructor parts (for the FakeWorker's
// createObjectURL below) while still being a REAL Blob otherwise --
// SYNCHAIN.decodePayload's own gzip decoding (gunzipDefault) does
// `new Blob([bytes]).stream()`, which needs a genuine Blob, not a stand-in.
function makeTrackedBlob(RealBlob) {
    return class TrackedBlob extends RealBlob {
        constructor(parts, opts) {
            super(parts, opts);
            this.__fakeParts = parts;
        }
    };
}

// node:vm's separate context is a genuinely separate V8 realm, same as a
// real Worker's own global scope -- but a real postMessage's structured
// clone actually MATERIALIZES fresh objects/typed arrays in the receiving
// realm, while just handing over the same cross-realm object references
// (the shortcut a naive stand-in would take) does not: V8's inline caches
// are per-realm, and running chain.js's hot loops over typed arrays that
// belong to a DIFFERENT realm was measured (scratch/c2/repro_realm.mjs, not
// committed) at ~3.5x slower than same-realm arrays of identical content --
// an artifact of this harness, not of the page. cloneIntoRealm below is a
// structured-clone stand-in restricted to what these messages actually
// carry (plain objects/arrays, strings, numbers, and typed arrays), so
// SYN.requestChain's own timing report reflects chain.js's real cost, not
// this harness's realm-crossing tax.
function cloneIntoRealm(value, ctors) {
    if (value === null || typeof value !== 'object') return value;
    if (ArrayBuffer.isView(value)) {
        const Ctor = ctors[value.constructor.name];
        return Ctor ? Ctor.from(value) : value.slice();
    }
    if (Array.isArray(value)) { return value.map((v) => cloneIntoRealm(v, ctors)); }
    const out = new ctors.Object();
    for (const k of Object.keys(value)) { out[k] = cloneIntoRealm(value[k], ctors); }
    return out;
}

// A fully synchronous Worker stand-in: postMessage on the "main thread" side
// clones its data into the worker's own realm (see cloneIntoRealm above),
// then calls straight into a second vm context running the worker source
// (chain.js + SYN.WORKER_SHIM_SRC); that context's own postMessage clones
// the reply back and calls straight into this instance's onmessage -- no
// timers, no real thread, so a chain request's whole round trip finishes
// inside one JS turn (this is also what makes SYN.requestChain's "at most
// one in flight" bookkeeping trivial to exercise/measure from the test below).
function makeFakeWorkerEnv(blobRegistry) {
    let nextId = 1;

    class FakeWorker {
        constructor(url) {
            const src = blobRegistry.get(url);
            if (src === undefined) throw new Error(`FakeWorker: unknown blob url ${url}`);
            const workerSandbox = { console, performance };
            workerSandbox.self = workerSandbox;
            vm.createContext(workerSandbox);
            const workerCtors = vm.runInContext('({Object, Array, Uint8Array, Uint16Array, Uint32Array, Float64Array})', workerSandbox);
            const mainCtors = { Object, Array, Uint8Array, Uint16Array, Uint32Array, Float64Array };
            workerSandbox.postMessage = (data) => {
                if (this.onmessage) { this.onmessage({ data: cloneIntoRealm(data, mainCtors) }); }
            };
            vm.runInContext(src, workerSandbox, { filename: 'worker.js' });
            this._workerSandbox = workerSandbox;
            this._workerCtors = workerCtors;
            this.onmessage = null;
        }
        postMessage(data) {
            this._workerSandbox.self.onmessage({ data: cloneIntoRealm(data, this._workerCtors) });
        }
        terminate() {}
    }

    return {
        FakeWorker,
        url: {
            createObjectURL(blob) {
                const id = `blob:fake-${nextId++}`;
                blobRegistry.set(id, blob.__fakeParts.join(''));
                return id;
            },
            revokeObjectURL(id) { blobRegistry.delete(id); },
        },
    };
}

function makeSandbox(chainJsRaw) {
    const sandbox = {};
    sandbox.window = sandbox; // window === globalThis, same as a real page
    sandbox.console = console;
    sandbox.TextDecoder = TextDecoder;
    sandbox.TextEncoder = TextEncoder;
    sandbox.DecompressionStream = DecompressionStream;
    sandbox.Response = Response;
    sandbox.atob = atob;
    sandbox.btoa = btoa;
    sandbox.performance = performance;
    sandbox.Blob = makeTrackedBlob(Blob);

    const blobRegistry = new Map();
    const { FakeWorker, url } = makeFakeWorkerEnv(blobRegistry);
    sandbox.Worker = FakeWorker;
    sandbox.URL = url;

    // SYN.startChainer's only DOM touch: reads the inert copy of chain.js
    // back out to build the Worker's source -- see build_page()'s own
    // comment on why that copy exists at all.
    sandbox.document = {
        getElementById(id) {
            return id === 'synchain-src' ? { textContent: chainJsRaw } : null;
        },
    };

    vm.createContext(sandbox);
    return sandbox;
}

// ---- stub Bokeh sources/widgets: only what SYN.init/SYN.applyChainResult
// actually touch (.data/.change.emit(), .value/.active, .title.text,
// .x_range/.y_range.start/end, .high) ----

function stubSource() {
    return { data: {}, selected: { indices: [] }, change: { emit() {} } };
}

function stubFig() {
    return { x_range: { start: 0, end: 0 }, y_range: { start: 0, end: 0 }, title: { text: '' } };
}

function buildStubUi() {
    return {
        querySource: stubSource(), subjectSource: stubSource(), labelSource: stubSource(),
        ribbonSource: stubSource(), gapSource: stubSource(),
        dpQuerySource: stubSource(), dpSubjectSource: stubSource(), dpGridSource: stubSource(),
        dpQueryLabelSource: stubSource(), dpSubjectLabelSource: stubSource(),
        dpSegmentSource: stubSource(), dpGapSource: stubSource(), dotplotFig: stubFig(),
        colorSpinner: { value: 10 },
        minBlockSpinner: { value: 3, high: 1000 },
        minIdentitySpinner: { value: 50 },
        maxGapSpinner: { value: 25 },
        hitRankSelect: { value: 'all' },
        selfLinksToggle: { active: false },
        hideSyntenyToggle: { active: false },
        chainStatusDiv: { text: '' },
        detailBarSource: stubSource(), detailRibbonSource: stubSource(), detailLabelSource: stubSource(),
        detailFig: stubFig(), detailGapSource: stubSource(),
    };
}

async function waitFor(cond, timeoutMs, label) {
    const t0 = Date.now();
    while (!cond()) {
        if (Date.now() - t0 > timeoutMs) throw new Error(`timeout waiting for: ${label}`);
        await new Promise((r) => setTimeout(r, 5));
    }
}

// ---- chain_blocks.mjs comparison ----

const REPO_ROOT = path.resolve(import.meta.dirname, '..', '..');
const CHAIN_BLOCKS = path.join(REPO_ROOT, 'bin', 'chain_blocks.mjs');

function runChainBlocks(targetHits, comparisonHits, flags, outPath) {
    execFileSync('node', [CHAIN_BLOCKS, '--query_hits', targetHits, '--subject_hits', comparisonHits,
        ...flags, '--out', outPath], { stdio: ['ignore', 'ignore', 'pipe'] });
    return readFileSync(outPath, 'utf8');
}

// ---- one param-set check: reset the page's UI to the given baseline (the
// auto/default values SYN.init itself settled on), apply this case's own
// overrides on top, re-chain, export, and diff against a fresh
// chain_blocks.mjs run with the equivalent flags -- each case is
// independent of whichever ran before it, matching a fresh `node
// bin/chain_blocks.mjs` invocation's own independence ----

async function checkParamSet(sandbox, ui, name, baseline, overrides, cliFlags, targetHits, comparisonHits, tmpPath) {
    const effective = Object.assign({}, baseline, overrides);
    ui.minIdentitySpinner.value = Math.round(effective.minPositive * 100);
    ui.maxGapSpinner.value = effective.maxGap;
    ui.hitRankSelect.value = { 1: 'best only', 2: '≤ 2', 3: '≤ 3' }[effective.maxHitRank] || 'all';
    ui.minBlockSpinner.value = effective.minBlock;

    const t0 = performance.now();
    sandbox.SYN.requestChain(sandbox.SYN.currentChainParams());
    const t1 = performance.now();

    const minScore = ui.minBlockSpinner.value;
    const filtered = sandbox.SYN.data.crossLinksFlat.filter((l) => l.score >= minScore);
    const pageTsv = sandbox.SYNCHAIN.linksToTsv(filtered);

    const cliTsv = runChainBlocks(targetHits, comparisonHits, cliFlags, tmpPath);

    if (pageTsv === cliTsv) {
        console.log(`  ${name}: OK (${filtered.length} block(s), page params->sources ${(t1 - t0).toFixed(1)} ms)`);
    } else {
        fail(`${name}: page export differs from chain_blocks.mjs (flags: ${cliFlags.join(' ') || '(none)'})`);
        const pageLines = pageTsv.split('\n'), cliLines = cliTsv.split('\n');
        console.error(`    page: ${pageLines.length} line(s), cli: ${cliLines.length} line(s)`);
        for (let i = 0; i < Math.max(pageLines.length, cliLines.length); i++) {
            if (pageLines[i] !== cliLines[i]) {
                console.error(`    first mismatch at line ${i}:\n      page: ${pageLines[i]}\n      cli:  ${cliLines[i]}`);
                break;
            }
        }
    }
    return t1 - t0;
}

async function main() {
    const [, , pagePath, targetHits, comparisonHits] = process.argv;
    if (!pagePath || !targetHits || !comparisonHits) {
        console.error('usage: node tests/js/page_chain_parity.mjs <page.html> <target.hits.tsv.gz> <comparison.hits.tsv.gz>');
        process.exit(2);
    }

    const { chainJsRaw, mainScript } = extractPage(pagePath);
    const sandbox = makeSandbox(chainJsRaw);
    vm.runInContext(mainScript, sandbox, { filename: 'page-main-script' });
    assert(typeof sandbox.SYNCHAIN === 'object', `${pagePath}: SYNCHAIN not defined after running the page's main script`);
    assert(typeof sandbox.SYN.hitsPayload === 'object', `${pagePath}: SYN.hitsPayload not set`);

    console.log(`Running ${pagePath}'s own init path in node:vm ...`);
    const ui = buildStubUi();
    sandbox.SYN.init(ui);
    // SYN.init's own chain kickoff is asynchronous (SYNCHAIN.decodePayload
    // uses real DecompressionStream machinery) -- wait for it exactly the
    // way a real page's first paint would, then everything after this is
    // synchronous (see makeFakeWorkerEnv's own comment).
    await waitFor(() => sandbox.SYN.chain && sandbox.SYN.chain.appliedSeq >= 0, 10000,
        'first SYN.applyChainResult after SYN.init');
    console.log(`  init settled: minIdentity=${ui.minIdentitySpinner.value}% maxGap=${ui.maxGapSpinner.value} `
        + `minBlock(auto)=${ui.minBlockSpinner.value} -- ${sandbox.SYN.data.crossLinksFlat.length} cross block(s)`);

    // the auto-tuned baseline SYN.init itself settled on -- every case below
    // starts from exactly this, then layers its own overrides on top, so
    // cases are independent of each other and of whatever order they run in
    const auto = sandbox.SYNCHAIN.autoParams(sandbox.SYN.chain.tables.target, sandbox.SYN.chain.tables.reference);
    const baseline = { minPositive: auto.minPositive, maxGap: 25, maxHitRank: 255, minBlock: auto.minBlock };

    const tmpPath = path.join(path.dirname(pagePath), '.page_chain_parity_tmp.tsv');
    const timings = [];
    try {
        console.log('- parameter sets (page export vs `node bin/chain_blocks.mjs`)');
        timings.push(await checkParamSet(sandbox, ui, 'defaults', baseline, {}, [],
            targetHits, comparisonHits, tmpPath));
        timings.push(await checkParamSet(sandbox, ui, '{minPositive: 0.6, maxGap: 10}', baseline,
            { minPositive: 0.6, maxGap: 10 }, ['--min_identity', '0.6', '--max_gap', '10'],
            targetHits, comparisonHits, tmpPath));
        timings.push(await checkParamSet(sandbox, ui, '{maxGap: 50, maxHitRank: 1}', baseline,
            { maxGap: 50, maxHitRank: 1 }, ['--max_gap', '50', '--max_hit_rank', '1'],
            targetHits, comparisonHits, tmpPath));
        timings.push(await checkParamSet(sandbox, ui, '{minBlock: 15}', baseline,
            { minBlock: 15 }, ['--min_block', '15'], targetHits, comparisonHits, tmpPath));
    } finally {
        try { unlinkSync(tmpPath); } catch { /* best effort */ }
    }

    // The plan's <=150 ms target is specifically for the axolotl page --
    // reported for every page, but only asserted (as a failure) there, since
    // a closer/more homeologous pair (e.g. an allotetraploid's two
    // subgenomes) can legitimately cost more in the two self-comparison
    // chains every request also runs (see SYN.runChainers).
    const maxMs = Math.max(...timings);
    const isAxolotl = /axolotl/i.test(pagePath);
    console.log(`\nPer-request time (params changed -> sources updated): ${timings.map((t) => t.toFixed(1)).join(', ')} ms `
        + `(max ${maxMs.toFixed(1)} ms${isAxolotl ? ', target <= 150 ms' : ' -- <=150ms target applies to the axolotl page only'})`);
    if (isAxolotl && maxMs > 150) {
        fail(`slowest request (${maxMs.toFixed(1)} ms) exceeds the 150 ms target`);
    }

    if (process.exitCode) {
        console.error('\nRESULT: FAIL');
    } else {
        console.log('\nRESULT: PASS');
    }
}

main().catch((e) => {
    console.error('FAIL: uncaught error:', e);
    process.exitCode = 1;
});
