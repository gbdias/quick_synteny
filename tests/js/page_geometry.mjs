#!/usr/bin/env node
// Compares the client-side ring/gap/dotplot geometry an old (pre-plan-B)
// synteny page baked into static Bokeh sources against what the new page's
// SHARED_JS computes for the same underlying data. Both pages are generated
// from the SAME input files (see the plan's Acceptance section), so any
// difference here is either a real regression or a genuine, deliberate
// change (e.g. the ribbon draw-order fix, or fewer gap-wedge vertices).
//
// Usage: node tests/js/page_geometry.mjs <old.html> <new.html>
//
// Both pages' embedded <script> (SHARED_JS followed by `SYN.data = {...}`)
// is pulled out with a plain string/brace scan (not a full HTML/JS parser --
// good enough since we control exactly how build_page() emits it) and run
// in a node:vm context with a `window` stub, so `window.SYN = window.SYN
// || {}` at the top of SHARED_JS bootstraps the same way it does in a real
// page (window IS the sandbox's own global object, exactly like a browser
// tab, so a bare top-level `SYN` reference resolves the same object).
//
// The old page's `format_link_label`/`format_gap_label` (Python) used single
// quotes in their HTML attributes; SHARED_JS's `SYN.formatLinkLabel`/
// `SYN.formatGapLabel` (both pages) use double quotes -- a pre-existing,
// harmless discrepancy that predates this plan. Labels are normalised
// (single quotes -> double, digit-group separators stripped) before
// comparison so that doesn't register as a mismatch, alongside Python's
// `{:,}` vs JS's locale-dependent `toLocaleString()`.

import { readFileSync } from 'node:fs';
import vm from 'node:vm';

const EPS = 1e-9;

function fail(msg) {
    console.error(`FAIL: ${msg}`);
    process.exitCode = 1;
}

function assert(cond, msg) {
    if (!cond) fail(msg);
}

// ---- extraction: pull SHARED_JS + SYN.data out of a generated page ----

function lastScriptBody(html, label) {
    const matches = [...html.matchAll(/<script[^>]*>([\s\S]*?)<\/script>/g)];
    if (matches.length === 0) throw new Error(`${label}: no <script> block found`);
    return matches[matches.length - 1][1];
}

// Scans forward from a leading '{' to find the matching '}', respecting
// (minimally escaped) JSON string contents -- enough to find the boundary
// of `json.dumps(syn_data)`'s output without a full JSON tokenizer.
function extractBalancedObject(text) {
    if (text[0] !== '{') throw new Error(`expected '{' at start of JSON, got: ${text.slice(0, 40)}`);
    let depth = 0, inString = false, escape = false;
    for (let i = 0; i < text.length; i++) {
        const c = text[i];
        if (inString) {
            if (escape) escape = false;
            else if (c === '\\') escape = true;
            else if (c === '"') inString = false;
            continue;
        }
        if (c === '"') { inString = true; continue; }
        if (c === '{') depth++;
        else if (c === '}') {
            depth--;
            if (depth === 0) return text.slice(0, i + 1);
        }
    }
    throw new Error('unbalanced JSON object in injected script');
}

// Returns {sharedJS, data} for a generated page: `data` is the parsed
// SYN.data object, `sharedJS` is everything before it (function/constant
// definitions only -- see build_page()'s own `injected` string).
function extractPage(path) {
    const html = readFileSync(path, 'utf8');
    const body = lastScriptBody(html, path);
    const marker = '\nSYN.data = ';
    const idx = body.indexOf(marker);
    if (idx === -1) throw new Error(`${path}: 'SYN.data = ' marker not found in injected script`);
    const sharedJS = body.slice(0, idx);
    const jsonText = extractBalancedObject(body.slice(idx + marker.length));
    return { sharedJS, data: JSON.parse(jsonText) };
}

// ---- sandbox: run SHARED_JS in node:vm with a stub `window` ----

function makeSandbox(sharedJS, data) {
    const sandbox = {};
    sandbox.window = sandbox; // window === globalThis, same as a real page
    sandbox.console = console;
    vm.createContext(sandbox);
    vm.runInContext(sharedJS, sandbox, { filename: 'SHARED_JS' });
    sandbox.SYN.data = data;
    return sandbox;
}

// Minimal stand-ins for the ColumnDataSource/figure args SYN.applyDotplotOrder
// writes into -- only `.data =` and `.change.emit()` (and fig.x_range/
// y_range.start/end) are ever touched, so a plain object with those shapes
// is enough; no real Bokeh runtime needed.
function stubSource() {
    return { data: {}, selected: { indices: [] }, change: { emit() {} } };
}

function stubDotplotSources() {
    return {
        query: stubSource(), subject: stubSource(), grid: stubSource(),
        queryLabel: stubSource(), subjectLabel: stubSource(),
        segment: stubSource(), gap: stubSource(),
        fig: { x_range: { start: 0, end: 0 }, y_range: { start: 0, end: 0 } },
    };
}

// ---- label/record normalisation ----

function normalizeLabel(s) {
    // '  -> "  (pre-existing Python/JS quoting difference, see module header)
    // strip digit-group separators (Python's `{:,}` vs JS's locale-dependent
    // toLocaleString()) -- safe here since these labels never contain a
    // legitimate comma otherwise (see format_link_label/formatLinkLabel)
    return String(s).replace(/'/g, '"').replace(/,/g, '');
}

// Splits a (normalised) label into its non-numeric template and the numbers
// embedded in it, in order.
function tokenizeLabel(s) {
    const numbers = [];
    const template = s.replace(/\d+\.\d+|\d+/g, (m) => { numbers.push(m); return '\u0000'; });
    return { template, numbers };
}

// Structural label equality: the non-numeric text must match exactly (so a
// wrong chromosome, orientation, or template change still fails), integer
// fields (coordinates, scores) must match exactly, but a one-decimal field
// (avg identity %, anchor density) is allowed to differ by up to 0.1 -- an
// exact tie (e.g. mean_identity=0.9625 -> 96.25%) renders as "96.2" under
// Python's f"{:.1f}" (round-half-to-even) but "96.3" under JS's
// Number.toFixed(1) (round-half-up per the ECMA-262 spec), a pre-existing,
// harmless Python/JS formatting difference that predates this plan (any
// client-side rebuild of a ribbon -- a reorder, a filter, a self-links
// toggle -- already re-rendered its label through SYN.formatLinkLabel, so
// this same tie-breaking difference already existed between the page's
// initial (Python) and post-interaction (JS) labels before this change).
function labelsEquivalent(a, b) {
    const ta = tokenizeLabel(a), tb = tokenizeLabel(b);
    if (ta.template !== tb.template) return false;
    if (ta.numbers.length !== tb.numbers.length) return false;
    for (let i = 0; i < ta.numbers.length; i++) {
        const na = ta.numbers[i], nb = tb.numbers[i];
        if (na === nb) continue;
        if (!na.includes('.') && !nb.includes('.')) return false; // integers must match exactly
        if (Math.abs(parseFloat(na) - parseFloat(nb)) > 0.1 + EPS) return false;
    }
    return true;
}

function ribbonSortKey(r) {
    return [r.is_homeolog ? 1 : 0, r.score, normalizeLabel(r.label)];
}

function compareKeys(a, b) {
    for (let i = 0; i < a.length; i++) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

function sortBy(list, keyFn) {
    return list.map((r) => [keyFn(r), r]).sort((a, b) => compareKeys(a[0], b[0])).map((p) => p[1]);
}

function numArraysClose(a, b, eps, path) {
    if (!Array.isArray(a) || !Array.isArray(b)) return `${path}: not both arrays`;
    if (a.length !== b.length) return `${path}: length ${a.length} != ${b.length}`;
    for (let i = 0; i < a.length; i++) {
        if (Math.abs(a[i] - b[i]) > eps) {
            return `${path}[${i}]: ${a[i]} vs ${b[i]} (diff ${Math.abs(a[i] - b[i])})`;
        }
    }
    return null;
}

// ---- the four checks ----

function checkRibbons(oldSandbox, newSandbox, label) {
    const oldRibbons = oldSandbox.SYN.data.ribbons;
    const k = newSandbox.SYN.data.palette.length;
    const newLayout = newSandbox.SYN.buildRingLayout(newSandbox.SYN.data.queryNames, newSandbox.SYN.data.subjectNames, k);
    const newRibbons = newLayout.ribbons;

    assert(Array.isArray(oldRibbons) && oldRibbons.length > 0, `${label}: old page has no SYN.data.ribbons`);
    assert(newRibbons.length === oldRibbons.length,
        `${label}: ribbon count mismatch: old=${oldRibbons.length} new=${newRibbons.length}`);

    const oldSorted = sortBy(oldRibbons, ribbonSortKey);
    const newSorted = sortBy(newRibbons, ribbonSortKey);
    let mismatches = 0;
    for (let i = 0; i < oldSorted.length && i < newSorted.length; i++) {
        const o = oldSorted[i], n = newSorted[i];
        const oLabel = normalizeLabel(o.label), nLabel = normalizeLabel(n.label);
        if (!labelsEquivalent(oLabel, nLabel)) { fail(`${label}: ribbon[${i}] label mismatch:\n  old=${oLabel}\n  new=${nLabel}`); mismatches++; continue; }
        if (o.score !== n.score) { fail(`${label}: ribbon[${i}] (${oLabel}) score ${o.score} != ${n.score}`); mismatches++; }
        if (Boolean(o.is_homeolog) !== Boolean(n.is_homeolog)) { fail(`${label}: ribbon[${i}] (${oLabel}) is_homeolog ${o.is_homeolog} != ${n.is_homeolog}`); mismatches++; }
        if (o.palette_index !== n.palette_index) { fail(`${label}: ribbon[${i}] (${oLabel}) palette_index ${o.palette_index} != ${n.palette_index}`); mismatches++; }
        if (Math.abs(o.alpha - n.alpha) > EPS) { fail(`${label}: ribbon[${i}] (${oLabel}) alpha ${o.alpha} != ${n.alpha}`); mismatches++; }
        const xErr = numArraysClose(o.xs, n.xs, EPS, `ribbon[${i}] (${oLabel}) xs`);
        if (xErr) { fail(`${label}: ${xErr}`); mismatches++; }
        const yErr = numArraysClose(o.ys, n.ys, EPS, `ribbon[${i}] (${oLabel}) ys`);
        if (yErr) { fail(`${label}: ${yErr}`); mismatches++; }
    }
    if (mismatches === 0) console.log(`  ribbons: ${oldRibbons.length} records match (old SYN.data.ribbons vs new SYN.buildRingLayout)`);

    // score non-decreasing within each of the three groups, in that group
    // order, in the NEW output (see B1's ribbon-sort fix)
    const crossCount = [].concat(...Object.values(newSandbox.SYN.data.linksByQuery)).length;
    const targetHomeologCount = newSandbox.SYN.data.targetHomeologLinks.length;
    const referenceHomeologCount = newSandbox.SYN.data.referenceHomeologLinks.length;
    const groups = [
        ['cross', newRibbons.slice(0, crossCount)],
        ['targetHomeolog', newRibbons.slice(crossCount, crossCount + targetHomeologCount)],
        ['referenceHomeolog', newRibbons.slice(crossCount + targetHomeologCount, crossCount + targetHomeologCount + referenceHomeologCount)],
    ];
    for (const [name, group] of groups) {
        for (let i = 1; i < group.length; i++) {
            if (group[i].score < group[i - 1].score) {
                fail(`${label}: new ribbons group '${name}' not non-decreasing at index ${i}: ${group[i - 1].score} -> ${group[i].score}`);
            }
        }
    }
    console.log(`  ribbon group order: cross(${crossCount}) + targetHomeolog(${targetHomeologCount}) + referenceHomeolog(${referenceHomeologCount}) all score-ascending in new output`);
}

function checkGapRecords(oldSandbox, newSandbox, label) {
    const oldGaps = oldSandbox.SYN.data.gapRecords || [];
    const groupGap = newSandbox.SYN.data.groupGap;
    const subjectOffsets = newSandbox.SYN.computeCircularOffsets(
        newSandbox.SYN.data.subjectNames, newSandbox.SYN.data.subjectSizes, Math.PI - groupGap, groupGap);
    const queryOffsets = newSandbox.SYN.computeCircularOffsets(
        newSandbox.SYN.data.queryNames, newSandbox.SYN.data.querySizes, Math.PI + groupGap, 2 * Math.PI - groupGap);
    const newGaps = newSandbox.SYN.buildGapRecords(queryOffsets, subjectOffsets, 48);

    if (oldGaps.length === 0 && newGaps.length === 0) {
        console.log('  gaps: no gap data in either page (no --target_gaps/--comparison_gaps given) -- skipped');
        return;
    }
    assert(newGaps.length === oldGaps.length,
        `${label}: gap record count mismatch: old=${oldGaps.length} new=${newGaps.length}`);

    const keyFn = (g) => [normalizeLabel(g.label)];
    const oldSorted = sortBy(oldGaps, keyFn);
    const newSorted = sortBy(newGaps, keyFn);
    let mismatches = 0;
    for (let i = 0; i < oldSorted.length && i < newSorted.length; i++) {
        const o = oldSorted[i], n = newSorted[i];
        const oLabel = normalizeLabel(o.label), nLabel = normalizeLabel(n.label);
        if (!labelsEquivalent(oLabel, nLabel)) { fail(`${label}: gap[${i}] label mismatch:\n  old=${oLabel}\n  new=${nLabel}`); mismatches++; continue; }
        const xErr = numArraysClose(o.xs, n.xs, EPS, `gap[${i}] (${oLabel}) xs`);
        if (xErr) { fail(`${label}: ${xErr}`); mismatches++; }
        const yErr = numArraysClose(o.ys, n.ys, EPS, `gap[${i}] (${oLabel}) ys`);
        if (yErr) { fail(`${label}: ${yErr}`); mismatches++; }
    }
    if (mismatches === 0) console.log(`  gaps: ${oldGaps.length} records match (old SYN.data.gapRecords vs new SYN.buildGapRecords(..., 48))`);
}

function checkDpCellAt(newSandbox, label) {
    const data = newSandbox.SYN.data;
    const k = data.palette.length;
    const queryOrder = data.queryNames, subjectOrder = data.subjectNames;
    const sources = stubDotplotSources();
    newSandbox.SYN.applyDotplotOrder(queryOrder, subjectOrder, k, data.mbaMin, sources);

    const state = newSandbox.SYN.state;
    let checked = 0;
    for (const qName of queryOrder) {
        const qx0 = state.dpQueryOffsets[qName], qsize = data.querySizes[qName];
        const x = qx0 + qsize / 2;
        for (const sName of subjectOrder) {
            const sy0 = state.dpSubjectOffsets[sName], ssize = data.subjectSizes[sName];
            const y = sy0 + ssize / 2;
            const hit = newSandbox.SYN.dpCellAt(x, y);
            if (!hit || hit.targetName !== qName || hit.subjectName !== sName) {
                fail(`${label}: dpCellAt(${x}, ${y}) (center of ${qName}x${sName}) returned ${JSON.stringify(hit)}, expected {${qName}, ${sName}}`);
            }
            checked++;
        }
    }
    console.log(`  dpCellAt: ${checked} (target, reference) cell centres hit exactly (${queryOrder.length} x ${subjectOrder.length})`);

    // null just outside each edge, and at negative x/y
    const midX = state.dpTotalX / 2, midY = state.dpTotalY / 2;
    const edgeCases = [
        ['x < 0', -1, midY],
        ['y < 0', midX, -1],
        ['just past right edge', state.dpTotalX + 1e-6, midY],
        ['just past top edge', midX, state.dpTotalY + 1e-6],
    ];
    for (const [name, x, y] of edgeCases) {
        const hit = newSandbox.SYN.dpCellAt(x, y);
        assert(hit === null, `${label}: dpCellAt(${x}, ${y}) [${name}] returned ${JSON.stringify(hit)}, expected null`);
    }
    console.log('  dpCellAt: null outside the grid on all four edges (x<0, y<0, x>=totalX, y>=totalY)');
}

function checkDotplotSegments(oldSandbox, newSandbox, label) {
    const oldData = oldSandbox.SYN.data, newData = newSandbox.SYN.data;
    const k = newData.palette.length;
    const minScore = newData.mbaMin;

    // "recompute it with the old page's JS on the old data" (see plan) --
    // the old page's own SYN.buildDotplotSegmentsForLayout is unchanged by
    // this stream, so this reproduces exactly what Python's build_dotplot_
    // sources() originally baked into dp_seg_src for the default order.
    const oldQueryOffsets = oldSandbox.SYN.computeOffsets(oldData.queryNames, oldData.querySizes);
    const oldSubjectOffsets = oldSandbox.SYN.computeOffsets(oldData.subjectNames, oldData.subjectSizes);
    const oldSeg = oldSandbox.SYN.buildDotplotSegmentsForLayout(oldQueryOffsets, oldSubjectOffsets, minScore, k);

    const newQueryOffsets = newSandbox.SYN.computeOffsets(newData.queryNames, newData.querySizes);
    const newSubjectOffsets = newSandbox.SYN.computeOffsets(newData.subjectNames, newData.subjectSizes);
    const newSeg = newSandbox.SYN.buildDotplotSegmentsForLayout(newQueryOffsets, newSubjectOffsets, minScore, k);

    const n = oldSeg.xs.length;
    assert(newSeg.xs.length === n, `${label}: dotplot segment count mismatch: old=${n} new=${newSeg.xs.length}`);

    const toRecords = (seg) => seg.xs.map((xs, i) => ({
        xs, ys: seg.ys[i], line_color: seg.line_color[i], alpha: seg.alpha[i],
        label: seg.label[i], palette_index: seg.palette_index[i],
    }));
    const keyFn = (r) => [r.palette_index, normalizeLabel(r.label)];
    const oldSorted = sortBy(toRecords(oldSeg), keyFn);
    const newSorted = sortBy(toRecords(newSeg), keyFn);
    let mismatches = 0;
    for (let i = 0; i < oldSorted.length && i < newSorted.length; i++) {
        const o = oldSorted[i], n2 = newSorted[i];
        const oLabel = normalizeLabel(o.label), nLabel = normalizeLabel(n2.label);
        if (!labelsEquivalent(oLabel, nLabel)) { fail(`${label}: dp segment[${i}] label mismatch:\n  old=${oLabel}\n  new=${nLabel}`); mismatches++; continue; }
        if (o.palette_index !== n2.palette_index) { fail(`${label}: dp segment[${i}] (${oLabel}) palette_index ${o.palette_index} != ${n2.palette_index}`); mismatches++; }
        if (o.line_color !== n2.line_color) { fail(`${label}: dp segment[${i}] (${oLabel}) line_color ${o.line_color} != ${n2.line_color}`); mismatches++; }
        if (Math.abs(o.alpha - n2.alpha) > EPS) { fail(`${label}: dp segment[${i}] (${oLabel}) alpha ${o.alpha} != ${n2.alpha}`); mismatches++; }
        const xErr = numArraysClose(o.xs, n2.xs, EPS, `dp segment[${i}] (${oLabel}) xs`);
        if (xErr) { fail(`${label}: ${xErr}`); mismatches++; }
        const yErr = numArraysClose(o.ys, n2.ys, EPS, `dp segment[${i}] (${oLabel}) ys`);
        if (yErr) { fail(`${label}: ${yErr}`); mismatches++; }
    }
    if (mismatches === 0) console.log(`  dotplot segments: ${n} records match (old page's own SYN.buildDotplotSegmentsForLayout on old data vs new page's on new data)`);
}

function main() {
    const [, , oldPath, newPath] = process.argv;
    if (!oldPath || !newPath) {
        console.error('usage: node tests/js/page_geometry.mjs <old.html> <new.html>');
        process.exit(2);
    }

    const oldPage = extractPage(oldPath);
    const newPage = extractPage(newPath);
    const oldSandbox = makeSandbox(oldPage.sharedJS, oldPage.data);
    const newSandbox = makeSandbox(newPage.sharedJS, newPage.data);

    console.log(`Comparing ${oldPath} (old) vs ${newPath} (new)`);

    console.log('- ribbons (SYN.buildRingLayout vs old SYN.data.ribbons)');
    checkRibbons(oldSandbox, newSandbox, oldPath);

    console.log('- gap records (SYN.buildGapRecords(..., 48) vs old SYN.data.gapRecords)');
    checkGapRecords(oldSandbox, newSandbox, oldPath);

    console.log('- dotplot cell hit-testing (SYN.dpCellAt)');
    checkDpCellAt(newSandbox, oldPath);

    console.log('- dotplot segments (SYN.buildDotplotSegmentsForLayout, old vs new)');
    checkDotplotSegments(oldSandbox, newSandbox, oldPath);

    if (process.exitCode) {
        console.error('\nRESULT: FAIL');
    } else {
        console.log('\nRESULT: PASS');
    }
}

main();
