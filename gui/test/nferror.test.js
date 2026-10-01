// lib/nferror.js on Nextflow's real console output (test/fixtures: captured
// from runs, with paths replaced)

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { ErrorCollector, withNetwork } = require('../lib/nferror');

function collect(fixture, code = 1) {
    const c = new ErrorCollector();
    fs.readFileSync(path.join(__dirname, 'fixtures', fixture), 'utf8').split('\n').forEach((l) => c.push(l));
    return c.result(code);
}

test('the pipeline rejecting its parameters: its own message, without the doubled ERROR:', () => {
    const e = collect('invalid-taxid.txt');
    assert.equal(e.title, 'The pipeline stopped before running');
    assert.equal(e.message, "--taxid must be an NCBI taxid (a positive whole number), got '0'");
    assert.equal(e.hint, undefined);
});

test('a failed step: its name, the error line, and its work folder', () => {
    const e = collect('unknown-taxid.txt');
    assert.equal(e.title, 'Step failed: RESOLVE_TAXONOMY:PARSE_LINEAGE (taxonomy.json)');
    assert.equal(e.message, '--taxid 999999999 is not a known NCBI taxid');
    assert.match(e.workdir, /^\/tmp\/run\/badtax\/work\/c5\//);
    assert.equal(e.hint, undefined);
});

test('the error line, not the boilerplate datasets prints after it', () => {
    const e = collect('download-stream-error.txt');
    assert.equal(e.message, 'Error: Download error: stream error: stream ID 5; INTERNAL_ERROR; received from peer');
});

test('NCBI unreachable: a plain message, the raw line kept as detail, and the network hint', () => {
    const e = collect('ncbi-offline.txt');
    assert.equal(e.message, "Couldn't reach NCBI (api.ncbi.nlm.nih.gov).");
    assert.match(e.detail, /no such host/);
    assert.match(e.hint, /network problem/);
});

test("a conda env that couldn't be built offline: named as such, retries collapsed", () => {
    const e = collect('conda-env-offline.txt');
    assert.equal(e.kind, 'conda-env');
    assert.match(e.title, /^Couldn't set up the tools for RENAME_SEQUENCES/);
    assert.match(e.message, /seems to be offline/);
    assert.match(e.detail, /Retrying in 2 seconds {2}\(×14\)/);
});

test('parameters nextflow_schema.json rejects: nf-schema\'s own list (stderr), not its generic ERROR line', () => {
    const e = collect('schema-invalid-taxid.txt');
    assert.equal(e.title, 'A parameter is not valid');
    // nf-schema prints a falsy value as "()"; dropped
    assert.equal(e.message, '--taxid: 0 is less than 1 (must be an NCBI taxid, a positive whole number)');
});

test('unrecognised parameters: listed although nf-schema prints no ERROR line, whatever the stream order', () => {
    const e = collect('schema-unrecognised.txt');
    assert.equal(e.title, 'A parameter is not valid');
    assert.equal(e.message, '--maxRank: genus');
});

test('no error block at all: the last lines of output', () => {
    const c = new ErrorCollector();
    ['Launching foo', 'Exception in thread main java.lang.OutOfMemoryError'].forEach((l) => c.push(l));
    const e = c.result(1);
    assert.equal(e.title, 'Nextflow exited with status 1');
    assert.equal(e.message, 'Exception in thread main java.lang.OutOfMemoryError');
});

test("withNetwork: offline turns a conda failure into the offline message, other failures get the hint", () => {
    const conda = { kind: 'conda-env', message: 'micromamba failed' };
    const step = { title: 'Step failed: X', message: 'x' };
    assert.match(withNetwork(conda, { online: false, conda: false }).message, /seems to be offline/);
    assert.equal(withNetwork(conda, { online: true, conda: true }).message, 'micromamba failed');
    assert.match(withNetwork(step, { online: false }).hint, /network problem/);
    assert.equal(withNetwork(step, { online: true, conda: true }).hint, undefined);
});
