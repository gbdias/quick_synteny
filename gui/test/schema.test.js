// nextflow_schema.json against everything else that names the pipeline's
// parameters: nextflow.config's params block, the form, and the runner. A
// parameter added, renamed or re-defaulted in one place and not the others
// fails here.

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { PARAMS } = require('../lib/runner');

const ROOT = path.join(__dirname, '..', '..');
const schema = JSON.parse(fs.readFileSync(path.join(ROOT, 'nextflow_schema.json'), 'utf8'));
const props = Object.assign({}, ...Object.values(schema.$defs).map((d) => d.properties));
const html = fs.readFileSync(path.join(__dirname, '..', 'renderer', 'index.html'), 'utf8');

// `name = value` lines of nextflow.config's params { } block
function configParams() {
    const config = fs.readFileSync(path.join(ROOT, 'nextflow.config'), 'utf8');
    const block = config.match(/^params \{([\s\S]*?)^\}/m)[1];
    const out = {};
    for (const m of block.matchAll(/^\s*(\w+)\s*=\s*([^\s/]+)/gm)) {
        const v = m[2];
        out[m[1]] = v === 'null' ? null : v === 'true' ? true : v === 'false' ? false
            : /^'.*'$/.test(v) ? v.slice(1, -1) : Number(v);
    }
    return out;
}

test('every schema group is in allOf, and every parameter has a description and a type', () => {
    const refs = schema.allOf.map((a) => a.$ref.replace('#/$defs/', ''));
    assert.deepEqual(refs.sort(), Object.keys(schema.$defs).sort());
    for (const [name, p] of Object.entries(props)) {
        assert.ok(p.description, `${name} has no description`);
        assert.ok(p.type, `${name} has no type`);
    }
});

test('nextflow.config and the schema list the same parameters, with the same defaults', () => {
    const config = configParams();
    assert.deepEqual(Object.keys(config).sort(), Object.keys(props).sort());
    for (const [name, value] of Object.entries(config)) {
        // false/null in nextflow.config mean "unset": the schema has no default
        const expected = value === null || value === false ? undefined : value;
        assert.deepEqual(props[name].default, expected, `default of ${name}`);
    }
});

test('every form field and every "?" names a schema parameter', () => {
    const notParams = new Set(['mode', 'engine', 'resume']);
    const fields = [...html.matchAll(/<(?:input|select)[^>]*\bname="(\w+)"/g)].map((m) => m[1]).filter((n) => !notParams.has(n));
    assert.ok(fields.length >= 9);
    for (const f of fields) assert.ok(props[f], `form field ${f} isn't in the schema`);
    const tips = [...html.matchAll(/class="help" data-param="(\w+)"/g)].map((m) => m[1]);
    assert.ok(tips.length >= 8);
    for (const t of tips) assert.ok(props[t] && (props[t].help_text || props[t].description), `tip ${t} has no text`);
});

test('no hidden parameter is in the form or passed by the runner', () => {
    // the chaining parameters: the result page re-chains with any values
    const hidden = Object.keys(props).filter((n) => props[n].hidden && n !== 'help');
    assert.deepEqual(hidden.sort(), ['max_gap', 'min_block', 'min_identity']);
    for (const n of hidden) {
        assert.ok(!html.includes(`name="${n}"`), `${n} is hidden but has a form field`);
        assert.ok(!PARAMS.includes(n), `${n} is hidden but the runner passes it`);
    }
});

test('every parameter the runner passes is in the schema', () => {
    for (const p of PARAMS) assert.ok(props[p], `the runner passes --${p}, which isn't in the schema`);
});

test("the app's run config pins every plugin: the weblog's and the pipeline's own", () => {
    // a plugins block in a -c config replaces nextflow.config's, so it must repeat them
    const runtime = require('../lib/runtime');
    runtime.configure({ pipelineDir: ROOT });
    const { guiConfigText } = require('../lib/runner');
    const text = guiConfigText({ engine: 'conda' }, 'http://127.0.0.1:1');
    const ids = [...text.matchAll(/id '([^']+)'/g)].map((m) => m[1]);
    assert.deepEqual(ids, [runtime.WEBLOG_PLUGIN, 'nf-schema@2.7.3']);
    assert.match(fs.readFileSync(path.join(ROOT, 'nextflow.config'), 'utf8'), /id 'nf-schema@2\.7\.3'/);
});
