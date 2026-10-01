// lib/inputs.js: FASTA sniffing, output-folder inspection, subfolder names

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const zlib = require('node:zlib');
const { sniffFasta, inspectOutdir, suggestRunName } = require('../lib/inputs');

const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'qs-inputs-'));
test.after(() => fs.rmSync(tmp, { recursive: true, force: true }));
const write = (name, data) => {
    const p = path.join(tmp, name);
    fs.writeFileSync(p, data);
    return p;
};

test('FASTA, gzipped FASTA, and anything else', () => {
    const fasta = '\n>chr1 test\nACGTACGT\n'.repeat(5000);
    assert.equal(sniffFasta(write('a.fa', fasta)), 'fasta');
    assert.equal(sniffFasta(write('a.fa.gz', zlib.gzipSync(fasta))), 'gzip-fasta');
    // only the first 64 KB of a large gzip is read: a truncated stream must still inflate
    assert.equal(sniffFasta(write('big.fa.gz', zlib.gzipSync(fasta.repeat(40)))), 'gzip-fasta');
    assert.equal(sniffFasta(write('notes.txt', 'just some notes\n')), 'other');
    assert.equal(sniffFasta(write('notes.txt.gz', zlib.gzipSync('just some notes\n'))), 'other');
    assert.equal(sniffFasta(write('broken.gz', Buffer.from([0x1f, 0x8b, 0, 1, 2, 3]))), 'other');
});

test('output folders: empty, a previous run, something else', () => {
    const empty = fs.mkdtempSync(path.join(tmp, 'empty-'));
    assert.deepEqual(inspectOutdir(empty), { path: empty, writable: true, empty: true, isRun: false });

    const run = fs.mkdtempSync(path.join(tmp, 'run-'));
    fs.mkdirSync(path.join(run, 'synteny'));
    assert.equal(inspectOutdir(run).isRun, true);

    const busy = fs.mkdtempSync(path.join(tmp, 'busy-'));
    fs.writeFileSync(path.join(busy, 'holiday.jpg'), '');
    fs.writeFileSync(path.join(busy, '.DS_Store'), '');
    assert.deepEqual(inspectOutdir(busy), { path: busy, writable: true, empty: false, isRun: false });

    const ds = fs.mkdtempSync(path.join(tmp, 'ds-'));
    fs.writeFileSync(path.join(ds, '.DS_Store'), '');
    assert.equal(inspectOutdir(ds).empty, true, 'a lone .DS_Store still counts as empty');
});

test('subfolder names: <target>_vs_<taxid>, filesystem-safe', () => {
    assert.equal(suggestRunName('/x/b_illinoisensis_ps_336_001.polished_assembly.fasta.gz', '4932'),
        'b_illinoisensis_ps_336_001.polished_assembly_vs_4932');
    assert.equal(suggestRunName('/x/my genome (v2).fa', ''), 'my_genome_v2_vs_reference');
    assert.equal(suggestRunName(null, '7227'), 'quick_synteny_vs_7227');
});
