// lib/record.js and lib/history.js: what a run looks like to the sidebar

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const { newRecord, applyEvent, finishRecord, fillFound, shortSpecies, fileStem } = require('../lib/record');
const { parseFound } = require('../lib/runner');
const { History } = require('../lib/history');

const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'qs-record-'));
test.after(() => fs.rmSync(tmp, { recursive: true, force: true }));

test('names: short species, and file names without their FASTA extension', () => {
    assert.equal(shortSpecies('Saccharomyces cerevisiae'), 'S. cerevisiae');
    assert.equal(shortSpecies('Saccharomyces cerevisiae x Saccharomyces paradoxus'), 'Saccharomyces cerevisiae x Saccharomyces paradoxus');
    assert.equal(fileStem('/x/GCF_000146045.2_R64_genomic.fna.gz'), 'GCF_000146045.2_R64_genomic');
    assert.equal(fileStem('/x/proteins.faa'), 'proteins');
});

test('a discovery run is titled by species, and gains the reference once it is found', () => {
    const r = newRecord({ assembly: '/x/R64.fna', taxid: '4932', outdir: tmp }, { targetSpecies: 'Saccharomyces cerevisiae' });
    assert.equal(r.discovery, true);
    assert.equal(r.title, 'S. cerevisiae vs a relative');
    applyEvent(r, new Map(), 'run:found', { kind: 'reference', species: 'Saccharomyces pastorianus', accession: 'GCA_1' });
    assert.equal(r.title, 'S. cerevisiae vs S. pastorianus');
});

test('a run with your own files is titled by its files, and has no discovery steps', () => {
    const r = newRecord({ assembly: '/x/R64.fna', reference: '/x/pastorianus.fa', proteome: '/x/p.faa', outdir: tmp });
    assert.equal(r.discovery, false);
    assert.equal(r.title, 'R64 vs pastorianus');
    assert.deepEqual(r.steps.map((s) => s.id), ['prepare', 'align', 'synteny', 'plot']);
});

test('finishing reads the summary from the run folder', () => {
    const out = fs.mkdtempSync(path.join(tmp, 'run-'));
    fs.mkdirSync(path.join(out, 'synteny'));
    fs.writeFileSync(path.join(out, 'synteny', 'target.reference.stats.json'), JSON.stringify({ proteome_total: 10760, query_aligned: 10697 }));
    fs.writeFileSync(path.join(out, 'synteny', 'target.reference.links.tsv'), 'header\na\nb\nc\n');
    fs.writeFileSync(path.join(out, 'synteny', 'target.reference.slider_links.tsv'), 'header\na\nb\nc\nd\ne\n');
    const r = newRecord({ assembly: '/x/a.fa', reference: '/x/b.fa', proteome: '/x/p.faa', outdir: out });
    const tasks = new Map();
    applyEvent(r, tasks, 'run:task', { id: 1, process: 'RENAME_SEQUENCES', status: 'COMPLETED' });
    finishRecord(r, tasks, { code: 0, cancelled: false, results: ['/x.html'] });
    assert.equal(r.status, 'finished');
    assert.equal(r.summary.blocks, 3);
    assert.equal(r.summary.proteome_total, 10760);
    assert.ok(r.steps.every((s) => s.state === 'done'));
});

test('history: newest first, and a run cut off by quitting comes back as cancelled', () => {
    const dir = fs.mkdtempSync(path.join(tmp, 'hist-'));
    const h = new History(dir);
    const a = newRecord({ assembly: '/x/a.fa', taxid: '1', outdir: tmp });
    const b = newRecord({ assembly: '/x/b.fa', taxid: '1', outdir: tmp });
    h.put(a);
    h.put(b);
    assert.deepEqual(h.list().map((r) => r.id), [b.id, a.id]);
    h.save(true);
    const again = new History(dir);
    assert.equal(again.get(a.id).status, 'cancelled');
    assert.equal(again.get(a.id).interrupted, true);
    again.remove(a.id);
    assert.equal(new History(dir).get(a.id), null);
});

test('the discovery log lines parse, including ranks of several words', () => {
    const genus = parseFound('quick_synteny: reference accession=GCA_056824455.1 (Saccharomyces pastorianus) rank=genus (Saccharomyces) from 47 candidate(s)');
    assert.deepEqual(genus, { kind: 'reference', accession: 'GCA_056824455.1', species: 'Saccharomyces pastorianus', rank: 'genus',
        rankName: 'Saccharomyces', candidates: 47, sameSpecies: false });
    // a real run's: the reference was found at a "species group"
    const group = parseFound('quick_synteny: proteome accession=GCF_003285875.2 (Drosophila novamexicana) rank=species group (virilis group) annotated=true from 2 candidate(s)');
    assert.equal(group.rank, 'species group');
    assert.equal(group.rankName, 'virilis group');
    assert.equal(group.candidates, 2);
    assert.equal(parseFound('quick_synteny: reference accession=GCA_1 (A b) rank=species (A b) [SAME SPECIES AS TARGET] from 3 candidate(s)').sameSpecies, true);
});

test('a finished run fills in what was found from pipeline_info/, when the log lines were missed', () => {
    const out = fs.mkdtempSync(path.join(tmp, 'found-'));
    fs.mkdirSync(path.join(out, 'pipeline_info'));
    const sel = (accession, organism) => JSON.stringify({ accession, rank: 'species group', name: 'virilis group', organism_name: organism,
        candidate_count: 4, same_species_as_target: false });
    fs.writeFileSync(path.join(out, 'pipeline_info', 'reference_selection.json'), sel('GCA_030788265.1', 'Drosophila americana'));
    fs.writeFileSync(path.join(out, 'pipeline_info', 'proteome_selection.json'), sel('GCF_003285875.2', 'Drosophila novamexicana'));
    const r = newRecord({ assembly: '/x/virilis.fna.gz', taxid: '7244', outdir: out }, { targetSpecies: 'Drosophila virilis' });
    assert.equal(r.title, 'D. virilis vs a relative');
    finishRecord(r, new Map(), { code: 0, cancelled: false, results: [] });
    assert.equal(r.found.reference.species, 'Drosophila americana');
    assert.equal(r.found.reference.rankName, 'virilis group');
    assert.equal(r.found.proteome.accession, 'GCF_003285875.2');
    assert.equal(r.title, 'D. virilis vs D. americana');
    // what the log lines did give is kept
    const r2 = newRecord({ assembly: '/x/v.fna', taxid: '7244', outdir: out });
    applyEvent(r2, new Map(), 'run:found', { kind: 'reference', species: 'From the log', accession: 'X' });
    assert.equal(fillFound(r2), true);
    assert.equal(r2.found.reference.species, 'From the log');
    assert.equal(r2.found.proteome.species, 'Drosophila novamexicana');
});
