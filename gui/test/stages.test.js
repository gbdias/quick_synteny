// lib/stages.js: Nextflow tasks -> the progress screen's steps

const test = require('node:test');
const assert = require('node:assert/strict');
const { stepOf, computeSteps, planOf } = require('../lib/stages');

const states = (steps) => Object.fromEntries(steps.map((s) => [s.id, s.state]));
// your own reference and proteome, or everything left to NCBI
const OWN = planOf({ reference: '/r.fa', proteome: '/p.faa' });
const ALL = planOf({ taxid: '4932' });

test('every process the pipeline runs belongs to a step', () => {
    const processes = {
        'RESOLVE_TAXONOMY:FETCH_TAXONOMY_JSON': 'lookup', 'RESOLVE_TAXONOMY:PARSE_LINEAGE': 'lookup',
        'FIND_REFERENCE_ASSEMBLY:QUERY_GENOME_LADDER_REFERENCE': 'reference', DOWNLOAD_GENOME: 'reference',
        'FIND_PROTEOME_ASSEMBLY:SELECT_PROTEOME_ASSEMBLY': 'proteome', DOWNLOAD_PROTEIN: 'proteome',
        RENAME_SEQUENCES: 'prepare', MINIPROT_INDEX: 'prepare', SPLIT_GENOME: 'prepare',
        MINIPROT_ALIGN: 'align', MINIPROT_ALIGN_CHUNK: 'align', MERGE_MINIPROT_GFF: 'align', RENAME_GFF: 'align',
        'RESOLVE_REFERENCE_TAXONOMY:FETCH_TAXONOMY_JSON': 'reference', FETCH_ACCESSION_SUMMARY: 'reference',
        DESCRIBE_REFERENCE_ACCESSION: 'reference',
        'BUILD_SYNTENY:CHAIN_CROSS': 'synteny', 'PLOT_SYNTENY:RENDER_SYNTENY_INTERACTIVE': 'plot',
    };
    for (const [p, step] of Object.entries(processes)) assert.equal(stepOf({ process: p }), step, p);
    // cached tasks come from the log, as names with a tag
    assert.equal(stepOf({ name: 'RESOLVE_TAXONOMY:PARSE_LINEAGE (taxonomy.json)' }), 'lookup');
});

test('a run with your own files has no discovery steps', () => {
    assert.deepEqual(computeSteps([], { plan: OWN }).map((s) => s.id), ['prepare', 'align', 'synteny', 'plot']);
    assert.equal(computeSteps([], { plan: ALL }).length, 7);
});

test('the steps before "prepare" follow what the form leaves to NCBI', () => {
    const ids = (form) => computeSteps([], { plan: planOf(form) }).map((s) => s.id).slice(0, -4);
    assert.deepEqual(ids({ taxid: '4932' }), ['lookup', 'reference', 'proteome']);
    assert.deepEqual(ids({ reference: '/r.fa', proteome: '/p.faa' }), []);
    // your reference, the proteome found from your organism
    assert.deepEqual(ids({ reference: '/r.fa', taxid: '4932' }), ['lookup', 'proteome']);
    // a chosen reference with your proteome: no organism to look up
    assert.deepEqual(ids({ reference_accession: 'GCA_1.1', proteome: '/p.faa' }), ['reference']);
    assert.deepEqual(ids({ reference_taxid: '27291', taxid: '4932' }), ['lookup', 'reference', 'proteome']);
});

test('mid-run: finished steps the run has moved past are done, running ones running, the rest pending', () => {
    const tasks = [
        { process: 'RESOLVE_TAXONOMY:FETCH_TAXONOMY_JSON', status: 'COMPLETED' },
        { process: 'FIND_REFERENCE_ASSEMBLY:SELECT_REFERENCE_ASSEMBLY', status: 'COMPLETED' },
        { process: 'DOWNLOAD_GENOME', status: 'RUNNING' },
        { process: 'RENAME_SEQUENCES', status: 'COMPLETED' },
    ];
    assert.deepEqual(states(computeSteps(tasks, { plan: ALL })),
        { lookup: 'done', reference: 'running', proteome: 'pending', prepare: 'running', align: 'pending', synteny: 'pending', plot: 'pending' });
});

test('a step whose tasks have finished is done once a later step has started', () => {
    const tasks = [{ process: 'RENAME_SEQUENCES', status: 'COMPLETED' }, { process: 'MINIPROT_INDEX', status: 'COMPLETED' },
        { process: 'MINIPROT_ALIGN', status: 'RUNNING' }];
    assert.deepEqual(states(computeSteps(tasks, { plan: OWN })), { prepare: 'done', align: 'running', synteny: 'pending', plot: 'pending' });
    // nothing later yet: it may still get tasks (the reference's index), so it's still going
    assert.equal(states(computeSteps(tasks.slice(0, 2), { plan: OWN })).prepare, 'running');
});

test('a failed task fails its step; a cancelled run stops the running ones', () => {
    const tasks = [{ process: 'RENAME_SEQUENCES', status: 'COMPLETED' }, { process: 'MINIPROT_ALIGN', status: 'FAILED' }];
    assert.equal(states(computeSteps(tasks, { plan: OWN, outcome: 'failed' })).align, 'failed');
    const running = [{ process: 'MINIPROT_ALIGN', status: 'RUNNING' }];
    assert.equal(states(computeSteps(running, { plan: OWN, outcome: 'cancelled' })).align, 'stopped');
});

test('a finished run has every step done, cached ones included', () => {
    const tasks = [{ name: 'RENAME_SEQUENCES (target)', status: 'CACHED' }, { process: 'PLOT_SYNTENY:RENDER', status: 'COMPLETED' }];
    assert.ok(computeSteps(tasks, { plan: OWN, outcome: 'finished' }).every((s) => s.state === 'done'));
});
