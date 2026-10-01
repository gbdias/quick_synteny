// lib/stages.js: Nextflow tasks -> the progress screen's steps

const test = require('node:test');
const assert = require('node:assert/strict');
const { stepOf, computeSteps } = require('../lib/stages');

const states = (steps) => Object.fromEntries(steps.map((s) => [s.id, s.state]));

test('every process the pipeline runs belongs to a step', () => {
    const processes = {
        'RESOLVE_TAXONOMY:FETCH_TAXONOMY_JSON': 'lookup', 'RESOLVE_TAXONOMY:PARSE_LINEAGE': 'lookup',
        'FIND_REFERENCE_ASSEMBLY:QUERY_GENOME_LADDER_REFERENCE': 'reference', DOWNLOAD_GENOME: 'reference',
        'FIND_PROTEOME_ASSEMBLY:SELECT_PROTEOME_ASSEMBLY': 'proteome', DOWNLOAD_PROTEIN: 'proteome',
        RENAME_SEQUENCES: 'prepare', MINIPROT_INDEX: 'prepare', SPLIT_GENOME: 'prepare',
        MINIPROT_ALIGN: 'align', MINIPROT_ALIGN_CHUNK: 'align', MERGE_MINIPROT_GFF: 'align', RENAME_GFF: 'align',
        'BUILD_SYNTENY:CHAIN_CROSS': 'synteny', 'PYGENOMEVIZ_PLOT:RENDER_SYNTENY_INTERACTIVE': 'plot',
    };
    for (const [p, step] of Object.entries(processes)) assert.equal(stepOf({ process: p }), step, p);
    // cached tasks come from the log, as names with a tag
    assert.equal(stepOf({ name: 'RESOLVE_TAXONOMY:PARSE_LINEAGE (taxonomy.json)' }), 'lookup');
});

test('a run with your own files has no discovery steps', () => {
    assert.deepEqual(computeSteps([], { discovery: false }).map((s) => s.id), ['prepare', 'align', 'synteny', 'plot']);
    assert.equal(computeSteps([], { discovery: true }).length, 7);
});

test('mid-run: finished steps the run has moved past are done, running ones running, the rest pending', () => {
    const tasks = [
        { process: 'RESOLVE_TAXONOMY:FETCH_TAXONOMY_JSON', status: 'COMPLETED' },
        { process: 'FIND_REFERENCE_ASSEMBLY:SELECT_REFERENCE_ASSEMBLY', status: 'COMPLETED' },
        { process: 'DOWNLOAD_GENOME', status: 'RUNNING' },
        { process: 'RENAME_SEQUENCES', status: 'COMPLETED' },
    ];
    assert.deepEqual(states(computeSteps(tasks, { discovery: true })),
        { lookup: 'done', reference: 'running', proteome: 'pending', prepare: 'running', align: 'pending', synteny: 'pending', plot: 'pending' });
});

test('a step whose tasks have finished is done once a later step has started', () => {
    const tasks = [{ process: 'RENAME_SEQUENCES', status: 'COMPLETED' }, { process: 'MINIPROT_INDEX', status: 'COMPLETED' },
        { process: 'MINIPROT_ALIGN', status: 'RUNNING' }];
    assert.deepEqual(states(computeSteps(tasks, { discovery: false })), { prepare: 'done', align: 'running', synteny: 'pending', plot: 'pending' });
    // nothing later yet: it may still get tasks (the reference's index), so it's still going
    assert.equal(states(computeSteps(tasks.slice(0, 2), { discovery: false })).prepare, 'running');
});

test('a failed task fails its step; a cancelled run stops the running ones', () => {
    const tasks = [{ process: 'RENAME_SEQUENCES', status: 'COMPLETED' }, { process: 'MINIPROT_ALIGN', status: 'FAILED' }];
    assert.equal(states(computeSteps(tasks, { discovery: false, outcome: 'failed' })).align, 'failed');
    const running = [{ process: 'MINIPROT_ALIGN', status: 'RUNNING' }];
    assert.equal(states(computeSteps(running, { discovery: false, outcome: 'cancelled' })).align, 'stopped');
});

test('a finished run has every step done, cached ones included', () => {
    const tasks = [{ name: 'RENAME_SEQUENCES (target)', status: 'CACHED' }, { process: 'PYGENOMEVIZ_PLOT:RENDER', status: 'COMPLETED' }];
    assert.ok(computeSteps(tasks, { discovery: false, outcome: 'finished' }).every((s) => s.state === 'done'));
});
