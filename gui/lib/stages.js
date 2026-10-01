// A run as the plain-language steps the progress screen lists, worked out
// from Nextflow's tasks (their process names and states). Discovery runs
// start with three steps of their own; a run with your own reference and
// proteome begins at "prepare".
//
// Steps can overlap: the target genome is prepared while NCBI is still
// being searched for a reference, so two steps may be running at once.

const STEPS = [
    { id: 'lookup', discovery: true, match: /^RESOLVE_TAXONOMY:/ },
    { id: 'reference', discovery: true, match: /^(FIND_REFERENCE_ASSEMBLY:.*|DOWNLOAD_GENOME)$/ },
    { id: 'proteome', discovery: true, match: /^(FIND_PROTEOME_ASSEMBLY:.*|DOWNLOAD_PROTEIN)$/ },
    { id: 'prepare', match: /^(RENAME_SEQUENCES|MINIPROT_INDEX|SPLIT_GENOME)$/ },
    { id: 'align', match: /^(MINIPROT_ALIGN|MINIPROT_ALIGN_CHUNK|MERGE_MINIPROT_GFF|RENAME_GFF)$/ },
    { id: 'synteny', match: /^BUILD_SYNTENY:/ },
    { id: 'plot', match: /^PYGENOMEVIZ_PLOT:/ },
];

const ACTIVE = new Set(['NEW', 'SUBMITTED', 'RUNNING']);

// "RESOLVE_TAXONOMY:PARSE_LINEAGE (taxonomy.json)" -> "RESOLVE_TAXONOMY:PARSE_LINEAGE"
const processOf = (task) => (task.process || task.name || '').replace(/\s*\(.*\)$/, '');

function stepOf(task) {
    const p = processOf(task);
    const s = STEPS.find((st) => st.match.test(p));
    return s ? s.id : null;
}

// tasks: iterable of { process|name, status }; discovery: whether the run
// searches NCBI; outcome: null while running, else 'finished'|'failed'|'cancelled'.
// Returns [{ id, state }], state one of pending|running|done|failed|stopped.
function computeSteps(tasks, { discovery, outcome = null }) {
    const steps = STEPS.filter((s) => discovery || !s.discovery);
    const byStep = new Map(steps.map((s) => [s.id, []]));
    for (const t of tasks) {
        const id = stepOf(t);
        if (byStep.has(id)) byStep.get(id).push(t.status);
    }
    const started = (id) => byStep.get(id).length > 0;
    return steps.map((s, i) => {
        const states = byStep.get(s.id);
        const later = steps.slice(i + 1).map((l) => l.id);
        let state = 'pending';
        if (states.some((st) => st === 'FAILED' || st === 'ABORTED')) state = 'failed';
        else if (states.some((st) => ACTIVE.has(st))) state = outcome ? 'stopped' : 'running';
        // done once its tasks have all finished and a later step has started:
        // a step's tasks don't all exist up front (the reference is indexed
        // only after it's downloaded), so "nothing running" alone isn't enough
        else if (states.length && (outcome === 'finished' || later.some(started))) state = 'done';
        else if (states.length) state = outcome ? 'stopped' : 'running';
        if (outcome === 'finished' && state === 'pending') state = 'done';
        return { id: s.id, state };
    });
}

module.exports = { STEPS, stepOf, computeSteps };
