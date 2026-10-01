// A run as the plain-language steps the progress screen lists, worked out
// from Nextflow's tasks (their process names and states). Up to three steps
// come before "prepare", each only when the run does it: looking up your
// genome's organism (when --taxid is given), getting the reference (unless
// it's your file: discovered, or the one you chose on NCBI), and finding
// the proteome (unless it's your file).
//
// Steps can overlap: the target genome is prepared while NCBI is still
// being searched for a reference, so two steps may be running at once.

const STEPS = [
    { id: 'lookup', discovery: true, match: /^RESOLVE_TAXONOMY:/ },
    { id: 'reference', discovery: true,
        match: /^(FIND_REFERENCE_ASSEMBLY:.*|DOWNLOAD_GENOME|RESOLVE_REFERENCE_TAXONOMY:.*|FETCH_ACCESSION_SUMMARY|DESCRIBE_REFERENCE_ACCESSION)$/ },
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

// Which of the first three steps a run has, from its form: { lookup,
// reference, proteome } (see above)
function planOf(form) {
    return { lookup: !!form.taxid, reference: !form.reference, proteome: !form.proteome };
}

// tasks: iterable of { process|name, status }; plan: planOf(form), or for
// records from before it existed `discovery` (all three steps, or none);
// outcome: null while running, else 'finished'|'failed'|'cancelled'.
// Returns [{ id, state }], state one of pending|running|done|failed|stopped.
function computeSteps(tasks, { plan, discovery, outcome = null }) {
    const has = plan || { lookup: !!discovery, reference: !!discovery, proteome: !!discovery };
    const steps = STEPS.filter((s) => !s.discovery || has[s.id]);
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

module.exports = { STEPS, stepOf, computeSteps, planOf };
