// Nextflow's tasks -> the progress screen's plain-language steps. The three
// before "prepare" appear only when the run does them: looking up the
// organism (--taxid given), getting the reference and finding the proteome
// (each unless it's the user's file).
//
// Steps can overlap: the target is prepared while NCBI is still searched.

const STEPS = [
    { id: 'lookup', discovery: true, match: /^RESOLVE_TAXONOMY:/ },
    { id: 'reference', discovery: true,
        match: /^(FIND_REFERENCE_ASSEMBLY:.*|DOWNLOAD_GENOME|RESOLVE_REFERENCE_TAXONOMY:.*|FETCH_ACCESSION_SUMMARY|DESCRIBE_REFERENCE_ACCESSION)$/ },
    { id: 'proteome', discovery: true, match: /^(FIND_PROTEOME_ASSEMBLY:.*|DOWNLOAD_PROTEIN)$/ },
    { id: 'prepare', match: /^(RENAME_SEQUENCES|MINIPROT_INDEX|SPLIT_GENOME)$/ },
    { id: 'align', match: /^(MINIPROT_ALIGN|MINIPROT_ALIGN_CHUNK|MERGE_MINIPROT_GFF|RENAME_GFF)$/ },
    { id: 'synteny', match: /^BUILD_SYNTENY:/ },
    { id: 'plot', match: /^PLOT_SYNTENY:/ },
];

const ACTIVE = new Set(['NEW', 'SUBMITTED', 'RUNNING']);

// "RESOLVE_TAXONOMY:PARSE_LINEAGE (taxonomy.json)" -> "RESOLVE_TAXONOMY:PARSE_LINEAGE"
const processOf = (task) => (task.process || task.name || '').replace(/\s*\(.*\)$/, '');

function stepOf(task) {
    const p = processOf(task);
    const s = STEPS.find((st) => st.match.test(p));
    return s ? s.id : null;
}

// which of the first three steps the form's run has
function planOf(form) {
    return { lookup: !!form.taxid, reference: !form.reference, proteome: !form.proteome };
}

// tasks: iterable of { process|name, status }; plan: planOf(form);
// outcome: null while running, else 'finished'|'failed'|'cancelled'.
// Returns [{ id, state }], state one of pending|running|done|failed|stopped.
function computeSteps(tasks, { plan, outcome = null }) {
    const steps = STEPS.filter((s) => !s.discovery || plan[s.id]);
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
        // done only once a later step has started too: a step's tasks don't
        // all exist up front (the reference is indexed after its download)
        else if (states.length && (outcome === 'finished' || later.some(started))) state = 'done';
        else if (states.length) state = outcome ? 'stopped' : 'running';
        if (outcome === 'finished' && state === 'pending') state = 'done';
        return { id: s.id, state };
    });
}

module.exports = { stepOf, computeSteps, planOf };
