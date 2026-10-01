// A run as the app remembers it: what the sidebar lists and the run screen
// draws. The main process builds it from the runner's events (lib/runner.js)
// and keeps it in the history (lib/history.js); the page only renders it.
//
//   { id, status: running|finished|failed|cancelled, startedAt, finishedAt,
//     discovery, form, outdir, title, inputs: { assembly, reference, proteome },
//     targetSpecies, found: { reference, proteome }, steps: [{ id, state }],
//     note, error, results: [html], summary, hasThumbnail }

const fs = require('node:fs');
const path = require('node:path');
const { computeSteps } = require('./stages');

// "Saccharomyces cerevisiae" -> "S. cerevisiae"; anything else unchanged
function shortSpecies(name) {
    const m = String(name || '').match(/^([A-Z])[a-z]+ ([a-z][\w-]*)$/);
    return m ? `${m[1]}. ${m[2]}` : name || '';
}

// "/x/GCF_000146045.2_R64_genomic.fna.gz" -> "GCF_000146045.2_R64_genomic"
function fileStem(p) {
    return p ? path.basename(p).replace(/\.gz$/, '').replace(/\.(fa|fasta|fna|fas|faa|fsa)$/i, '') : '';
}

function titleOf(r) {
    const target = shortSpecies(r.targetSpecies) || fileStem(r.form.assembly);
    const reference = r.discovery ? shortSpecies(r.found.reference?.species) || 'a relative' : fileStem(r.form.reference);
    return `${target} vs ${reference}`;
}

function newRecord(form, { targetSpecies = '' } = {}) {
    const discovery = !(form.reference && form.proteome);
    const r = {
        id: `${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 6)}`,
        status: 'running', startedAt: new Date().toISOString(), finishedAt: null,
        discovery, form, outdir: form.outdir, targetSpecies,
        inputs: { assembly: form.assembly, reference: form.reference || null, proteome: form.proteome || null },
        found: {}, steps: computeSteps([], { discovery }), note: '', error: null, results: [], summary: null,
        hasThumbnail: false,
    };
    r.title = titleOf(r);
    return r;
}

// One of the runner's events. `tasks` (id -> { process|name, status }) stays
// in memory with the run, since the record keeps only the steps. Returns
// whether the record changed.
function applyEvent(r, tasks, channel, payload) {
    switch (channel) {
    case 'run:task':
        tasks.set(payload.id, { process: payload.process, name: payload.name, status: payload.status });
        r.steps = computeSteps(tasks.values(), { discovery: r.discovery });
        return true;
    case 'run:found':
        r.found[payload.kind] = payload;
        r.title = titleOf(r);
        return true;
    case 'run:note':
        r.note = payload;
        return true;
    case 'run:error':
        r.error = payload;
        return true;
    default:
        return false;
    }
}

function finishRecord(r, tasks, { code, cancelled, results }) {
    r.status = cancelled ? 'cancelled' : code === 0 ? 'finished' : 'failed';
    r.finishedAt = new Date().toISOString();
    r.results = results || [];
    r.note = '';
    r.steps = computeSteps(tasks.values(), { discovery: r.discovery, outcome: r.status });
    if (r.status === 'finished') r.summary = readSummary(r.outdir);
    fillFound(r);
    r.title = titleOf(r);
}

// What discovery picked, from the selections the pipeline publishes to
// pipeline_info/ (FIND_REFERENCE_ASSEMBLY, FIND_PROTEOME_ASSEMBLY): the
// same facts as the log lines the runner parses while the run goes, for
// any of those it missed. Returns whether anything was added.
function fillFound(r) {
    if (!r.discovery) return false;
    let added = false;
    for (const kind of ['reference', 'proteome']) {
        if (r.found[kind]) continue;
        let sel;
        try { sel = JSON.parse(fs.readFileSync(path.join(r.outdir, 'pipeline_info', `${kind}_selection.json`), 'utf8')); } catch { continue; }
        if (!sel.accession) continue;
        r.found[kind] = { kind, accession: sel.accession, species: sel.organism_name || '', rank: sel.rank || '', rankName: sel.name || '',
            candidates: Number(sel.candidate_count) || 0, sameSpecies: !!sel.same_species_as_target };
        added = true;
    }
    if (added) r.title = titleOf(r);
    return added;
}

// the finished run's numbers: compute_alignment_stats.py's stats.json, and
// the blocks at the published min_block (links.tsv: a header, then a block a line)
function readSummary(outdir) {
    const dir = path.join(outdir, 'synteny');
    let names;
    try { names = fs.readdirSync(dir); } catch { return null; }
    const summary = {};
    const stats = names.find((n) => n.endsWith('.stats.json'));
    if (stats) {
        try { Object.assign(summary, JSON.parse(fs.readFileSync(path.join(dir, stats), 'utf8'))); } catch { /* no stats */ }
    }
    const links = names.find((n) => n.endsWith('.links.tsv') && !n.includes('slider'));
    if (links) {
        const lines = fs.readFileSync(path.join(dir, links), 'utf8').split('\n').filter((l) => l.trim());
        summary.blocks = Math.max(0, lines.length - 1);
    }
    return summary;
}

module.exports = { newRecord, applyEvent, finishRecord, fillFound, readSummary, shortSpecies, fileStem, titleOf };
