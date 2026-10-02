// The window: the sidebar of runs, and one view at a time beside it -- New
// run (form.js), a run (its steps in plain language, what was found, and
// once finished its plot's thumbnail and numbers), first-launch setup, and
// Settings. Runs are records the main process keeps (lib/record.js); this
// only draws them.

const fileStem = (p) => baseName(p).replace(/\.gz$/, '').replace(/\.(fa|fasta|fna|fas|faa|fsa)$/i, '');
const num = (n) => Number(n).toLocaleString('en-US');
const cleanIpcError = (e) => String(e?.message || e).replace(/^Error invoking remote method '[^']+': (Error: )?/, '');

// ---------- views ----------

const VIEWS = { new: '#form', run: '#view-run', setup: '#view-setup', settings: '#view-settings' };

function showView(name) {
    state.view = name;
    for (const [v, sel] of Object.entries(VIEWS)) $(sel).hidden = v !== name;
    if (name !== 'run') state.runId = null;
    if (name === 'new') refreshForm();
    if (name === 'settings') renderSettings();
    renderSidebar();
}

function showRun(id) {
    state.runId = id;
    state.view = 'run';
    for (const [v, sel] of Object.entries(VIEWS)) $(sel).hidden = v !== 'run';
    const r = state.runs.find((x) => x.id === id);
    if (r) renderRun(r);
    renderSidebar();
}

// a form "Run again" filled and that wasn't started is that run's, not a
// new one: New run clears it (a form filled by hand is kept)
$('#btn-new').addEventListener('click', () => {
    if (prefilled) resetForm();
    showView('new');
});
$('#btn-settings').addEventListener('click', () => showView('settings'));
$('#btn-open-result').addEventListener('click', () => window.qs.pickResult().catch((e) => alert(cleanIpcError(e))));
$('#tools-status').addEventListener('click', () => showView(state.setup === 'idle' ? 'settings' : 'setup'));

// ---------- sidebar ----------

const STATUS_ICON = { running: ['spinner', 'spin'], finished: ['check-circle', ''], failed: ['alert', ''], cancelled: ['minus', ''] };

function runMeta(r) {
    if (r.status === 'running') return 'Running…';
    if (r.status === 'failed') return `Stopped: ${r.error?.message || 'an error'}`;
    if (r.status === 'cancelled') return r.interrupted ? 'Interrupted' : 'Cancelled';
    return fmtWhen(r.finishedAt || r.startedAt);
}

function renderSidebar() {
    $('#run-list').replaceChildren(...state.runs.map((r) => {
        const [ic, cls] = STATUS_ICON[r.status] || STATUS_ICON.cancelled;
        return h('button', {
            type: 'button', class: `run-item ${r.status} ${state.view === 'run' && state.runId === r.id ? 'selected' : ''}`,
            onclick: () => showRun(r.id), title: r.title,
        }, icon(ic, cls), h('span', { class: 'run-text' }, h('span', { class: 'run-name', text: r.title }), h('span', { class: 'run-meta', text: runMeta(r) })));
    }));
    $('#runs-empty').hidden = state.runs.length > 0;
}

function renderTools() {
    const el = $('#tools-status');
    const nf = state.env?.nextflow;
    let cls = '';
    let text = 'Checking tools…';
    if (state.setup === 'running') [cls, text] = ['busy', 'Setting up the tools…'];
    else if (state.setup === 'failed') [cls, text] = ['bad', 'Setup failed: retry'];
    else if (nf?.ok) [cls, text] = ['ok', 'Tools ready'];
    else if (nf) [cls, text] = ['bad', 'Tools not working'];
    el.className = `tools-status ${cls}`;
    el.querySelector('.txt').textContent = text;
    el.title = nf?.ok ? `Nextflow ${nf.version}` : nf?.detail || '';
}

// ---------- a run ----------

// each step's [to do, doing, done] wording
const STEP_TEXT = {
    lookup: ['Look up your organism', 'Looking up your organism', 'Looked up your organism'],
    reference: ['Find a reference genome', 'Finding a reference genome', 'Found a reference genome'],
    proteome: ['Find a proteome', 'Finding a proteome', 'Found a proteome'],
    prepare: ['Prepare the genomes', 'Preparing the genomes', 'Prepared the genomes'],
    align: ['Align the proteins to both genomes', 'Aligning the proteins to both genomes', 'Aligned the proteins to both genomes'],
    synteny: ['Find synteny blocks', 'Finding synteny blocks', 'Found synteny blocks'],
    plot: ['Draw the interactive plot', 'Drawing the interactive plot', 'Drew the interactive plot'],
};

// the reference step, when it isn't a search from your species outwards
const REFERENCE_TEXT = {
    accession: ['Get the assembly you chose', 'Getting the assembly you chose', 'Got the assembly you chose'],
    taxid: ['Find the best genome of your choice', 'Finding the best genome of your choice', 'Found the best genome of your choice'],
};

// where the reference and proteome come from: a file, an exact assembly, a
// taxid, or discovery (as lib/record.js; runs from before it say less)
function sourcesOf(r) {
    if (r.sources) return r.sources;
    return r.discovery === false ? { reference: 'file', proteome: 'file' } : { reference: 'discovered', proteome: 'discovered' };
}

// "virilis group (species group)", "Saccharomyces (genus)"
const foundAt = (f) => (f.rankName ? `${f.rankName} (${f.rank})` : f.rank);
// ", chromosome level" or ", complete genome", when the selection says
const levelOf = (f) => {
    const level = (f?.level || '').toLowerCase();
    return !level ? '' : level === 'complete genome' ? ', complete genome' : `, ${level} level`;
};

function stepDetail(r, step, showNote) {
    if (showNote && r.note) return r.note;
    if (step.state === 'failed') return 'This is where it stopped: see above.';
    const s = r.summary || {};
    const ref = r.found.reference;
    const prot = r.found.proteome;
    const sources = sourcesOf(r);
    switch (step.id) {
    case 'lookup': return r.targetSpecies || `NCBI taxid ${r.form.taxid}`;
    case 'reference':
        if (sources.reference === 'accession') return ref ? `${ref.species} (${ref.accession})${levelOf(ref)}` : `${r.form.reference_accession}, from NCBI`;
        if (sources.reference === 'taxid') {
            return ref ? `${ref.species} (${ref.accession}): the best of ${ref.candidates}${levelOf(ref)}`
                : `The best chromosome-level assembly of NCBI taxid ${r.form.reference_taxid}`;
        }
        return ref ? `${ref.species} (${ref.accession}), found at ${foundAt(ref)}: the best of ${ref.candidates}`
            + `${ref.sameSpecies ? ', of your own species' : ''}` : 'Searching NCBI, from your species outwards';
    case 'proteome': return prot ? `${prot.species} (${prot.accession})${s.proteome_total ? `, ${num(s.proteome_total)} proteins` : ''}`
        : sources.reference === 'file' ? 'An annotated assembly, from your species outwards' : 'An annotated assembly, preferably of the reference species';
    case 'prepare': return 'Naming the sequences and indexing both genomes';
    case 'align': return 'Where each protein sits in both genomes: these matches are the synteny anchors';
    case 'synteny': return s.blocks !== undefined ? `${num(s.blocks)} blocks` : 'Runs of anchors in the same order on both genomes';
    case 'plot': return step.state === 'done' ? 'Ready: Open plot shows it' : 'Opens by itself when it\'s ready';
    default: return '';
    }
}

function renderSteps(r) {
    const firstRunning = r.steps.findIndex((s) => s.state === 'running');
    $('#run-steps').replaceChildren(...r.steps.map((step, i) => {
        const [todo, doing, done] = (step.id === 'reference' && REFERENCE_TEXT[sourcesOf(r).reference]) || STEP_TEXT[step.id];
        const title = step.state === 'running' ? doing : step.state === 'done' ? done : todo;
        const ic = { done: icon('check'), running: icon('spinner', 'spin'), failed: icon('bang') }[step.state] || null;
        return h('li', { class: `step ${step.state}` },
            h('span', { class: 'step-icon' }, ic),
            h('span', { class: 'step-text' },
                h('span', { class: 'step-title', text: title }),
                h('span', { class: 'step-detail', text: stepDetail(r, step, i === firstRunning) }),
                step.state === 'running' && i === firstRunning ? h('span', { class: 'bar' }) : null));
    }));
}

function foundItem(kind, main, sub, species = false) {
    return h('div', { class: 'found-item' }, h('span', { class: 'found-kind', text: kind }),
        h('span', { class: `found-main ${species ? 'species' : ''}`, text: main }), sub ? h('span', { class: 'found-sub', text: sub }) : null);
}

function renderFound(r) {
    const s = r.summary || {};
    const ref = r.found.reference;
    const prot = r.found.proteome;
    const proteins = s.proteome_total ? `${num(s.proteome_total)} proteins` : null;
    $('#found-label').textContent = r.discovery ? 'What we found' : 'What you gave';
    // not found yet: being searched for, still to come, or never reached
    const pending = (id) => {
        const st = r.steps.find((s) => s.id === id)?.state;
        if (r.status === 'running') return [st === 'running' ? 'Searching…' : 'Not yet', null];
        return ['—', 'The run stopped before this'];
    };
    const sources = sourcesOf(r);
    const refHow = {
        discovered: (f) => `found at ${foundAt(f)}, the best of ${f.candidates}`,
        taxid: (f) => `the best of ${f.candidates} for your choice${levelOf(f)}`,
        accession: (f) => `your choice${levelOf(f)}`,
    }[sources.reference];
    const reference = sources.reference === 'file'
        ? foundItem('Reference genome', fileStem(r.inputs.reference), baseName(r.inputs.reference))
        : foundItem('Reference genome', ...(ref ? [ref.species, `${ref.accession} · ${refHow(ref)}`] : pending('reference')), !!ref);
    const proteome = sources.proteome === 'file'
        ? foundItem('Proteome', proteins || fileStem(r.inputs.proteome), baseName(r.inputs.proteome))
        : foundItem('Proteome', ...(prot || proteins ? [proteins || prot.species, prot ? `${prot.accession}${proteins ? ` · ${prot.species}` : ''}` : null]
            : pending('proteome')), !proteins && !!prot);
    const genome = foundItem('Your genome', r.targetSpecies || fileStem(r.inputs.assembly), baseName(r.inputs.assembly), !!r.targetSpecies);
    // what came from NCBI first
    const items = r.discovery ? [reference, proteome, genome] : [genome, reference, proteome];
    items.push(foundItem('Saved in', baseName(r.outdir), tildify(r.outdir)));
    $('#run-found').replaceChildren(...items);
}

function renderHead(r) {
    $('#run-title').textContent = r.title;
    const chips = [];
    if (r.status === 'running') chips.push(h('span', { class: 'chip busy', text: 'Running' }), h('span', { class: 'chip', text: `Started ${fmtWhen(r.startedAt).replace(/^T/, 't').replace(/^Y/, 'y')}` }));
    else if (r.status === 'finished') chips.push(h('span', { class: 'chip ok', text: `Finished ${fmtWhen(r.finishedAt).replace(/^T/, 't').replace(/^Y/, 'y')}` }));
    else if (r.status === 'failed') chips.push(h('span', { class: 'chip bad', text: 'Stopped with an error' }));
    else chips.push(h('span', { class: 'chip', text: r.interrupted ? 'Interrupted: the app was closed' : 'Cancelled' }));
    if (r.finishedAt) chips.push(h('span', { class: 'chip', text: fmtDuration(new Date(r.finishedAt) - new Date(r.startedAt)) }));
    const sources = sourcesOf(r);
    const origin = sources.reference === 'file'
        ? sources.proteome === 'file' ? 'Your own reference and proteome' : 'Your reference, proteome from NCBI'
        : sources.reference === 'discovered' ? 'Reference found on NCBI' : 'Reference you chose on NCBI';
    chips.push(h('span', { class: 'chip', text: origin }));
    $('#run-chips').replaceChildren(...chips);

    const btn = (text, ic, onclick, cls = '') => h('button', { type: 'button', class: `btn ${cls}`, onclick }, ic ? icon(ic) : null, text);
    const remove = btn('Remove', null, async () => {
        if (!confirm('Remove this run from the list? Its folder and results stay where they are.')) return;
        await window.qs.removeRun(r.id).catch((e) => alert(cleanIpcError(e)));
        showView('new');
    }, 'ghost');
    // the same form either way: a finished run is run again, a stopped one
    // is mended first
    const again = (cls) => btn('Run again…', 'again', () => editRun(r), cls);
    const edit = (cls) => btn('Edit run…', 'edit', () => editRun(r), cls);
    const folder = btn('Show folder', 'folder', () => window.qs.reveal(r.outdir));
    let actions;
    if (r.status === 'running') actions = [btn('Cancel run', 'stop', () => window.qs.cancelRun())];
    else if (r.status === 'finished') actions = [remove, again(''), folder, btn('Open plot', 'open', () => r.results.forEach((f) => window.qs.openResult(f)), 'primary')];
    // with an error, its box has the primary Edit run… (renderError)
    else actions = [remove, folder, edit(r.error ? '' : 'primary')];
    $('#run-actions').replaceChildren(...actions);
}

function renderError(err, r) {
    const box = $('#run-error');
    box.hidden = !err;
    if (!err) return;
    box.querySelector('.err-hint').textContent = err.hint || '';
    box.querySelector('.err-hint').hidden = !err.hint;
    // a failed step is named as the steps list names it, not by its
    // Nextflow process ("Step failed: RESOLVE_TAXONOMY:PARSE_LINEAGE (…)"),
    // which stays as the title's tooltip
    const failed = r.steps.find((s) => s.state === 'failed');
    const title = err.process && failed && !/^Couldn't/.test(err.title)
        ? `Couldn't ${STEP_TEXT[failed.id][0].replace(/^./, (c) => c.toLowerCase())}` : err.title;
    box.querySelector('.title').textContent = title;
    box.querySelector('.title').title = err.process || '';
    box.querySelector('.msg').textContent = err.message || '';
    // a detail that only repeats the message (often with an "ERROR: " in
    // front) is left out; a traceback, or a tool's raw words behind a plainer
    // message, is worth showing
    const detail = (err.detail || '').split('\n').filter((l) => l.trim());
    const repeats = detail.length === 1 && detail[0].trim().endsWith((err.message || '').trim());
    box.querySelector('.detail').textContent = detail.length && !repeats ? detail.join('\n') : '';
    $('#err-workdir').hidden = !err.workdir;
    $('#err-workdir').onclick = () => window.qs.reveal(err.workdir);
    $('#err-log').hidden = !err.log;
    $('#err-log').onclick = () => window.qs.reveal(err.log);
    $('#err-edit').onclick = () => editRun(r);
}

// the New run form, filled as that run had it (form.js, prefillForm)
async function editRun(r) {
    showView('new');
    await prefillForm(r);
}

const thumbnails = new Map();
async function renderPreview(r) {
    const box = $('#run-preview');
    box.hidden = !(r.status === 'finished' && r.hasThumbnail);
    if (box.hidden) return;
    if (!thumbnails.has(r.id)) thumbnails.set(r.id, await window.qs.getThumbnail(r.id));
    if (state.runId !== r.id) return;
    $('#run-thumb').src = thumbnails.get(r.id) || '';
    box.hidden = !thumbnails.get(r.id);
    box.onclick = () => r.results.forEach((f) => window.qs.openResult(f));
}

function renderStats(r) {
    const s = r.summary;
    const box = $('#run-stats');
    box.hidden = !(r.status === 'finished' && s);
    if (box.hidden) return;
    const target = shortSpecies(r.targetSpecies) || fileStem(r.inputs.assembly);
    const sources = sourcesOf(r);
    const reference = sources.reference === 'file' ? fileStem(r.inputs.reference) : shortSpecies(r.found.reference?.species);
    const pct = (n) => (s.proteome_total ? `${((100 * n) / s.proteome_total).toFixed(1)}%` : '–');
    const stat = (value, label) => h('div', { class: 'stat' }, h('div', { class: 'stat-value', text: value }), h('div', { class: 'stat-label', text: label }));
    box.replaceChildren(
        stat(s.blocks !== undefined ? num(s.blocks) : '–', 'synteny blocks'),
        stat(pct(s.query_aligned), `of the proteins aligned to ${target}`),
        stat(pct(s.subject_aligned), `aligned to ${reference}`),
        sources.reference === 'discovered' && r.found.reference
            ? stat(r.found.reference.rank.replace(/^./, (c) => c.toUpperCase()), 'rank where the reference was found')
            : stat(s.proteome_total ? num(s.proteome_total) : '–', 'proteins in the proteome'));
}

function renderLog(r) {
    const log = $('#log');
    const text = state.logs.get(r.id);
    log.textContent = text || `This run's log isn't kept by the app: Nextflow's own is ${r.outdir}/.nextflow.log.`;
}

function renderRun(r) {
    renderHead(r);
    renderError(r.status === 'failed' ? r.error : null, r);
    renderPreview(r);
    renderStats(r);
    renderSteps(r);
    renderFound(r);
    renderLog(r);
}

function upsertRun(r) {
    const i = state.runs.findIndex((x) => x.id === r.id);
    if (i === -1) state.runs.unshift(r);
    else state.runs[i] = r;
}

window.qs.on('run:update', (r) => {
    upsertRun(r);
    renderSidebar();
    if (state.view === 'run' && state.runId === r.id) renderRun(r);
    if (state.view === 'new') refreshForm();
});
window.qs.on('runs:changed', (runs) => {
    state.runs = runs;
    renderSidebar();
});
window.qs.on('run:log', ({ id, text }) => {
    state.logs.set(id, (state.logs.get(id) || '') + text);
    if (state.view !== 'run' || state.runId !== id) return;
    const log = $('#log');
    const atBottom = log.scrollTop + log.clientHeight >= log.scrollHeight - 4;
    log.textContent = state.logs.get(id);
    if (atBottom) log.scrollTop = log.scrollHeight;
});

// ---------- tools: checking and first-launch setup ----------

async function checkEnv() {
    state.env = await window.qs.checkEnv();
    $('#runtime-root').textContent = tildify(state.env.runtimeRoot);
    $('#setup-root').textContent = tildify(state.env.runtimeRoot);
    renderTools();
    updatePreflight();
    refreshForm();
    if (state.view === 'settings') renderSettings();
}

// Nextflow, Java and the plugins into the runtime root, with micromamba
// (bundled, or fetched first when running from source). Starts by itself on
// a first launch; Retry setup after a failure.
async function setupRuntime() {
    state.setup = 'running';
    renderTools();
    showView('setup');
    $('#setup-progress').hidden = false;
    $('#setup-error').hidden = true;
    $('#setup-actions').hidden = true;
    $('#setup-log').textContent = '';
    try {
        await window.qs.setupRuntime();
        state.setup = 'idle';
        await checkEnv();
        showView('new');
    } catch (e) {
        state.setup = 'failed';
        $('#setup-log').textContent += `\n${cleanIpcError(e)}\n`;
        // micromamba's own offline message blames a corrupted cache
        const net = await window.qs.checkNetwork();
        const err = net.conda ? {
            title: 'Setup failed', message: cleanIpcError(e).split('\n').filter(Boolean).pop(),
            hint: 'The technical log below has micromamba\'s full output. Retry once the cause is fixed.',
        } : {
            title: 'Setup needs an internet connection',
            message: 'The first launch downloads Nextflow and Java from conda-forge and bioconda, which this computer can\'t reach right now. '
                + 'Connect to the internet and retry.',
            hint: 'This is needed once. After that the app works offline with your own reference genome and proteome; '
                + 'finding a relative on NCBI always needs a connection.',
        };
        const box = $('#setup-error');
        box.querySelector('.title').textContent = err.title;
        box.querySelector('.msg').textContent = err.message;
        box.querySelector('.err-hint').textContent = err.hint;
        box.querySelector('.err-hint').hidden = false;
        box.hidden = false;
        $('#setup-progress').hidden = true;
        $('#setup-actions').hidden = false;
        renderTools();
    }
}
$('#btn-setup').addEventListener('click', setupRuntime);
window.qs.on('setup:log', (text) => {
    const log = $('#setup-log');
    log.textContent += text;
    log.scrollTop = log.scrollHeight;
});

// ---------- settings ----------

async function loadSettings() {
    const s = await window.qs.getSettings();
    $('#set-nextflow').value = s.nextflowPath || '';
    $('#set-java').value = s.javaHome || '';
    $('#set-updates').checked = s.checkUpdates !== false;
}

function renderSettings() {
    const nf = state.env?.nextflow;
    $('#settings-tools').textContent = !nf ? 'Checking…' : nf.ok
        ? `Nextflow ${nf.version}, ${nf.source === 'settings' ? 'the one set below' : 'the app\'s own'}.`
        : `${nf.message}${nf.detail ? `: ${nf.detail}` : ''}`;
}

$('#btn-save-settings').addEventListener('click', async () => {
    await window.qs.setSettings({ nextflowPath: $('#set-nextflow').value.trim(), javaHome: $('#set-java').value.trim(),
        checkUpdates: $('#set-updates').checked });
    $('#btn-save-settings').textContent = 'Saved';
    setTimeout(() => { $('#btn-save-settings').textContent = 'Save'; }, 1500);
    await checkEnv();
});

// ---------- updates ----------

// a newer release, if there is one: never blocks or delays anything else
async function showUpdate() {
    const u = await window.qs.checkForUpdate().catch(() => null);
    if (!u) return;
    const banner = $('#update-banner');
    banner.querySelector('.msg').textContent = `quick_synteny ${u.version} is available (you have ${u.current}).`;
    $('#update-mac-hint').hidden = u.platform !== 'darwin';
    $('#update-download').onclick = () => window.qs.openUpdate(u.downloadUrl);
    $('#update-notes').onclick = () => window.qs.openUpdate(u.notesUrl);
    $('#update-dismiss').onclick = async () => {
        banner.hidden = true;
        await window.qs.dismissUpdate(u.version);
    };
    banner.hidden = false;
}

// ---------- start ----------

(async () => {
    document.body.classList.add(`platform-${window.qs.platform}`);
    loadSettings();
    applySchema();
    updateEngineHint();
    state.runs = await window.qs.listRuns();
    const running = state.runs.find((r) => r.status === 'running');
    if (running) showRun(running.id);
    else showView('new');
    await checkEnv();
    if (state.env.nextflow.needsSetup) setupRuntime();
    showUpdate();
})();
