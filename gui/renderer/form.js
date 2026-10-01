// The New run view: files (dropped or chosen), what to compare with (your
// own reference and proteome by default, or discovery on NCBI), where to
// save, and options; checked before a run starts so a bad input fails here,
// not minutes into a run.

const files = { assembly: null, reference: null, proteome: null };
let outdir = null;
let outdirInfo = null;
let outdirChosen = false;       // picked by the user: stop suggesting one

const mode = () => document.querySelector('input[name=mode]:checked').value;
const engine = () => document.querySelector('input[name=engine]:checked').value;

// ---------- files ----------

const FASTA = [{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna', 'faa', 'fas', 'fsa', 'gz'] }, { name: 'All files', extensions: ['*'] }];
const NOT_FASTA = 'This doesn\'t look like a FASTA file (it should start with a ">" line, optionally gzipped).';

function setFile(key, f) {
    files[key] = f;
    const meta = `${fmtBytes(f.bytes)}${f.kind === 'gzip-fasta' ? ' · gzipped FASTA' : f.kind === 'fasta' ? ' · FASTA' : ''}`;
    if (key === 'assembly') {
        $('#drop-assembly').hidden = true;
        const card = $('#card-assembly');
        card.hidden = false;
        card.querySelector('.file-name').textContent = baseName(f.path);
        card.querySelector('.file-name').title = f.path;
        card.querySelector('.file-meta').textContent = meta;
    } else {
        const slot = $(`#slot-${key}`);
        slot.classList.add('filled');
        slot.querySelector('.file-name').textContent = baseName(f.path);
        slot.querySelector('.file-name').title = f.path;
        slot.querySelector('.file-meta').textContent = meta;
        slot.querySelector('button').textContent = 'Change…';
    }
    setStatus(`#status-${key}`, f.kind === 'other' ? 'bad' : null, NOT_FASTA);
    updatePreflight();
    suggestOutdir();
    refreshForm();
}

document.querySelectorAll('[data-pick]').forEach((btn) => {
    btn.addEventListener('click', async () => {
        const f = await window.qs.pickFile({ title: btn.dataset.title, filters: FASTA }).catch((e) => {
            setStatus(`#status-${btn.dataset.pick}`, 'bad', e.message);
            return null;
        });
        if (f) setFile(btn.dataset.pick, f);
    });
});

// drag and drop onto the genome area or a reference/proteome slot; a file
// dropped anywhere else must not make the window navigate to it
window.addEventListener('dragover', (e) => e.preventDefault());
window.addEventListener('drop', (e) => e.preventDefault());
document.querySelectorAll('[data-drop]').forEach((zone) => {
    const key = zone.dataset.drop;
    zone.addEventListener('dragover', (e) => {
        e.preventDefault();
        zone.classList.add('over');
    });
    zone.addEventListener('dragleave', () => zone.classList.remove('over'));
    zone.addEventListener('drop', async (e) => {
        e.preventDefault();
        zone.classList.remove('over');
        const file = e.dataTransfer.files[0];
        if (!file) return;
        try {
            setFile(key, await window.qs.inspectFile(window.qs.pathForFile(file)));
        } catch (err) {
            setStatus(`#status-${key}`, 'bad', String(err.message).replace(/^Error invoking remote method '[^']+': (Error: )?/, ''));
        }
    });
});

// ---------- compare with ----------

document.querySelectorAll('input[name=mode]').forEach((r) => {
    r.addEventListener('change', () => {
        const manual = mode() === 'manual';
        $('#mode-auto').hidden = manual;
        $('#mode-manual').hidden = !manual;
        updatePreflight();
        suggestOutdir();
        refreshForm();
    });
});

// taxon name -> taxid lookup (lib/taxa.js), once typing pauses; an answer
// for text that has changed since is dropped
let suggestTimer = null;
$('#taxon-search').addEventListener('input', () => {
    clearTimeout(suggestTimer);
    const q = $('#taxon-search').value.trim();
    if (/^\d+$/.test(q)) {
        $('#taxid').value = q;
        $('#taxon-results').hidden = true;
        checkTaxid();
        return;
    }
    if (q.length < 3) { $('#taxon-results').hidden = true; return; }
    suggestTimer = setTimeout(async () => {
        let hits = [];
        try {
            hits = await window.qs.suggestTaxon(q);
            if ($('#taxon-search').value.trim() !== q) return;
            setStatus('#search-status', hits.length ? null : '', `NCBI knows no taxon called "${q}". Check the spelling, or type a taxid.`);
        } catch (e) {
            if ($('#taxon-search').value.trim() !== q) return;
            setStatus('#search-status', 'warn', /busy/.test(String(e.message))
                ? 'NCBI is busy right now: try again in a moment, or type a taxid.'
                : 'Can\'t reach NCBI to search by name. Check the internet connection, type a taxid instead, or use your own files.');
        }
        $('#taxon-results').replaceChildren(...hits.slice(0, 12).map((t) => h('li', {
            onclick: () => {
                $('#taxid').value = t.tax_id;
                $('#taxon-search').value = t.sci_name;
                $('#taxon-results').hidden = true;
                showTaxon({ status: 'found', taxid: t.tax_id, name: t.sci_name, rank: t.rank?.toLowerCase(), common: t.common_name });
            },
        }, `${t.sci_name} `, h('small', { text: [t.rank?.toLowerCase(), t.common_name, `taxid ${t.tax_id}`].filter(Boolean).join(' · ') }))));
        $('#taxon-results').hidden = hits.length === 0;
    }, 450);
});

// The taxid, checked against NCBI as it's typed, so a mistyped or unknown
// one is caught here rather than a minute into a run. NCBI being
// unreachable doesn't block a run -- the pipeline checks it again.
let taxon = { status: 'empty' };
let taxonCheck = Promise.resolve();
let taxidTimer = null;

function showTaxon(t) {
    taxon = t;
    const text = {
        empty: null,
        invalid: ['bad', 'Not a taxid: NCBI taxids are positive whole numbers (e.g. 4932).'],
        checking: ['', 'Checking with NCBI…'],
        found: ['ok', `✓ ${t.name}${t.rank ? ` · ${t.rank}` : ''}${t.common ? ` · ${t.common}` : ''}`],
        unknown: ['bad', `NCBI has no taxon ${t.taxid}.`],
        unchecked: ['', 'Couldn\'t reach NCBI to check this taxid; the run will check it.'],
    }[t.status];
    if (text) setStatus('#taxid-status', ...text);
    else setStatus('#taxid-status', null);
    if (t.status === 'found') suggestOutdir();
}

function checkTaxid() {
    clearTimeout(taxidTimer);
    const v = $('#taxid').value.trim();
    if (!v) return showTaxon({ status: 'empty' });
    if (!/^[1-9]\d*$/.test(v)) return showTaxon({ status: 'invalid', taxid: v });
    if (taxon.taxid === v && ['found', 'unknown'].includes(taxon.status)) return;
    showTaxon({ status: 'checking', taxid: v });
    taxonCheck = new Promise((resolve) => {
        taxidTimer = setTimeout(async () => {
            let t;
            try {
                const r = await window.qs.lookupTaxon(v);
                t = r.found ? { status: 'found', ...r } : { status: 'unknown', taxid: v };
            } catch {
                t = { status: 'unchecked', taxid: v };
            }
            // typed on since: a newer check owns the status
            if ($('#taxid').value.trim() === v) showTaxon(t);
            resolve();
        }, 400);
    });
}
$('#taxid').addEventListener('input', checkTaxid);

// ---------- where to save ----------

// Documents/quick_synteny/<target>_vs_<taxid or reference>, unless the user
// has picked a folder; a new one per run
let suggestTimerOut = null;
function suggestOutdir(force = false) {
    if (outdirChosen && !force) return;
    clearTimeout(suggestTimerOut);
    suggestTimerOut = setTimeout(async () => {
        if (!files.assembly) return;
        const d = await window.qs.defaultOutdir({
            assembly: files.assembly.path,
            taxid: mode() === 'auto' ? $('#taxid').value.trim() : '',
            reference: mode() === 'manual' ? files.reference?.path : '',
        });
        applyOutdir(d);
    }, 150);
}

$('#pick-outdir').addEventListener('click', async () => {
    const d = await window.qs.pickDir({ title: 'Save the run in…' });
    if (!d) return;
    outdirChosen = true;
    applyOutdir(d);
});

// what's in the folder decides the note under it; writing into one that
// holds other things (like ~/Downloads) would scatter results, work/ and
// .nextflow/ among them, so a subfolder is offered instead
function applyOutdir(d) {
    outdir = d.path;
    outdirInfo = d;
    $('#path-outdir').textContent = tildify(d.path);
    $('#path-outdir').title = d.path;
    $('#btn-subfolder').hidden = true;
    $('#resume-row').hidden = !d.isRun;
    if (!d.isRun) document.querySelector('[name=resume]').checked = false;
    if (d.isNew) {
        setStatus('#outdir-status', '', 'A new folder, made when the run starts.');
    } else if (d.unreadable) {
        setStatus('#outdir-status', 'bad', 'macOS won\'t let the app read this folder. Choose another, or allow it in '
            + 'System Settings → Privacy & Security → Files and Folders.');
    } else if (!d.writable) {
        setStatus('#outdir-status', 'bad', 'This folder can\'t be written to. Choose another one.');
    } else if (d.isRun) {
        setStatus('#outdir-status', '', 'An earlier run is in this folder.');
    } else if (!d.empty) {
        setStatus('#outdir-status', 'warn', 'This folder isn\'t empty: results, work/ and .nextflow/ would go straight into it.');
        $('#btn-subfolder').hidden = false;
        $('#btn-subfolder').textContent = 'Use a new subfolder instead';
    } else {
        setStatus('#outdir-status', null);
    }
    refreshForm();
}

$('#btn-subfolder').addEventListener('click', async () => {
    const name = await window.qs.suggestRunName({ assembly: files.assembly?.path, taxid: $('#taxid').value.trim() });
    applyOutdir({ path: `${outdirInfo.path.replace(/\/$/, '')}/${name}`, writable: true, empty: true, isRun: false, isNew: true });
});

// ---------- options ----------

function updateEngineHint() {
    const docker = state.docker;
    if (engine() === 'conda') {
        $('#engine-hint').textContent = 'Each tool runs in a conda environment the app builds from the pipeline\'s envs/*.yml: '
            + 'a few minutes and some downloads the first time, reused after that.';
        setStatus('#docker-status', null);
    } else {
        $('#engine-hint').textContent = 'Each tool runs in its pinned container image, the same ones the HPC profile uses.';
        if (docker === 'checking') setStatus('#docker-status', '', 'Checking Docker…');
        else if (docker?.ok) setStatus('#docker-status', 'ok', `✓ Docker ${docker.version}${docker.memBytes ? ` · ${(docker.memBytes / 2 ** 30).toFixed(1)} GB for containers` : ''}`);
        else if (docker) setStatus('#docker-status', 'bad', `${docker.message}: start Docker Desktop, or use Conda.`);
    }
}

async function refreshDocker() {
    state.docker = 'checking';
    updateEngineHint();
    state.docker = await window.qs.checkDocker();
    updateEngineHint();
    updatePreflight();
}
// Docker is checked afresh each time its engine is picked, never otherwise
document.querySelectorAll('input[name=engine]').forEach((r) => {
    r.addEventListener('change', () => {
        if (engine() === 'docker') refreshDocker();
        else updateEngineHint();
        updatePreflight();
    });
});

// miniprot's whole-genome index needs ~10 GB RAM per Gb of genome (README,
// "Large genomes"); warn before a run that can't fit -- in Docker's VM, or
// in this computer's RAM for the conda engine
function updatePreflight() {
    const box = $('#preflight');
    const inDocker = engine() === 'docker';
    const mem = inDocker ? state.docker?.memBytes : state.env?.host?.memBytes;
    const genomes = [files.assembly, mode() === 'manual' ? files.reference : null].filter(Boolean);
    if (!mem || !genomes.length) { box.hidden = true; return; }
    const gb = Math.max(...genomes.map((f) => (f.path.endsWith('.gz') ? 3.3 : 1) * f.bytes / 1e9));
    const memGb = mem / 2 ** 30;
    const chunk = parseFloat($('#chunk').value);
    const need = chunk > 0 ? Math.min(chunk, gb) * 11 + 4 : gb * 10 + 4;
    if (need <= memGb) { box.hidden = true; return; }
    const fit = Math.floor(((memGb - 4) / 11) * 10) / 10;
    box.hidden = false;
    box.textContent = `The largest genome is about ${gb.toFixed(2)} Gb: miniprot needs about ${need.toFixed(0)} GB of memory, `
        + `but ${inDocker ? 'Docker has' : 'this computer has'} ${memGb.toFixed(1)} GB. `
        + `${inDocker ? 'Raise Docker Desktop\'s memory limit, or set' : 'Set'} Options → Chunk size `
        + `to ${fit > 0 ? fit : 'a smaller value'} Gb or less (close to, but not identical to, a whole-genome run).`;
}
$('#chunk').addEventListener('input', updatePreflight);

// ---------- the pipeline's parameter schema ----------

// nextflow_schema.json, flattened to { name: property }: the "?" tips show
// each parameter's description, help_text and default (as --help <param>
// does), fields without their own placeholder show its default, the rank
// list is its enum, and numbers are checked against its limits before a run
let schema = {};

async function applySchema() {
    const s = await window.qs.getSchema();
    schema = Object.assign({}, ...Object.values(s.$defs || {}).map((d) => d.properties || {}));
    document.querySelectorAll('.help[data-param]').forEach((btn) => {
        const p = schema[btn.dataset.param] || {};
        const text = [p.description, p.help_text, p.default !== undefined && `Default: ${p.default}.`].filter(Boolean).join(' ');
        btn.dataset.help = text;
        btn.setAttribute('aria-label', text);
    });
    document.querySelectorAll('#form input[name]').forEach((input) => {
        const p = schema[input.name];
        if (p && p.default !== undefined && typeof p.default !== 'boolean' && !input.placeholder) input.placeholder = String(p.default);
    });
    const rank = schema.max_rank;
    $('select[name=max_rank]').replaceChildren(...rank.enum.map((r) => new Option(r, r, false, r === rank.default)));
}

// a typed value against its schema property: a message, or null if it's fine
function checkNumber(name, value) {
    const p = schema[name];
    if (!p) return null;
    const types = [].concat(p.type);
    const label = document.querySelector(`#form [name=${name}]`).closest('label').firstChild.textContent.trim();
    const n = Number(value);
    if (types.includes('integer') && !(/^-?\d+$/.test(value))) return `${label} must be a whole number.`;
    if ((types.includes('number') || types.includes('integer')) && !Number.isFinite(n)) return `${label} must be a number.`;
    if (p.minimum !== undefined && n < p.minimum) return `${label} must be at least ${p.minimum}.`;
    if (p.exclusiveMinimum !== undefined && n <= p.exclusiveMinimum) return `${label} must be more than ${p.exclusiveMinimum}.`;
    if (p.maximum !== undefined && n > p.maximum) return `${label} must be at most ${p.maximum}.`;
    return null;
}

// ---------- the run ----------

// the Options section's fields, each named after its pipeline parameter
const ADVANCED = ['min_seq_size', 'min_asm_gap', 'min_identity', 'max_gap', 'min_block', 'miniprot_m', 'miniprot_chunk_gb'];

function collectForm() {
    const fd = new FormData($('#form'));
    const form = { engine: engine(), outdir, assembly: files.assembly?.path, resume: fd.get('resume') === 'on' };
    for (const k of ADVANCED) {
        const v = String(fd.get(k) || '').trim();
        if (v) form[k] = v;
    }
    if (mode() === 'manual') {
        form.reference = files.reference?.path;
        form.proteome = files.proteome?.path;
    } else {
        form.taxid = String(fd.get('taxid') || '').trim();
        form.max_rank = fd.get('max_rank');
        form.exclude_target = fd.get('exclude_target') === 'on';
    }
    return form;
}

// what stops a run before it starts, as one sentence, or null
function validate(form) {
    if (!state.env?.nextflow?.ok) return state.setup === 'running' ? 'The tools are still being set up.' : 'The tools aren\'t ready: see Settings.';
    if (state.runs.some((r) => r.status === 'running')) return 'A run is in progress: wait for it, or cancel it.';
    if (form.engine === 'docker' && !state.docker?.ok) return `${state.docker?.message || 'Docker is not available'}: start it, or use Conda.`;
    if (!form.assembly) return 'Add your genome to start.';
    if (mode() === 'manual' && !form.reference) return 'Add a reference genome to compare with.';
    if (mode() === 'manual' && !form.proteome) return 'Add a proteome.';
    const LABELS = { assembly: 'Your genome', reference: 'The reference genome', proteome: 'The proteome' };
    const notFasta = ['assembly', ...(mode() === 'manual' ? ['reference', 'proteome'] : [])].filter((k) => files[k]?.kind === 'other');
    if (notFasta.length) return `${LABELS[notFasta[0]]} doesn't look like a FASTA file.`;
    if (mode() === 'auto') {
        if (!form.taxid) return 'Pick the organism your genome is from.';
        if (taxon.status === 'invalid') return `"${form.taxid}" is not a taxid: NCBI taxids are positive whole numbers.`;
        if (taxon.status === 'unknown') return `NCBI has no taxon ${form.taxid}. Search by name to find the right one.`;
    }
    for (const k of ADVANCED) {
        const bad = form[k] !== undefined && checkNumber(k, form[k]);
        if (bad) return bad;
    }
    if (!form.outdir) return 'Choose where to save the run.';
    if (outdirInfo?.unreadable || outdirInfo?.writable === false) return 'The folder to save to can\'t be used (see Save to).';
    return null;
}

// the footer's message and button, as things stand (not an error yet:
// errors show in red only after Start run is pressed)
function refreshForm(error = null) {
    const msg = $('#form-msg');
    const blocked = validate(collectForm());
    msg.textContent = error || blocked || 'Ready. A genome this size usually takes a few minutes.';
    msg.className = `form-msg ${error ? 'bad' : ''}`;
    $('#btn-run').disabled = !!blocked;
}

$('#form').addEventListener('submit', async (e) => {
    e.preventDefault();
    const form = collectForm();
    if (form.engine === 'docker') await refreshDocker();
    if (mode() === 'auto') {
        checkTaxid();
        await taxonCheck;
    }
    const err = validate(form);
    if (err) return refreshForm(err);
    // discovery fails minutes in without NCBI (datasets retries 11 times
    // first); say so now instead
    if (mode() === 'auto' && !(await window.qs.checkNetwork()).ncbi) {
        return refreshForm('Can\'t reach NCBI: finding a relative needs an internet connection. Connect, or compare with your own files.');
    }
    $('#btn-run').disabled = true;
    try {
        const id = await window.qs.startRun(form, { targetSpecies: mode() === 'auto' && taxon.status === 'found' ? taxon.name : '' });
        // the next run gets a folder of its own
        outdirChosen = false;
        suggestOutdir(true);
        showRun(id);
    } catch (ex) {
        refreshForm(String(ex.message).replace(/^Error invoking remote method '[^']+': (Error: )?/, ''));
    }
});

// "Run again": the form as that run had it, in a new folder
async function prefillForm(r) {
    const f = r.form;
    document.querySelector(`input[name=mode][value=${f.reference && f.proteome ? 'manual' : 'auto'}]`).click();
    for (const key of ['assembly', 'reference', 'proteome']) {
        if (!f[key]) continue;
        try {
            setFile(key, await window.qs.inspectFile(f[key]));
        } catch {
            setStatus(`#status-${key}`, 'bad', `${baseName(f[key])} isn't where it was any more: choose it again.`);
        }
    }
    if (f.taxid) {
        $('#taxid').value = f.taxid;
        $('#taxon-search').value = r.targetSpecies || f.taxid;
        checkTaxid();
    }
    if (f.max_rank) document.querySelector('[name=max_rank]').value = f.max_rank;
    document.querySelector('[name=exclude_target]').checked = !!f.exclude_target;
    for (const k of ADVANCED) document.querySelector(`[name=${k}]`).value = f[k] || '';
    document.querySelector(`input[name=engine][value=${f.engine || 'conda'}]`).click();
    if (ADVANCED.some((k) => f[k]) || f.engine === 'docker') $('#options').open = true;
    outdirChosen = false;
    suggestOutdir(true);
}

$('#form').addEventListener('input', () => refreshForm());
