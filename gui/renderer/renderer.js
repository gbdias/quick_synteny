const $ = (sel) => document.querySelector(sel);
const files = { assembly: null, reference: null, proteome: null };
let outdir = null;
let env = null;
// the last Docker check -- only made while the Docker engine is selected;
// 'checking' while one is in flight
let docker = null;
const tasks = new Map();

// ---------- environment ----------

function pill(text, cls) {
    const s = document.createElement('span');
    s.className = `pill ${cls}`;
    s.textContent = text;
    return s;
}

const gib = (bytes) => `${(bytes / 2 ** 30).toFixed(1)} GB`;
const engine = () => document.querySelector('input[name=engine]:checked').value;

async function checkEnv() {
    $('#env').replaceChildren(pill('checking environment…', ''));
    env = await window.qs.checkEnv();
    $('#btn-setup').hidden = !env.nextflow.needsSetup;
    $('#btn-setup').textContent = 'Retry setup';
    $('#runtime-root').textContent = env.runtimeRoot;
    renderEnv();
}

// header pills: the runtime's Nextflow, plus Docker only for the Docker engine
function renderEnv() {
    if (!env) return;
    const nf = env.nextflow;
    const nfText = nf.ok ? `Nextflow ${nf.version}${nf.source === 'settings' ? ' (from Settings)' : ''}` : nf.message;
    const pills = [pill(nfText, nf.ok ? 'ok' : nf.needsSetup ? 'warn' : 'bad')];
    if (nf.detail) pills[0].title = nf.detail;
    if (engine() === 'docker' && docker === 'checking') pills.push(pill('checking Docker…', ''));
    else if (engine() === 'docker' && docker) {
        pills.push(pill(docker.ok ? `Docker ${docker.version}${docker.memBytes ? ` · ${gib(docker.memBytes)}` : ''}` : docker.message,
            docker.ok ? 'ok' : 'bad'));
    }
    $('#env').replaceChildren(...pills);
    updateEngineHint();
    updatePreflight();
}

async function refreshDocker() {
    docker = 'checking';
    renderEnv();
    docker = await window.qs.checkDocker();
    renderEnv();
}

// first launch: Nextflow + Java into the runtime root, with micromamba
// (bundled, or fetched first when running from source). Starts on its own
// once; the button retries after a failure.
async function setupRuntime() {
    $('#btn-setup').hidden = true;
    $('#log').textContent = '';
    setRunning(true, 'setting up');
    $('#btn-cancel').hidden = true;
    $('#run-note').textContent = 'First launch: installing Nextflow and Java for the app (about a minute)…';
    $('#run-note').hidden = false;
    $('#run-error').hidden = true;
    let ok = false;
    try {
        await window.qs.setupRuntime();
        ok = true;
    } catch (e) {
        appendLog(`\n${e.message}\n`);
        // micromamba's own offline message blames a corrupted cache
        const net = await window.qs.checkNetwork();
        showError(net.conda ? {
            title: 'Setup failed',
            message: e.message.split('\n').filter(Boolean).pop(),
            hint: 'The Nextflow log below has micromamba\'s full output. Retry setup once the cause is fixed.',
        } : {
            title: 'Setup needs an internet connection',
            message: 'The first launch downloads Nextflow and Java from conda-forge and bioconda, '
                + 'which this computer can\'t reach right now. Connect to the internet and click Retry setup.',
            hint: 'This is needed once. After that the app works offline with your own reference genome and '
                + 'proteome; finding a reference on NCBI always needs a connection.',
        });
    }
    $('#run-note').hidden = true;
    setRunning(false, ok ? 'runtime ready' : 'setup failed');
    $('#run-status').className = `pill ${ok ? 'ok' : 'bad'}`;
    await checkEnv();
}
$('#btn-setup').onclick = setupRuntime;
window.qs.on('setup:log', (text) => appendLog(text));

function updateEngineHint() {
    const h = $('#engine-hint');
    if (engine() === 'conda') {
        h.textContent = 'Each tool runs in a conda environment built from the pipeline\'s envs/*.yml — '
            + 'a few minutes and some downloads on the first run, reused after that.';
    } else {
        h.textContent = docker?.ok ? 'Each tool runs in its pinned container image — the same ones the HPC profile uses.'
            : 'Needs Docker Desktop (or another Docker engine) running.';
    }
}
// Docker is checked afresh each time its engine is picked, never otherwise
document.querySelectorAll('input[name=engine]').forEach((r) => {
    r.onchange = () => (engine() === 'docker' ? refreshDocker() : renderEnv());
});

async function loadSettings() {
    const s = await window.qs.getSettings();
    $('#set-nextflow').value = s.nextflowPath || '';
    $('#set-java').value = s.javaHome || '';
    $('#set-updates').checked = s.checkUpdates !== false;
}

$('#btn-settings').onclick = () => { $('#settings').hidden = !$('#settings').hidden; };
$('#btn-save-settings').onclick = async () => {
    await window.qs.setSettings({ nextflowPath: $('#set-nextflow').value.trim(), javaHome: $('#set-java').value.trim(),
        checkUpdates: $('#set-updates').checked });
    await checkEnv();
};
$('#btn-open-result').onclick = () => window.qs.pickResult().catch((e) => alert(e.message));

// ---------- form ----------

const FASTA = [{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna', 'faa', 'fas', 'gz'] }, { name: 'All files', extensions: ['*'] }];

document.querySelectorAll('[data-pick]').forEach((btn) => {
    btn.onclick = async () => {
        const key = btn.dataset.pick;
        const f = await window.qs.pickFile({ title: btn.dataset.title, filters: FASTA });
        if (!f) return;
        files[key] = f;
        $(`#path-${key}`).textContent = `${f.path}  (${fmtBytes(f.bytes)})`;
        setStatus(`#status-${key}`, f.kind === 'other' ? 'bad' : null,
            'This doesn\'t look like a FASTA file (it should start with a ">" line, optionally gzipped).');
        updatePreflight();
    };
});

// a status line under a field: `cls` is ok/warn/bad/'' (neutral), or null to hide it
function setStatus(sel, cls, text) {
    const el = $(sel);
    el.hidden = cls === null;
    if (cls === null) return;
    el.className = `field-status ${cls}`;
    (el.querySelector('.msg') || el).textContent = text;
}

// what's in the chosen output folder decides the warning under it; writing
// into a folder that holds other things (like ~/Downloads) would scatter
// results, work/ and .nextflow/ among them, so offer a subfolder instead
let outdirInfo = null;
$('#pick-outdir').onclick = async () => {
    const d = await window.qs.pickDir({ title: 'Output folder' });
    if (d) applyOutdir(d);
};
function applyOutdir(d) {
    outdir = d.path;
    outdirInfo = d;
    $('#path-outdir').textContent = d.path;
    $('#btn-subfolder').hidden = true;
    if (d.unreadable) {
        setStatus('#outdir-status', 'bad', 'macOS won\'t let the app read this folder. Choose another, or allow it in '
            + 'System Settings → Privacy & Security → Files and Folders.');
    } else if (!d.writable) {
        setStatus('#outdir-status', 'bad', 'This folder can\'t be written to. Choose another one.');
    } else if (d.isRun) {
        setStatus('#outdir-status', '', 'A previous run is in this folder: tick Resume to reuse its finished steps.');
    } else if (!d.empty) {
        setStatus('#outdir-status', 'warn', 'This folder isn\'t empty: results, work/ and .nextflow/ would go straight into it.');
        $('#btn-subfolder').hidden = false;
        $('#btn-subfolder').textContent = 'Use a new subfolder instead';
    } else {
        setStatus('#outdir-status', null);
    }
}
$('#btn-subfolder').onclick = async () => {
    const name = await window.qs.suggestRunName({ assembly: files.assembly?.path, taxid: $('#taxid').value.trim() });
    outdir = `${outdirInfo.path.replace(/\/$/, '')}/${name}`;
    outdirInfo = { path: outdir, writable: true, empty: true, isRun: false };
    $('#path-outdir').textContent = outdir;
    $('#btn-subfolder').hidden = true;
    setStatus('#outdir-status', '', 'A new folder, created when the run starts.');
};

document.querySelectorAll('input[name=mode]').forEach((r) => {
    r.onchange = () => {
        const manual = mode() === 'manual';
        $('#mode-auto').hidden = manual;
        $('#mode-manual').hidden = !manual;
        updatePreflight();
    };
});
const mode = () => document.querySelector('input[name=mode]:checked').value;

// taxon name -> taxid lookup, debounced
let suggestTimer = null;
$('#taxon-search').oninput = () => {
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
            setStatus('#search-status', null);
        } catch {
            setStatus('#search-status', 'warn', 'Can\'t reach NCBI to search by name. Check the internet connection, '
                + 'type a taxid instead, or use your own files.');
        }
        const ul = $('#taxon-results');
        ul.replaceChildren(...hits.slice(0, 12).map((h) => {
            const li = document.createElement('li');
            li.textContent = `${h.sci_name} `;
            const small = document.createElement('small');
            small.textContent = [h.rank?.toLowerCase(), h.common_name, `taxid ${h.tax_id}`].filter(Boolean).join(' · ');
            li.append(small);
            li.onclick = () => {
                $('#taxid').value = h.tax_id;
                $('#taxon-search').value = h.sci_name;
                ul.hidden = true;
                showTaxon({ status: 'found', taxid: h.tax_id, name: h.sci_name, rank: h.rank?.toLowerCase(), common: h.common_name });
            };
            return li;
        }));
        ul.hidden = hits.length === 0;
    }, 300);
};

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
        invalid: [`bad`, 'Not a taxid: NCBI taxids are positive whole numbers (e.g. 4932).'],
        checking: ['', 'Checking with NCBI…'],
        found: ['ok', `✓ ${t.name}${t.rank ? ` · ${t.rank}` : ''}${t.common ? ` · ${t.common}` : ''}`],
        unknown: ['bad', `NCBI has no taxon ${t.taxid}.`],
        unchecked: ['', 'Couldn\'t reach NCBI to check this taxid; the run will check it.'],
    }[t.status];
    if (text) setStatus('#taxid-status', ...text);
    else setStatus('#taxid-status', null);
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
$('#taxid').oninput = checkTaxid;

// miniprot's whole-genome index needs ~10 GB RAM per Gb of genome (README,
// "Large genomes"); warn before a run that can't fit -- in Docker's VM, or
// in this computer's RAM for the conda engine
function updatePreflight() {
    const box = $('#preflight');
    const inDocker = engine() === 'docker';
    const mem = inDocker ? docker?.memBytes : env?.host?.memBytes;
    const genomes = [files.assembly, mode() === 'manual' ? files.reference : null].filter(Boolean);
    if (!mem || !genomes.length) { box.hidden = true; return; }
    const gb = Math.max(...genomes.map((f) => (f.path.endsWith('.gz') ? 3.3 : 1) * f.bytes / 1e9));
    const memGb = mem / 2 ** 30;
    const chunk = parseFloat($('#chunk').value);
    const need = chunk > 0 ? Math.min(chunk, gb) * 11 + 4 : gb * 10 + 4;
    if (need <= memGb) { box.hidden = true; return; }
    const fit = Math.floor(((memGb - 4) / 11) * 10) / 10;
    box.hidden = false;
    box.textContent = `The largest genome is ~${gb.toFixed(2)} Gb: miniprot needs ~${need.toFixed(0)} GB, `
        + `but ${inDocker ? 'Docker has' : 'this computer has'} ${memGb.toFixed(1)} GB. `
        + `${inDocker ? 'Raise Docker Desktop\'s memory limit, or set' : 'Set'} Advanced → `
        + `Chunk size ≤ ${fit > 0 ? fit : '—'} Gb (close to, but not identical to, a whole-genome run).`;
}
$('#chunk').oninput = updatePreflight;

function collectForm() {
    const fd = new FormData($('#form'));
    const form = { engine: engine(), outdir, assembly: files.assembly?.path, resume: fd.get('resume') === 'on' };
    for (const k of ['min_seq_size', 'min_asm_gap', 'min_identity', 'max_gap', 'min_block', 'miniprot_m', 'miniprot_chunk_gb']) {
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

function validate(form) {
    if (!env?.nextflow?.ok) return env?.nextflow?.needsSetup ? 'Set up the runtime first (button at the top).' : 'Nextflow is not available (see the header).';
    if (form.engine === 'docker' && !docker?.ok) return `${docker?.message || 'Docker is not available'} — start it, or choose Tools → Conda.`;
    if (!form.assembly) return 'Choose a target genome FASTA.';
    if (mode() === 'manual' && !(form.reference && form.proteome)) return 'Choose both a reference genome and a proteome.';
    const notFasta = ['assembly', ...(mode() === 'manual' ? ['reference', 'proteome'] : [])]
        .filter((k) => files[k]?.kind === 'other');
    const LABELS = { assembly: 'The target genome', reference: 'The reference genome', proteome: 'The proteome' };
    if (notFasta.length) return `${LABELS[notFasta[0]]} is not a FASTA file (see the note under it).`;
    if (mode() === 'auto') {
        if (!form.taxid) return 'Pick the target organism (a taxid).';
        if (taxon.status === 'invalid') return `"${form.taxid}" is not a taxid: NCBI taxids are positive whole numbers.`;
        if (taxon.status === 'unknown') return `NCBI has no taxon ${form.taxid}. Search by name to find the right one.`;
    }
    if (!form.outdir) return 'Choose an output folder.';
    if (outdirInfo?.unreadable || outdirInfo?.writable === false) return 'The output folder can\'t be used (see the note under it).';
    return null;
}

// ---------- run ----------

$('#form').onsubmit = async (e) => {
    e.preventDefault();
    const form = collectForm();
    if (form.engine === 'docker') await refreshDocker();
    if (mode() === 'auto') {
        checkTaxid();
        await taxonCheck;
    }
    const err = validate(form);
    if (err) { alert(err); return; }
    // discovery fails a minute in without NCBI (datasets retries 11 times
    // first); say so now instead
    if (mode() === 'auto' && !(await window.qs.checkNetwork()).ncbi) {
        alert('Can\'t reach NCBI. Finding a reference genome needs an internet connection: connect and try again, '
            + 'or choose "Use my own files" with a reference genome and proteome you already have.');
        return;
    }
    tasks.clear();
    $('#run-error').hidden = true;
    $('#tasks').replaceChildren();
    $('#results').replaceChildren();
    $('#log').textContent = '';
    setRunning(true);
    try {
        await window.qs.startRun(form);
    } catch (ex) {
        setRunning(false);
        alert(ex.message);
    }
};
$('#btn-cancel').onclick = () => window.qs.cancelRun();

function setRunning(on, label) {
    $('#btn-run').disabled = on;
    $('#btn-cancel').hidden = !on;
    const st = $('#run-status');
    st.textContent = label || (on ? 'running' : 'idle');
    st.className = `pill ${on ? 'warn' : ''}`;
}

window.qs.on('run:log', (text) => appendLog(text));

// what went wrong, from Nextflow's own report (lib/nferror.js), instead of
// leaving the user to find it in the log
window.qs.on('run:error', (err) => showError(err));

function showError(err) {
    const box = $('#run-error');
    box.querySelector('.err-hint').textContent = err.hint || '';
    box.querySelector('.err-hint').hidden = !err.hint;
    box.querySelector('.title').textContent = err.title;
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
    box.hidden = false;
}
window.qs.on('run:note', (text) => {
    $('#run-note').textContent = text;
    $('#run-note').hidden = !text;
});

function appendLog(text) {
    const log = $('#log');
    const atBottom = log.scrollTop + log.clientHeight >= log.scrollHeight - 4;
    log.textContent += text;
    if (atBottom) log.scrollTop = log.scrollHeight;
}

window.qs.on('run:task', (t) => {
    let li = tasks.get(t.id);
    if (!li) {
        li = document.createElement('li');
        li.innerHTML = '<span class="st"></span><span class="nm"></span><span class="meta"></span>';
        tasks.set(t.id, li);
        $('#tasks').append(li);
    }
    li.querySelector('.st').textContent = t.status;
    li.querySelector('.st').className = `st ${t.status}`;
    li.querySelector('.nm').textContent = t.name;
    const meta = [];
    if (t.duration) meta.push(`${(t.duration / 1000).toFixed(1)} s`);
    if (t.peakRss) meta.push(`${fmtBytes(t.peakRss)} RSS`);
    li.querySelector('.meta').textContent = meta.join(' · ');
});

window.qs.on('run:done', ({ code, cancelled, results }) => {
    $('#run-note').hidden = true;
    setRunning(false, cancelled ? 'cancelled' : code === 0 ? 'finished' : `failed (exit ${code})`);
    $('#run-status').className = `pill ${cancelled ? 'warn' : code === 0 ? 'ok' : 'bad'}`;
    // tasks killed on the way out don't always get a final weblog event
    for (const li of tasks.values()) {
        const st = li.querySelector('.st');
        if (['SUBMITTED', 'RUNNING'].includes(st.textContent)) {
            st.textContent = 'ABORTED';
            st.className = 'st FAILED';
        }
    }
    const box = $('#results');
    box.replaceChildren();
    for (const f of results) {
        const b = document.createElement('button');
        b.className = 'primary';
        b.textContent = `Open ${f.split('/').pop()}`;
        b.onclick = () => window.qs.openResult(f);
        box.append(b);
    }
    if (code === 0 && outdir) {
        const b = document.createElement('button');
        b.textContent = 'Show output folder';
        b.onclick = () => window.qs.reveal(outdir);
        box.append(b);
    }
    if (results.length === 1) window.qs.openResult(results[0]);
});

function fmtBytes(n) {
    const u = ['B', 'KB', 'MB', 'GB', 'TB'];
    let i = 0;
    while (n >= 1024 && i < u.length - 1) { n /= 1024; i++; }
    return `${n.toFixed(i ? 1 : 0)} ${u[i]}`;
}

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

loadSettings();
checkEnv().then(() => {
    if (env.nextflow.needsSetup) setupRuntime();
});
showUpdate();
