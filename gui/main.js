const { app, BrowserWindow, dialog, ipcMain, session, shell } = require('electron');
const fs = require('node:fs');
const path = require('node:path');
const env = require('./lib/env');
const { Run, findResults } = require('./lib/runner');
const runtime = require('./lib/runtime');
const inputs = require('./lib/inputs');
const updates = require('./lib/updates');
const network = require('./lib/network');
const taxa = require('./lib/taxa');
const { smokeTest } = require('./lib/smoke');
const { History } = require('./lib/history');
const record = require('./lib/record');

// the pipeline: this repo in development, a bundled copy once packaged
const PIPELINE_DIR = app.isPackaged ? path.join(process.resourcesPath, 'pipeline') : path.resolve(__dirname, '..');
// read by micromamba itself, so it can't live inside an asar archive
const RUNTIME_YML = app.isPackaged ? path.join(process.resourcesPath, 'runtime.yml') : path.join(__dirname, 'runtime.yml');
runtime.configure({ pipelineDir: PIPELINE_DIR });
// settings and the run history; QS_USER_DATA moves them, for testing
if (process.env.QS_USER_DATA) app.setPath('userData', process.env.QS_USER_DATA);
const SETTINGS_FILE = () => path.join(app.getPath('userData'), 'settings.json');

let mainWindow = null;
let currentRun = null;
let lastEnv = null;
let history = null;

function loadSettings() {
    try { return JSON.parse(fs.readFileSync(SETTINGS_FILE(), 'utf8')); } catch { return {}; }
}

function saveSettings(s) {
    fs.mkdirSync(path.dirname(SETTINGS_FILE()), { recursive: true });
    fs.writeFileSync(SETTINGS_FILE(), JSON.stringify(s, null, 2));
}

function send(channel, payload) {
    if (mainWindow && !mainWindow.isDestroyed()) mainWindow.webContents.send(channel, payload);
}

// macOS: the traffic lights sit inset over the sidebar, which doubles as the
// title bar (renderer/style.css marks the drag regions); elsewhere the
// window keeps its normal frame
function createMainWindow() {
    mainWindow = new BrowserWindow({
        width: 1200, height: 800, minWidth: 980, minHeight: 640, title: 'quick_synteny', backgroundColor: '#FFFFFF',
        ...(process.platform === 'darwin' ? { titleBarStyle: 'hiddenInset', trafficLightPosition: { x: 18, y: 18 } } : {}),
        webPreferences: { preload: path.join(__dirname, 'preload.js'), contextIsolation: true, sandbox: true },
    });
    mainWindow.loadFile(path.join(__dirname, 'renderer', 'index.html'));
}

// a result page is the pipeline's own self-contained HTML, opened as-is in
// its own sandboxed window with no access to the app's IPC
function openResult(file) {
    const win = new BrowserWindow({
        width: 1600, height: 1000, title: path.basename(file),
        webPreferences: { contextIsolation: true, sandbox: true },
    });
    win.webContents.setWindowOpenHandler(({ url }) => {
        shell.openExternal(url);
        return { action: 'deny' };
    });
    win.loadFile(file);
}

ipcMain.handle('env:check', async () => {
    lastEnv = await env.checkAll(loadSettings());
    return lastEnv;
});
ipcMain.handle('docker:check', () => env.checkDocker());
// the bundled pipeline's parameter schema: the form's tips, defaults and limits
ipcMain.handle('schema:get', () => JSON.parse(fs.readFileSync(path.join(PIPELINE_DIR, 'nextflow_schema.json'), 'utf8')));
ipcMain.handle('net:check', () => network.checkNetwork());
ipcMain.handle('runtime:setup', async () => {
    if (currentRun && !currentRun.done) throw new Error('a run is in progress');
    await runtime.setup(RUNTIME_YML, (text) => send('setup:log', text));
});
// once per launch, from the page: the newest GitHub release with a file for
// this machine, unless the user turned the check off or dismissed that
// version. QS_RELEASES_URL points it elsewhere, for testing.
ipcMain.handle('update:check', async () => {
    const s = loadSettings();
    if (s.checkUpdates === false) return null;
    const u = await updates.checkForUpdate({ current: app.getVersion(), url: process.env.QS_RELEASES_URL || undefined });
    return u && u.version !== s.dismissedUpdate ? { ...u, current: app.getVersion() } : null;
});
ipcMain.handle('update:dismiss', (_e, version) => saveSettings({ ...loadSettings(), dismissedUpdate: version }));
ipcMain.handle('update:open', (_e, url) => {
    if (!updates.isRepoUrl(url)) throw new Error('not a quick_synteny release link');
    return shell.openExternal(url);
});
ipcMain.handle('settings:get', () => loadSettings());
ipcMain.handle('settings:set', (_e, s) => saveSettings({ ...loadSettings(), ...s }));

// a picked or dropped file: its size, and whether it looks like FASTA
function inspectFile(file) {
    const st = fs.statSync(file);
    if (!st.isFile()) throw new Error(`${path.basename(file)} is a folder, not a file`);
    return { path: file, bytes: st.size, kind: inputs.sniffFasta(file) };
}
ipcMain.handle('dialog:file', async (_e, { title, filters }) => {
    const r = await dialog.showOpenDialog(mainWindow, { title, properties: ['openFile'], filters });
    return r.canceled ? null : inspectFile(r.filePaths[0]);
});
ipcMain.handle('file:inspect', (_e, file) => inspectFile(file));
ipcMain.handle('dialog:dir', async (_e, { title }) => {
    const r = await dialog.showOpenDialog(mainWindow, { title, properties: ['openDirectory', 'createDirectory'] });
    return r.canceled ? null : inputs.inspectOutdir(r.filePaths[0]);
});
ipcMain.handle('outdir:suggest', (_e, { assembly, taxid }) => inputs.suggestRunName(assembly, taxid));
// where a run goes unless the user picks a folder: Documents/quick_synteny/
// <target>_vs_<taxid or reference>, with -2, -3... once that's taken
ipcMain.handle('outdir:default', (_e, { assembly, taxid, reference }) => {
    const base = path.join(app.getPath('documents'), 'quick_synteny',
        inputs.suggestRunName(assembly, taxid || record.fileStem(reference)));
    let dir = base;
    for (let n = 2; fs.existsSync(dir); n++) dir = `${base}-${n}`;
    return { path: dir, writable: true, empty: true, isRun: false, isNew: true };
});

// scientific/common name -> taxa to pick from (lib/taxa.js)
ipcMain.handle('taxon:suggest', (_e, query) => taxa.searchTaxa(query));

// taxid -> its name and rank, or { found: false } when NCBI has no such
// taxon (a 200 carrying `errors`, not an HTTP error). Throws only when NCBI
// can't be reached, which the page treats as "not checked".
ipcMain.handle('taxon:lookup', async (_e, taxid) => {
    const url = `https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/${encodeURIComponent(taxid)}`;
    const r = await fetch(url, { signal: AbortSignal.timeout(10000) });
    if (!r.ok) throw new Error(`NCBI returned ${r.status}`);
    const node = (await r.json()).taxonomy_nodes?.[0];
    const t = node?.taxonomy;
    if (!t) return { found: false };
    return { found: true, taxid: String(t.tax_id), name: t.organism_name, rank: t.rank?.toLowerCase(), common: t.common_name };
});

// ---------- runs ----------

function publishRun(r) {
    history.put(r);
    send('run:update', r);
}

// a finished run's result page, rendered in a hidden window and captured for
// the run screen: the panels, below the page's title and controls (about
// the top 220 px at this width). Bokeh lays the page out after it loads,
// hence the wait.
async function captureThumbnail(html, out) {
    const win = new BrowserWindow({ show: false, width: 1600, height: 1000, webPreferences: { contextIsolation: true, sandbox: true } });
    try {
        await win.loadFile(html);
        await new Promise((r) => setTimeout(r, 3000));
        const image = await win.webContents.capturePage({ x: 0, y: 220, width: 1600, height: 780 });
        if (image.isEmpty()) return false;
        fs.mkdirSync(path.dirname(out), { recursive: true });
        fs.writeFileSync(out, image.resize({ width: 1200 }).toPNG());
        return true;
    } catch {
        return false;
    } finally {
        win.destroy();
    }
}

ipcMain.handle('runs:list', () => history.list());
ipcMain.handle('runs:remove', (_e, id) => {
    if (currentRun && !currentRun.done && currentRun.record.id === id) throw new Error('that run is still going');
    history.remove(id);
    send('runs:changed', history.list());
});
ipcMain.handle('thumb:get', (_e, id) => {
    try {
        return `data:image/png;base64,${fs.readFileSync(history.thumbnailPath(id)).toString('base64')}`;
    } catch {
        return null;
    }
});

ipcMain.handle('run:start', async (_e, form, meta = {}) => {
    if (currentRun && !currentRun.done) throw new Error('a run is already in progress');
    if (!lastEnv?.nextflow?.ok) throw new Error('Nextflow is not available');
    // checked now rather than trusting the page's last check: Docker
    // Desktop can have been quit or started since
    if (form.engine === 'docker') {
        const docker = await env.checkDocker();
        if (!docker.ok) throw new Error(`${docker.message} — start it, or choose Tools → Conda.`);
    }
    // offline, Nextflow shouldn't try to reach anything of its own (version
    // checks, plugins); what a run itself needs online is the page's call
    const online = (await network.checkNetwork()).online;
    const launch = online ? lastEnv.nextflow : { ...lastEnv.nextflow, env: { ...lastEnv.nextflow.env, NXF_OFFLINE: 'true' } };
    const r = record.newRecord(form, { targetSpecies: meta.targetSpecies });
    const tasks = new Map();
    // the runner's events become the record's; only the log goes to the
    // page as it is
    const onEvent = async (channel, payload) => {
        if (channel === 'run:log') return send('run:log', { id: r.id, text: payload });
        if (channel === 'run:done') {
            record.finishRecord(r, tasks, payload);
            publishRun(r);
            if (r.status === 'finished' && r.results[0]) {
                r.hasThumbnail = await captureThumbnail(r.results[0], history.thumbnailPath(r.id));
                publishRun(r);
            }
            if (r.results.length === 1) openResult(r.results[0]);
            return;
        }
        if (record.applyEvent(r, tasks, channel, payload)) publishRun(r);
    };
    currentRun = new Run({ pipelineDir: PIPELINE_DIR, launch, form, send: onEvent });
    currentRun.record = r;
    publishRun(r);
    try {
        await currentRun.start();
    } catch (e) {
        currentRun.done = true;
        r.status = 'failed';
        r.finishedAt = new Date().toISOString();
        r.error = { title: 'The run could not start', message: e.message };
        publishRun(r);
    }
    return r.id;
});
ipcMain.handle('run:cancel', () => currentRun?.cancel());

ipcMain.handle('result:open', (_e, file) => openResult(file));
ipcMain.handle('result:pick', async () => {
    const r = await dialog.showOpenDialog(mainWindow, {
        title: 'Open a quick_synteny result', properties: ['openFile', 'openDirectory'],
        filters: [{ name: 'Synteny page', extensions: ['html'] }],
    });
    if (r.canceled) return;
    const p = r.filePaths[0];
    // a run's outdir works too
    const files = fs.statSync(p).isDirectory() ? findResults(p) : [p];
    if (!files.length) throw new Error(`no *.synteny.interactive.html under ${p}/synteny`);
    files.forEach(openResult);
});
ipcMain.handle('path:reveal', (_e, p) => shell.openPath(p));

app.whenReady().then(async () => {
    // `--smoke-test`: check the app works here, print the result, exit
    // (lib/smoke.js; the CI workflow runs it on every build)
    if (process.argv.includes('--smoke-test')) {
        await env.fixPath();
        const ok = await smokeTest({ pipelineDir: PIPELINE_DIR, runtimeYml: RUNTIME_YML, log: (t) => process.stdout.write(t) });
        app.exit(ok ? 0 : 1);
        return;
    }
    // a packaged app gets its icon from the bundle; from source, the Dock
    // would show Electron's
    if (!app.isPackaged && process.platform === 'darwin') app.dock.setIcon(path.join(__dirname, 'build', 'icon-macos.png'));
    await env.fixPath();
    // a result page's exports (SVG/PNG/JPEG/TSV) are blob downloads: ask
    // where to save them, starting from Downloads under the page's own name
    // (QS_DOWNLOAD_DIR skips the dialog, for scripted testing)
    session.defaultSession.on('will-download', (_e, item) => {
        if (process.env.QS_DOWNLOAD_DIR) item.setSavePath(path.join(process.env.QS_DOWNLOAD_DIR, item.getFilename()));
        else item.setSaveDialogOptions({ defaultPath: path.join(app.getPath('downloads'), item.getFilename()) });
    });
    history = new History(app.getPath('userData'));
    // finished discovery runs from before the selections were read at the
    // end, whose log lines the runner couldn't parse (a rank of several
    // words): fill in what was found, and their titles
    for (const r of history.list()) if (r.status === 'finished') record.fillFound(r);
    history.save(true);
    createMainWindow();
});

// don't leave nextflow (and its containers) running behind a closed app
app.on('before-quit', () => currentRun?.cancel());
app.on('window-all-closed', () => app.quit());
