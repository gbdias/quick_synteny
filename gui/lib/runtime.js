// The app's own managed runtime, so running the pipeline needs nothing
// preinstalled: a pinned micromamba binary, and from it a conda env with
// Nextflow and Java (runtime.yml). The pipeline's tool environments
// (envs/*.yml, the conda profile) are then built by Nextflow with the same
// micromamba, on first use.
//
// Everything lives under ~/.quick_synteny ($QS_HOME overrides it) rather
// than Electron's userData ("~/Library/Application Support/..." on macOS):
// conda environments don't work reliably from a path with a space in it.
//
// A packaged app ships micromamba in its resources (see
// scripts/fetch-micromamba.js); running from source downloads it into the
// runtime root instead.

const { spawn, execFile } = require('node:child_process');
const crypto = require('node:crypto');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');

const ROOT = process.env.QS_HOME || path.join(os.homedir(), '.quick_synteny');
const MAMBA_ROOT = path.join(ROOT, 'mamba');             // micromamba's package cache
const RUNTIME_PREFIX = path.join(ROOT, 'runtime');       // nextflow + openjdk
const CONDA_CACHE = path.join(ROOT, 'conda_envs');       // Nextflow's per-envs/*.yml envs
const NXF_HOME = path.join(ROOT, 'nextflow');            // the app's Nextflow home: plugins, history
// The weblog the app reads progress from is a Nextflow plugin. Unpinned,
// every run asks registry.nextflow.io which version to use, and offline
// (NXF_OFFLINE) Nextflow refuses an unpinned plugin outright. So it's pinned
// here, set in each run's config (lib/runner.js), and installed by setup,
// leaving setup as the only step that needs a connection. Bump it together
// with runtime.yml's nextflow.
const WEBLOG_PLUGIN = 'nf-weblog@1.2.0';

// The pipeline's own plugins (nf-schema), as its nextflow.config pins them
// in `plugins { id '...' }`, read from there so their versions live in one
// place; setup installs them with the weblog, for the same reason.
let pipelineDir = null;
function configure(opts) {
    pipelineDir = opts.pipelineDir;
}

function pipelinePlugins() {
    if (!pipelineDir) return [];
    const config = fs.readFileSync(path.join(pipelineDir, 'nextflow.config'), 'utf8');
    const block = config.match(/^plugins\s*\{([^}]*)\}/m);
    return block ? [...block[1].matchAll(/id\s+['"]([\w-]+@[\w.-]+)['"]/g)].map((m) => m[1]) : [];
}

const PLUGINS = () => [WEBLOG_PLUGIN, ...pipelinePlugins()];
const pluginDir = (id) => path.join(NXF_HOME, 'plugins', id.replace('@', '-'));
const missingPlugins = () => PLUGINS().filter((id) => !fs.existsSync(pluginDir(id)));
const BUNDLED_MICROMAMBA = process.resourcesPath && path.join(process.resourcesPath, 'micromamba', 'micromamba');
const MICROMAMBA = BUNDLED_MICROMAMBA && fs.existsSync(BUNDLED_MICROMAMBA) ? BUNDLED_MICROMAMBA
    : path.join(ROOT, 'bin', 'micromamba');

// conda-forge's micromamba build per platform, with its published sha256
const MICROMAMBA_VERSION = '2.9.0';
const MICROMAMBA_BUILDS = {
    'darwin-arm64': ['osx-arm64', '500f5074feb8d02c4296ef9921c3650ed2874171805a9fbb8fbb53896433646b'],
    'darwin-x64': ['osx-64', '0426ecdc41636d369f57b8fe6acbf4385a69eca45b56d9ee7d3a840a9965d44f'],
    'linux-x64': ['linux-64', '8761c382127e6363bd9e0a2451aa3ef90d071a79133f736e2f759a3bf13040dd'],
    'linux-arm64': ['linux-aarch64', 'e705ffeed90ce0659eb546e4b1e1028c9eaf0bc9cc854867b19ac5ce0ba5852f'],
};

function run(cmd, args, opts = {}) {
    return new Promise((resolve) => {
        execFile(cmd, args, { timeout: 60000, ...opts }, (err, stdout, stderr) => {
            resolve({ ok: !err, stdout: String(stdout), stderr: String(stderr) });
        });
    });
}

// micromamba prints long solves and downloads; stream them to the UI
function stream(cmd, args, env, onLog) {
    return new Promise((resolve, reject) => {
        const child = spawn(cmd, args, { env });
        child.stdout.on('data', (d) => onLog(String(d)));
        child.stderr.on('data', (d) => onLog(String(d)));
        child.on('error', reject);
        child.on('exit', (code) => (code === 0 ? resolve() : reject(new Error(`${path.basename(cmd)} exited with ${code}`))));
    });
}

function javaHome() {
    // conda-forge's openjdk installs the JDK under lib/jvm on macOS and Linux
    const jvm = path.join(RUNTIME_PREFIX, 'lib', 'jvm');
    return fs.existsSync(path.join(jvm, 'bin', 'java')) ? jvm : RUNTIME_PREFIX;
}

// micromamba settings for every call, the app's and Nextflow's: packages
// cached under ROOT instead of ~/micromamba, and the user's own
// ~/.condarc/~/.mambarc ignored -- otherwise their channels (Anaconda's
// commercial `defaults`, typically) join each env's solve, and the envs
// would differ from one user to the next
const MAMBA_ENV = { MAMBA_ROOT_PREFIX: MAMBA_ROOT, MAMBA_NO_RC: 'true' };

// env for anything the runtime launches: micromamba and the runtime's
// bin first on PATH (Nextflow finds micromamba there)
function runtimeEnv(base = process.env) {
    return {
        ...base,
        PATH: [path.dirname(MICROMAMBA), path.join(RUNTIME_PREFIX, 'bin'), base.PATH].join(path.delimiter),
        ...MAMBA_ENV,
        JAVA_HOME: javaHome(),
        NXF_JAVA_HOME: javaHome(),
        NXF_HOME,
    };
}

async function status() {
    const hasMamba = fs.existsSync(MICROMAMBA);
    const nextflow = path.join(RUNTIME_PREFIX, 'bin', 'nextflow');
    if (!hasMamba || !fs.existsSync(nextflow) || missingPlugins().length) return { ok: false, root: ROOT, hasMamba };
    const r = await run(nextflow, ['-v'], { env: runtimeEnv() });
    const m = r.stdout.match(/nextflow version (\S+)/);
    return m ? { ok: true, root: ROOT, nextflow, version: m[1], condaCache: CONDA_CACHE }
        : { ok: false, root: ROOT, hasMamba, detail: (r.stdout + r.stderr).trim() };
}

// downloads conda-forge's micromamba for `platformArch` (e.g. darwin-arm64)
// to `dest`, checking its sha256
async function downloadMicromamba(platformArch, dest, onLog) {
    const build = MICROMAMBA_BUILDS[platformArch];
    if (!build) throw new Error(`no micromamba build for ${platformArch}`);
    const [subdir, sha256] = build;
    const url = `https://conda.anaconda.org/conda-forge/${subdir}/micromamba-${MICROMAMBA_VERSION}-0.tar.bz2`;
    onLog(`Downloading ${url}\n`);
    const r = await fetch(url);
    if (!r.ok) throw new Error(`micromamba download failed: HTTP ${r.status}`);
    const buf = Buffer.from(await r.arrayBuffer());
    const got = crypto.createHash('sha256').update(buf).digest('hex');
    if (got !== sha256) throw new Error(`micromamba checksum mismatch (got ${got})`);
    const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'micromamba-'));
    const tarball = path.join(tmp, 'micromamba.tar.bz2');
    fs.writeFileSync(tarball, buf);
    // node has no bzip2; macOS's and Linux's tar both read it
    const x = await run('tar', ['-xjf', tarball, '-C', tmp, 'bin/micromamba']);
    if (!x.ok) throw new Error(`could not unpack micromamba: ${x.stderr}`);
    fs.mkdirSync(path.dirname(dest), { recursive: true });
    fs.copyFileSync(path.join(tmp, 'bin', 'micromamba'), dest);
    fs.chmodSync(dest, 0o755);
    fs.rmSync(tmp, { recursive: true, force: true });
    onLog(`micromamba ${MICROMAMBA_VERSION} -> ${dest} (sha256 verified)\n`);
}

// does only what's missing: micromamba (from source), the runtime env, the
// plugins -- so a runtime set up before a plugin was needed just gains it
async function setup(runtimeYml, onLog) {
    fs.mkdirSync(ROOT, { recursive: true });
    if (MICROMAMBA === BUNDLED_MICROMAMBA) onLog(`Using the bundled micromamba ${MICROMAMBA_VERSION}\n`);
    else if (!fs.existsSync(MICROMAMBA)) await downloadMicromamba(`${process.platform}-${process.arch}`, MICROMAMBA, onLog);
    const nextflow = path.join(RUNTIME_PREFIX, 'bin', 'nextflow');
    if (!fs.existsSync(nextflow)) {
        onLog(`Creating ${RUNTIME_PREFIX} from ${path.basename(runtimeYml)}…\n`);
        await stream(MICROMAMBA, ['create', '--yes', '--prefix', RUNTIME_PREFIX, '--file', runtimeYml],
            { ...process.env, ...MAMBA_ENV }, onLog);
    }
    for (const id of missingPlugins()) {
        onLog(`Installing the Nextflow plugin ${id}…\n`);
        await stream(nextflow, ['plugin', 'install', id], runtimeEnv(), onLog);
    }
    const s = await status();
    if (!s.ok) throw new Error(`runtime installed, but nextflow does not start:\n${s.detail || ''}`);
    onLog(`Ready: Nextflow ${s.version}\n`);
    return s;
}

module.exports = { ROOT, CONDA_CACHE, WEBLOG_PLUGIN, MICROMAMBA_VERSION, configure, pipelinePlugins, status, setup, runtimeEnv,
    downloadMicromamba };
