// The Nextflow to launch (managed runtime, or one set in Settings) and,
// separately, whether Docker is running.
//
// An app started from Finder/the Dock gets launchd's minimal PATH, so the
// login shell's PATH is merged in first: that's how Docker Desktop or
// Homebrew installs are found.

const { execFile } = require('node:child_process');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const runtime = require('./runtime');

function run(cmd, args, opts = {}) {
    return new Promise((resolve) => {
        execFile(cmd, args, { timeout: 20000, maxBuffer: 16 << 20, ...opts }, (err, stdout, stderr) => {
            resolve({ ok: !err, stdout: String(stdout), stderr: String(stderr) });
        });
    });
}

async function loginShellPath() {
    const shell = process.env.SHELL || '/bin/zsh';
    const marker = '__QS_PATH__';
    const r = await run(shell, ['-ilc', `printf '${marker}%s${marker}' "$PATH"`], { timeout: 10000 });
    const m = r.stdout.match(new RegExp(`${marker}(.*)${marker}`));
    return m ? m[1] : null;
}

async function fixPath() {
    const shellPath = await loginShellPath();
    const parts = [...(shellPath ? shellPath.split(':') : []), ...(process.env.PATH || '').split(':'),
        '/opt/homebrew/bin', '/usr/local/bin'];
    process.env.PATH = [...new Set(parts.filter(Boolean))].join(':');
}

function which(name) {
    for (const dir of (process.env.PATH || '').split(path.delimiter)) {
        const p = path.join(dir, name);
        try {
            fs.accessSync(p, fs.constants.X_OK);
            return p;
        } catch { /* keep looking */ }
    }
    return null;
}

// settings' nextflowPath/javaHome win over the managed runtime. The env
// always carries the runtime's micromamba, which the conda engine needs
// whichever Nextflow runs.
async function nextflowLaunch(settings) {
    // the app pins its Nextflow, so the "newer version available" notice
    // only misleads (and offline, prints a curl error)
    const env = { ...runtime.runtimeEnv(), NXF_DISABLE_CHECK_LATEST: 'true' };
    if (settings.javaHome) {
        env.JAVA_HOME = settings.javaHome;
        env.NXF_JAVA_HOME = settings.javaHome;
    }
    if (settings.nextflowPath) {
        const r = await run(settings.nextflowPath, ['-v'], { env, timeout: 60000 });
        const m = r.stdout.match(/nextflow version (\S+)/);
        return m ? { ok: true, bin: settings.nextflowPath, version: m[1], source: 'settings', env }
            : { ok: false, bin: settings.nextflowPath, source: 'settings', message: 'the Nextflow set in Settings failed to start',
                detail: (r.stdout + r.stderr).trim().split('\n').slice(-4).join('\n') };
    }
    const s = await runtime.status();
    return s.ok ? { ok: true, bin: s.nextflow, version: s.version, source: 'managed', env }
        : { ok: false, source: 'managed', needsSetup: true, message: 'Runtime not set up', detail: s.detail };
}

async function checkDocker() {
    const bin = which('docker');
    if (!bin) return { ok: false, message: 'Docker not installed' };
    const r = await run(bin, ['info', '--format', '{{json .}}']);
    if (!r.ok) return { ok: false, message: 'Docker not running' };
    try {
        const info = JSON.parse(r.stdout);
        return { ok: true, version: info.ServerVersion, memBytes: info.MemTotal };
    } catch {
        return { ok: true, version: '?' };
    }
}

async function checkAll(settings) {
    const nextflow = await nextflowLaunch(settings);
    return { nextflow, host: { memBytes: os.totalmem() }, runtimeRoot: runtime.ROOT, home: os.homedir() };
}

module.exports = { fixPath, checkAll, checkDocker };
