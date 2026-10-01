// Runs the pipeline as a child `nextflow run` and reports progress.
//
// Progress comes from Nextflow's own weblog: it POSTs one
// JSON event per workflow/task state change to a URL we listen on, so we
// get structured task states without scraping the console log. The console
// log is still streamed for the user to read.

const { spawn } = require('node:child_process');
const http = require('node:http');
const fs = require('node:fs');
const path = require('node:path');
const { CONDA_CACHE, WEBLOG_PLUGIN, pipelinePlugins } = require('./runtime');
const { ErrorCollector, withNetwork } = require('./nferror');
const { checkNetwork } = require('./network');

// form field -> --param, in the order they appear in `--help`
const PARAMS = ['taxid', 'assembly', 'reference', 'proteome', 'max_rank', 'min_seq_size', 'exclude_target',
    'min_asm_gap', 'min_identity', 'max_gap', 'min_block', 'miniprot_m', 'miniprot_chunk_gb'];

// the tools' engine -> pipeline profile
const PROFILES = { conda: 'conda', docker: 'standard' };

function buildArgs(pipelineDir, form, guiConfig) {
    const args = ['run', path.join(pipelineDir, 'main.nf'), '-profile', PROFILES[form.engine]];
    for (const key of PARAMS) {
        const v = form[key];
        if (v === undefined || v === null || v === '' || v === false) continue;
        args.push(`--${key}`);
        if (v !== true) args.push(String(v));
    }
    args.push('--outdir', form.outdir, '-w', path.join(form.outdir, 'work'), '-c', guiConfig, '-ansi-log', 'false');
    if (form.resume) args.push('-resume');
    return args;
}

// config the app layers over the pipeline's for each run. The weblog is
// set here because -with-weblog is deprecated in favor of its config scope,
// with its plugin pinned so runs work offline (see runtime.js). A plugins
// block in a -c config replaces the pipeline's rather than adding to it, so
// it repeats the pipeline's own pins (nf-schema) too: otherwise Nextflow
// loads those unpinned, asking the registry for the latest version.
// The conda engine uses the runtime's micromamba (on PATH, see
// runtime.runtimeEnv) and keeps its environments under the runtime root,
// shared across runs.
function guiConfigText(form, weblogUrl) {
    const plugins = [WEBLOG_PLUGIN, ...pipelinePlugins()].map((id) => `    id '${id}'\n`).join('');
    let text = `plugins {\n${plugins}}\nweblog {\n    enabled = true\n    url = '${weblogUrl}'\n}\n`;
    if (form.engine === 'conda') {
        text += `conda {\n    useMicromamba = true\n    cacheDir = '${CONDA_CACHE}'\n}\n`;
    }
    return text;
}

function startWeblog(onEvent) {
    return new Promise((resolve) => {
        const server = http.createServer((req, res) => {
            let body = '';
            req.on('data', (c) => { body += c; });
            req.on('end', () => {
                res.end();
                try { onEvent(JSON.parse(body)); } catch { /* ignore non-JSON */ }
            });
        });
        server.listen(0, '127.0.0.1', () => resolve(server));
    });
}

function findResults(outdir) {
    const dir = path.join(outdir, 'synteny');
    try {
        return fs.readdirSync(dir).filter((f) => f.endsWith('.synteny.interactive.html')).map((f) => path.join(dir, f));
    } catch {
        return [];
    }
}

class Run {
    constructor({ pipelineDir, launch, form, send }) {
        Object.assign(this, { pipelineDir, launch, form, send });
        this.child = null;
        this.errors = new ErrorCollector();
        this.partial = { stdout: '', stderr: '' };
    }

    async start() {
        fs.mkdirSync(this.form.outdir, { recursive: true });
        this.server = await startWeblog((ev) => this.onWeblog(ev));
        const guiConfig = path.join(this.form.outdir, '.gui.config');
        fs.writeFileSync(guiConfig, guiConfigText(this.form, `http://127.0.0.1:${this.server.address().port}`));
        const args = buildArgs(this.pipelineDir, this.form, guiConfig);
        this.send('run:log', `$ nextflow ${args.join(' ')}\n`);
        // cwd = outdir, so .nextflow/ (the -resume cache) and .nextflow.log live with the run
        this.child = spawn(this.launch.bin, args, { cwd: this.form.outdir, env: this.launch.env });
        this.child.stdout.on('data', (d) => this.onOutput('stdout', String(d)));
        this.child.stderr.on('data', (d) => this.onOutput('stderr', String(d)));
        this.child.on('error', (err) => this.finish(-1, err.message));
        this.child.on('exit', (code, signal) => this.finish(code, signal));
    }

    // Both streams go to the log as-is, and line by line to the error
    // collector (lib/nferror.js). The weblog never reports tasks that -resume
    // skips, so take those from the log's "[4a/3f058f] Cached process > NAME"
    // lines. It doesn't report conda env builds either (minutes, on a first
    // run): "Creating env using micromamba: <yml> [cache <dir>]" starts one,
    // and the next submitted task means it's done.
    onOutput(stream, text) {
        this.send('run:log', text);
        const lines = (this.partial[stream] + text).split('\n');
        this.partial[stream] = lines.pop();
        for (const line of lines) this.onLine(line);
    }

    onLine(line) {
        this.errors.push(line);
        const m = line.match(/^\[(\w+\/\w+)\] Cached process > (.+)$/);
        if (m) this.send('run:task', { id: `cached:${m[1]}`, name: m[2], status: 'CACHED' });
        const env = line.match(/^Creating env using \w+: (\S+)/);
        if (env) this.send('run:note', `Building the ${path.basename(env[1])} environment (first run only)…`);
        else if (/ Submitted process > /.test(line)) this.send('run:note', '');
    }

    onWeblog(ev) {
        if (ev.trace) {
            const t = ev.trace;
            this.send('run:task', { id: t.task_id, name: t.name, process: t.process, status: t.status, exit: t.exit,
                duration: t.realtime, peakRss: t.peak_rss, workdir: t.workdir });
        } else {
            this.send('run:event', { event: ev.event, runName: ev.runName });
        }
    }

    async finish(code, detail) {
        if (this.done) return;
        this.done = true;
        this.server?.close();
        for (const rest of Object.values(this.partial)) if (rest) this.onLine(rest);
        if (code !== 0 && !this.cancelled) {
            const error = code === -1 ? { title: 'Nextflow could not be started', message: detail } : this.errors.result(code);
            this.send('run:error', { ...withNetwork(error, await checkNetwork()), log: path.join(this.form.outdir, '.nextflow.log') });
        }
        this.send('run:done', { code, detail, cancelled: !!this.cancelled, results: code === 0 ? findResults(this.form.outdir) : [] });
    }

    // SIGINT is what Ctrl-C sends: nextflow then kills its tasks (and their
    // containers) and exits on its own
    cancel() {
        if (!this.child || this.done) return;
        this.cancelled = true;
        this.child.kill('SIGINT');
    }
}

module.exports = { Run, buildArgs, findResults, guiConfigText, PARAMS };
