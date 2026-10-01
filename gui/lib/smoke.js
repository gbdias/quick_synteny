// `quick_synteny --smoke-test`: checks that a packaged app works on this
// machine, without a window -- the CI workflow runs it on every build.
// It goes through the app's own code paths:
//
//   1. first-launch setup (runtime.setup): the bundled micromamba, then
//      Nextflow, Java and the weblog plugin from runtime.yml, into QS_HOME
//      (set it to a scratch folder; ~/.quick_synteny otherwise);
//   2. the Nextflow it would launch (env.checkAll);
//   3. `nextflow run main.nf --help -profile conda` on the bundled pipeline,
//      which parses main.nf, nextflow.config and conf/conda.config and
//      loads nf-schema for the help. It runs with NXF_OFFLINE, so a plugin
//      setup didn't install fails here instead of being downloaded.
//
// It doesn't run the pipeline: that needs NCBI and tens of MB of genome.

const { spawn } = require('node:child_process');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const runtime = require('./runtime');
const env = require('./env');

function capture(cmd, args, opts) {
    return new Promise((resolve) => {
        const child = spawn(cmd, args, opts);
        let out = '';
        child.stdout.on('data', (d) => { out += d; });
        child.stderr.on('data', (d) => { out += d; });
        child.on('error', (err) => resolve({ code: -1, out: err.message }));
        child.on('exit', (code) => resolve({ code, out }));
    });
}

async function smokeTest({ pipelineDir, runtimeYml, log }) {
    const step = (s) => log(`\n== ${s}\n`);
    try {
        step(`platform ${process.platform}-${process.arch}, runtime root ${runtime.ROOT}`);
        for (const f of ['main.nf', 'nextflow.config', 'nextflow_schema.json', 'conf/conda.config', 'envs/miniprot.yml',
            'bin/chain_blocks.mjs']) {
            if (!fs.existsSync(path.join(pipelineDir, f))) throw new Error(`the bundled pipeline is missing ${f}`);
        }
        for (const f of ['bin/chain_blocks.mjs', 'bin/parse_lineage.py']) {
            fs.accessSync(path.join(pipelineDir, f), fs.constants.X_OK);
        }

        step('setup');
        await runtime.setup(runtimeYml, log);

        step('nextflow');
        const { nextflow } = await env.checkAll({});
        if (!nextflow.ok) throw new Error(`nextflow doesn't start: ${nextflow.message} ${nextflow.detail || ''}`);
        log(`Nextflow ${nextflow.version} (${nextflow.source})\n`);

        step('pipeline --help');
        const cwd = fs.mkdtempSync(path.join(os.tmpdir(), 'qs-smoke-'));
        const r = await capture(nextflow.bin, ['run', path.join(pipelineDir, 'main.nf'), '--help', '-profile', 'conda'],
            { cwd, env: { ...nextflow.env, NXF_OFFLINE: 'true' } });
        log(r.out.split('\n').slice(0, 8).join('\n') + '\n…\n');
        // nf-schema's help, built from nextflow_schema.json
        if (r.code !== 0 || !/taxonomy-guided synteny plotting/.test(r.out) || !/--max_rank/.test(r.out)) {
            throw new Error(`\`nextflow run main.nf --help\` exited ${r.code}:\n${r.out.slice(-2000)}`);
        }

        log('\nSMOKE TEST PASSED\n');
        return true;
    } catch (e) {
        log(`\nSMOKE TEST FAILED: ${e.message}\n`);
        return false;
    }
}

module.exports = { smokeTest };
