// Turns Nextflow's console error report into something to show the user.
//
// Nextflow prints one block per fatal error, starting at "ERROR ~ " and
// ending at " -- Check '.nextflow.log' file for details". It takes one of
// two forms:
//
//   ERROR ~ ERROR: --taxid must be an NCBI taxid ...     (the pipeline's own
//                                                         exit 1, "...")
//
//   ERROR ~ Error executing process > 'RESOLVE_TAXONOMY:PARSE_LINEAGE (taxonomy.json)'
//
//   Caused by:
//     Process `...` terminated with an error exit status (1)
//   Command executed:
//     ...
//   Command exit status:
//     1
//   Command output:
//     (empty)
//   Command error:
//     ERROR: --taxid 999999999 is not a known NCBI taxid
//   Work dir:
//     /.../work/c5/47d4a6...
//
//   Tip: you can replicate the issue by ...
//
// The weblog doesn't cover the first form at all (it runs before the
// workflow starts), so this reads the console for both.
//
// Parameters that fail nextflow_schema.json are reported by nf-schema, on
// stderr, as a list of its own; the ERROR line on stdout then only says
// "Validation of pipeline parameters failed!" (and for unrecognised
// parameters there's no ERROR line at all):
//
//   The following invalid input values have been detected:
//
//   * --taxid (): 0 is less than 1 (must be an NCBI taxid, a positive whole number)
//   * --max_rank (kingdom): Expected any of [species, genus, family, order, class, phylum]
//
//    -- Check script '.../main.nf' at line: 53 or see '.nextflow.log' file for more details

const START = /^ERROR ~ /;
const SCHEMA_START = /^The following invalid input values have been detected:/;
const SCHEMA_ITEM = /^\* (.+)$/;
// nf-schema colours its list, terminal or not
const ANSI = /\x1b\[[0-9;]*m/g;
const END = /^ -- Check '.*\.nextflow\.log' file for details/;
const END_SCRIPT = /^ -- Check script /;

// "Command error:" etc. -- flush-left and ending in a colon
const SECTION = /^([A-Z][A-Za-z ]+):$/;

function sections(lines) {
    const out = {};
    let current = null;
    for (const line of lines) {
        const m = line.match(SECTION);
        if (m) {
            current = m[1];
            out[current] = [];
        } else if (/^Tip: /.test(line)) {
            current = null;
        } else if (current) {
            out[current].push(line.replace(/^ {2}/, ''));
        }
    }
    for (const k of Object.keys(out)) {
        const body = out[k];
        while (body.length && !body[0].trim()) body.shift();
        while (body.length && !body[body.length - 1].trim()) body.pop();
        out[k] = body;
    }
    return out;
}

const meaningful = (body) => (body || []).filter((l) => l.trim() && l.trim() !== '(empty)');

// a run of identical lines as one, counted: micromamba says "Retrying in 2
// seconds" over and over while it can't connect
function collapse(lines) {
    const out = [];
    for (const l of lines) {
        const last = out[out.length - 1];
        if (last && last.line === l) last.n++;
        else out.push({ line: l, n: 1 });
    }
    return out.map(({ line, n }) => (n > 1 ? `${line}  (×${n})` : line));
}

// what the tools print when there's no connection: curl/micromamba ("Could
// not resolve host: conda.anaconda.org"), datasets' Go client ("dial tcp:
// lookup api.ncbi.nlm.nih.gov: no such host"), node
const NETWORK = /could not resolve host|no such host|network is unreachable|couldn't connect to server|connection refused|connection timed out|dial tcp|ENOTFOUND|EAI_AGAIN|temporary failure in name resolution/i;
const OFFLINE_HINT = 'This looks like a network problem: check the internet connection. Finding a reference on NCBI '
    + 'always needs one, and so does the first run of each step, which downloads its tools.';
const CONDA_OFFLINE = 'Downloading them failed: the computer seems to be offline, or can\'t reach conda-forge/bioconda. '
    + 'Connect and run again; once downloaded, they\'re kept for later runs.';

// strips the "ERROR: " our own messages carry, which Nextflow's
// "ERROR ~ " prefix would otherwise double up
const clean = (s) => s.replace(/^ERROR:\s*/, '').trim();

// the line that says what went wrong: the last one that announces an
// error, else the last line. Tools often print something after it --
// datasets ends every failure with "Use datasets ... --help for detailed
// help about a command."
function errorLine(body) {
    const announced = body.filter((l) => /^\s*(error\b|critical\b|\w*(error|exception)\b:)/i.test(l));
    return announced.length ? announced[announced.length - 1] : body[body.length - 1];
}

function parse(block) {
    const headline = block[0].replace(START, '');
    const proc = headline.match(/^Error executing process > '(.+)'$/);
    if (!proc) {
        const rest = block.slice(1).filter((l) => l.trim() && !END.test(l));
        return { title: 'The pipeline stopped before running', message: clean(headline), detail: rest.join('\n'),
            hint: NETWORK.test(block.join('\n')) ? OFFLINE_HINT : undefined };
    }
    const s = sections(block.slice(1));
    const err = meaningful(s['Command error']);
    const out = meaningful(s['Command output']);
    const exit = (s['Command exit status'] || [])[0];
    const network = NETWORK.test(block.join('\n'));
    // the step never ran: Nextflow couldn't build its conda env, and says
    // so under "Caused by:" alone (micromamba's own text below it, whose
    // "package cache may be corrupted" advice is wrong when offline)
    const caused = meaningful(s['Caused by']);
    if (/^Failed to create Conda environment/.test(caused[0] || '')) {
        return {
            kind: 'conda-env',
            title: `Couldn't set up the tools for ${proc[1]}`,
            process: proc[1],
            message: network ? CONDA_OFFLINE : 'micromamba failed to build its conda environment (details below).',
            detail: collapse(caused).slice(-40).join('\n'),
        };
    }
    // a traceback or a tool's complaint ends in the line that says what
    // went wrong; fall back to the output, then to what Nextflow says caused it
    const body = err.length ? err : out.length ? out : caused.filter((l) => !/^Process `.*` terminated/.test(l));
    // datasets' own wording is one long "Error: [gateway] Post ... giving up
    // after 11 attempt(s): ... no such host" line; say it plainly, and keep
    // that line as the detail
    const ncbiDown = network && /ncbi\.nlm\.nih\.gov/.test(body.join('\n'));
    return {
        title: `Step failed: ${proc[1]}`,
        process: proc[1],
        message: ncbiDown ? 'Couldn\'t reach NCBI (api.ncbi.nlm.nih.gov).'
            : body.length ? clean(errorLine(body)) : `exited with status ${exit ?? '?'}`,
        detail: body.slice(-40).join('\n'),
        exitStatus: exit,
        workdir: (s['Work dir'] || [])[0]?.trim(),
        hint: network ? OFFLINE_HINT : undefined,
    };
}

// Feed it console lines; it keeps the first error block, plus the last
// lines seen as a fallback for failures that print no block at all.
class ErrorCollector {
    constructor() {
        this.block = null;
        this.capturing = false;
        this.schemaItems = null;
        this.inSchemaList = false;
        this.tail = [];
    }

    push(raw) {
        const line = raw.replace(ANSI, '');
        this.tail.push(line);
        if (this.tail.length > 30) this.tail.shift();
        // nf-schema's list (stderr), wherever it falls among stdout's lines
        if (SCHEMA_START.test(line)) {
            this.schemaItems = this.schemaItems || [];
            this.inSchemaList = true;
            return;
        }
        if (this.inSchemaList) {
            const item = line.match(SCHEMA_ITEM);
            if (item) this.schemaItems.push(item[1].replace(/^(--\S+) \(\):/, '$1:'));
            else if (END_SCRIPT.test(line)) this.inSchemaList = false;
            if (item || !line.trim()) return;
        }
        if (!this.block && START.test(line)) {
            this.block = [line];
            this.capturing = true;
        } else if (this.capturing) {
            if (END.test(line)) this.capturing = false;
            else this.block.push(line);
        }
    }

    result(code) {
        if (this.schemaItems && this.schemaItems.length) {
            const n = this.schemaItems.length;
            return {
                title: n === 1 ? 'A parameter is not valid' : `${n} parameters are not valid`,
                message: n === 1 ? this.schemaItems[0] : this.schemaItems.join('; '),
                detail: this.schemaItems.join('\n'),
            };
        }
        if (this.block) return parse(this.block);
        const tail = this.tail.filter((l) => l.trim());
        return {
            title: `Nextflow exited with status ${code}`,
            message: tail.length ? clean(tail[tail.length - 1]) : 'no error message was printed',
            detail: tail.slice(-15).join('\n'),
            hint: NETWORK.test(tail.join('\n')) ? OFFLINE_HINT : undefined,
        };
    }
}

// The tools don't always say it's the network: micromamba with a cached
// index only reports "Could not load repodata.json ... after retry". So a
// failed run also gets a connection check (lib/network.js), and when the
// computer is offline the error says so whatever the tools printed.
function withNetwork(error, net) {
    if (net.online && (error.kind !== 'conda-env' || net.conda)) return error;
    if (error.kind === 'conda-env') return { ...error, message: CONDA_OFFLINE };
    return { ...error, hint: OFFLINE_HINT };
}

module.exports = { ErrorCollector, parse, withNetwork, NETWORK };
