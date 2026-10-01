// Cheap checks on the form's inputs, so an obviously wrong file or folder
// is caught before a run rather than minutes into one.

const fs = require('node:fs');
const path = require('node:path');
const zlib = require('node:zlib');

// First bytes of a FASTA: '>' after any leading whitespace, possibly under
// gzip (bgzip too -- same magic). Returns 'fasta', 'gzip-fasta' or 'other'.
function sniffFasta(file) {
    const fd = fs.openSync(file, 'r');
    const buf = Buffer.alloc(65536);
    const n = fs.readSync(fd, buf, 0, buf.length, 0);
    fs.closeSync(fd);
    let head = buf.subarray(0, n);
    let gz = false;
    if (head[0] === 0x1f && head[1] === 0x8b) {
        gz = true;
        try {
            // a truncated gzip stream: Z_SYNC_FLUSH returns what inflates
            // instead of throwing at the missing end
            head = zlib.gunzipSync(head, { finishFlush: zlib.constants.Z_SYNC_FLUSH });
        } catch {
            return 'other';
        }
    }
    const first = head.toString('latin1').trimStart()[0];
    if (first !== '>') return 'other';
    return gz ? 'gzip-fasta' : 'fasta';
}

// What running into `dir` would mean: whether it can be written, and
// whether it already holds something other than a previous run (whose
// synteny/, pipeline_info/ and .nextflow/ a rerun or -resume reuses)
function inspectOutdir(dir) {
    let writable = true;
    try { fs.accessSync(dir, fs.constants.W_OK); } catch { writable = false; }
    let entries;
    try {
        entries = fs.readdirSync(dir).filter((e) => e !== '.DS_Store');
    } catch {
        // macOS privacy settings can refuse a folder the permissions allow
        return { path: dir, writable: false, empty: false, isRun: false, unreadable: true };
    }
    const isRun = ['.nextflow', 'synteny', 'pipeline_info'].some((e) => entries.includes(e));
    return { path: dir, writable, empty: entries.length === 0, isRun };
}

// a folder name for a new run inside a chosen folder: <target>_vs_<taxid>
function suggestRunName(assembly, taxid) {
    const stem = assembly ? path.basename(assembly).replace(/\.gz$/, '').replace(/\.[^.]+$/, '') : 'quick_synteny';
    return `${stem}_vs_${taxid || 'reference'}`.replace(/[^\w.-]+/g, '_');
}

module.exports = { sniffFasta, inspectOutdir, suggestRunName };
