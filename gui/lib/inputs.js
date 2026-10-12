// Cheap checks that catch a wrong file or folder before a run starts.

const fs = require('node:fs');
const path = require('node:path');
const zlib = require('node:zlib');

// '>' after any leading whitespace, possibly gzipped (bgzip has the same
// magic). Returns 'fasta', 'gzip-fasta' or 'other'.
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
            // only the head was read: Z_SYNC_FLUSH inflates what's there
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

// whether `dir` is writable, and whether it holds a previous run (which a
// rerun or -resume reuses) or something else
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

// <target>_vs_<taxid>, filesystem-safe
function suggestRunName(assembly, taxid) {
    const safe = (s) => s.replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
    const stem = assembly ? safe(path.basename(assembly).replace(/\.gz$/, '').replace(/\.[^.]+$/, '')) : '';
    return `${stem || 'quick_synteny'}_vs_${safe(String(taxid || '')) || 'reference'}`;
}

module.exports = { sniffFasta, inspectOutdir, suggestRunName };
