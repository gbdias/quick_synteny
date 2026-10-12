// The sidebar's runs: userData/runs.json, newest first, plus thumbnails.
// Removing a run leaves its output folder alone.

const fs = require('node:fs');
const path = require('node:path');

class History {
    constructor(dir) {
        this.dir = dir;
        this.file = path.join(dir, 'runs.json');
        this.runs = [];
        try { this.runs = JSON.parse(fs.readFileSync(this.file, 'utf8')); } catch { /* none yet */ }
        // a run still "running" was cut off by the app quitting or crashing
        for (const r of this.runs) {
            if (r.status === 'running') {
                r.status = 'cancelled';
                r.note = '';
                r.steps = (r.steps || []).map((s) => (s.state === 'running' ? { ...s, state: 'stopped' } : s));
                r.interrupted = true;
            }
        }
        this.timer = null;
    }

    list() {
        return this.runs;
    }

    get(id) {
        return this.runs.find((r) => r.id === id) || null;
    }

    // add or replace; the write is debounced so a burst of task events costs one
    put(record) {
        const i = this.runs.findIndex((r) => r.id === record.id);
        if (i === -1) this.runs.unshift(record);
        else this.runs[i] = record;
        this.save();
    }

    remove(id) {
        this.runs = this.runs.filter((r) => r.id !== id);
        fs.rmSync(this.thumbnailPath(id), { force: true });
        this.save(true);
    }

    thumbnailPath(id) {
        return path.join(this.dir, 'thumbnails', `${id.replace(/[^\w-]/g, '')}.png`);
    }

    save(now = false) {
        clearTimeout(this.timer);
        const write = () => {
            fs.mkdirSync(this.dir, { recursive: true });
            fs.writeFileSync(this.file, JSON.stringify(this.runs, null, 1));
        };
        if (now) write();
        else this.timer = setTimeout(write, 300);
    }
}

module.exports = { History };
