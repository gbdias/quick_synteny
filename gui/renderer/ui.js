// Shared by form.js and app.js: element helpers, icons, formatting, and
// the page's state. Classic scripts, loaded in order (index.html): Chromium
// doesn't run ES modules from file:// pages.

const $ = (sel) => document.querySelector(sel);

const state = {
    env: null,          // window.qs.checkEnv(): the Nextflow to launch, host RAM
    docker: null,       // the last Docker check, only made for the Docker engine; 'checking' meanwhile
    runs: [],           // the history, newest first (lib/history.js records)
    view: 'new',        // new | run | setup | settings
    runId: null,        // the run the run view shows
    setup: 'idle',      // idle | running | failed
    logs: new Map(),    // run id -> this session's log text
};

// h('div', { class: 'x', text: 'y', onclick: f, title: 't' }, child, ...)
function h(tag, props = {}, ...children) {
    const el = document.createElement(tag);
    for (const [k, v] of Object.entries(props)) {
        if (v === undefined || v === null || v === false) continue;
        if (k === 'class') el.className = v;
        else if (k === 'text') el.textContent = v;
        else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
        else el.setAttribute(k, v === true ? '' : v);
    }
    for (const c of children.flat()) if (c !== null && c !== undefined && c !== false) el.append(c);
    return el;
}

// inline stroke icons (renderer/style.css .ico); static markup only
const ICONS = {
    check: '<path d="M5 12.5l4.5 4.5L19 7.5"></path>',
    'check-circle': '<circle cx="12" cy="12" r="9"></circle><path d="M8 12.5l3 3 5-6"></path>',
    alert: '<circle cx="12" cy="12" r="9"></circle><path d="M12 7.5v5.5"></path><path d="M12 16.5v.01"></path>',
    bang: '<path d="M12 6v8"></path><path d="M12 18v.01"></path>',
    minus: '<circle cx="12" cy="12" r="9"></circle><path d="M8 12h8"></path>',
    spinner: '<circle cx="12" cy="12" r="9" class="spin-track"></circle><path d="M21 12a9 9 0 0 0-9-9"></path>',
    open: '<path d="M14 4h6v6"></path><path d="M10 14L20 4"></path><path d="M19 13v5a2 2 0 0 1-2 2H6a2 2 0 0 1-2-2V7a2 2 0 0 1 2-2h5"></path>',
    folder: '<path d="M3 7a2 2 0 0 1 2-2h4l2 2h8a2 2 0 0 1 2 2v8a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2z"></path>',
    again: '<path d="M4 12a8 8 0 0 1 14-5.3L20 9"></path><path d="M20 4v5h-5"></path><path d="M20 12a8 8 0 0 1-14 5.3L4 15"></path><path d="M4 20v-5h5"></path>',
    stop: '<rect x="6" y="6" width="12" height="12" rx="2"></rect>',
};

function icon(name, cls = '') {
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('viewBox', '0 0 24 24');
    svg.setAttribute('class', `ico ${cls}`.trim());
    svg.setAttribute('aria-hidden', 'true');
    svg.innerHTML = ICONS[name];
    return svg;
}

// a status line under a field: `cls` is ok|warn|bad|'' (neutral), or null to hide it
function setStatus(sel, cls, text) {
    const el = $(sel);
    el.hidden = cls === null;
    if (cls === null) return;
    el.className = `field-status ${cls}`;
    (el.querySelector('.msg') || el).textContent = text;
}

function fmtBytes(n) {
    const u = ['B', 'KB', 'MB', 'GB', 'TB'];
    let i = 0;
    while (n >= 1024 && i < u.length - 1) { n /= 1024; i++; }
    return `${n.toFixed(i ? 1 : 0)} ${u[i]}`;
}

// "Today, 14:47" / "Yesterday, 09:12" / "23 Sep" / "23 Sep 2025"
function fmtWhen(iso) {
    if (!iso) return '';
    const d = new Date(iso);
    const now = new Date();
    const time = d.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' });
    const day = (x) => new Date(x.getFullYear(), x.getMonth(), x.getDate()).getTime();
    const days = Math.round((day(now) - day(d)) / 86400000);
    if (days === 0) return `Today, ${time}`;
    if (days === 1) return `Yesterday, ${time}`;
    return d.toLocaleDateString([], { day: 'numeric', month: 'short', ...(d.getFullYear() !== now.getFullYear() ? { year: 'numeric' } : {}) });
}

// 66000 -> "1 min 6 s"
function fmtDuration(ms) {
    const s = Math.round(ms / 1000);
    if (s < 60) return `${s} s`;
    const m = Math.floor(s / 60);
    if (m < 60) return `${m} min${s % 60 ? ` ${s % 60} s` : ''}`;
    return `${Math.floor(m / 60)} h ${m % 60} min`;
}

const baseName = (p) => (p ? String(p).split(/[\\/]/).pop() : '');

// "/Users/me/Documents/x" -> "~/Documents/x"
function tildify(p) {
    const home = state.env?.home;
    return home && p && (p === home || p.startsWith(`${home}/`)) ? `~${p.slice(home.length)}` : p || '';
}

// "Saccharomyces cerevisiae" -> "S. cerevisiae"
function shortSpecies(name) {
    const m = String(name || '').match(/^([A-Z])[a-z]+ ([a-z][\w-]*)$/);
    return m ? `${m[1]}. ${m[2]}` : name || '';
}
