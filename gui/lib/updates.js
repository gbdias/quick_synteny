// "A new version is available": the app isn't signed, so it can't update
// itself (electron-updater needs a signed app on macOS). It only looks for
// a newer release on GitHub and links to it.
//
// Only a release carrying an app file for this platform and architecture
// counts. The repo also has pipeline-only releases (v0.1.0-alpha), which
// have nothing to download. App files are named by package.json's
// build.artifactName: quick_synteny-<version>-<mac|linux>-<arch>.<ext>

const REPO = 'gbdias/quick_synteny';
const RELEASES = `https://api.github.com/repos/${REPO}/releases?per_page=30`;

// the file to offer, best first, per platform
const EXTENSIONS = { darwin: ['dmg'], linux: ['AppImage', 'deb', 'rpm'] };
const OS_NAMES = { darwin: 'mac', linux: 'linux' };

// semver order, enough for tags like 0.2.0 and 0.1.0-alpha: numbers first,
// then a release sorts after its own pre-releases
function parse(v) {
    const m = String(v).replace(/^v/, '').match(/^(\d+)\.(\d+)\.(\d+)(?:-([0-9A-Za-z.-]+))?/);
    return m ? { nums: m.slice(1, 4).map(Number), pre: m[4] ? m[4].split('.') : [] } : null;
}

function compare(a, b) {
    const x = parse(a);
    const y = parse(b);
    if (!x || !y) return 0;
    for (let i = 0; i < 3; i++) if (x.nums[i] !== y.nums[i]) return x.nums[i] - y.nums[i];
    if (!x.pre.length || !y.pre.length) return y.pre.length - x.pre.length;
    for (let i = 0; i < Math.max(x.pre.length, y.pre.length); i++) {
        const [p, q] = [x.pre[i], y.pre[i]];
        if (p === undefined) return -1;
        if (q === undefined) return 1;
        const [pn, qn] = [/^\d+$/.test(p), /^\d+$/.test(q)];
        if (pn && qn && Number(p) !== Number(q)) return Number(p) - Number(q);
        if (pn !== qn) return pn ? -1 : 1;
        if (p !== q) return p < q ? -1 : 1;
    }
    return 0;
}

// links the page may open: this repo's release pages and files only
function isRepoUrl(url) {
    return typeof url === 'string' && url.startsWith(`https://github.com/${REPO}/`);
}

// the newest release newer than `current` with a file for this machine,
// or null. Network trouble means "no news", never an error: checking is
// a courtesy, not something to fail over.
async function checkForUpdate({ current, platform = process.platform, arch = process.arch, url = RELEASES }) {
    const exts = EXTENSIONS[platform];
    if (!exts) return null;
    let releases;
    try {
        const r = await fetch(url, {
            headers: { Accept: 'application/vnd.github+json', 'User-Agent': 'quick_synteny-app' },
            signal: AbortSignal.timeout(8000),
        });
        if (!r.ok) return null;
        releases = await r.json();
    } catch {
        return null;
    }
    let best = null;
    for (const rel of Array.isArray(releases) ? releases : []) {
        if (rel.draft || rel.prerelease || compare(rel.tag_name, current) <= 0) continue;
        if (best && compare(rel.tag_name, best.version) <= 0) continue;
        const asset = exts.map((ext) => (rel.assets || []).find((a) => a.name.endsWith(`-${OS_NAMES[platform]}-${arch}.${ext}`)))
            .find(Boolean);
        if (!asset || !isRepoUrl(asset.browser_download_url) || !isRepoUrl(rel.html_url)) continue;
        best = {
            version: rel.tag_name.replace(/^v/, ''),
            name: rel.name || rel.tag_name,
            notesUrl: rel.html_url,
            downloadUrl: asset.browser_download_url,
            file: asset.name,
            platform,
        };
    }
    return best;
}

module.exports = { checkForUpdate, compare, isRepoUrl };
