// lib/updates.js: version order, and which release the notice offers, from
// a local mock of GitHub's releases API

const test = require('node:test');
const assert = require('node:assert/strict');
const http = require('node:http');
const { compare, checkForUpdate, isRepoUrl } = require('../lib/updates');

test('versions sort like semver, a release after its own pre-releases', () => {
    assert.ok(compare('0.1.0', '0.0.1') > 0);
    assert.ok(compare('0.10.0', '0.9.9') > 0);
    assert.ok(compare('0.1.0-alpha', '0.1.0') < 0);
    assert.ok(compare('0.1.0-alpha', '0.1.0-beta') < 0);
    assert.ok(compare('0.1.0-alpha.2', '0.1.0-alpha.10') < 0);
    assert.equal(compare('v1.2.3', '1.2.3'), 0);
});

test('only this repo\'s GitHub links may be opened', () => {
    assert.ok(isRepoUrl('https://github.com/gbdias/quick_synteny/releases/tag/v0.2.0'));
    assert.ok(!isRepoUrl('https://github.com/someone/else/releases'));
    assert.ok(!isRepoUrl('https://evil.example.com/quick_synteny.dmg'));
});

const gh = (p) => `https://github.com/gbdias/quick_synteny/${p}`;
const release = (tag, files) => ({
    tag_name: tag, name: tag, draft: false, prerelease: false, html_url: gh(`releases/tag/${tag}`),
    assets: files.map((n) => ({ name: n, browser_download_url: gh(`releases/download/${tag}/${n}`) })),
});
const RELEASES = [
    release('v0.3.0', ['quick_synteny-0.3.0-linux-x86_64.AppImage']), // newer, but nothing for macOS
    release('v0.2.0', ['quick_synteny-0.2.0-mac-arm64.dmg', 'quick_synteny-0.2.0-linux-x64.AppImage']),
    { ...release('v0.2.5', ['quick_synteny-0.2.5-mac-arm64.dmg']), prerelease: true },
    release('v0.1.5', []),                                          // pipeline-only
    release('v0.1.0-alpha', []),
];

async function withMockFeed(fn) {
    const server = http.createServer((_req, res) => {
        res.setHeader('content-type', 'application/json');
        res.end(JSON.stringify(RELEASES));
    });
    await new Promise((r) => server.listen(0, '127.0.0.1', r));
    try {
        return await fn(`http://127.0.0.1:${server.address().port}`);
    } finally {
        server.close();
    }
}

test('the newest release with a file for this platform and arch, skipping pre-releases and pipeline-only ones', () =>
    withMockFeed(async (url) => {
        const mac = await checkForUpdate({ current: '0.0.1', platform: 'darwin', arch: 'arm64', url });
        assert.equal(mac.version, '0.2.0');
        assert.equal(mac.file, 'quick_synteny-0.2.0-mac-arm64.dmg');
        assert.equal(mac.notesUrl, gh('releases/tag/v0.2.0'));
        // electron-builder names an x64 AppImage x86_64
        const linux = await checkForUpdate({ current: '0.0.1', platform: 'linux', arch: 'x64', url });
        assert.equal(linux.file, 'quick_synteny-0.3.0-linux-x86_64.AppImage');
    }));

test('nothing to offer: up to date, no file for this machine, unsupported platform', () =>
    withMockFeed(async (url) => {
        assert.equal(await checkForUpdate({ current: '0.2.0', platform: 'darwin', arch: 'arm64', url }), null);
        assert.equal(await checkForUpdate({ current: '0.0.1', platform: 'darwin', arch: 'x64', url }), null);
        assert.equal(await checkForUpdate({ current: '0.0.1', platform: 'win32', arch: 'x64', url }), null);
    }));

test('an unreachable feed is "no news", not an error', async () => {
    assert.equal(await checkForUpdate({ current: '0.0.1', url: 'http://127.0.0.1:1' }), null);
});
