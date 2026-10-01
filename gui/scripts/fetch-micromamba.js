#!/usr/bin/env node
// Puts conda-forge's pinned micromamba (lib/runtime.js) for each given
// platform-arch into vendor/micromamba/<platform-arch>/, where
// electron-builder picks it up as an extra resource (package.json).
//
//   node scripts/fetch-micromamba.js darwin-arm64 [darwin-x64 ...]

const fs = require('node:fs');
const path = require('node:path');
const { downloadMicromamba, MICROMAMBA_VERSION } = require('../lib/runtime');

(async () => {
    const targets = process.argv.slice(2);
    if (!targets.length) targets.push(`${process.platform}-${process.arch}`);
    for (const t of targets) {
        const dest = path.join(__dirname, '..', 'vendor', 'micromamba', t, 'micromamba');
        const stamp = `${dest}.version`;
        if (fs.existsSync(dest) && fs.existsSync(stamp) && fs.readFileSync(stamp, 'utf8') === MICROMAMBA_VERSION) {
            console.log(`${t}: micromamba ${MICROMAMBA_VERSION} already present`);
            continue;
        }
        await downloadMicromamba(t, dest, (s) => process.stdout.write(`${t}: ${s}`));
        fs.writeFileSync(stamp, MICROMAMBA_VERSION);
    }
})().catch((e) => {
    console.error(e.message);
    process.exit(1);
});
