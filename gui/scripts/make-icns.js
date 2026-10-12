// Builds build/icon.icns from build/icon-macos.svg: Chromium renders the
// vector art at each size, and iconutil packs them. macOS only; commit the
// result.
//
//   npm run icons
//
// electron-builder's PNG conversion writes the 16 and 32 px sizes (what
// Spotlight and Finder's list view draw) as noise, and scaling the PNG down
// with sips leaves a pale fringe on semi-transparent edges.

const { app, BrowserWindow } = require('electron');
const { execFileSync } = require('node:child_process');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');

const BUILD = path.join(__dirname, '..', 'build');
const SIZES = [16, 32, 128, 256, 512];

app.whenReady().then(async () => {
    const svg = fs.readFileSync(path.join(BUILD, 'icon-macos.svg'), 'utf8');
    const win = new BrowserWindow({ show: false, webPreferences: { offscreen: true } });
    await win.loadURL('data:text/html,<!doctype html><title>icons</title>');
    const iconset = path.join(fs.mkdtempSync(path.join(os.tmpdir(), 'icns-')), 'icon.iconset');
    fs.mkdirSync(iconset);
    for (const size of SIZES) {
        for (const scale of [1, 2]) {
            const px = size * scale;
            // rasterized at px, not shrunk from a larger render
            const dataUrl = await win.webContents.executeJavaScript(`new Promise((resolve, reject) => {
                const img = new Image();
                img.onload = () => {
                    const c = document.createElement('canvas');
                    c.width = c.height = ${px};
                    const g = c.getContext('2d');
                    g.imageSmoothingQuality = 'high';
                    g.drawImage(img, 0, 0, ${px}, ${px});
                    resolve(c.toDataURL('image/png'));
                };
                img.onerror = () => reject(new Error('the SVG did not load'));
                img.src = 'data:image/svg+xml;base64,' + ${JSON.stringify(Buffer.from(svg).toString('base64'))};
            })`);
            const name = `icon_${size}x${size}${scale === 2 ? '@2x' : ''}.png`;
            fs.writeFileSync(path.join(iconset, name), Buffer.from(dataUrl.split(',')[1], 'base64'));
        }
    }
    execFileSync('iconutil', ['-c', 'icns', '-o', path.join(BUILD, 'icon.icns'), iconset]);
    fs.rmSync(path.dirname(iconset), { recursive: true, force: true });
    console.log(`build/icon.icns: ${(fs.statSync(path.join(BUILD, 'icon.icns')).size / 1024).toFixed(0)} KB`);
    app.quit();
}).catch((e) => {
    console.error(e.message);
    app.exit(1);
});
