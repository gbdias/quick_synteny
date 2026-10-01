const { contextBridge, ipcRenderer, webUtils } = require('electron');

contextBridge.exposeInMainWorld('qs', {
    platform: process.platform,
    checkEnv: () => ipcRenderer.invoke('env:check'),
    setupRuntime: () => ipcRenderer.invoke('runtime:setup'),
    checkDocker: () => ipcRenderer.invoke('docker:check'),
    checkNetwork: () => ipcRenderer.invoke('net:check'),
    getSchema: () => ipcRenderer.invoke('schema:get'),
    checkForUpdate: () => ipcRenderer.invoke('update:check'),
    dismissUpdate: (version) => ipcRenderer.invoke('update:dismiss', version),
    openUpdate: (url) => ipcRenderer.invoke('update:open', url),
    getSettings: () => ipcRenderer.invoke('settings:get'),
    setSettings: (s) => ipcRenderer.invoke('settings:set', s),
    pickFile: (opts) => ipcRenderer.invoke('dialog:file', opts),
    pickDir: (opts) => ipcRenderer.invoke('dialog:dir', opts),
    // a dropped file's path: File.path is gone from Electron's renderer
    pathForFile: (file) => webUtils.getPathForFile(file),
    inspectFile: (p) => ipcRenderer.invoke('file:inspect', p),
    defaultOutdir: (opts) => ipcRenderer.invoke('outdir:default', opts),
    suggestTaxon: (q) => ipcRenderer.invoke('taxon:suggest', q),
    lookupTaxon: (taxid) => ipcRenderer.invoke('taxon:lookup', taxid),
    lookupAssembly: (accession) => ipcRenderer.invoke('assembly:lookup', accession),
    suggestRunName: (opts) => ipcRenderer.invoke('outdir:suggest', opts),
    startRun: (form, meta) => ipcRenderer.invoke('run:start', form, meta),
    cancelRun: () => ipcRenderer.invoke('run:cancel'),
    listRuns: () => ipcRenderer.invoke('runs:list'),
    removeRun: (id) => ipcRenderer.invoke('runs:remove', id),
    getThumbnail: (id) => ipcRenderer.invoke('thumb:get', id),
    openResult: (file) => ipcRenderer.invoke('result:open', file),
    pickResult: () => ipcRenderer.invoke('result:pick'),
    reveal: (p) => ipcRenderer.invoke('path:reveal', p),
    on: (channel, fn) => {
        if (!['run:log', 'run:update', 'runs:changed', 'setup:log'].includes(channel)) return;
        ipcRenderer.on(channel, (_e, payload) => fn(payload));
    },
});
