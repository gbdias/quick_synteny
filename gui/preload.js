const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('qs', {
    checkEnv: () => ipcRenderer.invoke('env:check'),
    setupRuntime: () => ipcRenderer.invoke('runtime:setup'),
    checkDocker: () => ipcRenderer.invoke('docker:check'),
    checkNetwork: () => ipcRenderer.invoke('net:check'),
    checkForUpdate: () => ipcRenderer.invoke('update:check'),
    dismissUpdate: (version) => ipcRenderer.invoke('update:dismiss', version),
    openUpdate: (url) => ipcRenderer.invoke('update:open', url),
    getSettings: () => ipcRenderer.invoke('settings:get'),
    setSettings: (s) => ipcRenderer.invoke('settings:set', s),
    pickFile: (opts) => ipcRenderer.invoke('dialog:file', opts),
    pickDir: (opts) => ipcRenderer.invoke('dialog:dir', opts),
    suggestTaxon: (q) => ipcRenderer.invoke('taxon:suggest', q),
    lookupTaxon: (taxid) => ipcRenderer.invoke('taxon:lookup', taxid),
    suggestRunName: (opts) => ipcRenderer.invoke('outdir:suggest', opts),
    startRun: (form) => ipcRenderer.invoke('run:start', form),
    cancelRun: () => ipcRenderer.invoke('run:cancel'),
    openResult: (file) => ipcRenderer.invoke('result:open', file),
    pickResult: () => ipcRenderer.invoke('result:pick'),
    reveal: (p) => ipcRenderer.invoke('path:reveal', p),
    on: (channel, fn) => {
        if (!['run:log', 'run:task', 'run:event', 'run:note', 'run:error', 'run:done', 'setup:log'].includes(channel)) return;
        ipcRenderer.on(channel, (_e, payload) => fn(payload));
    },
});
