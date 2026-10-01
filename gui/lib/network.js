// Is there a connection, to the places the app needs one for? navigator.onLine
// and Electron's net.isOnline() only know whether a network interface is
// up, which says nothing about a captive portal, a firewall or a dead
// Wi-Fi, so this asks the two services themselves. Any HTTP answer, even
// an error status, means reachable.

const PROBES = {
    ncbi: 'https://api.ncbi.nlm.nih.gov/datasets/v2/',        // reference discovery, taxid checks
    conda: 'https://conda.anaconda.org/conda-forge/',         // runtime setup, tool environments
};

async function reachable(url) {
    try {
        await fetch(url, { method: 'HEAD', signal: AbortSignal.timeout(5000) });
        return true;
    } catch {
        return false;
    }
}

async function checkNetwork() {
    const [ncbi, conda] = await Promise.all([reachable(PROBES.ncbi), reachable(PROBES.conda)]);
    return { ncbi, conda, online: ncbi || conda };
}

module.exports = { checkNetwork };
