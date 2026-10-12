// Can the app reach NCBI and conda-forge? navigator.onLine and
// net.isOnline() only know whether an interface is up (not captive portals,
// firewalls, dead Wi-Fi), so this asks the services. Any HTTP answer, even
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
