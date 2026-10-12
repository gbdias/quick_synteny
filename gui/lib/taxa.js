// The New run form's species search.
//
// taxon_suggest alone misses some species even by exact name ("Babesia
// ovis" gets other Babesia, not itself). Filtered to taxa with genomes it
// lists it, but the user's genome may be its species' first; the taxon
// endpoint takes only exact scientific names. So all three are asked,
// merged without repeats, closest first: the exact name, then prefix
// matches (species first, then shorter names), then the rest.
//
// Without an API key NCBI answers bursts with 429: retried once, a second later.

const API = 'https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy';

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

async function getJson(fetchImpl, url, retryAfter = 1000) {
    let r = await fetchImpl(url, { signal: AbortSignal.timeout(10000) });
    if (r.status === 429) {
        await sleep(retryAfter);
        r = await fetchImpl(url, { signal: AbortSignal.timeout(10000) });
    }
    if (!r.ok) {
        const e = new Error(r.status === 429 ? 'NCBI is busy' : `NCBI returned ${r.status}`);
        e.busy = r.status === 429;
        throw e;
    }
    return r.json();
}

// -> [{ sci_name, tax_id, rank, common_name }], best first. Throws only when
// every request failed, i.e. NCBI can't be reached.
async function searchTaxa(query, fetchImpl = fetch, { retryAfter = 1000 } = {}) {
    const q = query.trim();
    const enc = encodeURIComponent(q);
    const [all, genomes, exact] = await Promise.allSettled([
        getJson(fetchImpl, `${API}/taxon_suggest/${enc}`, retryAfter),
        getJson(fetchImpl, `${API}/taxon_suggest/${enc}?taxon_resource_filter=TAXON_RESOURCE_FILTER_GENOME`, retryAfter),
        getJson(fetchImpl, `${API}/taxon/${enc}`, retryAfter),
    ]);
    if ([all, genomes, exact].every((p) => p.status === 'rejected')) throw all.reason;
    const hits = [];
    const t = exact.status === 'fulfilled' ? exact.value.taxonomy_nodes?.[0]?.taxonomy : null;
    if (t) hits.push({ sci_name: t.organism_name, tax_id: String(t.tax_id), rank: t.rank, common_name: t.common_name });
    for (const p of [all, genomes]) {
        if (p.status === 'fulfilled') hits.push(...(p.value.sci_name_and_ids || []).map((h) => ({ ...h, tax_id: String(h.tax_id) })));
    }
    const seen = new Set();
    const lower = q.toLowerCase();
    const closeness = (h) => {
        const name = (h.sci_name || '').toLowerCase();
        return name === lower ? 0 : name.startsWith(lower) ? 1 : 2;
    };
    const species = (h) => (String(h.rank).toUpperCase() === 'SPECIES' ? 0 : 1);
    return hits.filter((h) => !seen.has(h.tax_id) && seen.add(h.tax_id))
        .map((h, i) => ({ h, i }))
        .sort((a, b) => closeness(a.h) - closeness(b.h)
            || (closeness(a.h) === 1 && (species(a.h) - species(b.h) || a.h.sci_name.length - b.h.sci_name.length))
            || a.i - b.i)
        .map(({ h }) => h);
}

module.exports = { searchTaxa };
