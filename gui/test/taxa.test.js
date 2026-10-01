// lib/taxa.js: the species search, against canned NCBI answers (from real
// ones: taxon_suggest leaves Babesia ovis out of its own name's results)

const test = require('node:test');
const assert = require('node:assert/strict');
const { searchTaxa } = require('../lib/taxa');

const reply = (body) => Promise.resolve({ ok: true, json: () => Promise.resolve(body) });
const suggest = (...names) => ({ sci_name_and_ids: names.map(([sci_name, tax_id, rank]) => ({ sci_name, tax_id, rank })) });

function fakeNcbi({ all, genomes, exact, down = false }) {
    return (url) => {
        if (down) return Promise.reject(new Error('fetch failed'));
        if (url.includes('TAXON_RESOURCE_FILTER_GENOME')) return reply(genomes);
        if (url.includes('/taxon_suggest/')) return reply(all);
        return reply(exact);
    };
}

test('a species taxon_suggest leaves out is found by its exact name, and listed first', async () => {
    const hits = await searchTaxa('Babesia ovis', fakeNcbi({
        all: suggest(['Babesia ovata', '189622', 'SPECIES'], ['Babesia bigemina', '5866', 'SPECIES']),
        genomes: suggest(['Babesia bigemina', '5866', 'SPECIES'], ['Babesia ovis', '5869', 'SPECIES']),
        exact: { taxonomy_nodes: [{ taxonomy: { tax_id: 5869, organism_name: 'Babesia ovis', rank: 'SPECIES' } }] },
    }));
    assert.deepEqual(hits.map((h) => h.sci_name), ['Babesia ovis', 'Babesia ovata', 'Babesia bigemina']);
    assert.equal(hits[0].tax_id, '5869');
});

test('names that start with the query come first, the plain species before its hybrids and strains; nothing repeats', async () => {
    const hits = await searchTaxa('saccharomyces cer', fakeNcbi({
        all: suggest(['Saccharomyces cerevisiae x Saccharomyces kudriavzevii', '332112', 'SPECIES'], ['Saccharomyces cerevisiae virus L-A', '11008'],
            ['Brewer yeast thing', '1'], ['Saccharomyces cerevisiae S288C', '559292', 'STRAIN']),
        genomes: suggest(['Saccharomyces cerevisiae', '4932', 'SPECIES']),
        exact: { taxonomy_nodes: [{ errors: [{ reason: 'not a recognized NCBI Taxonomy name' }] }] },
    }));
    assert.deepEqual(hits.map((h) => h.tax_id), ['4932', '332112', '559292', '11008', '1']);
});

test('NCBI saying "too many requests" is retried once; still busy, the error says so', async () => {
    let calls = 0;
    const flaky = (url) => {
        calls++;
        if (calls <= 3) return Promise.resolve({ ok: false, status: 429 });
        return fakeNcbi({ all: suggest(['Babesia ovis', '5869', 'SPECIES']), genomes: suggest(), exact: { taxonomy_nodes: [{}] } })(url);
    };
    const hits = await searchTaxa('Babesia ovis', flaky, { retryAfter: 1 });
    assert.equal(hits[0].tax_id, '5869');
    const busy = () => Promise.resolve({ ok: false, status: 429 });
    await assert.rejects(searchTaxa('Babesia ovis', busy, { retryAfter: 1 }), (e) => e.busy === true);
});

test('NCBI unreachable: an error, so the form can say so', async () => {
    await assert.rejects(searchTaxa('Babesia ovis', fakeNcbi({ down: true })));
});
