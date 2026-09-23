// Recall and purity against planted ground truth: 25k genes on 20
// chromosomes, genome B = A cut into 50-500-gene segments, 30 % inverted,
// shuffled across chromosomes; isoforms and random low-identity secondary
// hits as noise (see helpers.mjs:syntheticPair).
import test from 'node:test';
import assert from 'node:assert/strict';
import {SYNCHAIN, P, syntheticPair} from './helpers.mjs';

const S = syntheticPair({});

function geneOfLocus(t, primaryProt) {
    const protGene = new Map(primaryProt.map((p, g) => [p, g]));
    const out = new Int32Array(t.nLoci).fill(-1);
    for (let i = 0; i < t.n; i++) {
        if (t.rank[i] !== 1) { continue; }
        const name = `p${String(t.prot[i]).padStart(6, '0')}`;
        // prot indices index the sorted accession list, which is exactly
        // p000000.. in order, since every accession is zero-padded
        const g = protGene.get(name);
        if (g !== undefined) { out[t.locus[i]] = g; }
    }
    return out;
}

test('recovers planted synteny with inversions, and blocks stay within one planted segment', () => {
    const gA = geneOfLocus(S.A, S.primaryProt), gB = geneOfLocus(S.B, S.primaryProt);
    const t = performance.now();
    const r = SYNCHAIN.chain(S.A, S.B, P({minPositive: 0.6, minBlock: 5}));
    const ms = performance.now() - t;

    // A noise anchor can only join a chain where strict monotonicity leaves
    // room -- beyond a true block's last gene -- so false members may extend
    // a block's ends by a few positions but must never sit inside it.
    const END_ZONE = 2;
    const found = new Uint8Array(S.nGenes);
    let impure = 0, falseMembers = 0, falseInside = 0, members = 0;
    for (const b of r.blocks) {
        const segs = new Set();
        b.members.forEach((m, k) => {
            members++;
            const ga = gA[S.A.locus[r.anchors.ha[m]]], gb = gB[S.B.locus[r.anchors.hb[m]]];
            if (ga < 0 || ga !== gb) {
                falseMembers++;
                if (Math.min(k, b.members.length - 1 - k) > END_ZONE) { falseInside++; }
                return;
            }
            found[ga] = 1;
            segs.add(S.segOfGene[ga]);
        });
        if (segs.size > 1) { impure++; }
    }
    const recall = found.reduce((s, v) => s + v, 0) / S.nGenes;
    const inv = r.blocks.filter((b) => b.orientation === '-').length;
    console.log(`  synthetic: ${S.A.n}+${S.B.n} hits, ${r.stats.nAnchors} anchors -> ${r.blocks.length} blocks `
        + `(${inv} '-'), recall ${(recall * 100).toFixed(2)} %, impure blocks ${impure}, `
        + `false members ${falseMembers}/${members} (${falseInside} inside a block), chain ${ms.toFixed(0)} ms`);
    assert.ok(recall >= 0.98, `recall ${recall}`);
    assert.ok(impure <= Math.ceil(0.01 * r.blocks.length), `impure ${impure}`);
    assert.equal(falseInside, 0);
    assert.ok(falseMembers <= 0.01 * members, `false members ${falseMembers}`);
    assert.ok(inv > 0);
    assert.ok(ms < 1000, `chain took ${ms} ms`);
});
