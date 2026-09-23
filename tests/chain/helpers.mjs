// Shared builders for the chainer tests: small hand-made genomes and a
// larger synthetic pair with planted rearrangements.
import {createRequire} from 'node:module';

const require = createRequire(import.meta.url);
export const SYNCHAIN = require('../../bin/chain.js');

// hits: [{chrom, pos, prot, strand = '+', positive = 0.95, rank = 1, locus?}]
// -> spec-ordered TSV-style rows. Loci default to one per distinct
// (chrom, pos), numbered in genomic order.
export function rowsFrom(hits, chromNames) {
    const ci = new Map(chromNames.map((c, i) => [c, i]));
    const sorted = hits.slice().sort((a, b) => (ci.get(a.chrom) - ci.get(b.chrom)) || (a.pos - b.pos)
        || (a.prot < b.prot ? -1 : a.prot > b.prot ? 1 : 0) || ((a.rank || 1) - (b.rank || 1)));
    const locusOf = new Map();
    const lociInOrder = [];
    for (const h of sorted) {
        const key = h.locus !== undefined ? `L${h.locus}` : `${h.chrom}:${h.pos}`;
        if (!locusOf.has(key)) { locusOf.set(key, lociInOrder.length); lociInOrder.push(key); }
    }
    return sorted.map((h) => [h.chrom, String(h.pos - 500), String(h.pos + 500), h.strand || '+',
        (h.positive ?? 0.95).toFixed(4), (h.positive ?? 0.95).toFixed(4), '1000', String(h.rank || 1), h.prot,
        String(locusOf.get(h.locus !== undefined ? `L${h.locus}` : `${h.chrom}:${h.pos}`))]);
}

export function tablesFrom(hitsA, hitsB, chromsA = ['A1'], chromsB = ['B1']) {
    const genomes = [{name: 'target', rows: rowsFrom(hitsA, chromsA), chromNames: chromsA}];
    if (hitsB) { genomes.push({name: 'reference', rows: rowsFrom(hitsB, chromsB), chromNames: chromsB}); }
    return SYNCHAIN.makeTables(genomes).tables;
}

export const P = (over = {}) => ({minPositive: 0.5, maxGap: 25, maxLookback: 50, gapPenalty: 0, minBlock: 5,
                                  maxHitRank: 255, ...over});

function rng(seed) {
    let s = seed >>> 0;
    return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296);
}

// A: nGenes genes over nChrom chromosomes; B: A cut into 50-500-gene
// segments, 30 % inverted, shuffled, redistributed. Every gene gets 1-3
// isoform proteins (same locus) and up to `secMax` random secondary hits per
// genome at other genes' loci with lower identity. Returns tables plus the
// ground truth needed to score recall and block purity.
export function syntheticPair({seed = 42, nChrom = 20, nGenes = 25000, secMax = 4, spacing = 1_000_000} = {}) {
    const R = rng(seed);
    const perChrom = Math.ceil(nGenes / nChrom);
    const segOfGene = new Int32Array(nGenes);
    const segs = [];
    for (let i = 0; i < nGenes;) {
        const L = 50 + Math.floor(R() * 450);
        const s = [];
        for (let g = i; g < Math.min(i + L, nGenes); g++) { s.push(g); segOfGene[g] = segs.length; }
        segs.push({genes: s, inv: R() < 0.3});
        i += L;
    }
    const order = segs.map((_, i) => i);
    for (let i = order.length - 1; i > 0; i--) { const j = Math.floor(R() * (i + 1)); [order[i], order[j]] = [order[j], order[i]]; }
    const layA = [], layB = new Array(nGenes);
    for (let g = 0; g < nGenes; g++) { layA.push({chrom: Math.floor(g / perChrom), idx: g % perChrom}); }
    let k = 0;
    for (const si of order) {
        const gs = segs[si].inv ? segs[si].genes.slice().reverse() : segs[si].genes;
        for (const g of gs) { layB[g] = {chrom: Math.floor(k / perChrom), idx: k % perChrom}; k++; }
    }
    const chromsA = Array.from({length: nChrom}, (_, i) => `a${i + 1}`);
    const chromsB = Array.from({length: nChrom}, (_, i) => `b${i + 1}`);
    const hitsA = [], hitsB = [];
    const posOf = (l) => l.idx * spacing + 5000;
    let pid = 0;
    const primaryProt = [];
    for (let g = 0; g < nGenes; g++) {
        const nIso = 1 + Math.floor(R() * 3);
        for (let iso = 0; iso < nIso; iso++, pid++) {
            const prot = `p${String(pid).padStart(6, '0')}`;
            if (iso === 0) { primaryProt.push(prot); }
            hitsA.push({chrom: chromsA[layA[g].chrom], pos: posOf(layA[g]), prot, positive: 0.85 + R() * 0.14, locus: `A${g}`});
            hitsB.push({chrom: chromsB[layB[g].chrom], pos: posOf(layB[g]), prot, positive: 0.85 + R() * 0.14, locus: `B${g}`});
            for (const [hits, lay, chroms, tag] of [[hitsA, layA, chromsA, 'A'], [hitsB, layB, chromsB, 'B']]) {
                const nSec = Math.floor(R() * (secMax + 1));
                for (let s = 0; s < nSec; s++) {
                    const og = Math.floor(R() * nGenes);
                    hits.push({chrom: chroms[lay[og].chrom], pos: posOf(lay[og]), prot, positive: 0.3 + R() * 0.5,
                               rank: 2 + s, locus: `${tag}${og}`});
                }
            }
        }
    }
    // loci must be numbered in genomic order: rowsFrom numbers them by first
    // appearance in sorted order, and each locus has a single position here
    const [A, B] = tablesFrom(hitsA, hitsB, chromsA, chromsB);
    return {A, B, nGenes, segOfGene, primaryProt, layA, layB, chromsA, chromsB};
}

// union / summed length of block query spans, for coverage & depth
export function coverage(blocks, nChrom) {
    const by = Array.from({length: nChrom}, () => []);
    let sum = 0;
    for (const b of blocks) { by[b.qChrom].push([b.qStart, b.qEnd]); sum += b.qEnd - b.qStart; }
    let union = 0;
    for (const iv of by) {
        iv.sort((x, y) => x[0] - y[0]);
        let cs = -1, ce = -1;
        for (const [s, e] of iv) {
            if (s > ce) { if (ce > cs) { union += ce - cs; } cs = s; ce = e; } else if (e > ce) { ce = e; }
        }
        if (ce > cs) { union += ce - cs; }
    }
    return {union, sum, depth: union ? sum / union : 0};
}
