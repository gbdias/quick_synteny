// Synteny chainer over per-genome hit tables -- the only implementation,
// shared by the pipeline (node, via bin/chain_blocks.mjs) and the
// interactive page (inlined into the HTML and run in a Web Worker). The
// contract -- table layout, parameters, and the exact chaining semantics --
// is docs/specs/hit_table.md; everything here follows its section numbers.
//
// Dependency-free and DOM-free on purpose: the same source text must run
// as a CommonJS module, a classic <script>, and a Worker body.
(function (root) {
    'use strict';

    const VERSION = 1;
    const MIN_CHAIN_LENGTH = 2;   // a single anchor is never a block, whatever minBlock says

    const DEFAULTS = {
        maxHitRank: 255,
        maxGap: 25,
        maxLookback: 50,
        gapPenalty: 0,
    };

    // ---------------------------------------------------------------- tables

    function parseHitsTsv(text) {
        const lines = text.split('\n');
        const header = lines[0].replace(/\r$/, '').split('\t');
        const expect = ['chrom', 'start', 'end', 'strand', 'positive', 'identity', 'score', 'rank', 'protein', 'locus'];
        if (header.join('\t') !== expect.join('\t')) {
            throw new Error(`hits TSV: unexpected header "${header.join('\\t')}"`);
        }
        const rows = [];
        for (let i = 1; i < lines.length; i++) {
            const line = lines[i];
            if (!line) { continue; }
            rows.push(line.replace(/\r$/, '').split('\t'));
        }
        return rows;
    }

    // Builds HitTables (spec 2) from parsed TSV rows for one comparison --
    // two genomes, or one for a self-comparison. The protein list is the
    // sorted union of every genome's accessions, so prot indices are shared.
    // chromNames defaults to first appearance order, which is chrom.sizes
    // order restricted to sequences with hits (rows are sorted that way).
    function makeTables(genomes) {
        const accessions = new Set();
        for (const g of genomes) { for (const r of g.rows) { accessions.add(r[8]); } }
        const proteins = Array.from(accessions).sort();
        const protIndex = new Map(proteins.map((p, i) => [p, i]));
        const tables = genomes.map((g) => {
            const chromNames = g.chromNames ? g.chromNames.slice() : [];
            const chromIndex = new Map(chromNames.map((c, i) => [c, i]));
            const n = g.rows.length;
            const t = {
                name: g.name, chromNames, chromSizes: g.chromSizes ? g.chromSizes.slice() : null, n,
                chrom: new Uint16Array(n), start: new Float64Array(n), end: new Float64Array(n),
                strand: new Uint8Array(n), positive: new Uint16Array(n), rank: new Uint8Array(n),
                prot: new Uint32Array(n), locus: new Uint32Array(n), nLoci: 0,
            };
            for (let i = 0; i < n; i++) {
                const r = g.rows[i];
                let ci = chromIndex.get(r[0]);
                if (ci === undefined) {
                    if (g.chromNames) { throw new Error(`hits: chromosome "${r[0]}" not in chromNames`); }
                    ci = chromNames.length; chromNames.push(r[0]); chromIndex.set(r[0], ci);
                }
                t.chrom[i] = ci;
                t.start[i] = +r[1];
                t.end[i] = +r[2];
                t.strand[i] = r[3] === '+' ? 1 : 0;
                t.positive[i] = Math.round(+r[4] * 10000);
                t.rank[i] = Math.min(+r[7], 255);
                t.prot[i] = protIndex.get(r[8]);
                t.locus[i] = +r[9];
                if (t.locus[i] + 1 > t.nLoci) { t.nLoci = t.locus[i] + 1; }
            }
            return t;
        });
        return {proteins, tables};
    }

    // Derived, parameter-independent indexes: hits grouped by protein
    // (CSR), and each locus's chromosome.
    function prepareTable(t, nProteins) {
        if (t._idx && t._idx.nProteins >= nProteins) { return t; }
        const protOff = new Uint32Array(nProteins + 1);
        for (let i = 0; i < t.n; i++) { protOff[t.prot[i] + 1]++; }
        for (let p = 0; p < nProteins; p++) { protOff[p + 1] += protOff[p]; }
        const fill = protOff.slice(0, nProteins), byProt = new Uint32Array(t.n);
        for (let i = 0; i < t.n; i++) { byProt[fill[t.prot[i]]++] = i; }
        const locusChrom = new Uint16Array(t.nLoci);
        for (let i = 0; i < t.n; i++) { locusChrom[t.locus[i]] = t.chrom[i]; }
        t._idx = {nProteins, protOff, byProt, locusChrom};
        return t;
    }

    // ---------------------------------------------------------------- params

    // mean over proteins of the rank-1 positive (falling back to the best
    // positive for a protein without a rank-1 hit) -- spec 3
    function meanBestPositive(t) {
        const r1 = new Map(), best = new Map();
        for (let i = 0; i < t.n; i++) {
            const p = t.prot[i], v = t.positive[i];
            if (t.rank[i] === 1) { r1.set(p, v); }
            if (!best.has(p) || v > best.get(p)) { best.set(p, v); }
        }
        if (best.size === 0) { return 0; }
        let sum = 0;
        for (const [p, v] of best) { sum += r1.has(p) ? r1.get(p) : v; }
        return sum / best.size / 10000;
    }

    function autoParams(A, B) {
        const w = Math.min(meanBestPositive(A), B && B !== A ? meanBestPositive(B) : Infinity);
        return {
            minPositive: Math.max(0.3, Math.min(0.9, w)),
            minBlock: w >= 0.8 ? 15 : 5,
            weakerMeanPositive: w,
            ...DEFAULTS,
        };
    }

    function resolveParams(p, selfMode) {
        const q = Object.assign({}, DEFAULTS, p || {});
        if (typeof q.minPositive !== 'number' || !(q.minPositive >= 0 && q.minPositive <= 1)) {
            throw new Error('chain: minPositive must be a number in [0, 1]');
        }
        for (const k of ['maxHitRank', 'maxGap', 'maxLookback']) {
            if (!Number.isInteger(q[k]) || q[k] < 1) { throw new Error(`chain: ${k} must be a positive integer`); }
        }
        if (!(q.gapPenalty >= 0)) { throw new Error('chain: gapPenalty must be >= 0'); }
        q.minBlock = Math.max(MIN_CHAIN_LENGTH, Math.floor(q.minBlock || MIN_CHAIN_LENGTH));
        q.selfMode = !!selfMode;
        // integer form of the identity floor, matching the table's x10000 storage
        q._minPositiveInt = Math.ceil(q.minPositive * 10000 - 1e-6);
        return q;
    }

    // ---------------------------------------------------------------- sorting

    // stable counting sort of `idx` by key[idx[k]] in [0, range)
    function countingSort(idx, key, range) {
        const count = new Uint32Array(range + 1);
        for (let k = 0; k < idx.length; k++) { count[key[idx[k]] + 1]++; }
        for (let r = 0; r < range; r++) { count[r + 1] += count[r]; }
        const out = new Uint32Array(idx.length);
        for (let k = 0; k < idx.length; k++) { out[count[key[idx[k]]]++] = idx[k]; }
        return out;
    }

    // ---------------------------------------------------------------- stage 1: anchors

    function locusRanks(t, minPosInt, maxHitRank) {
        const pass = new Uint8Array(t.nLoci);
        for (let i = 0; i < t.n; i++) {
            if (t.positive[i] >= minPosInt && t.rank[i] <= maxHitRank) { pass[t.locus[i]] = 1; }
        }
        const lc = t._idx.locusChrom, r = new Int32Array(t.nLoci);
        let cur = -1, k = 0, maxRank = 0;
        for (let l = 0; l < t.nLoci; l++) {
            if (lc[l] !== cur) { cur = lc[l]; k = 0; }
            if (pass[l]) { r[l] = ++k; if (k > maxRank) { maxRank = k; } } else { r[l] = -1; }
        }
        return {rank: r, maxRank};
    }

    // spec 4.1-4.3: one anchor per (locusA, locusB), ordered by (chrom pair, rankA, rankB)
    function buildAnchors(A, B, q) {
        const thr = q._minPositiveInt, mr = q.maxHitRank, self = q.selfMode;
        const LA = locusRanks(A, thr, mr), LB = self ? LA : locusRanks(B, thr, mr);
        const rkA = LA.rank, rkB = LB.rank;
        const pass = (t, i) => t.positive[i] >= thr && t.rank[i] <= mr;
        const offB = B._idx.protOff, byB = B._idx.byProt;

        // best candidate per (locusA, locusB); tie-breaks per spec 4.2
        const bestOf = new Map();
        const cHa = [], cHb = [], cQ = [], cR = [];
        for (let a = 0; a < A.n; a++) {
            if (!pass(A, a)) { continue; }
            const la = A.locus[a], p = A.prot[a];
            for (let y = offB[p]; y < offB[p + 1]; y++) {
                const b = byB[y];
                if (!pass(B, b)) { continue; }
                const lb = B.locus[b];
                if (self) {
                    if (la >= lb) { continue; }
                    if (A._idx.locusChrom[la] === B._idx.locusChrom[lb] && Math.abs(rkA[la] - rkB[lb]) <= q.maxGap) { continue; }
                }
                const qual = Math.min(A.positive[a], B.positive[b]), rs = A.rank[a] + B.rank[b];
                const key = la * B.nLoci + lb;
                const prev = bestOf.get(key);
                if (prev !== undefined) {
                    const better = qual > cQ[prev] || (qual === cQ[prev] && (rs < cR[prev] ||
                        (rs === cR[prev] && (a < cHa[prev] || (a === cHa[prev] && b < cHb[prev])))));
                    if (!better) { continue; }
                    cHa[prev] = a; cHb[prev] = b; cQ[prev] = qual; cR[prev] = rs;
                } else {
                    bestOf.set(key, cHa.length);
                    cHa.push(a); cHb.push(b); cQ.push(qual); cR.push(rs);
                }
            }
        }

        const n = cHa.length, nCB = B.chromNames.length;
        const pair = new Uint32Array(n), rA = new Uint32Array(n), rB = new Uint32Array(n);
        let maxPair = 0;
        for (let k = 0; k < n; k++) {
            pair[k] = A.chrom[cHa[k]] * nCB + B.chrom[cHb[k]];
            rA[k] = rkA[A.locus[cHa[k]]];
            rB[k] = rkB[B.locus[cHb[k]]];
            if (pair[k] > maxPair) { maxPair = pair[k]; }
        }
        // LSD radix: rankB, then rankA, then pair -- (pair, rankA, rankB) is
        // unique per anchor, so this order is total
        let ord = new Uint32Array(n);
        for (let k = 0; k < n; k++) { ord[k] = k; }
        ord = countingSort(ord, rB, LB.maxRank + 1);
        ord = countingSort(ord, rA, LA.maxRank + 1);
        ord = countingSort(ord, pair, maxPair + 1);

        const anchors = {
            n, ha: new Uint32Array(n), hb: new Uint32Array(n), qual: new Uint16Array(n),
            pair: new Uint32Array(n), rA: new Uint32Array(n), rB: new Uint32Array(n),
        };
        for (let k = 0; k < n; k++) {
            const c = ord[k];
            anchors.ha[k] = cHa[c]; anchors.hb[k] = cHb[c]; anchors.qual[k] = cQ[c];
            anchors.pair[k] = pair[c]; anchors.rA[k] = rA[c]; anchors.rB[k] = rB[c];
        }
        return anchors;
    }

    // ---------------------------------------------------------------- stage 2: DP + extraction

    function chainAnchors(an, q) {
        const n = an.n, G = q.maxGap, H = q.maxLookback, pen = q.gapPenalty;
        const rA = an.rA, rB = an.rB;
        const f = [new Float64Array(n), new Float64Array(n)];
        const par = [new Int32Array(n), new Int32Array(n)];
        for (let s = 0; s < n;) {
            let e = s;
            const pr = an.pair[s];
            while (e < n && an.pair[e] === pr) { e++; }
            for (let o = 0; o < 2; o++) {
                const F = f[o], P = par[o];
                for (let i = s; i < e; i++) {
                    let best = 1, bj = -1, valid = 0;
                    for (let j = i - 1; j >= s; j--) {
                        const dA = rA[i] - rA[j];
                        if (dA > G) { break; }
                        if (dA === 0) { continue; }
                        const dB = o === 0 ? rB[i] - rB[j] : rB[j] - rB[i];
                        if (dB < 1 || dB > G) { continue; }
                        const sc = F[j] + 1 - pen * (Math.max(dA, dB) - 1);
                        if (sc > best) { best = sc; bj = j; }
                        if (++valid >= H) { break; }
                    }
                    F[i] = best; P[i] = bj;
                }
            }
            s = e;
        }

        // spec 4.5: candidates c = o*n + i, by f descending, ties by c
        // ascending ('+' first, then lower anchor index)
        let cand;
        if (pen === 0) {
            // f is an integer chain length in [1, n] -- counting sort by
            // (maxF - f) keeps the tie order for free, since it's stable
            let maxF = 1;
            for (let o = 0; o < 2; o++) { for (let i = 0; i < n; i++) { if (f[o][i] > maxF) { maxF = f[o][i]; } } }
            const key = new Uint32Array(2 * n), idx = new Uint32Array(2 * n);
            for (let c = 0; c < 2 * n; c++) { idx[c] = c; key[c] = maxF - f[c < n ? 0 : 1][c < n ? c : c - n]; }
            cand = countingSort(idx, key, maxF);
        } else {
            cand = new Uint32Array(2 * n);
            for (let c = 0; c < 2 * n; c++) { cand[c] = c; }
            const fv = (c) => (c < n ? f[0][c] : f[1][c - n]);
            cand.sort((x, y) => (fv(y) - fv(x)) || (x - y));
        }

        const used = new Uint8Array(n);
        const members = new Uint32Array(n);   // every anchor lands in exactly one chain
        const chainOff = [0], chainOri = [];
        let m = 0;
        for (let k = 0; k < cand.length; k++) {
            const c = cand[k], o = c < n ? 0 : 1;
            let i = o === 0 ? c : c - n;
            if (used[i]) { continue; }
            const from = m;
            while (i >= 0 && !used[i]) { used[i] = 1; members[m++] = i; i = par[o][i]; }
            // walked end -> start; store in chain order (ascending rankA)
            for (let a = from, b = m - 1; a < b; a++, b--) { const t = members[a]; members[a] = members[b]; members[b] = t; }
            chainOff.push(m); chainOri.push(o);
        }
        return {members, chainOff: Uint32Array.from(chainOff), chainOri: Uint8Array.from(chainOri)};
    }

    // ---------------------------------------------------------------- stage 3: blocks

    function makeBlock(A, B, an, ch, c) {
        const lo = ch.chainOff[c], hi = ch.chainOff[c + 1];
        let qs = Infinity, qe = -Infinity, ss = Infinity, se = -Infinity, qsum = 0;
        for (let k = lo; k < hi; k++) {
            const x = ch.members[k], ha = an.ha[x], hb = an.hb[x];
            if (A.start[ha] < qs) { qs = A.start[ha]; }
            if (A.end[ha] > qe) { qe = A.end[ha]; }
            if (B.start[hb] < ss) { ss = B.start[hb]; }
            if (B.end[hb] > se) { se = B.end[hb]; }
            qsum += an.qual[x];
        }
        const nA = hi - lo, x0 = ch.members[lo];
        return {
            qChrom: A.chrom[an.ha[x0]], sChrom: B.chrom[an.hb[x0]],
            qStart: qs, qEnd: qe, sStart: ss, sEnd: se,
            orientation: ch.chainOri[c] === 0 ? '+' : '-',
            nAnchors: nA,
            meanPositive: qsum / nA / 10000,
            anchorDensity: nA / (Math.max(qe - qs, 1) / 1e6),
            members: ch.members.subarray(lo, hi),
        };
    }

    function compareBlocks(a, b) {
        return (a.qChrom - b.qChrom) || (a.qStart - b.qStart) || (a.sChrom - b.sChrom) || (a.sStart - b.sStart)
            || (a.orientation < b.orientation ? -1 : a.orientation > b.orientation ? 1 : 0)
            || (a.qEnd - b.qEnd) || (a.sEnd - b.sEnd) || (a.members[0] - b.members[0]);
    }

    // ---------------------------------------------------------------- driver

    // Caches each stage by the parameters it depends on (spec 4): anchors
    // by (minPositive, maxHitRank, and maxGap in self mode), chains by the
    // anchor key plus (maxGap, maxLookback, gapPenalty); minBlock is a pure
    // filter over the cached chains, so changing it re-chains nothing.
    function createChainer(A, B, opts) {
        const selfMode = !!(opts && opts.selfMode);
        if (selfMode && B && B !== A) { throw new Error('chain: selfMode needs the same table on both sides'); }
        B = selfMode ? A : B;
        const nProteins = (opts && opts.nProteins) || Math.max(maxPlusOne(A.prot), maxPlusOne(B.prot));
        prepareTable(A, nProteins); prepareTable(B, nProteins);
        const cache = {anchorKey: null, anchors: null, chainKey: null, chains: null, blocksByChain: null};

        function run(params) {
            const q = resolveParams(params, selfMode);
            const t0 = now();
            const anchorKey = [q._minPositiveInt, q.maxHitRank, selfMode ? q.maxGap : 0].join(',');
            if (cache.anchorKey !== anchorKey) {
                cache.anchors = buildAnchors(A, B, q);
                cache.anchorKey = anchorKey; cache.chainKey = null;
            }
            const t1 = now();
            const chainKey = [anchorKey, q.maxGap, q.maxLookback, q.gapPenalty].join(',');
            if (cache.chainKey !== chainKey) {
                cache.chains = chainAnchors(cache.anchors, q);
                cache.chainKey = chainKey; cache.blocksByChain = new Map();
            }
            const t2 = now();
            const ch = cache.chains, blocks = [];
            for (let c = 0; c + 1 < ch.chainOff.length; c++) {
                if (ch.chainOff[c + 1] - ch.chainOff[c] < q.minBlock) { continue; }
                let b = cache.blocksByChain.get(c);
                if (!b) { b = makeBlock(A, B, cache.anchors, ch, c); cache.blocksByChain.set(c, b); }
                blocks.push(b);
            }
            blocks.sort(compareBlocks);
            const t3 = now();
            return {
                version: VERSION, params: publicParams(q), A, B, anchors: cache.anchors, blocks,
                stats: {nAnchors: cache.anchors.n, nChains: ch.chainOff.length - 1, nBlocks: blocks.length,
                        ms: {anchors: t1 - t0, chains: t2 - t1, blocks: t3 - t2, total: t3 - t0}},
            };
        }
        return {run, selfMode};
    }

    function chain(A, B, params) {
        return createChainer(A, B, {selfMode: params && params.selfMode}).run(params);
    }

    function publicParams(q) {
        const out = {};
        for (const k of Object.keys(q)) { if (k[0] !== '_') { out[k] = q[k]; } }
        return out;
    }

    function maxPlusOne(arr) { let m = -1; for (let i = 0; i < arr.length; i++) { if (arr[i] > m) { m = arr[i]; } } return m + 1; }

    function now() {
        return (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
    }

    // ---------------------------------------------------------------- outputs

    // link objects in the shape the page already renders (spec 5)
    function blocksToLinks(result) {
        const A = result.A, B = result.B;
        return result.blocks.map((b) => ({
            q_chrom: A.chromNames[b.qChrom], q_start: b.qStart, q_end: b.qEnd,
            s_chrom: B.chromNames[b.sChrom], s_start: b.sStart, s_end: b.sEnd,
            score: b.nAnchors, orientation: b.orientation,
            mean_identity: b.meanPositive, anchor_density: b.anchorDensity,
        }));
    }

    const TSV_HEADER = 'query_chrom\tquery_start\tquery_end\tsubject_chrom\tsubject_start\tsubject_end\t'
        + 'score\torientation\tmean_identity\tanchor_density\n';

    function linksToTsv(links) {
        let s = TSV_HEADER;
        for (const l of links) {
            s += `${l.q_chrom}\t${l.q_start}\t${l.q_end}\t${l.s_chrom}\t${l.s_start}\t${l.s_end}\t`
                + `${l.score}\t${l.orientation}\t${l.mean_identity.toFixed(4)}\t${l.anchor_density.toFixed(2)}\n`;
        }
        return s;
    }

    function blocksToTsv(result) { return linksToTsv(blocksToLinks(result)); }

    // ---------------------------------------------------------------- embedded payload (spec 6)

    function b64ToBytes(b64) {
        if (typeof atob === 'function') {
            const bin = atob(b64), out = new Uint8Array(bin.length);
            for (let i = 0; i < bin.length; i++) { out[i] = bin.charCodeAt(i); }
            return out;
        }
        return new Uint8Array(Buffer.from(b64, 'base64'));
    }

    async function gunzipDefault(bytes) {
        const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
        return new Uint8Array(await new Response(stream).arrayBuffer());
    }

    const DTYPES = {u8: Uint8Array, u16: Uint16Array, u32: Uint32Array, f64: Float64Array};

    async function decodeColumn(col, gunzip) {
        const bytes = await gunzip(b64ToBytes(col.data));
        const T = DTYPES[col.dtype];
        if (!T) { throw new Error(`payload: unknown dtype ${col.dtype}`); }
        // copy into an aligned buffer: the decompressed bytes may not start
        // on an element boundary
        const buf = new ArrayBuffer(bytes.byteLength);
        new Uint8Array(buf).set(bytes);
        return new T(buf);
    }

    async function decodeGenome(g, chroms, gunzip) {
        const c = {};
        for (const name of Object.keys(g.cols)) { c[name] = await decodeColumn(g.cols[name], gunzip); }
        const n = g.n;
        const t = {
            name: chroms.name, chromNames: chroms.names.slice(), chromSizes: chroms.sizes.slice(), n,
            chrom: Uint16Array.from(c.chrom), start: new Float64Array(n), end: new Float64Array(n),
            strand: Uint8Array.from(c.strand), positive: Uint16Array.from(c.positive), rank: Uint8Array.from(c.rank),
            prot: Uint32Array.from(c.prot), locus: Uint32Array.from(c.locus), nLoci: g.nLoci,
        };
        let prevChrom = -1, prev = 0;
        for (let i = 0; i < n; i++) {
            if (t.chrom[i] !== prevChrom) { prevChrom = t.chrom[i]; prev = 0; }
            prev += c.startDelta[i];
            t.start[i] = prev;
            t.end[i] = prev + c.length[i];
        }
        return t;
    }

    // -> {proteins, target, reference}; chroms = {target: {names, sizes},
    // reference: {names, sizes}}; gunzip is injectable for environments
    // without DecompressionStream
    async function decodePayload(payload, chroms, gunzip) {
        if (payload.version !== VERSION) { throw new Error(`payload: version ${payload.version}, expected ${VERSION}`); }
        gunzip = gunzip || gunzipDefault;
        const text = new TextDecoder().decode(await gunzip(b64ToBytes(payload.proteins)));
        const proteins = text.length ? text.split('\n') : [];
        const target = await decodeGenome(payload.genomes.target, {name: 'target', ...chroms.target}, gunzip);
        const reference = payload.genomes.reference === 'same_as_target'
            ? Object.assign({}, target, {name: 'reference', chromNames: chroms.reference.names.slice(),
                                          chromSizes: chroms.reference.sizes.slice()})
            : await decodeGenome(payload.genomes.reference, {name: 'reference', ...chroms.reference}, gunzip);
        return {proteins, target, reference};
    }

    const api = {
        VERSION, DEFAULTS, MIN_CHAIN_LENGTH,
        parseHitsTsv, makeTables, prepareTable, meanBestPositive, autoParams,
        createChainer, chain, blocksToLinks, linksToTsv, blocksToTsv, decodePayload,
    };
    if (typeof module !== 'undefined' && module.exports) { module.exports = api; } else { root.SYNCHAIN = api; }
})(typeof globalThis !== 'undefined' ? globalThis : this);
