// Synteny chainer over per-genome hit tables -- the only implementation,
// shared by the pipeline (node, via bin/chain_blocks.mjs) and the
// interactive page (inlined into the HTML and run in a Web Worker). This
// header is the contract every producer and consumer follows; the on-disk
// hit-table format is specified in bin/extract_hits.py's docstring.
// Changing any rule below means bumping VERSION.
//
// Dependency-free and DOM-free on purpose: the same source text must run
// as a CommonJS module, a classic <script>, and a Worker body.
//
// HIT TABLE (in memory; built by makeTables from the TSV, or by
// decodePayload from the page's embedded copy -- both must give identical
// arrays):
//   {name, chromNames (chrom.sizes order), chromSizes, n,
//    chrom: Uint16 (index into chromNames), start/end: Float64 (exact bp,
//    safe beyond 4.29 Gb), strand: Uint8 (1 '+', 0 '-'),
//    positive: Uint16 (round(Positive * 10000)), rank: Uint8 (min(Rank, 255)),
//    prot: Uint32 (index into the shared proteins list), locus: Uint32, nLoci}
//   proteins: the union of both genomes' accessions, sorted, shared by both
//   tables of a comparison. Rows keep the TSV order. prepareTable attaches
//   derived indexes under _idx; nothing outside this file may rely on them.
//
// PARAMETERS:
//   minPositive  0-1          a hit passes only if positive >= minPositive
//                             (default 0.7)
//   maxHitRank   int >= 1     ... and rank <= maxHitRank (default 255: all)
//   maxGap       int (loci)   max rank step between consecutive chain
//                             members, on both genomes (default 25)
//   maxLookback  int          max valid predecessors examined per anchor (50)
//   gapPenalty   float >= 0   score cost per skipped locus (default 0)
//   minBlock     int (loci)   minimum chain length reported as a block
//                             (default 15)
//   selfMode     bool         both sides are the same table (homeolog scan)
//   Nothing is tuned to the data: the defaults are where the page's controls
//   start. gapPenalty stays 0 because, while 0.1 trims the end noise of
//   consequence (b) below, it also splits real, gappy ancient-duplication
//   blocks. The page runs one parameter set for all three of its chains
//   (cross and both self). meanBestPositive (mean over proteins of their
//   rank-1 positive, else their best) only serves chain_blocks.mjs's sanity
//   check.
//
// CHAINING RULES (normative):
//   1. A locus passes if any of its hits passes. A passing locus's rank
//      is its 1-based position among passing loci of its chromosome, in
//      locus-id order; other loci have no rank.
//   2. Anchors: for each protein, every (passing hit a in A, passing hit b in
//      B), collapsed to one per (locusA, locusB) keeping the highest
//      min(positive_a, positive_b); ties -> lower rank_a + rank_b, then lower
//      a, then lower b. Then MCScanX's repeat collapse (its -w,
//      COLLAPSE_WINDOW = 5), per (chromA, chromB): anchors sharing rankA,
//      each with rankB within 5 of the previous one's, form a run, of which
//      only the highest min(positive) stays (the lowest rankB on ties); then
//      the same along B over what's left (shared rankB, rankA within 5,
//      lowest rankA on ties). Tandem copies would otherwise chain in parallel.
//   3. Self mode keeps only locusA < locusB and drops same-chromosome pairs
//      with |rankA - rankB| <= maxGap (tandem neighbours).
//   4. DP per (chromA, chromB) and orientation o in {+, -}, anchors sorted by
//      (rankA, rankB): for anchor i, scan earlier anchors j backwards,
//      stopping at the first with rankA_i - rankA_j > maxGap or after
//      maxLookback valid predecessors. j is valid iff rankA_i - rankA_j >= 1
//      and 1 <= dB <= maxGap, dB = rankB_i - rankB_j for '+' and
//      rankB_j - rankB_i for '-'. f_i = max(1, max_j f_j + 1 -
//      gapPenalty * (max(dA, dB) - 1)); the first j found wins ties.
//   5. Extraction, per (chromA, chromB), in rounds. Each round runs rule 4
//      over the anchors still free, then takes chain ends that have a
//      parent (chains of two or more) by f descending, ties '+' first, then
//      lower anchor index. From each end not yet taken this round, walk the
//      parent links: a walk that reaches an anchor taken this round ends the
//      round untaken; otherwise the whole walk is one chain. The anchors
//      taken leave, and the next round reruns the DP; a round with no such
//      walk is the last. Anchors left over are in no chain. This is
//      MCScanX's dagchainer without its score cutoff: a walk is never cut
//      short, and a lone anchor never takes another chain's anchor.
//   6. Blocks are chains with length >= minBlock. Since extraction ignores
//      minBlock, the blocks at threshold T are exactly the blocks at any
//      T' < T filtered to length >= T -- a min-block control is a pure filter.
//   7. Strict monotonicity on both genomes makes every chain member a distinct
//      locus on both sides, so block size = number of independent loci.
//   Consequences:
//   (a) Anchor deserts don't break chains: distance counts passing loci, so two
//       loci with nothing between them are adjacent however many bp apart.
//       That is what makes chaining genome-size invariant.
//   (b) Noise can extend a block's ends, never its interior: monotonicity
//       only leaves room beyond a block's last true anchor.
//   (c) A block can bridge a shorter inversion: the '+' chain jumps from one
//       flank to the inversion's middle anchor and on to the other flank
//       (any inversion under about 2 * maxGap loci), while the rest of the
//       inversion is its own '-' block, so the two spans overlap, as in
//       MCScanX. Rule 5's rounds are what keep that '-' block: in a single
//       pass, the inverted anchors' high '+' scores (carried on from the
//       flank) would take them first, one anchor at a time.
//
// OUTPUTS: Block = {qChrom, sChrom (chrom indices), qStart, qEnd, sStart,
//   sEnd (min start / max end over member hits), orientation, nAnchors
//   (chain length), meanPositive (mean of min(positive_a, positive_b)),
//   anchorDensity (nAnchors per Mb of query span, span floored at 1 bp),
//   members (anchor indices in chain order)}. blocksToLinks gives the link
//   objects the page renders ({q_chrom, ..., score: nAnchors, orientation,
//   mean_identity, anchor_density}); links.tsv (blocksToTsv, the node CLI and
//   the page's "blocks TSV" download) has the columns target_chrom,
//   target_start, target_end, reference_chrom, reference_start, reference_end,
//   score, orientation, mean_identity (4 dp), anchor_density (2 dp), rows
//   ordered by (target chrom index, target start, reference chrom index,
//   reference start, orientation). Given the parameters, linksToTsv puts them
//   on a first line ahead of the header, named after chain_blocks.mjs's
//   flags: `# min_identity=0.7 max_gap=25 min_block=5 max_hit_rank=255
//   max_lookback=50 gap_penalty=0`, then ` self=true` for a self chain and
//   any further filters the caller names (the page's min_seq_size, names).
//
// EMBEDDED PAYLOAD (SYN.hitsPayload, written by plot_synteny_interactive.py,
// read by decodePayload):
//   {version, proteins: b64(gzip('\n'-joined accessions)),
//    genomes: {target: GenomeColumns, reference: GenomeColumns | 'same_as_target'}}
//   GenomeColumns = {n, nLoci, cols: {name: {dtype: u8|u16|u32|f64,
//   data: b64(gzip(little-endian bytes))}}}, with columns chrom (u16),
//   startDelta (u32: start minus the previous row's start on the same chrom,
//   first row of a chrom relative to 0; f64 if any value >= 2^32), length
//   (u32 end - start, f64 if needed), strand (u8), positive (u16), rank (u8),
//   prot (u32), locus (u32). 'same_as_target' only when both genomes' hit
//   tables and chrom.sizes orders are identical. Chrom names/sizes come from
//   the page's own size tables.
(function (root) {
    'use strict';

    const VERSION = 2;
    const MIN_CHAIN_LENGTH = 2;   // a single anchor is never a block, whatever minBlock says
    const COLLAPSE_WINDOW = 5;    // CHAINING RULES 2: MCScanX's default -w

    const DEFAULTS = {
        minPositive: 0.7,
        minBlock: 15,
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

    // HIT TABLEs for one comparison (two genomes, or one for self) from parsed
    // TSV rows. chromNames is first-appearance order, i.e. chrom.sizes order
    // restricted to sequences with hits (rows are sorted that way).
    function makeTables(genomes) {
        const accessions = new Set();
        for (const g of genomes) { for (const r of g.rows) { accessions.add(r[8]); } }
        const proteins = Array.from(accessions).sort();
        const protIndex = new Map(proteins.map((p, i) => [p, i]));
        const tables = genomes.map((g) => {
            const chromNames = [], chromIndex = new Map();
            const n = g.rows.length;
            const t = {
                name: g.name, chromNames, chromSizes: null, n,
                chrom: new Uint16Array(n), start: new Float64Array(n), end: new Float64Array(n),
                strand: new Uint8Array(n), positive: new Uint16Array(n), rank: new Uint8Array(n),
                prot: new Uint32Array(n), locus: new Uint32Array(n), nLoci: 0,
            };
            for (let i = 0; i < n; i++) {
                const r = g.rows[i];
                let ci = chromIndex.get(r[0]);
                if (ci === undefined) { ci = chromNames.length; chromNames.push(r[0]); chromIndex.set(r[0], ci); }
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

    // see PARAMETERS above
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

    // CHAINING RULES 1-3: one anchor per (locusA, locusB), ordered by (chrom pair, rankA, rankB)
    function buildAnchors(A, B, q) {
        const thr = q._minPositiveInt, mr = q.maxHitRank, self = q.selfMode;
        const LA = locusRanks(A, thr, mr), LB = self ? LA : locusRanks(B, thr, mr);
        const rkA = LA.rank, rkB = LB.rank;
        const pass = (t, i) => t.positive[i] >= thr && t.rank[i] <= mr;
        const offB = B._idx.protOff, byB = B._idx.byProt;

        // best candidate per (locusA, locusB); tie-breaks per rule 2
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

        const keep = collapseRepeats(ord, pair, rA, rB, cQ);
        let nk = 0;
        for (let k = 0; k < n; k++) { if (keep[k]) { nk++; } }
        const anchors = {
            n: nk, ha: new Uint32Array(nk), hb: new Uint32Array(nk), qual: new Uint16Array(nk),
            pair: new Uint32Array(nk), rA: new Uint32Array(nk), rB: new Uint32Array(nk),
        };
        for (let k = 0, x = 0; k < n; k++) {
            if (!keep[k]) { continue; }
            const c = ord[k];
            anchors.ha[x] = cHa[c]; anchors.hb[x] = cHb[c]; anchors.qual[x] = cQ[c];
            anchors.pair[x] = pair[c]; anchors.rA[x] = rA[c]; anchors.rB[x] = rB[c];
            x++;
        }
        return anchors;
    }

    // CHAINING RULES 2, second half: keep flags by position in `ord`
    // (candidate indices sorted by pair, rankA, rankB)
    function collapseRepeats(ord, pair, rA, rB, qual) {
        const n = ord.length, keep = new Uint8Array(n).fill(1);
        const pass = (pos, same, step) => {
            for (let s = 0; s < pos.length;) {
                let best = pos[s], e = s + 1;
                while (e < pos.length && same(pos[e - 1], pos[e]) && step(pos[e - 1], pos[e]) <= COLLAPSE_WINDOW) {
                    if (qual[ord[pos[e]]] > qual[ord[best]]) { best = pos[e]; }
                    e++;
                }
                for (let k = s; k < e; k++) { if (pos[k] !== best) { keep[pos[k]] = 0; } }
                s = e;
            }
        };
        const at = (k) => ord[k];
        // along A: `ord` is already in (pair, rankA, rankB) order
        const byA = Array.from({length: n}, (_, k) => k);
        pass(byA, (x, y) => pair[at(x)] === pair[at(y)] && rA[at(x)] === rA[at(y)], (x, y) => rB[at(y)] - rB[at(x)]);
        // along B: the survivors in (pair, rankB, rankA) order
        const byB = byA.filter((k) => keep[k]);
        byB.sort((x, y) => (pair[at(x)] - pair[at(y)]) || (rB[at(x)] - rB[at(y)]) || (rA[at(x)] - rA[at(y)]));
        pass(byB, (x, y) => pair[at(x)] === pair[at(y)] && rB[at(x)] === rB[at(y)], (x, y) => rA[at(y)] - rA[at(x)]);
        return keep;
    }

    // ---------------------------------------------------------------- stage 2: DP + extraction

    // CHAINING RULES 4-5. A round's first candidate is always taken whole, so
    // every round makes progress.
    function chainAnchors(an, q) {
        const n = an.n, G = q.maxGap, H = q.maxLookback, pen = q.gapPenalty;
        const rA = an.rA, rB = an.rB;
        const free = new Uint8Array(n).fill(1);
        const members = new Uint32Array(n);
        const chainOff = [0], chainOri = [];
        let m = 0;
        for (let s = 0; s < n;) {
            let e = s;
            const pr = an.pair[s];
            while (e < n && an.pair[e] === pr) { e++; }
            for (;;) {
                const idx = [];
                for (let i = s; i < e; i++) { if (free[i]) { idx.push(i); } }
                const k = idx.length;
                if (k < MIN_CHAIN_LENGTH) { break; }
                const f = [new Float64Array(k), new Float64Array(k)];
                const par = [new Int32Array(k), new Int32Array(k)];
                for (let o = 0; o < 2; o++) {
                    const F = f[o], P = par[o];
                    for (let x = 0; x < k; x++) {
                        const i = idx[x];
                        let best = 1, bj = -1, valid = 0;
                        for (let y = x - 1; y >= 0; y--) {
                            const j = idx[y];
                            const dA = rA[i] - rA[j];
                            if (dA > G) { break; }
                            if (dA === 0) { continue; }
                            const dB = o === 0 ? rB[i] - rB[j] : rB[j] - rB[i];
                            if (dB < 1 || dB > G) { continue; }
                            const sc = F[y] + 1 - pen * (Math.max(dA, dB) - 1);
                            if (sc > best) { best = sc; bj = y; }
                            if (++valid >= H) { break; }
                        }
                        F[x] = best; P[x] = bj;
                    }
                }

                // rule 5: chain ends c = o*k + x (two or more members), by f
                // descending, ties by c ascending ('+' first, then lower anchor)
                const cand = [];
                for (let c = 0; c < 2 * k; c++) { if (par[c < k ? 0 : 1][c < k ? c : c - k] >= 0) { cand.push(c); } }
                if (!cand.length) { break; }
                const fv = (c) => (c < k ? f[0][c] : f[1][c - k]);
                cand.sort((x, y) => (fv(y) - fv(x)) || (x - y));

                const used = new Uint8Array(k);
                let collided = false;
                for (const c of cand) {
                    const o = c < k ? 0 : 1, P = par[o];
                    const end = c < k ? c : c - k;
                    if (used[end]) { continue; }
                    let len = 0, y = end;
                    while (y >= 0 && !used[y]) { len++; y = P[y]; }
                    if (y >= 0) { collided = true; break; }
                    // walked end -> start; store in chain order (ascending rankA)
                    m += len;
                    for (let w = m - 1, z = end; z >= 0; w--, z = P[z]) { used[z] = 1; members[w] = idx[z]; }
                    chainOff.push(m); chainOri.push(o);
                }
                for (let x = 0; x < k; x++) { if (used[x]) { free[idx[x]] = 0; } }
                if (!collided) { break; }
            }
            s = e;
        }
        return {members: members.slice(0, m), chainOff: Uint32Array.from(chainOff), chainOri: Uint8Array.from(chainOri)};
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

    // Caches each stage by the parameters it depends on: anchors by
    // (minPositive, maxHitRank, and maxGap in self mode), chains by those plus
    // (maxGap, maxLookback, gapPenalty); minBlock only filters cached chains.
    function createChainer(A, B, opts) {
        const selfMode = !!(opts && opts.selfMode);
        if (selfMode && B && B !== A) { throw new Error('chain: selfMode needs the same table on both sides'); }
        B = selfMode ? A : B;
        const nProteins = Math.max(maxPlusOne(A.prot), maxPlusOne(B.prot));
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
            const chainKey = [anchorKey, q.maxGap, q.maxLookback, q.gapPenalty].join(',');
            if (cache.chainKey !== chainKey) {
                cache.chains = chainAnchors(cache.anchors, q);
                cache.chainKey = chainKey; cache.blocksByChain = new Map();
            }
            const ch = cache.chains, blocks = [];
            for (let c = 0; c + 1 < ch.chainOff.length; c++) {
                if (ch.chainOff[c + 1] - ch.chainOff[c] < q.minBlock) { continue; }
                let b = cache.blocksByChain.get(c);
                if (!b) { b = makeBlock(A, B, cache.anchors, ch, c); cache.blocksByChain.set(c, b); }
                blocks.push(b);
            }
            blocks.sort(compareBlocks);
            return {
                params: publicParams(q), A, B, anchors: cache.anchors, blocks,
                stats: {nAnchors: cache.anchors.n, nChains: ch.chainOff.length - 1, nBlocks: blocks.length,
                        ms: now() - t0},
            };
        }
        return {run};
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

    // OUTPUTS above
    function blocksToLinks(result) {
        const A = result.A, B = result.B;
        return result.blocks.map((b) => ({
            q_chrom: A.chromNames[b.qChrom], q_start: b.qStart, q_end: b.qEnd,
            s_chrom: B.chromNames[b.sChrom], s_start: b.sStart, s_end: b.sEnd,
            score: b.nAnchors, orientation: b.orientation,
            mean_identity: b.meanPositive, anchor_density: b.anchorDensity,
        }));
    }

    const TSV_HEADER = 'target_chrom\ttarget_start\ttarget_end\treference_chrom\treference_start\treference_end\t'
        + 'score\torientation\tmean_identity\tanchor_density\n';

    // OUTPUTS above; extra: further [key, value] pairs, appended as given
    function settingsLine(params, extra) {
        const p = Object.assign({}, DEFAULTS, params);
        const pairs = [['min_identity', +p.minPositive.toFixed(4)], ['max_gap', p.maxGap], ['min_block', p.minBlock],
            ['max_hit_rank', p.maxHitRank], ['max_lookback', p.maxLookback], ['gap_penalty', p.gapPenalty]];
        if (p.selfMode) { pairs.push(['self', true]); }
        return '# ' + pairs.concat(extra || []).map(([k, v]) => `${k}=${v}`).join(' ') + '\n';
    }

    function linksToTsv(links, params, extra) {
        let s = (params ? settingsLine(params, extra) : '') + TSV_HEADER;
        for (const l of links) {
            s += `${l.q_chrom}\t${l.q_start}\t${l.q_end}\t${l.s_chrom}\t${l.s_start}\t${l.s_end}\t`
                + `${l.score}\t${l.orientation}\t${l.mean_identity.toFixed(4)}\t${l.anchor_density.toFixed(2)}\n`;
        }
        return s;
    }

    function blocksToTsv(result) { return linksToTsv(blocksToLinks(result), result.params); }

    // ---------------------------------------------------------------- embedded payload (EMBEDDED PAYLOAD above)

    function b64ToBytes(b64) {
        const bin = atob(b64), out = new Uint8Array(bin.length);
        for (let i = 0; i < bin.length; i++) { out[i] = bin.charCodeAt(i); }
        return out;
    }

    async function gunzip(bytes) {
        const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
        return new Uint8Array(await new Response(stream).arrayBuffer());
    }

    const DTYPES = {u8: Uint8Array, u16: Uint16Array, u32: Uint32Array, f64: Float64Array};

    async function decodeColumn(col) {
        const bytes = await gunzip(b64ToBytes(col.data));
        const T = DTYPES[col.dtype];
        if (!T) { throw new Error(`payload: unknown dtype ${col.dtype}`); }
        // the decompressed bytes may not start on an element boundary
        const buf = new ArrayBuffer(bytes.byteLength);
        new Uint8Array(buf).set(bytes);
        return new T(buf);
    }

    async function decodeGenome(g, chroms) {
        const c = {};
        for (const name of Object.keys(g.cols)) { c[name] = await decodeColumn(g.cols[name]); }
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
    // reference: {names, sizes}}
    async function decodePayload(payload, chroms) {
        if (payload.version !== VERSION) { throw new Error(`payload: version ${payload.version}, expected ${VERSION}`); }
        const text = new TextDecoder().decode(await gunzip(b64ToBytes(payload.proteins)));
        const proteins = text.length ? text.split('\n') : [];
        const target = await decodeGenome(payload.genomes.target, {name: 'target', ...chroms.target});
        const reference = payload.genomes.reference === 'same_as_target'
            ? Object.assign({}, target, {name: 'reference', chromNames: chroms.reference.names.slice(),
                                          chromSizes: chroms.reference.sizes.slice()})
            : await decodeGenome(payload.genomes.reference, {name: 'reference', ...chroms.reference});
        return {proteins, target, reference};
    }

    const api = {
        VERSION, DEFAULTS,
        parseHitsTsv, makeTables, meanBestPositive,
        createChainer, chain, blocksToLinks, settingsLine, linksToTsv, blocksToTsv, decodePayload,
    };
    if (typeof module !== 'undefined' && module.exports) { module.exports = api; } else { root.SYNCHAIN = api; }
})(typeof globalThis !== 'undefined' ? globalThis : this);
