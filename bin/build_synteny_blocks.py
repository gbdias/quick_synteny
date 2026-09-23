#!/usr/bin/env python3
"""Build synteny blocks directly from two miniprot GFFs -- no separate
ortholog-finding aligner needed.

The same proteome was already aligned against both genomes to build gene
models (see miniprot_align.nf). Each protein's hit positions in genome A and
genome B *are* the raw synteny anchors: for every protein, every (hit in A,
hit in B) pair is a candidate anchor. Running the same pairing within one
genome's own hits (a protein hitting two different loci in the SAME genome)
gives homeolog anchors for free, from data already produced -- no self-
comparison run needed.

Guardrails (kept deliberately simple, "quick and dirty" but not naive):
  - drop anchors below --min_identity (miniprot's own Positive= per hit --
    see parse_gff() for why Positive, not the stricter Identity), auto-tuned
    from the data by default instead of a fixed cutoff -- see
    auto_min_identity()
  - collapse anchors landing within a small window on both sides (tandem
    duplicates), keeping the higher-identity/lower-rank one
  - chain anchors per chromosome pair via a longest-monotonic-subsequence
    (increasing for '+' blocks, decreasing for '-'/inverted blocks), with a
    small allowed regression to tolerate local noise, same idea MCScanX
    uses without its full machinery
  - drop chains shorter than --min_block_anchors DISTINCT PROTEINS (not raw
    anchor count -- see find_blocks()), also auto-tuned by default -- see
    auto_min_block_anchors()

Output is the same links.tsv format plot_synteny_interactive.py already reads
(query_chrom, query_start, query_end, subject_chrom, subject_start,
subject_end, score, orientation), plus two summary columns computed from the
anchors before they're discarded: mean_identity (mean of each anchor's
min(query hit identity, subject hit identity), i.e. the same per-anchor
"quality" the chainer already scored on) and anchor_density (score per Mb of
query span -- how tightly packed the block's supporting genes are). score
itself is the number of DISTINCT PROTEINS backing the block, not the raw
count of (hit-in-A, hit-in-B) anchor pairs -- see find_blocks()'s docstring
for why raw count can wildly overstate independent evidence when a gene has
several hits on one or both sides. No anchors/BED intermediate needed, since
this works in genomic coordinates from the start.
"""
import argparse
import re
import sys
from collections import defaultdict

MRNA_ATTR_RE = re.compile(r'(\w+)=([^;]+)')
TANDEM_WINDOW = 5000  # bp; anchors within this on both sides are treated as one


def parse_gff(path):
    """protein_id -> list of hits: (chrom, start, end, strand, identity, rank)

    "identity" here is actually miniprot's Positive= score (fraction of
    aligned residues that are identical OR a positive-scoring substitution
    under its protein scoring matrix -- the standard BLAST-style "similarity"
    /"positives" metric), not its stricter Identity= (exact-residue-match
    fraction only). Positive is the more appropriate divergence signal for
    real orthologs: conservative amino acid substitutions (chemically similar
    residues, e.g. Leu<->Ile) accumulate well before non-conservative ones do
    as two genomes diverge, so a genuine ortholog's Positive score stays
    meaningfully higher than its Identity score long before Identity alone
    would make it look like noise -- exactly the gap that made --min_identity
    need auto-tuning down so far for divergent pairs in the first place (see
    auto_min_identity()). Falls back to Identity= for the rare hit missing
    Positive=, though miniprot always reports both in practice."""
    hits = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'mRNA':
                continue
            chrom, start, end, strand, attrs_str = fields[0], fields[3], fields[4], fields[6], fields[8]
            attrs = dict(MRNA_ATTR_RE.findall(attrs_str))
            target = attrs.get('Target', '').split()[0]
            if not target:
                continue
            identity = float(attrs.get('Positive', attrs.get('Identity', 0)))
            rank = int(attrs.get('Rank', 1))
            hits[target].append((chrom, int(start) - 1, int(end), strand, identity, rank))
    return hits


def mean_best_identity(hits):
    """Average, across all aligned proteins, of each protein's best hit
    identity (its rank-1 hit, falling back to its highest-identity hit for
    the rare protein that somehow has none) -- same "how well did this
    genome's alignment go overall" signal compute_alignment_stats.py reports
    for the plots' stats inset, computed here too so --min_identity and
    --min_block_anchors can both auto-scale off it without an extra
    file/process dependency."""
    rank1, best = {}, {}
    for protein, protein_hits in hits.items():
        for h in protein_hits:
            identity = h[4]
            if h[5] == 1:
                rank1[protein] = identity
            if protein not in best or identity > best[protein]:
                best[protein] = identity
    combined = {**best, **rank1}
    return sum(combined.values()) / len(combined) if combined else 0.0


AUTO_MIN_IDENTITY_FLOOR = 0.3    # never auto-drop below this, however divergent the pair
AUTO_MIN_IDENTITY_CEIL = 0.9     # never auto-exceed the old fixed default


def auto_min_identity(weaker_mean_identity):
    """min_identity, scaled to the data instead of hard-coded at 0.9.

    0.9 was tuned against the close-relative pairs this pipeline was built
    and tested on (A. thaliana/A. suecica, sister Drosophila species), where
    a genuine ortholog's protein alignment should score very high on both
    sides -- so a hit scoring below 90% was almost certainly spurious noise,
    not a real anchor. That assumption breaks for genuinely divergent pairs
    (different genera, say): real orthologs there may only align at ~60-75%
    identity (see parse_gff()'s use of miniprot's Positive= score), and a
    fixed 90% floor silently discards nearly all of them, producing a
    near-empty result that looks like "no synteny" when the actual issue is
    just an overly strict guardrail.

    No margin is subtracted below the observed mean -- the threshold is
    exactly the weaker genome's own mean best-hit identity, clamped to
    [0.3, 0.9]. An earlier version subtracted a fixed margin first, but any
    fixed margin is itself another hard-coded number that has to be re-tuned
    per pair anyway (a margin picked to keep close relatives clean turned
    out to still let a divergent pair's noise floor through, and a margin
    picked for the divergent pair reopened the close-relative noise problem
    this whole guardrail exists to prevent -- see the results/parameter_tuning
    sweep). Using the mean directly needs no such tuning: it's a bar that's
    relative to this specific pair's own alignment quality by construction
    -- always passing the above-average (and so more likely genuine) hits
    and rejecting the below-average (and so more likely spurious) ones,
    proportionally, rather than by an arbitrary fixed amount that happens to
    suit one pair's distribution and not another's.

    Takes the weaker genome's own observed mean best-hit identity (the same
    metric compute_alignment_stats.py reports), not the average of both --
    computed once by the caller and shared with auto_min_block_anchors()
    rather than recomputed here -- because an anchor's own quality is already
    min(hit_a_identity, hit_b_identity): it's bounded by whichever genome
    aligned worse, so that's the more relevant number to scale against.

    The floor exists because a low mean identity is ambiguous: it may mean
    "these are just distant relatives" (auto-relaxing is the right call) or
    "this alignment/species pairing is bad" (auto-relaxing just lets junk
    through as spurious blocks) -- there's no way to tell those apart from
    the identity distribution alone, so this stops short of disabling the
    guardrail entirely no matter how low the observed identity gets. The
    ceiling keeps behavior identical to the old fixed default for any pair
    at or above 90% mean identity."""
    return max(AUTO_MIN_IDENTITY_FLOOR, min(AUTO_MIN_IDENTITY_CEIL, weaker_mean_identity))


AUTO_MIN_BLOCK_ANCHORS_THRESHOLD = 0.8  # weaker-genome mean identity cutoff between the two tiers below
AUTO_MIN_BLOCK_ANCHORS_HIGH = 15        # close relatives: real synteny already forms long chains,
                                         # so requiring more anchors per block costs little and
                                         # keeps out more noise
AUTO_MIN_BLOCK_ANCHORS_LOW = 5          # divergent pairs: per the results/parameter_tuning sweep,
                                         # real homologous-chromosome signal there is carried by many
                                         # SHORT locally-collinear chains (chromosome-wide order isn't
                                         # conserved even though short local runs still are), so
                                         # requiring long chains throws away most of the real signal


def auto_min_block_anchors(weaker_mean_identity):
    """min_block_anchors, scaled to the data instead of hard-coded at one
    value -- see the module docstring's chaining guardrail and
    auto_min_identity() just above for the same close-relative-vs-divergent
    reasoning. A close-relative pair (weaker mean identity >= 0.8) can afford
    a stricter, higher bar since real synteny blocks there are naturally
    long; a divergent pair needs a lower bar to catch short chains, per the
    empirical sweep in results/parameter_tuning (mba 5-15 tested against
    both a close-relative and a divergent real pair)."""
    return AUTO_MIN_BLOCK_ANCHORS_HIGH if weaker_mean_identity >= AUTO_MIN_BLOCK_ANCHORS_THRESHOLD \
        else AUTO_MIN_BLOCK_ANCHORS_LOW


def build_raw_anchors(hits_a, hits_b, min_identity):
    """(protein, chrom_a, pos_a, chrom_b, pos_b, quality, hit_a, hit_b) per pair.
    Pass the same dict for both to get within-genome (homeolog) anchors --
    self-pairs and mirrored duplicate pairs are skipped automatically.

    Every hit on side A pairs with every hit on side B (bounded already by
    miniprot's own -N/--outs caps, so this is small in practice). An earlier
    version restricted side A to its single primary (rank 1) hit to avoid
    over-pairing weak secondary hits on both sides at once -- but that's a
    query/target-vs-subject/comparison-genome asymmetry, not a biological
    one: when the polyploid genome is passed as side A (e.g. as
    --query_gff, the "target" role), restricting side A to rank 1 throws
    away exactly the homeolog copies the whole point of --show_homeologs is
    to find, and the cross-genome plot silently degrades to near-1:1
    whenever the polyploid side happens to be the target rather than the
    comparison genome. Symmetric pairing plus the identity/gap/block-size
    guardrails below controls the noise just as well without that blind
    spot."""
    self_mode = hits_a is hits_b
    anchors = []
    # sorted(), not the bare set intersection: Python's hash randomization
    # (on by default, a fresh seed per process) makes plain set iteration
    # order vary run to run even given byte-identical inputs, which -- via
    # collapse_tandem()'s and find_blocks()'s tie-breaking on ties in
    # position/quality -- can make two runs of this exact script on the exact
    # same GFFs pick a different (but same-size) set of anchors/blocks. That
    # briefly surfaced as a real inconsistency: the interactive plot's
    # min-block-anchors slider and the static plots are now built from two
    # separate invocations of this script (see build_synteny.nf's *_FOR_SLIDER
    # processes) that need to agree exactly at the same threshold.
    proteins = sorted(set(hits_a) & set(hits_b))
    for protein in proteins:
        list_a, list_b = hits_a[protein], hits_b[protein]
        for i, ha in enumerate(list_a):
            for j, hb in enumerate(list_b):
                if self_mode and i >= j:
                    continue
                if ha[4] < min_identity or hb[4] < min_identity:
                    continue
                pos_a = (ha[1] + ha[2]) // 2
                pos_b = (hb[1] + hb[2]) // 2
                quality = min(ha[4], hb[4])
                anchors.append({
                    'protein': protein, 'chrom_a': ha[0], 'pos_a': pos_a, 'hit_a': ha,
                    'chrom_b': hb[0], 'pos_b': pos_b, 'hit_b': hb, 'quality': quality,
                })
    return anchors


def collapse_tandem(anchors):
    """Within each chromosome-pair group, collapse anchors that land within
    TANDEM_WINDOW of an already-kept anchor on both sides, keeping the
    better (higher identity, then lower/primary rank) one."""
    by_pair = defaultdict(list)
    for a in anchors:
        by_pair[(a['chrom_a'], a['chrom_b'])].append(a)

    kept = []
    for group in by_pair.values():
        group.sort(key=lambda a: (-a['quality'], a['hit_a'][5] + a['hit_b'][5]))
        accepted = []
        for a in group:
            if any(abs(a['pos_a'] - k['pos_a']) < TANDEM_WINDOW and
                   abs(a['pos_b'] - k['pos_b']) < TANDEM_WINDOW for k in accepted):
                continue
            accepted.append(a)
        kept.extend(accepted)
    return kept


def longest_monotonic_run(anchors_sorted_by_a, increasing, max_dist):
    """DP longest-subsequence of pos_b monotonic in pos_a order. Consecutive
    anchors in the chain must be within max_dist (a proxy for intergenic
    distance) on BOTH axes -- bounding only the pos_b regression (as an
    earlier version of this did) lets two anchors hundreds of Mb apart on
    pos_a still "chain" as long as pos_b happens to be ordered, which isn't
    collinearity, it's coincidence. Real synteny chaining (DAGchainer,
    MCScanX) always bounds the distance on both sides; this is that same
    idea, simplified. O(n^2); fine per chromosome-pair group, not
    genome-wide."""
    n = len(anchors_sorted_by_a)
    if n == 0:
        return []
    key = (lambda a: a['pos_b']) if increasing else (lambda a: -a['pos_b'])
    dp = [1] * n
    parent = [-1] * n
    for i in range(n):
        for j in range(i):
            a_dist = anchors_sorted_by_a[i]['pos_a'] - anchors_sorted_by_a[j]['pos_a']
            b_dist = key(anchors_sorted_by_a[i]) - key(anchors_sorted_by_a[j])
            if 0 <= a_dist <= max_dist and -max_dist <= b_dist <= max_dist and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                parent[i] = j
    best_i = max(range(n), key=lambda i: dp[i])
    chain, i = [], best_i
    while i != -1:
        chain.append(anchors_sorted_by_a[i])
        i = parent[i]
    return list(reversed(chain))


def find_blocks(anchors, min_block_anchors, max_dist):
    """min_block_anchors filters on DISTINCT SUPPORTING PROTEINS, not raw
    anchor-pair count -- see build_raw_anchors()'s docstring: a single
    protein with several hits on one or both sides (miniprot secondary
    alignments) contributes one raw anchor per hit combination, so a chain
    can rack up many anchors from just a handful of genes whenever one of
    them has a locally tandem-duplicated relative (a real, common case: a
    small gene family expanded in one lineage, all copies clustered within
    a few kb, each hitting the same single ancestral locus in the other
    genome). That inflates raw anchor count without adding independent
    evidence for collinearity -- a block should mean "N different genes
    agree on this order," not "one gene's repeat-driven secondary hits
    happened to land N times." Counting distinct proteins instead catches
    this without disturbing genuine multi-copy signal (e.g. a polyploid's
    duplicated subgenomes): each subgenome copy still forms its own
    separate chain/block, each backed by many distinct genes, since the
    duplication there is genome-wide and collinear, not one gene's local
    tandem cluster."""
    by_pair = defaultdict(list)
    for a in anchors:
        by_pair[(a['chrom_a'], a['chrom_b'])].append(a)

    blocks = []
    for (chrom_a, chrom_b), group in by_pair.items():
        remaining = sorted(group, key=lambda a: a['pos_a'])
        while True:
            fwd = longest_monotonic_run(remaining, True, max_dist)
            rev = longest_monotonic_run(remaining, False, max_dist)
            chosen, orientation = (fwd, '+') if len(fwd) >= len(rev) else (rev, '-')
            # a lower bound on distinct proteins, cheap to check before the
            # (slightly more expensive) set() below -- if there aren't even
            # this many raw anchors, there certainly aren't this many
            # distinct proteins among them
            if len(chosen) < min_block_anchors:
                break
            if len({a['protein'] for a in chosen}) < min_block_anchors:
                break
            blocks.append((chrom_a, chrom_b, chosen, orientation))
            chosen_ids = {id(a) for a in chosen}
            remaining = [a for a in remaining if id(a) not in chosen_ids]
            if len(remaining) < min_block_anchors:
                break
    return blocks


def write_links(blocks, out_path):
    with open(out_path, 'w') as f:
        f.write("query_chrom\tquery_start\tquery_end\t"
                "subject_chrom\tsubject_start\tsubject_end\tscore\torientation\t"
                "mean_identity\tanchor_density\n")
        for chrom_a, chrom_b, anchors, orientation in blocks:
            starts_a = [a['hit_a'][1] for a in anchors]
            ends_a = [a['hit_a'][2] for a in anchors]
            starts_b = [a['hit_b'][1] for a in anchors]
            ends_b = [a['hit_b'][2] for a in anchors]
            q_start, q_end = min(starts_a), max(ends_a)
            mean_identity = sum(a['quality'] for a in anchors) / len(anchors)
            # distinct proteins, not raw anchor count -- see find_blocks()
            distinct_proteins = len({a['protein'] for a in anchors})
            span_mb = max(q_end - q_start, 1) / 1_000_000
            anchor_density = distinct_proteins / span_mb
            f.write(f"{chrom_a}\t{q_start}\t{q_end}\t"
                    f"{chrom_b}\t{min(starts_b)}\t{max(ends_b)}\t"
                    f"{distinct_proteins}\t{orientation}\t"
                    f"{mean_identity:.4f}\t{anchor_density:.2f}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--query_gff', required=True)
    parser.add_argument('--subject_gff', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--min_identity', type=float, default=None,
                         help='drop anchors below this per-hit identity (0-1). Default: '
                              'auto-tuned to the observed mean best-hit identity of the '
                              'weaker genome, clamped to [0.3, 0.9] (see auto_min_identity()) '
                              'instead of a fixed 0.9, so a divergent species pair is not '
                              'forced through a threshold tuned for close relatives -- pass '
                              'a value explicitly to disable auto-tuning and pin it yourself')
    parser.add_argument('--min_block_anchors', type=int, default=None,
                         help='drop chains shorter than this many anchors. Default: '
                              'auto-tuned from the same weaker-genome mean identity '
                              '--min_identity uses (see auto_min_block_anchors()) -- 15 for '
                              'a close-relative pair (mean identity >= 0.8), 5 for a '
                              'divergent one, since divergent pairs carry real signal in '
                              'many short locally-collinear chains rather than long ones. '
                              'Lower values surface more, shorter blocks at the cost of '
                              'more short chains that are coincidental rather than real; '
                              'cross-reference against how many other blocks support the '
                              'same chromosome pair before trusting any single short one. '
                              'Pass a value explicitly to disable auto-tuning and pin it '
                              'yourself')
    parser.add_argument('--max_dist', type=int, default=300000,
                         help='max bp between consecutive anchors in a chain, on both '
                              'sides -- a proxy for intergenic distance, bounding '
                              'collinearity locally so anchors that are individually real '
                              'but scattered across huge distances '
                              "can't chain together into one spurious genome-spanning block")
    parser.add_argument('--self', action='store_true',
                         help='query and subject are the same genome (homeolog scan)')
    args = parser.parse_args()

    query_hits = parse_gff(args.query_gff)
    subject_hits = parse_gff(args.subject_gff) if not args.self else query_hits

    # both auto-tunes key off the same weaker-genome signal, computed once
    # and shared rather than each function reparsing/recomputing it
    query_mean = mean_best_identity(query_hits)
    subject_mean = mean_best_identity(subject_hits)
    weaker_mean_identity = min(query_mean, subject_mean)

    if args.min_identity is None:
        min_identity = auto_min_identity(weaker_mean_identity)
        print(f"[build_synteny_blocks] auto-tuned --min_identity to {min_identity:.2f} "
              f"(mean best-hit identity: query={query_mean:.2f}, subject={subject_mean:.2f})",
              file=sys.stderr)
    else:
        min_identity = args.min_identity

    if args.min_block_anchors is None:
        min_block_anchors = auto_min_block_anchors(weaker_mean_identity)
        print(f"[build_synteny_blocks] auto-tuned --min_block_anchors to {min_block_anchors} "
              f"(weaker-genome mean identity: {weaker_mean_identity:.2f})", file=sys.stderr)
    else:
        min_block_anchors = args.min_block_anchors

    raw = build_raw_anchors(query_hits, subject_hits, min_identity)
    collapsed = collapse_tandem(raw)
    blocks = find_blocks(collapsed, min_block_anchors, args.max_dist)
    write_links(blocks, args.out)

    print(f"[build_synteny_blocks] {len(raw)} raw anchor(s) -> "
          f"{len(collapsed)} after tandem collapse -> {len(blocks)} block(s) "
          f"written to {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
