#!/usr/bin/env python3
"""Summarize the proteome-vs-genome miniprot alignments for the interactive
HTML's stats panel (see plot_synteny_interactive.py --stats).

Reuses the same miniprot GFF Target=/Positive=/Rank= attributes
build_synteny_blocks.py parses (including its choice of Positive= over the
stricter Identity=, see that file's parse_gff() for why), but only needs
per-protein counts and each protein's best-hit score here, not a full
anchors structure -- kept as a separate small script rather than importing
that module, since these are read-only summary numbers with no chaining
involved. Using the same metric here as the chainer uses matters: this is
what the auto-tuned --min_identity/--min_block_anchors values are computed
from, so the stats panel would otherwise show a different "identity" number
than the one actually driving the guardrails.
"""
import argparse
import json
import re
import sys

MRNA_ATTR_RE = re.compile(r'(\w+)=([^;]+)')


def count_proteome(path):
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                n += 1
    return n


def gff_stats(path):
    """(n_aligned_proteins, mean identity of each aligned protein's best hit).

    "identity" here is miniprot's Positive= score, not its stricter Identity=
    -- see build_synteny_blocks.py's parse_gff() docstring."""
    rank1_identity = {}
    max_identity = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'mRNA':
                continue
            attrs = dict(MRNA_ATTR_RE.findall(fields[8]))
            target = attrs.get('Target', '').split()[0]
            if not target:
                continue
            identity = float(attrs.get('Positive', attrs.get('Identity', 0)))
            rank = int(attrs.get('Rank', 1))
            if rank == 1:
                rank1_identity[target] = identity
            if target not in max_identity or identity > max_identity[target]:
                max_identity[target] = identity

    # every aligned protein has a rank-1 hit in practice, but fall back to its
    # best-seen hit for the rare protein that somehow doesn't
    best_identity = {**max_identity, **rank1_identity}
    n = len(best_identity)
    mean_identity = sum(best_identity.values()) / n if n else 0.0
    return n, mean_identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proteome', required=True)
    parser.add_argument('--query_gff', required=True)
    parser.add_argument('--subject_gff', required=True)
    parser.add_argument('--query_name', required=True)
    parser.add_argument('--subject_name', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    proteome_total = count_proteome(args.proteome)
    query_aligned, query_identity = gff_stats(args.query_gff)
    subject_aligned, subject_identity = gff_stats(args.subject_gff)

    stats = {
        'proteome_total': proteome_total,
        'query_name': args.query_name,
        'query_aligned': query_aligned,
        'query_mean_identity': query_identity,
        'subject_name': args.subject_name,
        'subject_aligned': subject_aligned,
        'subject_mean_identity': subject_identity,
    }
    with open(args.out, 'w') as f:
        json.dump(stats, f, indent=2)

    print(f"[compute_alignment_stats] proteome={proteome_total} "
          f"{args.query_name}_aligned={query_aligned} ({query_identity * 100:.1f}% avg identity) "
          f"{args.subject_name}_aligned={subject_aligned} ({subject_identity * 100:.1f}% avg identity) "
          f"-> {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
