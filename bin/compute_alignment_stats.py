#!/usr/bin/env python3
"""Summarize the proteome-vs-genome miniprot alignments for the interactive
HTML's stats panel (see plot_synteny_interactive.py --stats).

Reads the same miniprot GFF Target=/Positive=/Rank= attributes
extract_hits.py puts in the hit table, but only needs per-protein counts and
each protein's best-hit score. Uses the same metric as the chainer on
purpose: Positive= (identical or positively-scoring residues) is what the
auto-tuned --min_identity/--min_block values in bin/chain.js are computed
from, so the stats panel shows the same "identity" number that drives
those defaults.
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

    "identity" here is miniprot's Positive= score, not its stricter Identity=:
    conservative substitutions accumulate long before radical ones as two
    genomes diverge, so Positive stays informative for real orthologs well
    after Identity alone would make them look like noise."""
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
    # How the page's stats panel names each input: the file name if the user
    # supplied it, else the NCBI accession it was auto-discovered from; plus
    # its species when known (from --taxid's lineage, or the discovered
    # assembly's own record). Empty means unknown -- nothing is shown for it.
    for role in ('query', 'subject', 'proteome'):
        parser.add_argument(f'--{role}_source', default='',
                             help=f'{role}: input file name if user-supplied, else its NCBI accession')
        parser.add_argument(f'--{role}_species', default='',
                             help=f'{role}: species name when known, else empty')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    proteome_total = count_proteome(args.proteome)
    query_aligned, query_identity = gff_stats(args.query_gff)
    subject_aligned, subject_identity = gff_stats(args.subject_gff)

    stats = {
        'proteome_total': proteome_total,
        'proteome_source': args.proteome_source or None,
        'proteome_species': args.proteome_species or None,
        'query_name': args.query_name,
        'query_source': args.query_source or None,
        'query_species': args.query_species or None,
        'query_aligned': query_aligned,
        'query_mean_identity': query_identity,
        'subject_name': args.subject_name,
        'subject_source': args.subject_source or None,
        'subject_species': args.subject_species or None,
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
