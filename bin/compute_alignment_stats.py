#!/usr/bin/env python3
"""Summarize the proteome-vs-genome miniprot alignments for the interactive
HTML's stats panel (see plot_synteny_interactive.py --stats).

"identity" here is miniprot's Positive= (identical or positively-scoring
residues), the metric --min_identity filters on, so the panel's number is
the one a viewer compares the page's Min identity control against.
Positive= also stays informative for real orthologs long after the stricter
Identity= would make diverged ones look like noise.
"""
import argparse
import gzip
import json
import re
import sys

MRNA_ATTR_RE = re.compile(r'(\w+)=([^;]+)')


def count_proteome(path):
    """Number of '>' lines, plain or gzipped (detected by magic bytes, so
    any file name works, as for miniprot)."""
    with open(path, 'rb') as f:
        gzipped = f.read(2) == b'\x1f\x8b'
    n = 0
    with (gzip.open if gzipped else open)(path, 'rb') as f:
        for line in f:
            if line.startswith(b'>'):
                n += 1
    return n


def gff_stats(path):
    """(n_aligned_proteins, mean identity of each aligned protein's best hit)."""
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
            target = attrs.get('Target', '').split()
            if not target:
                continue
            target = target[0]
            identity = float(attrs.get('Positive', attrs.get('Identity', 0)))
            rank = int(attrs.get('Rank', 1))
            if rank == 1:
                rank1_identity[target] = identity
            if target not in max_identity or identity > max_identity[target]:
                max_identity[target] = identity

    # a protein without a rank-1 hit falls back to its best hit
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
    # How the stats panel labels each input; empty means unknown and shows nothing.
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
