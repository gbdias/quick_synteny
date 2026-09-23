#!/usr/bin/env python3
"""Rewrite a GFF's column 1 (seqid) from original sequence IDs to the short,
renamed IDs rename_sequences.py assigned -- streaming, so a multi-GB GFF
never needs to fit in memory. Comment/directive lines (and any other line
with no tab in it) pass through untouched; every other line is byte-
identical to the input except for its seqid field, so source, type,
coordinates, score, strand, phase and attributes can't drift even by
whitespace.

miniprot now runs on the ORIGINAL genome fasta (see miniprot_align.nf), not
a renamed copy, so its GFF still carries the genome's original sequence
names. This is the one place those get swapped for the plot-friendly
chr1/chr2/... labels everything downstream (build_synteny.nf, the
interactive plot) expects.
"""
import argparse
import sys


def load_lookup(path):
    """rename_sequences.py --out_lookup: original_id\\tnew_id\\tlength\\ttype,
    one header line."""
    lookup = {}
    with open(path) as f:
        next(f, None)  # header
        for line in f:
            if not line.strip():
                continue
            original_id, new_id = line.rstrip('\n').split('\t')[:2]
            lookup[original_id] = new_id
    return lookup


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gff', required=True)
    parser.add_argument('--lookup', required=True, help='rename_sequences.py --out_lookup output')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    lookup = load_lookup(args.lookup)

    with open(args.gff, newline='') as fin, open(args.out, 'w', newline='') as fout:
        for line in fin:
            if line.startswith('#') or '\t' not in line:
                fout.write(line)
                continue
            seqid, rest = line.split('\t', 1)
            if seqid not in lookup:
                sys.exit(f"ERROR: sequence ID '{seqid}' in {args.gff} has no entry in {args.lookup}")
            fout.write(lookup[seqid])
            fout.write('\t')
            fout.write(rest)


if __name__ == '__main__':
    main()
