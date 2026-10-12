#!/usr/bin/env python3
"""Rewrite a GFF's column 1 (seqid) from original sequence IDs to the
renamed IDs rename_sequences.py assigned, streaming. Lines without a tab
(comments, directives) pass through; every other line is byte-identical
apart from its seqid.

miniprot aligns against the ORIGINAL genome fasta, so this is where its
seqids become the chr1/chr2/... labels everything downstream expects.
"""
import argparse
import sys


def load_lookup(path):
    """rename_sequences.py --out_lookup: original_id\\tnew_id\\tlength\\ttype,
    one header line."""
    lookup = {}
    with open(path) as f:
        next(f, None)
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
