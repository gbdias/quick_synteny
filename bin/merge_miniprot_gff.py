#!/usr/bin/env python3
"""Merge per-chunk miniprot GFFs (--miniprot_chunk_gb) back into one
GFF equivalent to a single whole-genome miniprot run.

Why this is exact: for each query
protein, miniprot keeps alignments with score >= --outs * bestScore, up to
-N secondaries. A hit in the GLOBAL top N+1 is also in its own chunk's top
N+1 (a chunk's best score for that protein is <= the whole genome's best),
so the per-chunk outputs -- each already filtered to ITS OWN local best --
are together a superset of the whole-genome output. Re-applying both
filters (--outs and -N) here, using the score across ALL chunks as "best",
recovers the whole-genome result exactly.

Per input GFF, a record is one mRNA line plus its child lines (CDS,
stop_codon, ... -- anything with Parent=<the mRNA's ID>), read as a
contiguous block the way miniprot itself emits them. ##PAF lines are
dropped (nothing downstream reads them) and so is each chunk's own
##gff-version header (one is written once for the merged output).

IDs collide across chunks (miniprot numbers them MP000001, MP000002, ...
from 1 within EACH chunk), so every ID=/Parent= is prefixed with
"c<chunk index>_" first, chunk index being the input files' position after
sorting by filename (chunk_000.*, chunk_001.*, ... sort in chunk order).

Per protein (first token of a record's Target= attribute): best = the
highest mRNA score (GFF column 6) across every chunk's record for that
protein. Records with score >= --outs * best survive, ranked by score
descending (ties: lower chunk index, then the record's position within
that chunk's file), truncated to 1 + --max_secondary. Rank= is then
rewritten 1..k, on the mRNA line and every one of its children, to match.

Output: one "##gff-version 3" line, then every surviving record (mRNA line
followed by its children, in their original relative order), sorted by
(protein accession, new rank).
"""
import argparse
import os
import sys
from collections import defaultdict


def parse_attrs(attr_str):
    """'ID=x;Rank=1;Target=y 1 2' -> [['ID', 'x'], ['Rank', '1'], ...],
    order preserved (mutated in place by set_attr, then rejoined)."""
    pairs = []
    for token in attr_str.split(';'):
        if not token:
            continue
        key, _sep, value = token.partition('=')
        pairs.append([key, value])
    return pairs


def attrs_to_str(pairs):
    return ';'.join(f"{k}={v}" for k, v in pairs)


def get_attr(pairs, key):
    for k, v in pairs:
        if k == key:
            return v
    return None


def set_attr(pairs, key, value):
    for pair in pairs:
        if pair[0] == key:
            pair[1] = value
            return
    pairs.append([key, value])


def parse_gff(path):
    """-> [record, ...] in file order. record: {'mrna': (fields, pairs),
    'children': [(fields, pairs), ...], 'protein': str, 'score': float,
    'orig_id': str (the ID= value BEFORE chunk-prefixing, used only to
    match up children)}."""
    records = []
    current = None
    with open(path) as f:
        for lineno, line in enumerate(f, start=1):
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) != 9:
                sys.exit(f"ERROR: {path}:{lineno}: expected 9 GFF columns, got {len(fields)}: {line}")
            pairs = parse_attrs(fields[8])
            kind = fields[2]
            if kind == 'mRNA':
                orig_id = get_attr(pairs, 'ID')
                target = get_attr(pairs, 'Target')
                protein = target.split()[0] if target else None
                if orig_id is None or protein is None:
                    sys.exit(f"ERROR: {path}:{lineno}: mRNA line missing ID= or Target=: {line}")
                current = {
                    'mrna': (fields, pairs),
                    'children': [],
                    'protein': protein,
                    'orig_id': orig_id,
                    'score': float(fields[5]),
                }
                records.append(current)
            else:
                if current is None:
                    sys.exit(f"ERROR: {path}:{lineno}: '{kind}' line before any mRNA: {line}")
                parent = get_attr(pairs, 'Parent')
                if parent != current['orig_id']:
                    sys.exit(f"ERROR: {path}:{lineno}: '{kind}' Parent={parent} doesn't match "
                             f"the preceding mRNA's ID={current['orig_id']}: {line}")
                current['children'].append((fields, pairs))
    return records


def prefix_ids(records, chunk_idx):
    """Rewrite ID=/Parent= in place so they're unique across chunks --
    chunk index becomes part of every record's identity from here on."""
    for rec in records:
        _fields, pairs = rec['mrna']
        set_attr(pairs, 'ID', f"c{chunk_idx}_{rec['orig_id']}")
        for _cfields, cpairs in rec['children']:
            set_attr(cpairs, 'Parent', f"c{chunk_idx}_{rec['orig_id']}")


def merge(gff_paths, outs, max_secondary):
    """-> (n_input, [record, ...]) -- n_input is the total candidate count
    across all chunks before filtering; the records are the survivors, each
    with 'new_rank' set, NOT yet sorted for output."""
    all_records = []
    # sort by filename so "c<chunk index>_" and the tie-break "lower chunk
    # index" both mean the same thing chunk_NNN.* naming implies
    for chunk_idx, path in enumerate(sorted(gff_paths, key=os.path.basename)):
        records = parse_gff(path)
        prefix_ids(records, chunk_idx)
        for orig_order, rec in enumerate(records):
            rec['chunk_idx'] = chunk_idx
            rec['orig_order'] = orig_order
        all_records.extend(records)

    by_protein = defaultdict(list)
    for rec in all_records:
        by_protein[rec['protein']].append(rec)

    survivors = []
    for group in by_protein.values():
        best = max(rec['score'] for rec in group)
        keep = [rec for rec in group if rec['score'] >= outs * best]
        keep.sort(key=lambda rec: (-rec['score'], rec['chunk_idx'], rec['orig_order']))
        keep = keep[:1 + max_secondary]
        for new_rank, rec in enumerate(keep, start=1):
            rec['new_rank'] = new_rank
            _fields, pairs = rec['mrna']
            set_attr(pairs, 'Rank', str(new_rank))
            for _cfields, cpairs in rec['children']:
                set_attr(cpairs, 'Rank', str(new_rank))
        survivors.extend(keep)
    return len(all_records), survivors


def write_gff(records, out_path):
    records = sorted(records, key=lambda rec: (rec['protein'], rec['new_rank']))
    with open(out_path, 'w') as f:
        f.write("##gff-version 3\n")
        for rec in records:
            fields, pairs = rec['mrna']
            fields[8] = attrs_to_str(pairs)
            f.write('\t'.join(fields) + '\n')
            for cfields, cpairs in rec['children']:
                cfields[8] = attrs_to_str(cpairs)
                f.write('\t'.join(cfields) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gff', required=True, nargs='+', help='per-chunk miniprot --gff outputs')
    parser.add_argument('--out', required=True)
    parser.add_argument('--outs', type=float, default=0.7,
                         help='keep alignments scoring >= outs * that protein\'s best (default: 0.7, '
                              'matching MINIPROT_ALIGN\'s --outs)')
    parser.add_argument('--max_secondary', type=int, default=5,
                         help='keep at most 1 + this many alignments per protein (default: 5, '
                              'matching MINIPROT_ALIGN\'s -N)')
    args = parser.parse_args()

    n_input, survivors = merge(args.gff, args.outs, args.max_secondary)
    write_gff(survivors, args.out)

    n_proteins = len({rec['protein'] for rec in survivors})
    print(f"[merge_miniprot_gff] {len(args.gff)} chunk(s), {n_input} candidate alignment(s) -> "
          f"{len(survivors)} kept for {n_proteins} protein(s) -> {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
