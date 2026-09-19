#!/usr/bin/env python3
"""Pick the closest available genome assembly from pre-fetched ladder results.

The datasets CLI and this script deliberately run in separate containers (the
biocontainers ncbi-datasets-cli image has no Python), so the ladder-walk is
split in two: a process using the datasets container queries every rank up to
--max_rank and writes one JSON-lines file per rank (<rank>.jsonl, empty if no
hits); this script then picks the first rank (species-first order) with a
non-empty file and ranks its candidates by assembly quality.

Writes, in the current directory:
  <outprefix>_selection.json   chosen accession + rank + candidate count
  <outprefix>_candidates.tsv   every candidate considered at the winning rank,
                                best first
  <outprefix>_search_log.tsv   rank, taxid, name, hit count for every rank
                                actually searched (including zero-hit ranks)

Exits non-zero with a clear message if no rank up to --max_rank has a hit.
"""
import argparse
import json
import os
import sys

RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum']

LEVEL_RANK = {
    'chromosome': 3,
    'complete genome': 3,
    'scaffold': 1,
    'contig': 0,
}


def level_score(record):
    level = record.get('assembly_info', {}).get('assembly_level', '').lower()
    return LEVEL_RANK.get(level, -1)


def sort_key(record):
    stats = record.get('assembly_stats', {})
    return (
        level_score(record),
        stats.get('contig_n50', 0) or 0,
        stats.get('scaffold_n50', 0) or 0,
    )


def read_lineage(path):
    lineage = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            rank, taxid, name = line.split('\t', 2)
            lineage[rank] = (taxid, name)
    return lineage


def read_jsonl(path):
    if not os.path.exists(path):
        return []
    records = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lineage', required=True)
    parser.add_argument('--max_rank', required=True, choices=RANKS)
    parser.add_argument('--jsonl_dir', default='.',
                         help='directory containing <rank>.jsonl files from the ladder query')
    parser.add_argument('--outprefix', required=True)
    args = parser.parse_args()

    lineage = read_lineage(args.lineage)
    ladder = RANKS[:RANKS.index(args.max_rank) + 1]

    search_log = []
    chosen_rank = None
    chosen_taxid = None
    chosen_name = None
    candidates = []

    for rank in ladder:
        if rank not in lineage:
            continue
        taxid, name = lineage[rank]
        records = read_jsonl(os.path.join(args.jsonl_dir, f"{rank}.jsonl"))
        search_log.append((rank, taxid, name, len(records)))
        if records and not candidates:
            chosen_rank, chosen_taxid, chosen_name = rank, taxid, name
            candidates = sorted(records, key=sort_key, reverse=True)

    with open(f"{args.outprefix}_search_log.tsv", 'w') as f:
        f.write("rank\ttaxid\tname\thits\n")
        for rank, taxid, name, hits in search_log:
            f.write(f"{rank}\t{taxid}\t{name}\t{hits}\n")

    if not candidates:
        sys.exit(
            f"ERROR: no genome assembly found for {args.outprefix} up to "
            f"rank '{args.max_rank}' (searched: "
            f"{','.join(r for r, *_ in search_log)}). "
            f"Re-run with a higher --max_rank to climb further."
        )

    best = candidates[0]
    with open(f"{args.outprefix}_candidates.tsv", 'w') as f:
        f.write("accession\tassembly_level\tcontig_n50\tscaffold_n50\tannotated\n")
        for c in candidates:
            stats = c.get('assembly_stats', {})
            f.write(
                f"{c['accession']}\t"
                f"{c.get('assembly_info', {}).get('assembly_level', '')}\t"
                f"{stats.get('contig_n50', '')}\t"
                f"{stats.get('scaffold_n50', '')}\t"
                f"{'annotation_info' in c}\n"
            )

    selection = {
        'accession': best['accession'],
        'rank': chosen_rank,
        'taxid': chosen_taxid,
        'name': chosen_name,
        'candidate_count': len(candidates),
        'annotated': 'annotation_info' in best,
    }
    with open(f"{args.outprefix}_selection.json", 'w') as f:
        json.dump(selection, f, indent=2)

    print(f"[find_closest_assembly] chosen accession={best['accession']} "
          f"rank={chosen_rank} ({chosen_name}) "
          f"from {len(candidates)} candidate(s)", file=sys.stderr)


if __name__ == '__main__':
    main()
