#!/usr/bin/env python3
"""Parse `datasets summary taxonomy taxon <taxid>` JSON into a lineage TSV.

Usage: parse_lineage.py taxonomy.json > lineage.tsv

Output columns: rank, taxid, name -- ordered species -> genus -> family ->
order -> class -> phylum (ranks absent from NCBI's classification for this
taxon are skipped, not padded).
"""
import json
import sys

RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum']


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} taxonomy.json")

    with open(sys.argv[1]) as f:
        data = json.load(f)

    reports = data.get('reports') or []
    if not reports:
        sys.exit(f"ERROR: no taxonomy report returned (input: {sys.argv[1]})")

    classification = reports[0].get('taxonomy', {}).get('classification', {})
    if not classification:
        sys.exit("ERROR: taxonomy report has no 'classification' block -- "
                  "the datasets CLI JSON schema may have changed")

    rows = []
    for rank in RANKS:
        entry = classification.get(rank)
        if entry is None:
            continue
        rows.append((rank, str(entry['id']), entry['name']))

    if not rows:
        sys.exit("ERROR: none of the expected ranks "
                  f"({','.join(RANKS)}) were found in the classification block")

    for rank, taxid, name in rows:
        print(f"{rank}\t{taxid}\t{name}")


if __name__ == '__main__':
    main()
