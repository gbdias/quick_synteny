#!/usr/bin/env python3
"""--reference_accession's assembly as a reference_selection.json.

Usage: describe_assembly.py --summary summary.jsonl --accession GCA_... --out reference_selection.json

summary.jsonl is `datasets summary genome accession <acc> --as-json-lines`.
The selection has the fields find_closest_assembly.py writes for a searched
reference, so main.nf, the download and the page treat both alike; `chosen`
says it was picked by accession, and `rank`/`name` say so too ("assembly",
the accession), since no taxon was climbed to.
"""
import argparse
import json
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', required=True)
    parser.add_argument('--accession', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    with open(args.summary) as f:
        records = [json.loads(line) for line in f if line.strip()]
    # an accession NCBI doesn't have isn't an error to datasets: no records
    if not records:
        sys.exit(f"ERROR: --reference_accession {args.accession}: NCBI has no assembly with that accession")
    r = records[0]
    organism = r.get('organism', {})
    selection = {
        'accession': r.get('accession', args.accession),
        'rank': 'assembly',
        'taxid': '',
        'name': r.get('accession', args.accession),
        'organism_name': organism.get('organism_name', ''),
        'organism_taxid': str(organism.get('tax_id', '')),
        'assembly_level': r.get('assembly_info', {}).get('assembly_level', ''),
        'preferred_species': False,
        'candidate_count': 1,
        'annotated': 'annotation_info' in r,
        'same_species_as_target': False,
        'chosen': 'accession',
    }
    with open(args.out, 'w') as f:
        json.dump(selection, f, indent=2)


if __name__ == '__main__':
    main()
