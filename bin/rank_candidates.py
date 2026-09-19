#!/usr/bin/env python3
"""Sort `datasets summary genome` JSON-lines records by assembly quality.

Sort key: assembly level (chromosome > complete genome > scaffold > contig),
then contig N50 descending, then scaffold N50 descending. Standalone CLI use:

    datasets summary genome taxon ... --as-json-lines | rank_candidates.py

reads JSON-lines from stdin and writes them back out, best first.
"""
import json
import sys

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


def sort_candidates(records):
    return sorted(records, key=sort_key, reverse=True)


def main():
    records = [json.loads(line) for line in sys.stdin if line.strip()]
    for record in sort_candidates(records):
        print(json.dumps(record))


if __name__ == '__main__':
    main()
