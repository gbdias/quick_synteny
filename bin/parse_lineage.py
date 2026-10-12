#!/usr/bin/env python3
"""Parse `datasets summary taxonomy taxon <taxid> --parents` JSON into a lineage TSV.

Usage: parse_lineage.py taxonomy.json > lineage.tsv

Output columns: rank, taxid, name -- one row per taxon from species outward
to phylum, nearest first, INCLUDING every intermediate taxon (tribe,
superfamily, unranked clades, ...): e.g. for Trox cadaverinus, skipping
superfamily Scarabaeoidea would jump from a family with no other assemblies
straight to all of Coleoptera. Rank is NCBI's, lowercased; 'clade' and 'no
rank' can repeat, so key on taxid. Taxa below species and above phylum are
dropped.

Without --parents, falls back to the report's classification block (major
ranks only).
"""
import json
import sys

# the major ranks --max_rank chooses from; the lineage ends at the outermost
# of these present (phylum, when NCBI has one)
RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum']


def full_lineage(reports, query):
    """Every taxon from `query` outward, nearest first, from the --parents
    reports; None when they're missing."""
    by_id = {r['taxonomy']['tax_id']: r['taxonomy'] for r in reports if 'taxonomy' in r}
    # the root (taxid 1) itself comes back without a report of its own
    parents = [p for p in query.get('parents') or [] if p in by_id]
    if not parents:
        return None
    rows = []
    for t in [query] + [by_id[p] for p in reversed(parents)]:
        rows.append(((t.get('rank') or 'no rank').lower().replace('_', ' '),
                     str(t['tax_id']),
                     t.get('current_scientific_name', {}).get('name', '')))
    return rows


def classification_lineage(query):
    classification = query.get('classification', {})
    if not classification:
        sys.exit("ERROR: taxonomy report has no 'classification' block -- "
                  "the datasets CLI JSON schema may have changed")
    return [(rank, str(classification[rank]['id']), classification[rank]['name'])
            for rank in RANKS if rank in classification]


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} taxonomy.json")

    with open(sys.argv[1]) as f:
        data = json.load(f)

    reports = data.get('reports') or []
    if not reports:
        sys.exit(f"ERROR: no taxonomy report returned (input: {sys.argv[1]})")
    # an unknown taxid isn't a datasets failure: it exits 0 with a report
    # carrying only `errors` (reason: "... is not a recognized NCBI Taxonomy
    # name.") and no `taxonomy` block
    if not any('taxonomy' in r for r in reports):
        unknown = sorted({i for r in reports for e in r.get('errors') or [] for i in e.get('invalid_identifiers') or []})
        reasons = sorted({e['reason'] for r in reports for e in r.get('errors') or [] if e.get('reason')})
        if unknown:
            sys.exit(f"ERROR: --taxid {', '.join(unknown)} is not a known NCBI taxid")
        sys.exit(f"ERROR: NCBI returned no taxonomy for this taxid: {'; '.join(reasons) or 'no reason given'}")

    # --parents returns the ancestors (and immediate children) in no
    # particular order, every report echoing the same `query`: the queried
    # taxon is the one whose own taxid is that query
    asked = {str(q) for r in reports for q in r.get('query') or []}
    query = next((r['taxonomy'] for r in reports
                  if str(r.get('taxonomy', {}).get('tax_id')) in asked),
                 reports[0]['taxonomy'])

    rows = full_lineage(reports, query) or classification_lineage(query)

    # species (not a subspecies query) through the outermost major rank
    ranks = [rank for rank, _, _ in rows]
    if 'species' in ranks:
        rows = rows[ranks.index('species'):]
    major = [i for i, (rank, _, _) in enumerate(rows) if rank in RANKS]
    if not major:
        sys.exit("ERROR: none of the expected ranks "
                  f"({','.join(RANKS)}) were found in the lineage")
    rows = rows[:major[-1] + 1]

    for rank, taxid, name in rows:
        print(f"{rank}\t{taxid}\t{name}")


if __name__ == '__main__':
    main()
