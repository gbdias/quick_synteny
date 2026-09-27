#!/usr/bin/env python3
"""Pick the closest available genome assembly from pre-fetched ladder results.

The datasets CLI and this script deliberately run in separate containers (the
datasets image, envs/ncbi_datasets.yml, has no Python), so the ladder-walk is
split in two: query_genome_ladder.sh, in the datasets container, queries the
lineage's taxa from species outward (intermediate ranks and clades included)
and writes one JSON-lines file per taxon (<taxid>.jsonl, empty if no hits)
plus query_log.tsv listing the taxa it queried, in order; this script then
picks the first of those taxa with a usable hit and ranks its candidates by
assembly quality.

Writes, in the current directory:
  <outprefix>_selection.json   chosen accession + rank + candidate count
  <outprefix>_candidates.tsv   every candidate considered at the winning rank,
                                best first
  <outprefix>_search_log.tsv   rank, taxid, name, usable hit count and how the
                                query was cut (see query_genome_ladder.sh) for
                                every taxon actually queried (including zero-
                                hit ones)

Exits non-zero with a clear message if no taxon up to --max_rank has a hit.
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
    """(rank, taxid, name) rows, nearest first -- rank isn't unique ('clade')."""
    lineage = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            rank, taxid, name = line.split('\t', 2)
            lineage.append((rank, taxid, name))
    return lineage


def read_query_log(path):
    """(rank, taxid, name, query mode) for every taxon queried, in order."""
    with open(path) as f:
        next(f)
        return [tuple(line.rstrip('\n').split('\t')[i] for i in (0, 1, 2, 4))
                for line in f if line.strip()]


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


def is_chromosome_level(record):
    return record.get('assembly_info', {}).get('assembly_level', '').lower() in \
        ('chromosome', 'complete genome')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lineage', required=True)
    parser.add_argument('--max_rank', required=True, choices=RANKS)
    parser.add_argument('--ladder_dir', default='ladder',
                         help='directory query_genome_ladder.sh wrote <taxid>.jsonl, '
                              'preferred.jsonl and query_log.tsv into')
    parser.add_argument('--outprefix', required=True)
    parser.add_argument('--exclude_target', action='store_true',
                         help="exclude every candidate of the target's own species outright, "
                              "regardless of assembly quality. Off by default: a same-species "
                              "candidate is a legitimate result (e.g. a second, independently "
                              "submitted assembly of the same organism) and is kept as long as "
                              "it's chromosome-level or better -- only a same-species candidate "
                              "below that quality bar is dropped either way, same-species or not.")
    parser.add_argument('--require_chromosome_level', action='store_true',
                         help="require EVERY candidate, same-species or not, to be chromosome- or "
                              "complete-genome-level -- used for reference-genome discovery, where "
                              "a low-quality genome makes a poor synteny comparison regardless of "
                              "species. Off for proteome discovery, which only needs a candidate's "
                              "protein sequences and is deliberately more permissive (a scaffold-"
                              "level annotated assembly is fine there). When this flag is off, a "
                              "same-species candidate is still individually held to chromosome "
                              "level -- see --exclude_target -- but a DIFFERENT-species one is not.")
    parser.add_argument('--prefer_taxid', default='',
                         help="species taxid to prefer: if it has any kept candidate "
                              "(preferred.jsonl, queried on its own), the best of those wins over "
                              "the usual ranking. Proteome discovery passes the chosen reference "
                              "genome's species, so the proteome comes from the reference species "
                              "whenever it has an annotated assembly -- the ranking below knows "
                              "assembly quality, not relatedness.")
    parser.add_argument('--prefer_rank_taxid', default='',
                         help="the lineage taxon --prefer_taxid's species joins the target's "
                              "lineage at (the reference selection's own taxid), reported as the "
                              "rank of a preferred pick")
    args = parser.parse_args()

    lineage = read_lineage(args.lineage)
    ladder = read_query_log(os.path.join(args.ladder_dir, 'query_log.tsv'))
    # 'species' in the lineage file is, by construction, the taxon the whole
    # search started from -- i.e. the target's own species -- not merely "a"
    # species. A same-species candidate is real, useful data in general (a
    # different assembly of the same organism, e.g. for validating synteny
    # within a species, or the only chromosome-level resource available for
    # a sparsely-sequenced group -- this pipeline was pointed at exactly
    # that case once already), so it's not excluded outright by default --
    # only gated to chromosome/complete level (see --require_chromosome_level
    # below for why that gate applies to a same-species candidate even when
    # it's off for everyone else). --exclude_target restores the stricter
    # behavior for a user who wants a same-species result excluded
    # regardless of quality. Applied at every rank a same-species record
    # could in principle appear at, not just species -- there's no
    # structural guarantee it can't.
    target_taxid = next((taxid for rank, taxid, _ in lineage if rank == 'species'), None)
    # ...and every taxon under it: assemblies are often filed under a
    # strain's own taxid (S. cerevisiae S288C is 559292, not 4932), which is
    # still the target's species. query_genome_ladder.sh lists them.
    target_taxids = {target_taxid} if target_taxid else set()
    taxids_file = os.path.join(args.ladder_dir, 'target_taxids.txt')
    if target_taxid and os.path.exists(taxids_file):
        with open(taxids_file) as f:
            target_taxids |= {line.strip() for line in f if line.strip()}

    def is_same_species(r):
        return str(r.get('organism', {}).get('tax_id', '')) in target_taxids

    def keep(r):
        # explicit and self-contained here, rather than relying solely on
        # find_reference_assembly.nf's own upstream --assembly-level
        # chromosome,complete query filter to make this true incidentally
        if args.require_chromosome_level and not is_chromosome_level(r):
            return False
        if not is_same_species(r):
            return True
        if args.exclude_target:
            return False
        return is_chromosome_level(r)

    search_log = []
    chosen_rank = None
    chosen_taxid = None
    chosen_name = None
    candidates = []
    kept_by_rank = []

    for rank, taxid, name, mode in ladder:
        records = read_jsonl(os.path.join(args.ladder_dir, f"{taxid}.jsonl"))
        records = [r for r in records if keep(r)]
        search_log.append((rank, taxid, name, len(records), mode))
        kept_by_rank.append((taxid, records))
        if records and not candidates:
            chosen_rank, chosen_taxid, chosen_name = rank, taxid, name
            candidates = sorted(records, key=sort_key, reverse=True)

    # a candidate of the preferred species beats the ranking above -- listed
    # first, then the rest of the rank it joins the lineage at, if queried
    preferred = False
    if args.prefer_taxid:
        match = [r for r in read_jsonl(os.path.join(args.ladder_dir, 'preferred.jsonl'))
                 if keep(r) and str(r.get('organism', {}).get('tax_id', '')) == str(args.prefer_taxid)]
        if match:
            match = sorted(match, key=sort_key, reverse=True)
            chosen_rank, chosen_taxid, chosen_name = next(
                (row for row in lineage if row[1] == str(args.prefer_rank_taxid)), ('', '', ''))
            match_ids = {r['accession'] for r in match}
            rest = next((records for taxid, records in kept_by_rank if taxid == chosen_taxid), [])
            rest = [r for r in sorted(rest, key=sort_key, reverse=True) if r['accession'] not in match_ids]
            candidates = match + rest
            preferred = True

    with open(f"{args.outprefix}_search_log.tsv", 'w') as f:
        f.write("rank\ttaxid\tname\thits\tquery\n")
        for rank, taxid, name, hits, mode in search_log:
            f.write(f"{rank}\t{taxid}\t{name}\t{hits}\t{mode}\n")

    if not candidates:
        sys.exit(
            f"ERROR: no genome assembly found for {args.outprefix} up to "
            f"rank '{args.max_rank}' (searched: "
            f"{', '.join(f'{r} {n}' for r, _, n, *_ in search_log)}). "
            f"Re-run with a higher --max_rank to climb further."
        )

    best = candidates[0]
    with open(f"{args.outprefix}_candidates.tsv", 'w') as f:
        f.write("accession\tassembly_level\tcontig_n50\tscaffold_n50\tannotated\tsame_species_as_target\n")
        for c in candidates:
            stats = c.get('assembly_stats', {})
            f.write(
                f"{c['accession']}\t"
                f"{c.get('assembly_info', {}).get('assembly_level', '')}\t"
                f"{stats.get('contig_n50', '')}\t"
                f"{stats.get('scaffold_n50', '')}\t"
                f"{'annotation_info' in c}\t"
                f"{is_same_species(c)}\n"
            )

    same_species = is_same_species(best)
    selection = {
        'accession': best['accession'],
        'rank': chosen_rank,
        'taxid': chosen_taxid,
        'name': chosen_name,
        # the chosen assembly's OWN species -- `name` above is the taxon of
        # the rank the ladder stopped at, which is a genus/family name
        # whenever the pick came from above species rank
        'organism_name': best.get('organism', {}).get('organism_name', ''),
        'organism_taxid': str(best.get('organism', {}).get('tax_id', '')),
        # true when --prefer_taxid decided the pick (see above)
        'preferred_species': preferred,
        'candidate_count': len(candidates),
        'annotated': 'annotation_info' in best,
        # true when the chosen candidate is of the target's own species --
        # a legitimate result by default (see --exclude_target above), but
        # a meaningfully different kind of comparison worth surfacing
        # rather than leaving indistinguishable from a cross-species pick
        'same_species_as_target': same_species,
    }
    with open(f"{args.outprefix}_selection.json", 'w') as f:
        json.dump(selection, f, indent=2)

    same_species_note = ' [SAME SPECIES AS TARGET]' if same_species else ''
    print(f"[find_closest_assembly] chosen accession={best['accession']} "
          f"({selection['organism_name']}) rank={chosen_rank} ({chosen_name}){same_species_note} "
          f"from {len(candidates)} candidate(s)", file=sys.stderr)


if __name__ == '__main__':
    main()
