#!/usr/bin/env bash
# Query NCBI for genome assemblies at each taxon of a lineage, nearest first.
#
# Usage: query_genome_ladder.sh <lineage.tsv> <max_rank> <prefer_taxid|''> [datasets filter flags...]
#
# Walks lineage.tsv (from parse_lineage.py: every taxon from species outward,
# intermediate ranks and clades included) up to <max_rank>, running
# `datasets summary genome taxon <taxid> <filters>` for each, and stops after
# the first taxon that returns an assembly of any species other than the
# target's own. Every taxon above that one is a strict superset of it, and
# find_closest_assembly.py only ever picks from the lowest taxon with a
# usable hit, so querying further would only fetch farther relatives. The
# stop is deliberately conservative: a same-species hit alone never stops
# the climb, since whether it's usable (--exclude_target, the chromosome-
# level gate on same-species candidates) is decided in the selection step.
#
# At most LADDER_CAP (default 1000) records are kept per taxon. A taxon with
# more than that is re-queried with --reference (NCBI's designated reference
# genome, one per species) instead of keeping an arbitrary first
# LADDER_CAP of whatever order NCBI returns -- a principled one-per-species
# cut rather than a truncated sample. Only if even that exceeds the cap is it
# truncated, and query_log.tsv says so. Records average ~5 kB of JSON each,
# so the cap also bounds the selection step's memory, which is what the old
# fixed --limit 200 was for.
#
# <prefer_taxid>, when given, is also queried on its own (preferred.jsonl):
# proteome discovery prefers the chosen reference species' own annotation,
# and that species may join the lineage above where this climb stops.
#
# Writes into ./ladder/: <taxid>.jsonl per taxon queried, preferred.jsonl,
# and query_log.tsv (rank, taxid, name, records, query) listing every taxon
# queried in order -- the selection step reads the ladder from that log.
set -u

lineage=$1
max_rank=$2
prefer=$3
shift 3
filters=("$@")
cap=${LADDER_CAP:-1000}

mkdir -p ladder
target=$(awk -F'\t' '$1=="species"{print $2; exit}' "$lineage")

# query <taxid> <out.jsonl>: sets $mode to how the records were cut
query() {
    local taxid=$1 out=$2
    datasets summary genome taxon "$taxid" "${filters[@]}" --as-json-lines --limit $((cap + 1)) > "$out" || true
    mode=all
    if [ "$(grep -c . "$out")" -gt "$cap" ]; then
        datasets summary genome taxon "$taxid" "${filters[@]}" --reference --as-json-lines --limit $((cap + 1)) > "$out" || true
        mode=reference_only
        if [ "$(grep -c . "$out")" -gt "$cap" ]; then
            head -n "$cap" "$out" > "$out.tmp" && mv "$out.tmp" "$out"
            mode=reference_only_truncated
        fi
    fi
}

# the lineage rows up to max_rank -- or the first major rank above it, when
# NCBI has no taxon at max_rank itself for this lineage
awk -F'\t' -v max="$max_rank" '
    BEGIN { n = split("species genus family order class phylum", major, " ")
            for (i = 1; i <= n; i++) idx[major[i]] = i }
    { print; if (($1 in idx) && idx[$1] >= idx[max]) exit }
' "$lineage" > ladder/ranks.tsv

printf 'rank\ttaxid\tname\trecords\tquery\n' > ladder/query_log.tsv
while IFS=$'\t' read -r rank taxid name; do
    out="ladder/${taxid}.jsonl"
    query "$taxid" "$out"
    printf '%s\t%s\t%s\t%s\t%s\n' "$rank" "$taxid" "$name" "$(grep -c . "$out")" "$mode" >> ladder/query_log.tsv
    others=0
    if [ -s "$out" ]; then
        others=$(dataformat tsv genome --inputfile "$out" --fields organism-tax-id --elide-header | grep -cvx "$target")
    fi
    if [ "$others" -gt 0 ]; then
        break
    fi
done < ladder/ranks.tsv
rm ladder/ranks.tsv

if [ -n "$prefer" ]; then
    query "$prefer" ladder/preferred.jsonl
fi
