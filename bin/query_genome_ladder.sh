#!/usr/bin/env bash
# Query NCBI for genome assemblies at each taxon of a lineage, nearest first.
#
# Usage: query_genome_ladder.sh <lineage.tsv> <max_rank> <prefer_taxid|''> [datasets filter flags...]
#
# Walks lineage.tsv (parse_lineage.py) up to <max_rank>, running
# `datasets summary genome taxon <taxid> <filters>` for each, and stops after
# the first taxon with an assembly of a species other than the target's:
# every taxon above is a superset, and find_closest_assembly.py picks from
# the lowest taxon with a usable hit. A same-species hit alone never stops
# the climb, since whether it's usable is decided in the selection step.
#
# At most LADDER_CAP (default 1000) records are kept per taxon. A larger
# taxon is re-queried with --reference (one designated genome per species)
# rather than keeping an arbitrary first LADDER_CAP; only if that still
# exceeds the cap is it truncated, as query_log.tsv records. At ~5 kB of
# JSON per record, the cap also bounds the selection step's memory.
#
# <prefer_taxid>, when given, is also queried on its own (preferred.jsonl):
# that species may join the lineage above where the climb stops.
#
# Writes into ./ladder/: <taxid>.jsonl per taxon queried, preferred.jsonl,
# target_taxids.txt (the target species and every taxon under it), and
# query_log.tsv (rank, taxid, name, records, query), in query order.
set -u

lineage=$1
max_rank=$2
prefer=$3
shift 3
filters=("$@")
cap=${LADDER_CAP:-1000}

mkdir -p ladder
target=$(awk -F'\t' '$1=="species"{print $2; exit}' "$lineage")

# the target species and every strain/subspecies under it: assemblies are
# often filed under a strain's taxid (S. cerevisiae S288C is 559292, not 4932)
: > ladder/target_taxids.txt
if [ -n "$target" ]; then
    echo "$target" > ladder/target_taxids.txt
    datasets summary taxonomy taxon "$target" --children --report ids_only --as-json-lines \
        | grep -o '"tax_id":[0-9]*' | cut -d: -f2 >> ladder/target_taxids.txt || true
fi

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
        others=$(dataformat tsv genome --inputfile "$out" --fields organism-tax-id --elide-header | grep -cvxF -f ladder/target_taxids.txt)
    fi
    if [ "$others" -gt 0 ]; then
        break
    fi
done < ladder/ranks.tsv
rm ladder/ranks.tsv

if [ -n "$prefer" ]; then
    query "$prefer" ladder/preferred.jsonl
fi
