#!/usr/bin/env bash
# Preflight: verify tools, proteome and every accession in genomes.tsv before
# burning cluster time. Takes a few seconds. Run this first.
#
#   ./preflight.sh [genomes.tsv] [proteome.faa]

set -uo pipefail

TABLE=${1:-genomes.tsv}
PROTEOME=${2:-}
fail=0

say()  { printf '%-14s %s\n' "$1" "$2"; }
ok()   { say "  [ok]"   "$1"; }
bad()  { say "  [FAIL]" "$1"; fail=1; }
warn() { say "  [warn]" "$1"; }

# MINIPROT/DATASETS may be bare commands or full container invocations
# (see setup_containers.sh). Actually run them rather than checking PATH.
MINIPROT=${MINIPROT:-miniprot}
DATASETS=${DATASETS:-datasets}
read -r -a MP <<< "$MINIPROT"
read -r -a DS <<< "$DATASETS"

echo "== tools =="
for pair in "miniprot|${MP[*]}" "datasets|${DS[*]}"; do
    name=${pair%%|*}; cmd=${pair##*|}
    read -r -a arr <<< "$cmd"
    if ver=$("${arr[@]}" --version 2>&1 | head -1) && [[ -n "$ver" ]]; then
        ok "$name -> $cmd  ($ver)"
    else
        bad "$name not runnable as: $cmd"
        [[ "$cmd" == "$name" ]] && \
            echo "         Nothing is installed? Run ./setup_containers.sh first —" && \
            echo "         it pulls the same images the pipeline uses, no install needed."
    fi
done

# GNU time, not the shell builtin and not BSD time. This is what measures RSS.
if /usr/bin/time -v true >/dev/null 2>&1; then
    ok "/usr/bin/time -v (GNU time)"
else
    bad "/usr/bin/time -v unavailable — install GNU time or load a module"
fi

for t in samtools seqkit; do
    command -v "$t" >/dev/null 2>&1 && ok "$t (optional, speeds up sizing)" \
                                    || warn "$t not found (optional)"
done

command -v sbatch >/dev/null 2>&1 && ok "sbatch" || bad "sbatch not on PATH"

echo
echo "== proteome =="
if [[ -z "$PROTEOME" ]]; then
    bad "no proteome given. Pass one: ./preflight.sh $TABLE /path/to/proteome.faa"
    echo "         Any proteome works — memory is index-dominated, so its identity"
    echo "         barely matters. The script subsamples it to 2000 sequences."
elif [[ -s "$PROTEOME" ]]; then
    n=$(grep -c '^>' "$PROTEOME" 2>/dev/null || echo 0)
    if [[ "$n" -ge 2000 ]]; then ok "$PROTEOME ($n sequences)"
    elif [[ "$n" -gt 0 ]]; then warn "$PROTEOME has only $n sequences (<2000, fine but noted)"
    else bad "$PROTEOME contains no FASTA records"; fi
else
    bad "$PROTEOME not found or empty"
fi

echo
echo "== disk =="
avail=$(df -BG --output=avail . 2>/dev/null | tail -1 | tr -dc '0-9')
if [[ -n "$avail" ]]; then
    # axolotl alone is ~32 Gb of sequence; gz + unpacked needs headroom
    if [[ "$avail" -ge 250 ]]; then ok "${avail}G available here"
    else warn "${avail}G available — the big genomes need ~250G unpacked. Run from scratch/project storage."; fi
fi

echo
echo "== accessions =="
[[ -r "$TABLE" ]] || { bad "cannot read $TABLE"; exit 1; }

while IFS=$'\t' read -r label source approx mem wall; do
    [[ -z "${label:-}" || "$label" == \#* ]] && continue

    if [[ "$source" == /* ]]; then
        [[ -s "$source" ]] && ok "$label -> local file $source" \
                           || bad "$label -> local file missing: $source"
        continue
    fi

    # `datasets summary` is cheap and tells us the real assembly size, which is
    # a better sanity check than trusting the approx_gb column.
    json=$("${DS[@]}" summary genome accession "$source" --as-json-lines 2>/dev/null)
    if [[ -z "$json" ]]; then
        bad "$label -> $source NOT FOUND on NCBI (replace this accession)"
        continue
    fi

    real_gb=$(printf '%s' "$json" \
        | grep -o '"total_sequence_length":"[0-9]*"' | head -1 \
        | grep -o '[0-9]*' \
        | awk '{printf "%.1f", $1/1e9}')
    name=$(printf '%s' "$json" \
        | grep -o '"organism_name":"[^"]*"' | head -1 | cut -d'"' -f4)

    if [[ -z "$real_gb" || "$real_gb" == "0.0" ]]; then
        warn "$label -> $source found ($name) but size unreported"
    else
        delta=$(awk -v a="$approx" -v b="$real_gb" \
                'BEGIN{d=(a>b?a-b:b-a); print (a>0 && d/a>0.25) ? "off" : "ok"}')
        if [[ "$delta" == "off" ]]; then
            warn "$label -> $source is ${real_gb}Gb ($name), table says ${approx}Gb — update approx_gb/mem_gb"
        else
            ok "$label -> $source  ${real_gb}Gb  $name"
        fi
    fi
done < "$TABLE"

echo
if [[ "$fail" -eq 0 ]]; then
    echo "Preflight passed. Submit with:  ./submit.sh $TABLE ${PROTEOME:-<proteome.faa>}"
else
    echo "Preflight FAILED — fix the items above before submitting."
    exit 1
fi
