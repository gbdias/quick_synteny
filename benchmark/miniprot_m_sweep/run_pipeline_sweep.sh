#!/usr/bin/env bash
# Run the REAL quick_synteny pipeline at several miniprot -M values on one
# already-validated dataset, then score how much synteny-calling accuracy
# each value costs via compare_links.py.
#
# This is the accuracy side of the -M investigation; sweep_m.sh (this
# directory's sibling) is the pure RAM/runtime side on synthetic index-only
# runs. Both are needed before recommending a value: sweep_m.sh alone can't
# tell you what a RAM saving costs in missed synteny.
#
#   ./run_pipeline_sweep.sh <assembly.fa> <comparison.fa> <proteome.faa> \
#       <outdir_prefix> [nextflow_profile]
#
# Env:  MVALS (default "1 2 3 4")   NEXTFLOW (default "nextflow")
#       NF_EXTRA (extra args passed through to every `nextflow run`, e.g.
#                 '--show_homeologs both')
#
# M=1 is included deliberately as the first/baseline value, not omitted --
# miniprot's own manual states its stock default IS M=1 ("[1]"), so this run
# reproduces default behavior exactly while keeping every invocation in the
# sweep structurally identical (no special-cased "no flag" run).
#
# Runs are SEQUENTIAL, not parallel: this is a runtime/RSS measurement too,
# and concurrent Nextflow runs on the same node would contaminate both (see
# BENCHMARK_PROTOCOL.md's note on this from the earlier concurrent-run
# incident).

set -euo pipefail

ASSEMBLY=${1:?assembly.fa}
COMPARISON=${2:?comparison.fa}
PROTEOME=${3:?proteome.faa}
OUTPREFIX=${4:?outdir_prefix, e.g. msweep_suecica}
PROFILE=${5:-standard}

MVALS=${MVALS:-"1 2 3 4"}
NEXTFLOW=${NEXTFLOW:-nextflow}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$HERE/../.." && pwd)

declare -a M_LIST=($MVALS)
BASELINE_M=${M_LIST[0]}

SUMMARY="${OUTPREFIX}_summary.tsv"
printf 'm\twall_sec\ttarget_peak_rss_gb\tcomparison_peak_rss_gb\n' > "$SUMMARY"

echo "assembly:   $ASSEMBLY"
echo "comparison: $COMPARISON"
echo "proteome:   $PROTEOME"
echo "M values:   $MVALS  (baseline = M=$BASELINE_M)"
echo

for m in "${M_LIST[@]}"; do
    outdir="${OUTPREFIX}_M${m}"
    echo "== M=$m -> $outdir =="
    t0=$(date +%s)
    "$NEXTFLOW" run "$REPO_ROOT/main.nf" -profile "$PROFILE" \
        --assembly "$ASSEMBLY" --comparison "$COMPARISON" --proteome "$PROTEOME" \
        --miniprot_m "$m" --outdir "$outdir" ${NF_EXTRA:-}
    t1=$(date +%s)

    trace="$outdir/pipeline_info/execution_trace.txt"
    tgt_rss="" ; cmp_rss=""
    if [[ -s "$trace" ]]; then
        # Header-driven column lookup -- Nextflow's default trace field set
        # can differ by version, so don't hardcode column positions.
        read -r -a hdr <<< "$(head -1 "$trace")"
        name_col=0; rss_col=0
        for i in "${!hdr[@]}"; do
            [[ "${hdr[$i]}" == "name" ]] && name_col=$((i+1))
            [[ "${hdr[$i]}" == "peak_rss" ]] && rss_col=$((i+1))
        done
        if [[ $name_col -gt 0 && $rss_col -gt 0 ]]; then
            tgt_rss=$(awk -F'\t' -v nc="$name_col" -v rc="$rss_col" \
                '$nc ~ /MINIPROT_ALIGN/ && $nc ~ /\(target\)/ {print $rc}' "$trace" | tail -1)
            cmp_rss=$(awk -F'\t' -v nc="$name_col" -v rc="$rss_col" \
                '$nc ~ /MINIPROT_ALIGN/ && $nc ~ /\(comparison\)/ {print $rc}' "$trace" | tail -1)
        fi
    fi
    printf '%s\t%s\t%s\t%s\n' "$m" "$((t1-t0))" "${tgt_rss:-NA}" "${cmp_rss:-NA}" >> "$SUMMARY"
    echo "  wall=${t1-t0}s  target_peak_rss=${tgt_rss:-NA}  comparison_peak_rss=${cmp_rss:-NA}"
    echo
done

echo "== $SUMMARY =="
column -t "$SUMMARY" 2>/dev/null || cat "$SUMMARY"

# --------------------------------------------------- concordance vs baseline
echo
echo "== concordance vs M=$BASELINE_M baseline =="
compare() {  # compare <glob under synteny/> <label>
    local pattern=$1 label=$2
    local base_file="${OUTPREFIX}_M${BASELINE_M}/synteny/$pattern"
    [[ -s "$base_file" ]] || { echo "  ($label: no $pattern in baseline run, skipping)"; return; }

    local args=(--baseline "M${BASELINE_M}=$base_file")
    for m in "${M_LIST[@]}"; do
        [[ "$m" == "$BASELINE_M" ]] && continue
        f="${OUTPREFIX}_M${m}/synteny/$pattern"
        [[ -s "$f" ]] && args+=(--other "M${m}=$f")
    done
    echo "-- $label ($pattern) --"
    python3 "$HERE/compare_links.py" "${args[@]}" \
        --out "${OUTPREFIX}_${label}_concordance.tsv"
    echo
}

compare "target.comparison.slider_links.tsv"          cross
compare "target.slider_homeolog_links.tsv"             homeolog_target
compare "comparison.slider_homeolog_links.tsv"          homeolog_comparison

echo "Done. Per-M runtime/RSS: $SUMMARY"
echo "      Concordance tables: ${OUTPREFIX}_*_concordance.tsv"
