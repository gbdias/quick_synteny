#!/usr/bin/env bash
# Gather per-genome rows into one table, recover MaxRSS for jobs that were
# OOM-killed (their own row never got written), fit RAM vs genome size using
# only the default-M baseline rows, and report where the ceiling falls.
# Separately summarizes any -M sweep rows (see sweep_m.sh) if present.
#
#   ./collect.sh [outdir] [node_ram_gb]

set -uo pipefail

OUTDIR=${1:-$PWD/memprofile_out}
NODE_RAM=${2:-128}
OUT="$OUTDIR/memory_profile.tsv"

[[ -d "$OUTDIR/results" ]] || { echo "no results dir in $OUTDIR" >&2; exit 1; }

HDR='label\tm\tgenome_gb\tthreads\tindex_rss_gb\tindex_sec\talign_rss_gb\talign_sec\tstatus\tjobid\tnode'
printf "$HDR\n" > "$OUT"

for f in "$OUTDIR"/results/*.tsv; do
    [[ -e "$f" ]] || continue
    row=$(tail -n +2 "$f" | head -1)
    [[ -z "$row" ]] && continue

    status=$(cut -f9 <<<"$row")
    jobid=$(cut -f10 <<<"$row")

    # A PENDING row means the cgroup was killed before the script could finish.
    # sacct still knows the peak and the reason.
    if [[ "$status" == PENDING && "$jobid" != local ]]; then
        state=$(sacct -j "$jobid" --noheader --format=State%20 2>/dev/null | head -1 | xargs)
        maxrss=$(sacct -j "$jobid" --noheader --format=MaxRSS --units=G 2>/dev/null \
                 | tr -d ' G' | grep -E '^[0-9.]+$' | sort -g | tail -1)
        case "$state" in
            OUT_OF_MEMORY*) status="OOM" ;;
            TIMEOUT*)       status="TIMEOUT" ;;
            RUNNING*|PENDING*) status="STILL_RUNNING" ;;
            *)              status="KILLED_${state// /_}" ;;
        esac
        # field 7 = align_rss_gb
        row=$(awk -F'\t' -v OFS='\t' -v s="$status" -v m="${maxrss:-}" \
              '{ $9=s; if (m!="" && $7=="") $7=m; print }' <<<"$row")
    fi
    printf '%s\n' "$row" >> "$OUT"
done

echo "== $OUT =="
# Not `column -t -s$'\t'`: many `column` builds (BSD, older util-linux without
# -n) silently collapse consecutive/empty tab fields, which shifts columns
# and looks like data corruption when a field is legitimately blank (e.g. an
# OOM row with no align_rss_gb). This awk printer preserves every field.
awk -F'\t' '
NR==1 { for (i=1;i<=NF;i++) w[i]=length($i) }
{ for (i=1;i<=NF;i++) if (length($i)>w[i]) w[i]=length($i); rows[NR]=$0 }
END {
    for (r=1;r<=NR;r++) {
        n=split(rows[r], f, "\t")
        line=""
        for (i=1;i<=n;i++) line = line sprintf("%-*s  ", w[i], f[i])
        print line
    }
}' "$OUT"
echo

# --------------------------------------------------------- the scaling fit
# Only default-M rows go into the fit -- an -M sweep row is a deliberately
# different configuration and would distort a RAM-vs-genome-size curve meant
# to characterize miniprot's stock behavior. Least squares forced through the
# origin: the index dominates and scales with sequence length, so a
# proportional model is the right one and extrapolates more honestly than an
# unconstrained line.
awk -F'\t' -v cap="$NODE_RAM" '
NR>1 && $2=="default" && $9=="OK" && $3+0>0 && $7+0>0 {
    n++; sxy+=$3*$7; sxx+=$3*$3; if ($3+0>maxg) maxg=$3+0
}
END {
    if (n < 2) { print "Not enough default-M successful points to fit (need >=2, have " n ")."; exit }
    k = sxy/sxx
    printf "points fitted      : %d (largest %.1f Gb)\n", n, maxg
    printf "slope              : %.2f GB RAM per Gb of genome\n", k
    printf "predicted ceiling  : %.1f Gb genome hits %d GB RAM\n", cap/k, cap
    printf "\nprediction 8 said 15-20 Gb -> %s\n", \
        (cap/k >= 15 && cap/k <= 20) ? "CONFIRMED" : "REFUTED, update the protocol"
}' "$OUT"

echo
echo "OOM / failed rows (these are data, not errors):"
awk -F'\t' 'NR>1 && $9!="OK" {printf "  %-16s m=%-8s %-8s %s\n", $1, $2, $3, $9}' "$OUT"

# ------------------------------------------------------------- the M sweep
if awk -F'\t' 'NR>1 && $2!="default"{f=1} END{exit !f}' "$OUT"; then
    echo
    echo "== -M sweep rows (see sweep_m.sh) =="
    awk -F'\t' 'NR==1 || $2!="default"' "$OUT" \
        | awk -F'\t' '{printf "  %-14s m=%-6s %-8s idx=%-8s aln=%-8s %s\n", $1,$2,$3,$5,$7,$9}'
    echo
    echo "  Reading this: compare align_rss_gb (aln=) across m= for the same"
    echo "  label. A genome that OOM'd at m=default succeeding at a higher m"
    echo "  is the headline result -- pair it with compare_links.py's"
    echo "  concordance numbers (benchmark/miniprot_m_sweep/) before citing a"
    echo "  RAM saving, since this table alone doesn't show the sensitivity cost."
fi
