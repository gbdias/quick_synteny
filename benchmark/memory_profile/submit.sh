#!/usr/bin/env bash
# Submit one SLURM job per genome, each with its own --mem/--time from the
# table. Separate jobs rather than an array so the small genomes can schedule
# on busy nodes instead of all waiting for a 128 GB slot.
#
#   ./submit.sh [genomes.tsv] <proteome.faa> [outdir]
#
# Env overrides:  PARTITION  ACCOUNT  CPUS  MINIPROT  DATASETS  NPROT  INDEX_OUT
#                 NODELIST  EXCLUSIVE  MEM_OVERRIDE
#
# NODELIST=node-c01,node-c02,node-c04   pin jobs to specific (e.g. high-RAM) nodes
# EXCLUSIVE=1                           whole-node allocation -- IMPORTANT if this
#                                        cluster doesn't hard-enforce --mem per job
#                                        (see collect.sh's sacct diagnosis note):
#                                        without it, an innocent small job can be
#                                        collaterally OOM-killed by a co-scheduled
#                                        job's overshoot, even while itself using a
#                                        fraction of what it requested.
# MEM_OVERRIDE=0                        replace every row's mem_gb with this value.
#                                        0 is a real sbatch idiom meaning "no per-job
#                                        cap, use whatever the node has" -- pair with
#                                        EXCLUSIVE=1 so that's the whole node's RAM.

set -euo pipefail

TABLE=${1:-genomes.tsv}
PROTEOME=${2:?usage: ./submit.sh <genomes.tsv> <proteome.faa> [outdir]}
OUTDIR=${3:-$PWD/memprofile_out}

PARTITION=${PARTITION:-}
ACCOUNT=${ACCOUNT:-}
NODELIST=${NODELIST:-}
EXCLUSIVE=${EXCLUSIVE:-}
MEM_OVERRIDE=${MEM_OVERRIDE:-}
CPUS=${CPUS:-8}            # keep FIXED across genomes -- thread count affects
                           # per-thread buffers, so varying it dirties the curve

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROTEOME=$(readlink -f "$PROTEOME")
mkdir -p "$OUTDIR"/{results,logs,genomes}

extra=()
[[ -n "$PARTITION" ]] && extra+=(--partition="$PARTITION")
[[ -n "$ACCOUNT"   ]] && extra+=(--account="$ACCOUNT")
[[ -n "$NODELIST"  ]] && extra+=(--nodelist="$NODELIST")
[[ -n "$EXCLUSIVE" ]] && extra+=(--exclusive)

echo "outdir:   $OUTDIR"
echo "proteome: $PROTEOME"
echo "cpus:     $CPUS (fixed across all genomes)"
echo "nodelist: ${NODELIST:-(any)}"
if [[ -n "$EXCLUSIVE" ]]; then
    echo "exclusive: yes"
else
    echo "exclusive: no  <-- if jobs keep failing well under their mem request (see"
    echo "                   collect.sh's sacct diagnosis note), a co-scheduled job"
    echo "                   on a shared node is probably why. Set EXCLUSIVE=1."
fi
echo "mem_override: ${MEM_OVERRIDE:-(use genomes.tsv per-row values)}"
echo "miniprot: ${MINIPROT:-miniprot}"
echo "datasets: ${DATASETS:-datasets}"
echo

n=0
while IFS=$'\t' read -r label source approx mem wall; do
    [[ -z "${label:-}" || "$label" == \#* ]] && continue
    [[ -n "$MEM_OVERRIDE" ]] && mem="$MEM_OVERRIDE"

    # sbatch's --mem=0 is a real SLURM idiom meaning "all memory on the node"
    # -- unlike every other value it must NOT get a G suffix, or it becomes a
    # literal (and useless) 0-byte request instead.
    mem_arg="${mem}G"; [[ "$mem" == "0" ]] && mem_arg="0"

    jid=$(sbatch --parsable \
        --job-name="mp_${label}" \
        --cpus-per-task="$CPUS" \
        --mem="$mem_arg" \
        --time="$wall" \
        --output="$OUTDIR/logs/${label}.slurm.log" \
        "${extra[@]}" \
        --wrap "MINIPROT='${MINIPROT:-miniprot}' DATASETS='${DATASETS:-datasets}' \
                NPROT='${NPROT:-2000}' INDEX_OUT='${INDEX_OUT:-/dev/null}' \
                bash '$HERE/profile_genome.sh' '$label' '$source' '$OUTDIR' '$PROTEOME'")

    printf '  %-14s %5s Gb  mem=%-5s time=%-9s job=%s\n' \
        "$label" "$approx" "${mem}G" "$wall" "$jid"
    echo "$jid" >> "$OUTDIR/.jobids"
    n=$((n+1))
done < "$TABLE"

echo
echo "submitted $n jobs. Watch with:   squeue -u \$USER -n \$(sed 's/^/mp_/' <<<'')"
echo "                   or simply:    squeue -u \$USER"
echo
echo "When they finish:   ./collect.sh $OUTDIR"
