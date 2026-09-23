#!/usr/bin/env bash
# Sweep miniprot's -M k-mer sampling exponent for ONE genome, to see whether
# it rescues a genome that OOMs at the default (effectively M=1) setting, and
# by how much RAM it buys back. Pure RAM/runtime measurement -- reuses
# ../memory_profile/profile_genome.sh directly rather than duplicating its
# measurement logic; does not touch synteny quality (see compare_links.py /
# run_pipeline_sweep.sh in this directory for the accuracy side).
#
# Intended targets: the genomes that failed the default-M memory ladder
# (pleurodeles, axolotl per BENCHMARK_PROTOCOL.md's rung 3) -- does raising M
# bring them under your node's RAM ceiling at all?
#
#   ./sweep_m.sh <label> <accession|path> <proteome.faa> <mem_gb> <wall> [outdir]
#
# <mem_gb> of "0" is a real sbatch idiom meaning "no per-job cap, use whatever
# the node has" -- pair it with EXCLUSIVE=1 and NODELIST so that's the whole
# node's RAM, not a shared one. See memory_profile/submit.sh's comment on why
# this matters: on a cluster that doesn't hard-enforce --mem per job, a job
# can be collaterally OOM-killed by a co-scheduled neighbor's overshoot even
# while using a fraction of what it requested.
#
# Env:  MVALS (default "1 2 3 4")  PARTITION  ACCOUNT  NODELIST  EXCLUSIVE
#       CPUS  MINIPROT  DATASETS

set -euo pipefail

LABEL=${1:?label, e.g. axolotl}
SOURCE=${2:?accession or path}
PROTEOME=${3:?proteome.faa}
MEM=${4:?mem_gb, e.g. 128}
WALL=${5:?walltime, e.g. 24:00:00}
OUTDIR=${6:-$PWD/msweep_out}

MVALS=${MVALS:-"1 2 3 4"}
PARTITION=${PARTITION:-}
ACCOUNT=${ACCOUNT:-}
NODELIST=${NODELIST:-}
EXCLUSIVE=${EXCLUSIVE:-}
CPUS=${CPUS:-8}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MEMPROFILE_DIR=$(cd "$HERE/../memory_profile" && pwd)
PROTEOME=$(readlink -f "$PROTEOME")
mkdir -p "$OUTDIR"/{results,logs,genomes}

extra=()
[[ -n "$PARTITION" ]] && extra+=(--partition="$PARTITION")
[[ -n "$ACCOUNT"   ]] && extra+=(--account="$ACCOUNT")
[[ -n "$NODELIST"  ]] && extra+=(--nodelist="$NODELIST")
[[ -n "$EXCLUSIVE" ]] && extra+=(--exclusive)

# See note above: --mem=0 must NOT get a G suffix, or it becomes a literal
# (and useless) 0-byte request instead of "all memory on the node".
mem_arg="${MEM}G"; [[ "$MEM" == "0" ]] && mem_arg="0"

echo "genome:    $LABEL ($SOURCE)"
echo "M values:  $MVALS"
echo "mem/wall:  ${mem_arg} / $WALL  (same for every M -- if a high M needs less,"
echo "           that's exactly the result being measured; this just avoids"
echo "           under-provisioning a value that doesn't help)"
echo "nodelist:  ${NODELIST:-(any)}"
echo "exclusive: ${EXCLUSIVE:+yes}${EXCLUSIVE:-no}"
echo "outdir:    $OUTDIR"
echo

for m in $MVALS; do
    jid=$(sbatch --parsable \
        --job-name="msweep_${LABEL}_M${m}" \
        --cpus-per-task="$CPUS" \
        --mem="$mem_arg" \
        --time="$WALL" \
        --output="$OUTDIR/logs/${LABEL}_M${m}.slurm.log" \
        "${extra[@]}" \
        --wrap "MINIPROT='${MINIPROT:-miniprot}' DATASETS='${DATASETS:-datasets}' \
                MINIPROT_M='$m' \
                bash '$MEMPROFILE_DIR/profile_genome.sh' '$LABEL' '$SOURCE' '$OUTDIR' '$PROTEOME' '$m'")
    printf '  M=%-3s job=%s\n' "$m" "$jid"
    echo "$jid" >> "$OUTDIR/.jobids"
done

echo
echo "When they finish:"
echo "  $MEMPROFILE_DIR/collect.sh $OUTDIR $MEM"
echo "(collect.sh's own -M-sweep section activates automatically since these"
echo " rows carry m != default)"
