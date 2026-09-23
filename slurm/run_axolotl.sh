#!/usr/bin/env bash
#SBATCH --job-name=quick_synteny_axolotl
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --output=quick_synteny_axolotl_%j.log

set -euo pipefail

########################################################################
# quick_synteny -- axolotl (Ambystoma mexicanum, taxid 8296) test run,
# full auto-discovery: only --taxid and --assembly are given, so the
# pipeline picks the comparison genome AND the proteome itself from NCBI
# (FIND_COMPARISON_ASSEMBLY / FIND_PROTEOME_ASSEMBLY climbing the lineage
# ladder from species up to --max_rank, default 'order') -- see main.nf
# --help. Nothing here pins a specific comparison/proteome accession.
#
# Why this needs more than `nextflow run main.nf --taxid 8296 ...
# -profile slurm`: axolotl is enormous (~29 Gb assembled), and miniprot's
# peak RAM is dominated by its genome index, ~10-12 GB per Gb of genome.
# A single whole-genome alignment peaked at 355 GB (target) / 340 GB
# (comparison) in the 2026-09-22 run, which only one node here can hold.
#
# So by default this script aligns in CHUNKS (--miniprot_chunk_gb, see
# main.nf --help and modules/local/miniprot_align.nf): the genome is split
# into ~MINIPROT_CHUNK_GB-sized pieces of whole sequences, each aligned as
# its own SLURM task with memory sized to that chunk (~11 GB per Gb of
# chunk + 4 GB), then merged into a result equivalent to one whole-genome
# run. At the default 2 Gb every chunk task asks for at most ~26 GB, so no
# special node is needed and the chunks run in parallel wherever they fit.
# Whatever comparison genome auto-discovery lands on (likely another large
# Caudata genome) is chunked the same way.
#
# The old single whole-genome job is still available: set MINIPROT_CHUNK_GB
# to an empty value (MINIPROT_CHUNK_GB= sbatch ...). Only then do the
# MINIPROT_NODE / MINIPROT_CPUS / MINIPROT_MEM_GB / MINIPROT_EXCLUSIVE
# settings below apply, pinning MINIPROT_ALIGN to the high-memory node.
#
# Chunked mode hasn't run at this scale yet: after the first run, check
# peak_rss for the MINIPROT_ALIGN_CHUNK tasks in
# $OUTDIR/pipeline_info/execution_trace.txt, and adjust MINIPROT_CHUNK_GB or
# the pipeline's --miniprot_gb_per_gb (default 11) if they're far off.
#
# !! EDIT THE "cluster-specific" SECTION BELOW BEFORE SUBMITTING. !!
########################################################################

# ---- cluster-specific: EDIT THESE -------------------------------------
REPO_DIR="${REPO_DIR:-$HOME/quick_synteny}"                        # this checkout (has main.nf)
ASSEMBLY_FASTA="${ASSEMBLY_FASTA:-/active/guilherme/axolotl/GCF_040938575.1_UKY_AmexF1_1_genomic.fna.gz}"    # target FASTA, plain or gzipped -- must already exist, see below
ACCOUNT="${ACCOUNT:-}"                                             # leave empty if your cluster needs none
SINGULARITY_CACHE="${SINGULARITY_CACHE:-/scratch/$USER/quick_synteny_singularity}"
OUTDIR="${OUTDIR:-$PWD/results_axolotl}"

# -- chunked alignment (default) --
MINIPROT_CHUNK_GB="${MINIPROT_CHUNK_GB-2}"                         # Gb of genome per miniprot task. Axolotl's largest
                                                                     # chromosome is 1.74 Gb, so at 2 each chunk holds
                                                                     # one big chromosome (or several smaller ones):
                                                                     # <= ~26 GB per task. No ':' in the default on
                                                                     # purpose: MINIPROT_CHUNK_GB= (set but empty)
                                                                     # switches to the whole-genome job below, which
                                                                     # ':-' would make impossible to tell from unset
MINIPROT_CHUNK_CPUS="${MINIPROT_CHUNK_CPUS:-8}"                    # threads per chunk task
MINIPROT_CHUNK_TIME="${MINIPROT_CHUNK_TIME:-12.h}"                 # Nextflow duration syntax (NOT SLURM's
                                                                     # HH:MM:SS -- process.time doesn't accept that)

# -- whole-genome alignment (only when MINIPROT_CHUNK_GB is empty) --
MINIPROT_NODE="${MINIPROT_NODE:-node-c01}"                        # the one high-mem node MINIPROT_ALIGN must land
                                                                     # on -- everything else is left unconstrained
MINIPROT_CPUS="${MINIPROT_CPUS:-32}"                               # threads for MINIPROT_ALIGN; confirmed to fit
                                                                     # node-c01's actual core count -- 2026-09-22
MINIPROT_MEM_GB="${MINIPROT_MEM_GB:-320}"                          # confirmed to fit node-c01 -- 2026-09-22; that
                                                                     # run peaked at 355 GB, above this request (the
                                                                     # cluster doesn't hard-enforce --mem)
MINIPROT_TIME="${MINIPROT_TIME:-48.h}"
MINIPROT_EXCLUSIVE="${MINIPROT_EXCLUSIVE-1}"                       # whole-node allocation -- on a cluster that
                                                                     # doesn't hard-enforce --mem, a job this size
                                                                     # can otherwise be OOM-killed by a co-scheduled
                                                                     # neighbour's overshoot even while under its
                                                                     # own cap (see benchmark/memory_profile/submit.sh);
                                                                     # MINIPROT_EXCLUSIVE= (set but empty) disables it

# -M: miniprot's k-mer sampling exponent, trades sensitivity for less RAM
# (see main.nf --miniprot_m / modules/local/miniprot_align.nf). Chunking
# caps RAM without that sensitivity cost, so leave this empty unless you
# have a specific reason; benchmark/miniprot_m_sweep/ measures the tradeoff.
MINIPROT_M="${MINIPROT_M:-}"
# --exclude_target: never pick a same-species comparison genome. The first
# run's discovery picked the target's own assembly (GCF_040938575.1) as the
# comparison, i.e. axolotl vs itself; set EXCLUDE_TARGET=1 for a
# cross-species comparison instead.
EXCLUDE_TARGET="${EXCLUDE_TARGET:-}"

# module loads -- EDIT to your cluster's actual module names, or drop this
# entirely if nextflow/apptainer are already on PATH
module load nextflow apptainer 2>/dev/null || true
# ------------------------------------------------------------------------

# Ambystoma mexicanum (axolotl)
TAXID=8296

if [[ ! -f "$ASSEMBLY_FASTA" ]]; then
    echo "ERROR: axolotl assembly not found at ASSEMBLY_FASTA=$ASSEMBLY_FASTA" >&2
    echo "--assembly is always user-supplied (never auto-discovered, only the" >&2
    echo "comparison genome and proteome are) -- fetch it first, e.g. from a" >&2
    echo "login/data-transfer node rather than spending this job's allocation" >&2
    echo "on a download:" >&2
    echo "  datasets download genome accession GCA_002915635.3 --include genome --filename axolotl.zip" >&2
    echo "  unzip axolotl.zip" >&2
    echo "  mv ncbi_dataset/data/GCA_002915635.3/*_genomic.fna \"$ASSEMBLY_FASTA\"" >&2
    echo "(GCA_002915635.3 is the accession this repo's own memory benchmark used" >&2
    echo "for axolotl -- see benchmark/memory_profile/genomes.tsv -- swap in a" >&2
    echo "different one if you specifically want another axolotl assembly.)" >&2
    exit 1
fi
# A gzipped --assembly is passed through as is: rename_sequences.py,
# split_genome.py and miniprot all read gzip directly.

mkdir -p "$SINGULARITY_CACHE" "$OUTDIR"

# per-run resource override, layered on top of -profile slurm (conf/
# slurm.config) via -c rather than edited in place, so that shipped
# profile stays correct for every genome that ISN'T this size. withName
# beats withLabel and the process's own directives regardless of file/CLI
# order.
OVERRIDE_CONFIG="$(mktemp)"
trap 'rm -f "$OVERRIDE_CONFIG"' EXIT

account_opt=""
[[ -n "$ACCOUNT" ]] && account_opt="--account=${ACCOUNT}"

if [[ -n "$MINIPROT_CHUNK_GB" ]]; then
    # memory is left to MINIPROT_ALIGN_CHUNK's own per-chunk formula (setting
    # it here would override that); no node pin -- chunks land anywhere
    miniprot_block="    withName: 'MINIPROT_ALIGN_CHUNK' {
        cpus           = ${MINIPROT_CHUNK_CPUS}
        time           = '${MINIPROT_CHUNK_TIME}'
        clusterOptions = '${account_opt}'
    }"
    miniprot_summary="chunked, ~${MINIPROT_CHUNK_GB} Gb per task (${MINIPROT_CHUNK_CPUS} cpus, ${MINIPROT_CHUNK_TIME}, memory sized per chunk), any node"
else
    cluster_opts=("--nodelist=${MINIPROT_NODE}")
    [[ -n "$account_opt" ]] && cluster_opts+=("$account_opt")
    [[ -n "$MINIPROT_EXCLUSIVE" ]] && cluster_opts+=("--exclusive")
    # pinned to the one node with enough RAM -- see MINIPROT_NODE above
    miniprot_block="    withName: 'MINIPROT_ALIGN' {
        cpus           = ${MINIPROT_CPUS}
        memory         = '${MINIPROT_MEM_GB}.GB'
        time           = '${MINIPROT_TIME}'
        clusterOptions = '${cluster_opts[*]}'
    }"
    miniprot_summary="whole genome on ${MINIPROT_NODE} (${MINIPROT_CPUS} cpus, ${MINIPROT_MEM_GB}GB, ${MINIPROT_TIME}, exclusive=${MINIPROT_EXCLUSIVE:-no})"
fi

cat > "$OVERRIDE_CONFIG" <<EOF
process {
    // this cluster has no partitions -- conf/slurm.config's global
    // \`queue = params.slurm_queue ?: 'normal'\` would otherwise still
    // attach a nonexistent \`--partition=normal\` to every job, since
    // params.slurm_queue defaults to the literal string 'normal' in
    // nextflow.config even when this script passes nothing for it.
    // Nulling it here (a later, top-level assignment beats the profile's
    // earlier one) makes every process submit with no --partition at all.
    queue = null

${miniprot_block}
    // RENAME_SEQUENCES streams the FASTA once (lengths + gaps in the same
    // pass), so RAM stays fine at process_low's default -- but one pass over
    // a ~30 GB FASTA can still run past its 30m default time cap. Time only;
    // memory is untouched.
    withLabel: 'process_low' {
        time = '6.h'
    }
}
EOF

echo "assembly:           $ASSEMBLY_FASTA"
echo "taxid:              $TAXID (Ambystoma mexicanum)"
echo "outdir:             $OUTDIR"
echo "miniprot:           $miniprot_summary"
echo "everything else:    unconstrained (no partitions on this cluster)"
echo "miniprot_m:         ${MINIPROT_M:-(pipeline default)}"
echo "exclude_target:     ${EXCLUDE_TARGET:-(off)}"
echo "singularity cache:  $SINGULARITY_CACHE"
echo

cd "$REPO_DIR"

nextflow run main.nf \
    --taxid "$TAXID" \
    --assembly "$ASSEMBLY_FASTA" \
    --outdir "$OUTDIR" \
    --singularity_cache_dir "$SINGULARITY_CACHE" \
    ${MINIPROT_M:+--miniprot_m "$MINIPROT_M"} \
    ${MINIPROT_CHUNK_GB:+--miniprot_chunk_gb "$MINIPROT_CHUNK_GB"} \
    ${EXCLUDE_TARGET:+--exclude_target} \
    -profile slurm \
    -c "$OVERRIDE_CONFIG" \
    -resume
