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
# -profile slurm`: axolotl is enormous (~28 Gb assembled) and this repo's
# OWN prior memory benchmark (APPLICATION_NOTE_PLAN.md section 2, measured
# 2026-09-20) put miniprot's peak RSS for it at ~260.6 GB at default
# settings -- far past conf/slurm.config's 'process_high' label (32 GB),
# which is sized for ordinary genomes, not this one. Rather than edit that
# shipped profile (which would then be wrong for every other genome), this
# script layers a run-specific override on top of it with -c, targeting
# just the one process (MINIPROT_ALIGN) that needs it. Whatever comparison
# genome auto-discovery lands on is likely to be another large, chromosome-
# level Caudata genome (the same benchmark's own axolotl/Pleurodeles pairing
# note) and to need close to as much RAM, so the override applies to BOTH
# of MINIPROT_ALIGN's two calls (target + comparison), not just axolotl's.
#
# This has NOT been run end-to-end against the real pipeline before --
# the 260 GB figure comes from raw miniprot index/align profiling
# (benchmark/memory_profile/), not from quick_synteny's own MINIPROT_ALIGN
# process, which uses the same miniprot commands but has never been
# measured at this genome size. Treat the numbers below as an informed
# starting point, not a guarantee -- watch the first run and adjust.
#
# !! EDIT THE "cluster-specific" SECTION BELOW BEFORE SUBMITTING. !!
########################################################################

# ---- cluster-specific: EDIT THESE -------------------------------------
REPO_DIR="${REPO_DIR:-$HOME/quick_synteny}"                        # this checkout (has main.nf)
ASSEMBLY_FASTA="${ASSEMBLY_FASTA:-/active/guilherme/axolotl/GCF_040938575.1_UKY_AmexF1_1_genomic.fna.gz}"    # target FASTA -- must already exist, see below
MINIPROT_NODE="${MINIPROT_NODE:-node-c01}"                        # the one high-mem node MINIPROT_ALIGN must land
                                                                     # on -- everything else is left unconstrained
                                                                     # (no `queue`/partition set anywhere below: this
                                                                     # cluster has none, and nextflow.config's own
                                                                     # default --slurm_queue='normal' would otherwise
                                                                     # still attach a nonexistent `--partition=normal`
                                                                     # to every job -- see the override config below)
ACCOUNT="${ACCOUNT:-}"                                             # leave empty if your cluster needs none
SINGULARITY_CACHE="${SINGULARITY_CACHE:-/scratch/$USER/quick_synteny_singularity}"
OUTDIR="${OUTDIR:-$PWD/results_axolotl}"
MINIPROT_CPUS="${MINIPROT_CPUS:-32}"                               # threads for MINIPROT_ALIGN; confirmed to fit
                                                                     # node-c01's actual core count -- 2026-09-22
MINIPROT_MEM_GB="${MINIPROT_MEM_GB:-320}"                          # ~260 GB measured + headroom, see header --
                                                                     # measured on GCA_002915635.3, not the
                                                                     # GCF_040938575.1 assembly below, so treat this
                                                                     # as ballpark, not exact, for this specific file.
                                                                     # Confirmed to fit node-c01 -- 2026-09-22
MINIPROT_TIME="${MINIPROT_TIME:-48.h}"                             # Nextflow duration syntax (NOT SLURM's
                                                                     # HH:MM:SS -- process.time doesn't accept that)
MINIPROT_EXCLUSIVE="${MINIPROT_EXCLUSIVE-1}"                       # whole-node allocation -- on a cluster that
                                                                     # doesn't hard-enforce --mem, a job this size
                                                                     # can otherwise be OOM-killed by a co-scheduled
                                                                     # neighbour's overshoot even while under its
                                                                     # own cap (see benchmark/memory_profile/submit.sh);
                                                                     # pass MINIPROT_EXCLUSIVE= (set but empty) to
                                                                     # disable -- note the missing ':', deliberate:
                                                                     # with ':-' an explicit empty value can never be
                                                                     # told apart from "unset", so it'd always fall
                                                                     # back to the default and the toggle would be
                                                                     # impossible to actually turn off
# -M: miniprot's k-mer sampling exponent, trades sensitivity for less RAM
# (see main.nf --miniprot_m / modules/local/miniprot_align.nf). Leave empty
# to use miniprot's own default (what the 260 GB figure above was measured
# at). benchmark/miniprot_m_sweep/ exists to find a value that brings this
# genome under a smaller RAM ceiling, but per APPLICATION_NOTE_PLAN.md
# section 2 that sweep hasn't actually been run yet -- so there is no
# validated M value to default to here. Set one yourself only if
# MINIPROT_MEM_GB above isn't available on your cluster, and treat the
# resulting synteny as unverified against the M=default result.
MINIPROT_M="${MINIPROT_M:-}"

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

# --assembly must be a plain (uncompressed) FASTA: RENAME_SEQUENCES --
# the very first process it hits -- runs bin/rename_sequences.py, which
# does a plain `open(fasta_path)` text read (see its parse_fasta_headers);
# fed a .gz it doesn't crash cleanly, it reads gzip's binary bytes as text
# and either throws a UnicodeDecodeError or, worse, silently produces a
# garbage rename mapping. Decompressed once here (skipped on a rerun/
# -resume if already done) rather than inside the pipeline itself, so this
# stays a plain shell concern, not a pipeline one.
if [[ "$ASSEMBLY_FASTA" == *.gz ]]; then
    DECOMPRESSED_DIR="/scratch/$USER/quick_synteny_axolotl"
    DECOMPRESSED_FASTA="$DECOMPRESSED_DIR/$(basename "${ASSEMBLY_FASTA%.gz}")"
    if [[ ! -f "$DECOMPRESSED_FASTA" ]]; then
        echo "ASSEMBLY_FASTA is gzipped ($ASSEMBLY_FASTA) -- decompressing once to $DECOMPRESSED_FASTA"
        mkdir -p "$DECOMPRESSED_DIR"
        gunzip -c "$ASSEMBLY_FASTA" > "$DECOMPRESSED_FASTA"
    else
        echo "Using already-decompressed $DECOMPRESSED_FASTA (from a previous run)"
    fi
    ASSEMBLY_FASTA="$DECOMPRESSED_FASTA"
fi

mkdir -p "$SINGULARITY_CACHE" "$OUTDIR"

# per-run resource override, layered on top of -profile slurm (conf/
# slurm.config) via -c rather than edited in place, so that shipped
# profile stays correct for every genome that ISN'T this size. withName
# beats withLabel regardless of file/CLI order, so this cleanly overrides
# just MINIPROT_ALIGN's inherited 'process_high' values below.
OVERRIDE_CONFIG="$(mktemp)"
trap 'rm -f "$OVERRIDE_CONFIG"' EXIT

cluster_opts=("--nodelist=${MINIPROT_NODE}")
[[ -n "$ACCOUNT" ]] && cluster_opts+=("--account=${ACCOUNT}")
[[ -n "$MINIPROT_EXCLUSIVE" ]] && cluster_opts+=("--exclusive")
cluster_opts_str="${cluster_opts[*]}"

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

    withName: 'MINIPROT_ALIGN' {
        cpus           = ${MINIPROT_CPUS}
        memory         = '${MINIPROT_MEM_GB}.GB'
        time           = '${MINIPROT_TIME}'
        // pinned to the one node with enough RAM -- see MINIPROT_NODE
        // above -- everything else is left free to land anywhere
        clusterOptions = '${cluster_opts_str}'
    }
    // RENAME_SEQUENCES/FIND_ASSEMBLY_GAPS/CHROM_SIZES all stream the FASTA
    // line-by-line rather than loading a whole sequence into memory (see
    // bin/find_assembly_gaps.py's docstring), so RAM stays fine at
    // process_low's default -- but a pure-Python byte-by-byte pass over a
    // multi-Gb chromosome can still run well past its 30m default time cap.
    // Time only; memory is untouched.
    withLabel: 'process_low' {
        time = '6.h'
    }
}
EOF

echo "assembly:          $ASSEMBLY_FASTA"
echo "taxid:              $TAXID (Ambystoma mexicanum)"
echo "outdir:             $OUTDIR"
echo "miniprot node:      $MINIPROT_NODE (MINIPROT_ALIGN: ${MINIPROT_CPUS} cpus, ${MINIPROT_MEM_GB}GB, ${MINIPROT_TIME}, exclusive=${MINIPROT_EXCLUSIVE:-no})"
echo "everything else:    unconstrained (no partitions on this cluster)"
echo "miniprot_m:         ${MINIPROT_M:-(pipeline default)}"
echo "singularity cache:  $SINGULARITY_CACHE"
echo

cd "$REPO_DIR"

nextflow run main.nf \
    --taxid "$TAXID" \
    --assembly "$ASSEMBLY_FASTA" \
    --outdir "$OUTDIR" \
    --singularity_cache_dir "$SINGULARITY_CACHE" \
    ${MINIPROT_M:+--miniprot_m "$MINIPROT_M"} \
    -profile slurm \
    -c "$OVERRIDE_CONFIG" \
    -resume
