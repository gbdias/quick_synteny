#!/usr/bin/env bash
# Worker: measure miniprot peak RSS for ONE genome, optionally at a given
# -M k-mer sampling exponent (see modules/local/miniprot_align.nf and
# benchmark/miniprot_m_sweep/ for the full RAM-vs-concordance sweep).
#
# Measures two numbers:
#   index_rss_gb  -- `miniprot -d` alone. The floor: index construction only.
#                    -M is applied HERE -- it's an indexing-time parameter,
#                    baked into the .mpi file, not redone per query.
#   align_rss_gb  -- `miniprot genome proteins`. The number that actually
#                    matters, = index + per-thread alignment buffers.
#
# The proteome is subsampled to 2000 sequences: alignment memory is dominated
# by the index and thread buffers, not by protein count, so a small query
# approximates the full-run peak at a fraction of the runtime.
#
# Usage (normally invoked by submit.sh, but runs standalone too):
#   ./profile_genome.sh <label> <accession|path> <outdir> <proteome.faa> [M]
#
# M (optional 5th arg, or env MINIPROT_M): miniprot's -M sampling exponent.
# Unset/empty uses miniprot's own default. Passing it also renames the row's
# label to "<label>_M<M>" so multiple M values for the same genome can share
# one results/ directory without overwriting each other -- see sweep_m.sh.

set -uo pipefail

LABEL=${1:?label}
SOURCE=${2:?accession or path}
OUTDIR=${3:?outdir}
PROTEOME=${4:?proteome}
MINIPROT_M=${5:-${MINIPROT_M:-}}   # optional; '' = miniprot's own default

THREADS=${SLURM_CPUS_PER_TASK:-8}
NPROT=${NPROT:-2000}
INDEX_OUT=${INDEX_OUT:-/dev/null}   # avoids writing a >50GB index for axolotl

# The genome file itself doesn't depend on M -- only the row/log naming does,
# so several M values for one genome share the downloaded file but don't
# clobber each other's results.
ROW_LABEL="$LABEL"
M_FLAG=()
if [[ -n "$MINIPROT_M" ]]; then
    ROW_LABEL="${LABEL}_M${MINIPROT_M}"
    M_FLAG=(-M "$MINIPROT_M")
fi

# MINIPROT / DATASETS may be a bare command ("miniprot") or a full container
# invocation ("apptainer exec --bind /scratch img.sif miniprot"). Split into
# arrays so multi-word values work. See setup_containers.sh.
MINIPROT=${MINIPROT:-miniprot}
DATASETS=${DATASETS:-datasets}
read -r -a MP  <<< "$MINIPROT"
read -r -a DS  <<< "$DATASETS"

GDIR="$OUTDIR/genomes"
RDIR="$OUTDIR/results"
LDIR="$OUTDIR/logs"
mkdir -p "$GDIR" "$RDIR" "$LDIR"

ROW="$RDIR/${ROW_LABEL}.tsv"
HDR="label\tm\tgenome_gb\tthreads\tindex_rss_gb\tindex_sec\talign_rss_gb\talign_sec\tstatus\tjobid\tnode"
JOBID=${SLURM_JOB_ID:-local}
NODE=$(hostname -s)

# Write a PENDING row up front. If SLURM OOM-kills the whole cgroup, nothing
# below gets to run -- but this row survives, and collect.sh reconciles it
# against sacct to recover the real MaxRSS and mark it OOM.
printf "$HDR\n%s\t%s\t\t%s\t\t\t\t\tPENDING\t%s\t%s\n" \
    "$LABEL" "${MINIPROT_M:-default}" "$THREADS" "$JOBID" "$NODE" > "$ROW"

emit() {  # emit <genome_gb> <idx_rss> <idx_sec> <aln_rss> <aln_sec> <status>
    printf "$HDR\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$LABEL" "${MINIPROT_M:-default}" "$1" "$THREADS" "$2" "$3" "$4" "$5" "$6" "$JOBID" "$NODE" > "$ROW"
}

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# ---------------------------------------------------------------- get genome
if [[ "$SOURCE" == /* ]]; then
    GENOME="$SOURCE"
    log "using local genome $GENOME"
else
    GENOME=$(ls "$GDIR/$LABEL"/*.fna "$GDIR/$LABEL"/*.fa 2>/dev/null | head -1)
    if [[ -z "${GENOME:-}" || ! -s "$GENOME" ]]; then
        log "downloading $SOURCE"
        tmpzip="$GDIR/${LABEL}.zip"
        if ! "${DS[@]}" download genome accession "$SOURCE" \
                --include genome --filename "$tmpzip" >"$LDIR/${LABEL}.download.log" 2>&1; then
            log "download FAILED (see $LDIR/${LABEL}.download.log)"
            emit "" "" "" "" "" "DOWNLOAD_FAILED"; exit 1
        fi
        rm -rf "$GDIR/$LABEL"; mkdir -p "$GDIR/$LABEL"
        unzip -q -o "$tmpzip" -d "$GDIR/$LABEL" && rm -f "$tmpzip"
        GENOME=$(find "$GDIR/$LABEL" -name '*.fna' -o -name '*.fa' | head -1)
    fi
fi
[[ -s "${GENOME:-}" ]] || { log "no genome file"; emit "" "" "" "" "" "NO_GENOME"; exit 1; }

# decompress if needed -- miniprot reads .gz, but sizing is simpler uncompressed
if [[ "$GENOME" == *.gz ]]; then
    log "decompressing"
    gunzip -kf "$GENOME" && GENOME="${GENOME%.gz}"
fi

# ------------------------------------------------------------- measure size
log "sizing $GENOME"
if command -v samtools >/dev/null 2>&1; then
    samtools faidx "$GENOME" 2>/dev/null \
        && BP=$(awk '{n+=$2} END{print n}' "${GENOME}.fai")
fi
[[ -n "${BP:-}" ]] || BP=$(awk '/^>/{next}{n+=length($0)} END{print n}' "$GENOME")
GB=$(awk -v b="$BP" 'BEGIN{printf "%.3f", b/1e9}')
log "genome = ${GB} Gb"

# --------------------------------------------------------------- subsample
QUERY="$OUTDIR/query_${NPROT}.faa"
if [[ ! -s "$QUERY" ]]; then
    awk -v n="$NPROT" '/^>/{c++} c>n{exit} {print}' "$PROTEOME" > "$QUERY"
    log "subsampled proteome -> $QUERY ($(grep -c '^>' "$QUERY") seqs)"
fi

# ------------------------------------------------------- run under GNU time
# GNU time reports ru_maxrss for the child, in kbytes. Log/tmp files are
# keyed on ROW_LABEL (not LABEL) so multiple -M values for the same genome
# don't clobber each other's .time/.err files.
run_timed() {  # run_timed <tag> <cmd...>  -> echoes "<rss_gb> <seconds>"
    local tag=$1; shift
    local tf="$LDIR/${ROW_LABEL}.${tag}.time"
    local t0 t1
    t0=$(date +%s)
    /usr/bin/time -v -o "$tf" "$@" >/dev/null 2>"$LDIR/${ROW_LABEL}.${tag}.err"
    local rc=$?
    t1=$(date +%s)
    local kb
    kb=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tf" 2>/dev/null)
    if [[ $rc -ne 0 || -z "$kb" ]]; then echo "ERR $((t1-t0))"; return 1; fi
    awk -v k="$kb" -v s="$((t1-t0))" 'BEGIN{printf "%.2f %d", k/1048576, s}'
}

# NOTE: `if read -r A B < <(cmd)` is a bash trap -- `if` only sees read's own
# exit status, never cmd's, because the command substitution runs in an
# async subshell whose $? is discarded. If run_timed fails and echoes
# "ERR <sec>", read still parses two tokens and returns 0, so the success
# branch ran unconditionally regardless of whether run_timed actually
# failed. Fixed by capturing into a variable first, which preserves $?.

# Non-fatal: `-d /dev/null` may be rejected on some builds, and the index
# figure is only the floor. The align measurement below is the one that
# matters, so a failure here must not abort the job.
log "index-only build (-d $INDEX_OUT)${MINIPROT_M:+, -M $MINIPROT_M}"
_out=$(run_timed index "${MP[@]}" -t"$THREADS" "${M_FLAG[@]}" -d "$INDEX_OUT" "$GENOME"); _rc=$?
if [[ $_rc -eq 0 ]]; then
    read -r IDX_RSS IDX_SEC <<< "$_out"
    log "index peak = ${IDX_RSS} GB in ${IDX_SEC}s"
else
    read -r _ IDX_SEC <<< "$_out"   # discard the literal "ERR" token
    IDX_RSS=""
    log "index step failed/OOM after ${IDX_SEC}s (continuing to align -- set INDEX_OUT=<path> if it was /dev/null)"
fi

# NOTE: this step passes the raw genome fasta (not a saved .mpi index) since
# INDEX_OUT defaults to /dev/null above -- miniprot re-indexes in memory on
# the fly here, so M_FLAG must be repeated on this line too, unlike the real
# pipeline (miniprot_align.nf), which saves and reuses one index and so only
# needs -M once, on its -d line.
log "align $NPROT proteins (index + alignment buffers)${MINIPROT_M:+, -M $MINIPROT_M}"
_out=$(run_timed align "${MP[@]}" -t"$THREADS" "${M_FLAG[@]}" "$GENOME" "$QUERY"); _rc=$?
if [[ $_rc -eq 0 ]]; then
    read -r ALN_RSS ALN_SEC <<< "$_out"
    log "align peak = ${ALN_RSS} GB in ${ALN_SEC}s"
else
    read -r _ ALN_SEC <<< "$_out"
    log "alignment FAILED/OOM after ${ALN_SEC}s"
    emit "$GB" "$IDX_RSS" "$IDX_SEC" "" "${ALN_SEC:-}" "ALIGN_OOM_OR_FAIL"; exit 137
fi

emit "$GB" "$IDX_RSS" "$IDX_SEC" "$ALN_RSS" "$ALN_SEC" "OK"
log "done: $ROW_LABEL  ${GB}Gb  index ${IDX_RSS}GB  align ${ALN_RSS}GB"
