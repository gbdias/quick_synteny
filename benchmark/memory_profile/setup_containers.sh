#!/usr/bin/env bash
# Pull the same container images the pipeline uses, so nothing needs
# installing and the profiled miniprot is byte-identical to the one
# quick_synteny runs.
#
#   ./setup_containers.sh [image_dir] [bind_path]
#
# Then copy-paste the two exports it prints, and run preflight/submit in that
# same shell.
#
# bind_path: any filesystem the container must see that Apptainer does not
# auto-bind. It auto-binds $HOME, $PWD and /tmp; if your genomes or proteome
# live on /scratch, /proj, /crex etc., pass that here. Multiple paths can be
# comma-separated.

set -euo pipefail

IMGDIR=${1:-$PWD/containers}
BIND=${2:-}

MINIPROT_IMG_URI="docker://quay.io/biocontainers/miniprot:0.18--h577a1d6_0"
DATASETS_IMG_URI="docker://quay.io/staphb/ncbi-datasets:18.35.0"

# Match the versions pinned in modules/local/{miniprot_align,download_assembly}.nf
MINIPROT_SIF="$IMGDIR/miniprot_0.18.sif"
DATASETS_SIF="$IMGDIR/ncbi-datasets_18.35.0.sif"

command -v apptainer >/dev/null 2>&1 && RUNTIME=apptainer \
    || { command -v singularity >/dev/null 2>&1 && RUNTIME=singularity \
         || { echo "ERROR: neither apptainer nor singularity on PATH." >&2
              echo "       Try: module avail apptainer" >&2; exit 1; }; }

mkdir -p "$IMGDIR"
echo "runtime: $RUNTIME"
echo "images:  $IMGDIR"
echo

for pair in "$MINIPROT_SIF|$MINIPROT_IMG_URI" "$DATASETS_SIF|$DATASETS_IMG_URI"; do
    sif=${pair%%|*}; uri=${pair##*|}
    if [[ -s "$sif" ]]; then
        echo "  [have] $(basename "$sif")"
    else
        echo "  [pull] $(basename "$sif") <- $uri"
        "$RUNTIME" pull --force "$sif" "$uri"
    fi
done

bindarg=""
[[ -n "$BIND" ]] && bindarg="--bind $BIND "

cat <<EOF

Done. Export these, then run preflight in the same shell:

  export MINIPROT="$RUNTIME exec ${bindarg}$MINIPROT_SIF miniprot"
  export DATASETS="$RUNTIME exec ${bindarg}$DATASETS_SIF datasets"

  ./preflight.sh genomes.tsv /path/to/proteome.faa

submit.sh forwards both variables to every SLURM job automatically.
EOF

[[ -z "$BIND" ]] && cat <<'EOF'

NOTE: no bind path given. Apptainer auto-binds $HOME, $PWD and /tmp only.
      If your genomes, proteome or output live elsewhere (/scratch, /proj,
      /crex ...), re-run as:  ./setup_containers.sh containers /scratch
EOF
