// plot_synteny_interactive.py embeds each genome's hit table (from
// BUILD_SYNTENY's EXTRACT_HITS) directly -- see bin/chain.js's EMBEDDED
// PAYLOAD header -- and chains them client-side (bin/chain.js, inlined into the page), so
// this module ships raw hits, not precomputed links/blocks. build_synteny.nf
// still emits links.tsv/slider_links.tsv for other consumers, just not this one.
//
// This pipeline's one plot output is the interactive HTML below: a Circos-
// style ring, a linear zoom panel, and a whole-genome dotplot, all built
// with Bokeh/CustomJS and all individually exportable as SVG once the
// viewer has picked a palette/chromosome order/zoom they like. There used
// to also be a pair of static matplotlib renders (plot_synteny.py linear +
// plot_synteny_circos.py circular) feeding fixed PDF/SVG/PNG files, but
// those couldn't be personalized at all and duplicated the ring/zoom
// geometry the interactive version already draws -- retired in favor of
// putting all further effort into the one interactive view.

// plot_synteny_interactive.py draws a Circos-style ring, a linear detail
// panel, and a whole-genome dotplot as one self-contained interactive HTML
// (Bokeh, CustomJS only, BokehJS inlined -- no server, works offline).
// It needs Bokeh, which isn't on
// bioconda (so no biocontainers image exists for it, and the plain-Python
// image the non-plotting steps use -- quay.io/biocontainers/python:3.13.7 --
// is missing libstdc++.so.6 entirely, with no package manager available in
// the image to add it, so pip-installing Bokeh's numpy dependency there
// fails regardless of version) -- rather than repurposing some unrelated
// bioinformatics tool's container for its incidental working numpy (this
// process used to run on the pygenomeviz biocontainer for exactly that
// reason, even though pygenomeviz itself is never actually used anywhere in
// this repo), this runs on a container built directly from conda-forge's
// real `bokeh` package via Seqera Containers (Wave; docs.seqera.io/wave) --
// community.wave.seqera.io builds and permanently hosts images from
// Bioconda/conda-forge/PyPI on demand, pullable like any other registry
// image, no Nextflow wave plugin needed. Bokeh ships in the image already,
// so there's no runtime pip install left to fail or to need network access
// for. numpy (used to encode the embedded hit-table payload -- see
// bin/plot_synteny_interactive.py's build_hits_payload) ships alongside
// Bokeh in the same conda-forge image, so it needs no separate handling
// here. Pinned to bokeh=3.10.0 (Seqera's own build hash, stable for years per
// their docs) and to linux/amd64, matching docker.runOptions =
// '--platform=linux/amd64' in nextflow.config's standard profile -- every
// container in this pipeline runs as amd64 there (Apple Silicon dev
// machines go through Rosetta, per the README). linux/arm64 builds of the
// same bokeh=3.10.0 also exist (community.wave.seqera.io/library/bokeh:
// 3.10.0--35eff5f46e2379fb for Docker, --77f10b7a44d9daf0 for Singularity)
// if a native-arm64 profile is ever added.
//
// Two different image references, not one: Seqera Containers builds a
// genuinely different artifact per target engine -- the Singularity one is
// a native SIF (single sylabs.sif.layer blob, no Docker-style layered
// filesystem), which plain `docker pull` cannot consume at all, so the same
// SIF this pipeline's -profile slurm (Apptainer) needs would break -profile
// standard (Docker) outright, not just run unoptimally there. Singularity
// gets that SIF as a direct HTTPS download of its registry blob (the SIF
// build of oras://community.wave.seqera.io/library/bokeh:3.10.0--
// 3fdfec626f703f33), the convention nf-core modules use for Seqera images:
// some Singularity/Apptainer versions' ORAS clients reject the manifest type
// Seqera serves those under (seen on the cluster, 2026-09-23, for the
// nodejs image built the same way). workflow.containerEngine (set by whichever profile is active --
// see nextflow.config/conf/slurm.config) picks the right one per run,
// restoring the same "one container line, works under either profile"
// property every other process in this pipeline already has.
process RENDER_SYNTENY_INTERACTIVE {
    tag "${target_name} vs ${reference_name}"
    label 'process_low'
    container { workflow.containerEngine == 'docker'
        ? 'community.wave.seqera.io/library/bokeh:3.10.0--daa8ae8c4a0b7001'
        : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d1d2c3a4140cc15a434b599754b57408b4c840e56724c5f8ba5e8e2213a91dd/data' }
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(target_name), path(target_hits)
    tuple val(reference_name), path(reference_hits)
    path target_chrom_sizes
    path reference_chrom_sizes
    val min_identity        // '' auto-tunes client-side (SYNCHAIN.autoParams); otherwise an explicit 0-1 floor
    val max_gap
    val min_block           // '' auto-tunes client-side; otherwise the initial Min block size control value
    path stats               // compute_alignment_stats.py output, for the page's stats panel -- always a real file (COMPUTE_ALIGNMENT_STATS is unconditional)
    val target_subtitle      // input file name/accession, shown under the target label
    val reference_subtitle  // ditto, for the reference label
    path target_gaps         // RENAME_SEQUENCES gaps output -- always a real file (unconditional)
    path reference_gaps     // ditto

    output:
    path "${target_name}.${reference_name}.synteny.interactive.html"

    script:
    def idFlag = min_identity ? "--min_identity ${min_identity}" : ''
    def blockFlag = min_block ? "--min_block ${min_block}" : ''
    """
    plot_synteny_interactive.py \\
        --query_name ${target_name} --subject_name ${reference_name} \\
        --query_chrom_sizes ${target_chrom_sizes} \\
        --subject_chrom_sizes ${reference_chrom_sizes} \\
        --target_hits ${target_hits} --reference_hits ${reference_hits} \\
        ${idFlag} --max_gap ${max_gap} ${blockFlag} \\
        --stats ${stats} \\
        --query_subtitle "${target_subtitle}" --subject_subtitle "${reference_subtitle}" \\
        --target_gaps ${target_gaps} --reference_gaps ${reference_gaps} \\
        --out_prefix ${target_name}.${reference_name}.synteny
    """
}

workflow PYGENOMEVIZ_PLOT {
    take:
    hits                   // tuple(name, path hits.tsv.gz) -- BUILD_SYNTENY's EXTRACT_HITS output, target and reference mixed together
    target_chrom_sizes     // tuple(name, path chrom.sizes)
    reference_chrom_sizes
    min_identity           // '' auto-tunes; otherwise an explicit 0-1 floor (the page's initial control value)
    max_gap                // max rank step between consecutive chain members (the page's initial control value)
    min_block              // '' auto-tunes; otherwise the page's initial Min block size control value
    stats                  // path -- compute_alignment_stats.py output
    target_subtitle        // val -- target's input file name, for the page's subtitle
    reference_subtitle    // val -- reference genome's input file name/accession, ditto
    target_gaps            // tuple(name, path gaps.tsv) -- RENAME_SEQUENCES gaps output
    reference_gaps

    main:
    hits_by_role = hits.branch {
        target: it[0] == 'target'
        reference: it[0] == 'reference'
    }

    RENDER_SYNTENY_INTERACTIVE(
        hits_by_role.target,
        hits_by_role.reference,
        target_chrom_sizes.map { it[1] },
        reference_chrom_sizes.map { it[1] },
        min_identity, max_gap, min_block,
        stats,
        target_subtitle, reference_subtitle,
        target_gaps.map { it[1] },
        reference_gaps.map { it[1] },
    )
}
