// build_synteny.nf already emits links.tsv directly in the format
// plot_synteny_interactive.py needs (genomic coordinates from the start), so
// there's no anchors-to-BED join step here anymore -- this just needs the
// chrom sizes for track layout and whichever links files were produced.
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

// plot_synteny_interactive.py draws a Circos-style ring, a linear zoom
// panel, and a whole-genome dotplot as one self-contained interactive HTML
// (Bokeh, CustomJS only -- no server). It needs Bokeh, which isn't on
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
// for. Pinned to bokeh=3.10.0 (Seqera's own build hash, stable for years per
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
// oras:// reference this pipeline's -profile slurm (Apptainer) needs would
// break -profile standard (Docker) outright, not just run unoptimally
// there. workflow.containerEngine (set by whichever profile is active --
// see nextflow.config/conf/slurm.config) picks the right one per run,
// restoring the same "one container line, works under either profile"
// property every other process in this pipeline already has.
process RENDER_SYNTENY_INTERACTIVE {
    tag "${target_name} vs ${comparison_name}"
    label 'process_low'
    container { workflow.containerEngine == 'docker'
        ? 'community.wave.seqera.io/library/bokeh:3.10.0--daa8ae8c4a0b7001'
        : 'oras://community.wave.seqera.io/library/bokeh:3.10.0--3fdfec626f703f33' }
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(target_name), val(comparison_name), path(links)
    path target_chrom_sizes
    path comparison_chrom_sizes
    path target_homeolog_links     // optional: a real file, or NO_HOMEOLOGS placeholder
    path comparison_homeolog_links // ditto
    path stats              // compute_alignment_stats.py output, for the page's stats panel -- always a real file (COMPUTE_ALIGNMENT_STATS is unconditional)
    val target_subtitle     // input file name/accession, shown under the target label
    val comparison_subtitle // ditto, for the comparison label
    path target_gaps        // RENAME_SEQUENCES gaps output -- always a real file (unconditional, unlike the homeolog links above)
    path comparison_gaps    // ditto

    output:
    path "${target_name}.${comparison_name}.synteny.interactive.html"

    script:
    def targetHomeologFlag = (target_homeolog_links.name != 'NO_HOMEOLOGS') ? "--target_homeolog_links ${target_homeolog_links}" : ''
    def comparisonHomeologFlag = (comparison_homeolog_links.name != 'NO_HOMEOLOGS') ? "--comparison_homeolog_links ${comparison_homeolog_links}" : ''
    """
    plot_synteny_interactive.py \\
        --query_name ${target_name} --subject_name ${comparison_name} \\
        --query_chrom_sizes ${target_chrom_sizes} \\
        --subject_chrom_sizes ${comparison_chrom_sizes} \\
        --links ${links} \\
        ${targetHomeologFlag} \\
        ${comparisonHomeologFlag} \\
        --stats ${stats} \\
        --query_subtitle "${target_subtitle}" --subject_subtitle "${comparison_subtitle}" \\
        --target_gaps ${target_gaps} --comparison_gaps ${comparison_gaps} \\
        --out_prefix ${target_name}.${comparison_name}.synteny
    """
}

workflow PYGENOMEVIZ_PLOT {
    take:
    links_slider           // tuple(target_name, comparison_name, path links.tsv) -- always computed at mba=5, feeds the interactive plot's slider
    target_chrom_sizes     // tuple(name, path chrom.sizes)
    comparison_chrom_sizes
    homeolog_links_slider  // tuple(name, path links.tsv), 0-2 items (one per requested side), always at mba=5
    stats                  // path -- compute_alignment_stats.py output
    target_subtitle        // val -- target's input file name, for the page's subtitle
    comparison_subtitle    // val -- comparison genome's input file name/accession, ditto
    target_gaps            // tuple(name, path gaps.tsv) -- RENAME_SEQUENCES gaps output
    comparison_gaps

    main:
    // homeolog_links_slider can carry a 'target' tuple, a 'comparison' tuple,
    // both, or neither (see BUILD_SYNTENY's find_homeologs) -- split back out
    // by role rather than passing it straight through, since
    // RENDER_SYNTENY_INTERACTIVE needs each side wired to its own CLI flag
    // (a single process input can't carry a variable 0-2 item list in lockstep
    // with the other per-run inputs).
    homeolog_branches = homeolog_links_slider.branch {
        target: it[0] == 'target'
        comparison: it[0] == 'comparison'
    }
    target_homeolog_file = homeolog_branches.target.map { it[1] }
        .ifEmpty(file("${projectDir}/assets/NO_HOMEOLOGS"))
    comparison_homeolog_file = homeolog_branches.comparison.map { it[1] }
        .ifEmpty(file("${projectDir}/assets/NO_HOMEOLOGS"))

    RENDER_SYNTENY_INTERACTIVE(
        links_slider,
        target_chrom_sizes.map { it[1] },
        comparison_chrom_sizes.map { it[1] },
        target_homeolog_file,
        comparison_homeolog_file,
        stats,
        target_subtitle, comparison_subtitle,
        target_gaps.map { it[1] },
        comparison_gaps.map { it[1] },
    )
}
