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
// bioconda (so no biocontainers image exists for it) -- rather than hand-
// rolling a custom image for one package, this pip-installs a pinned Bokeh
// at runtime instead.
//
// That install has to land in the pygenomeviz container above, not the
// plain-Python container the non-plotting steps use: bokeh pulls in numpy,
// and the plain-Python image (quay.io/biocontainers/python:3.13.7) is
// missing libstdc++.so.6 entirely, with no package manager available in the
// image to add it -- numpy's C extensions fail to import there regardless
// of which numpy version pip resolves. pygenomeviz already has a working
// numpy (it's a matplotlib dependency), so bokeh installs cleanly on top of
// it instead of introducing a new base image just to fix a missing system
// library.
process RENDER_SYNTENY_INTERACTIVE {
    tag "${target_name} vs ${comparison_name}"
    label 'process_low'
    container 'quay.io/biocontainers/pygenomeviz:0.4.4--pyhdfd78af_0'
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

    output:
    path "${target_name}.${comparison_name}.synteny.interactive.html"

    script:
    def targetHomeologFlag = (target_homeolog_links.name != 'NO_HOMEOLOGS') ? "--target_homeolog_links ${target_homeolog_links}" : ''
    def comparisonHomeologFlag = (comparison_homeolog_links.name != 'NO_HOMEOLOGS') ? "--comparison_homeolog_links ${comparison_homeolog_links}" : ''
    """
    pip install --quiet bokeh==3.10.0

    plot_synteny_interactive.py \\
        --query_name ${target_name} --subject_name ${comparison_name} \\
        --query_chrom_sizes ${target_chrom_sizes} \\
        --subject_chrom_sizes ${comparison_chrom_sizes} \\
        --links ${links} \\
        ${targetHomeologFlag} \\
        ${comparisonHomeologFlag} \\
        --stats ${stats} \\
        --query_subtitle "${target_subtitle}" --subject_subtitle "${comparison_subtitle}" \\
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
    )
}
