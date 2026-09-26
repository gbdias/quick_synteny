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
// It runs on envs/bokeh.yml, conda-forge's bokeh package, whose numpy
// dependency also encodes the embedded hit-table payload (see
// bin/plot_synteny_interactive.py's build_hits_payload). Bokeh has no
// bioconda recipe, and the plain-Python image once used for the other steps
// had no libstdc++ to pip-install numpy into, hence its own environment.
process RENDER_SYNTENY_INTERACTIVE {
    tag "${target_name} vs ${reference_name}"
    label 'process_low'
    label 'env_bokeh'
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
    path target_lookup       // RENAME_SEQUENCES lookup -- original sequence IDs for the ideogram hovers
    path reference_lookup   // ditto

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
        --target_lookup ${target_lookup} --reference_lookup ${reference_lookup} \\
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
    target_lookup          // tuple(name, path rename_lookup.tsv) -- RENAME_SEQUENCES lookup output
    reference_lookup

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
        target_lookup.map { it[1] },
        reference_lookup.map { it[1] },
    )
}
