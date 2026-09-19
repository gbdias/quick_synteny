// Replaces the old extract_cds.nf + jcvi_ortholog.nf + jcvi_synteny_screen.nf
// chain entirely: bin/build_synteny_blocks.py works directly off the two
// miniprot GFFs (already produced by MINIPROT_ALIGN with loosened --outs/-N)
// with no separate ortholog aligner needed -- see that script's docstring.
// Runs in a plain Python container; no jcvi, no LAST, no gffread.

// Feeds the interactive plot's min-block-anchors slider (see
// plot_synteny_interactive.py's MBA_SLIDER_MIN). That slider works by
// filtering already-computed blocks by their anchor count client-side,
// which only reproduces what a real re-run at a higher threshold would have
// found because find_blocks() in build_synteny_blocks.py discovers blocks
// largest-first per chromosome pair and simply stops once the next chain is
// too small -- so a run's output at any threshold is exactly the "score >=
// threshold" subset of a run at any lower threshold. Forcing 5 (the
// slider's minimum) here means every block the slider could ever want to
// show already exists in this one file; nothing above 5 needs its own run,
// and this is the only synteny-chaining run this pipeline does per
// comparison -- there used to be a second, auto-tuned-threshold run feeding
// a since-removed static plot, but the interactive HTML is the only
// consumer now and it always wants the permissive floor.
process BUILD_CROSS_SYNTENY_FOR_SLIDER {
    tag "${target_name} vs ${comparison_name} (slider)"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(target_name), path(target_gff)
    tuple val(comparison_name), path(comparison_gff)
    val min_identity  // '' auto-tunes (see build_synteny_blocks.py); otherwise an explicit 0-1 override

    output:
    tuple val(target_name), val(comparison_name), path("${target_name}.${comparison_name}.slider_links.tsv"), emit: links

    script:
    def minIdentityFlag = min_identity ? "--min_identity ${min_identity}" : ''
    """
    build_synteny_blocks.py \\
        --query_gff ${target_gff} --subject_gff ${comparison_gff} \\
        ${minIdentityFlag} \\
        --min_block_anchors 5 \\
        --out ${target_name}.${comparison_name}.slider_links.tsv
    """
}

process BUILD_HOMEOLOG_SYNTENY_FOR_SLIDER {
    tag "${name} self (slider)"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(name), path(gff)
    val min_identity  // '' auto-tunes (see build_synteny_blocks.py); otherwise an explicit 0-1 override

    output:
    tuple val(name), path("${name}.slider_homeolog_links.tsv"), emit: links

    script:
    def minIdentityFlag = min_identity ? "--min_identity ${min_identity}" : ''
    """
    build_synteny_blocks.py \\
        --query_gff ${gff} --subject_gff ${gff} --self \\
        ${minIdentityFlag} \\
        --min_block_anchors 5 \\
        --out ${name}.slider_homeolog_links.tsv
    """
}

// Small standalone summary for the interactive HTML's stats panel (see
// pygenomeviz_plot.nf and plot_synteny_interactive.py's --stats) -- proteome
// size plus each genome's aligned-protein count/mean identity, straight from
// the same miniprot GFFs BUILD_CROSS_SYNTENY_FOR_SLIDER already consumes.
// Kept as its own process rather than folded into that one since it needs
// the proteome fasta too, which that process has no other use for.
process COMPUTE_ALIGNMENT_STATS {
    tag "${target_name} vs ${comparison_name}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(target_name), path(target_gff)
    tuple val(comparison_name), path(comparison_gff)
    path proteome

    output:
    path "${target_name}.${comparison_name}.stats.json", emit: stats

    script:
    """
    compute_alignment_stats.py \\
        --proteome ${proteome} \\
        --query_gff ${target_gff} --subject_gff ${comparison_gff} \\
        --query_name ${target_name} --subject_name ${comparison_name} \\
        --out ${target_name}.${comparison_name}.stats.json
    """
}

workflow BUILD_SYNTENY {
    take:
    target_gff       // tuple(name, path gff) from MINIPROT_ALIGN
    comparison_gff
    proteome         // path -- the same proteome fasta MINIPROT_ALIGN aligned against both genomes
    find_homeologs   // '', 'target', 'comparison', or 'both' -- which genome(s) to self-compare
    min_identity     // '' auto-tunes (see build_synteny_blocks.py); otherwise an explicit 0-1 override

    main:
    cross_slider = BUILD_CROSS_SYNTENY_FOR_SLIDER(target_gff, comparison_gff, min_identity).links
    stats        = COMPUTE_ALIGNMENT_STATS(target_gff, comparison_gff, proteome).stats

    // target_gff/comparison_gff tuples already carry their own role name as
    // element 0 (from main.nf's .branch{} split), so one process call over
    // whichever of them was asked for handles 'target', 'comparison', and
    // 'both' alike -- and an empty input channel (find_homeologs == '')
    // naturally yields an empty output channel, so there's nothing to
    // special-case there either.
    homeolog_inputs = Channel.empty()
    if (find_homeologs == 'target' || find_homeologs == 'both') {
        homeolog_inputs = homeolog_inputs.mix(target_gff)
    }
    if (find_homeologs == 'comparison' || find_homeologs == 'both') {
        homeolog_inputs = homeolog_inputs.mix(comparison_gff)
    }
    homeolog_slider = BUILD_HOMEOLOG_SYNTENY_FOR_SLIDER(homeolog_inputs, min_identity).links

    emit:
    cross_slider    = cross_slider
    homeolog_slider = homeolog_slider
    stats           = stats
}
