// Synteny blocks straight from the two miniprot GFFs -- no separate
// ortholog aligner: the proteome was already aligned to both genomes, and
// every protein's hit positions in the two genomes ARE the synteny anchors.
// bin/extract_hits.py turns each GFF into a per-genome hit table, and
// bin/chain.js (via bin/chain_blocks.mjs) chains them -- the same chainer the
// interactive page runs client-side. Contract: bin/chain.js's header.

// CHAIN_* run conda-forge's nodejs 26.8.2 from Seqera Containers, the same
// way RENDER_SYNTENY_INTERACTIVE gets bokeh (see pygenomeviz_plot.nf for why
// there are two references: Docker and Singularity get different image
// formats). chain.js is plain, dependency-free JS, so any node >= 18 would
// do; what matters for the image is that it ships bash and procps, both of
// which Nextflow needs in every task container (node:*-slim lacks `ps`).
// Singularity gets the SIF as a direct HTTPS download of its registry blob
// (the SIF build of oras://community.wave.seqera.io/library/nodejs:26.8.2--
// 79cbd548ac9675ad), the convention nf-core modules use for Seqera images,
// rather than the oras:// reference itself -- the cluster's Singularity
// refused that with "could not get image manifest, received mediaType:
// application/vnd.docker.distribution.manifest.v2+json" (2026-09-23). Both
// are linux/amd64 builds, matching docker.runOptions in the standard
// profile; linux/arm64 builds of the same package also exist
// (nodejs:26.8.2--6646c230528b0eae for Docker, --864372a79e757566 for
// Singularity) if a native-arm64 profile is ever added.
//
// The interactive page doesn't read any links file -- it embeds the hit
// tables and re-chains client-side (see pygenomeviz_plot.nf). links.tsv is
// the published result at the run's parameters; the *.slider_*links.tsv
// files are chained with --min_block 5 for downstream comparisons:
// extraction in chain.js is
// independent of the minimum block size (bin/chain.js's CHAINING RULES header, rule 6), so the blocks
// at any higher threshold are exactly that file's blocks with score >= it.

process EXTRACT_HITS {
    tag "${name}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(name), path(gff), path(chrom_sizes)

    output:
    tuple val(name), path("${name}.hits.tsv.gz"), emit: hits

    script:
    """
    extract_hits.py --gff ${gff} --chrom_sizes ${chrom_sizes} --out ${name}.hits.tsv.gz
    """
}

process CHAIN_CROSS {
    tag "${target_name} vs ${reference_name}"
    label 'process_low'
    container { workflow.containerEngine == 'docker'
        ? 'community.wave.seqera.io/library/nodejs:26.8.2--e0ce3f03c0e9c0f0'
        : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d56dd8ae13cc6187cea921d2c32cc89f5ca81a6277dbb65cd0906fbddcf1da5a/data' }
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(target_name), path(target_hits)
    tuple val(reference_name), path(reference_hits)
    val min_identity  // '' auto-tunes (see chain.js autoParams); otherwise an explicit 0-1 floor
    val max_gap
    val min_block     // '' auto-tunes; the slider file always uses 5 (see top of file)

    output:
    tuple val(target_name), val(reference_name), path("${target_name}.${reference_name}.slider_links.tsv"), emit: slider
    tuple val(target_name), val(reference_name), path("${target_name}.${reference_name}.links.tsv"), emit: links

    script:
    def idFlag = min_identity ? "--min_identity ${min_identity}" : ''
    def blockFlag = min_block ? "--min_block ${min_block}" : ''
    """
    chain_blocks.mjs --query_hits ${target_hits} --subject_hits ${reference_hits} \\
        ${idFlag} --max_gap ${max_gap} ${blockFlag} \\
        --out ${target_name}.${reference_name}.links.tsv
    chain_blocks.mjs --query_hits ${target_hits} --subject_hits ${reference_hits} \\
        ${idFlag} --max_gap ${max_gap} --min_block 5 \\
        --out ${target_name}.${reference_name}.slider_links.tsv
    """
}

process CHAIN_SELF {
    tag "${name} self"
    label 'process_low'
    container { workflow.containerEngine == 'docker'
        ? 'community.wave.seqera.io/library/nodejs:26.8.2--e0ce3f03c0e9c0f0'
        : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d56dd8ae13cc6187cea921d2c32cc89f5ca81a6277dbb65cd0906fbddcf1da5a/data' }
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(name), path(hits)
    val min_identity
    val max_gap

    output:
    tuple val(name), path("${name}.slider_homeolog_links.tsv"), emit: slider

    script:
    def idFlag = min_identity ? "--min_identity ${min_identity}" : ''
    """
    chain_blocks.mjs --query_hits ${hits} --self \\
        ${idFlag} --max_gap ${max_gap} --min_block 5 \\
        --out ${name}.slider_homeolog_links.tsv
    """
}

// Small standalone summary for the interactive HTML's stats panel (see
// pygenomeviz_plot.nf and plot_synteny_interactive.py's --stats) -- proteome
// size plus each genome's aligned-protein count/mean identity, straight from
// the miniprot GFFs. Kept as its own process since it needs the proteome
// fasta too, which nothing else in this module reads.
process COMPUTE_ALIGNMENT_STATS {
    tag "${target_name} vs ${reference_name}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(target_name), path(target_gff)
    tuple val(reference_name), path(reference_gff)
    path proteome
    // what the page names each input by: its filename if user-supplied, else
    // the NCBI accession it was discovered from; plus its species when known
    // ('' otherwise) -- see main.nf's input_sources
    tuple val(query_source), val(query_species), val(subject_source), val(subject_species),
          val(proteome_source), val(proteome_species)

    output:
    path "${target_name}.${reference_name}.stats.json", emit: stats

    script:
    // single-quoted for the shell, so file names/species with spaces or quotes survive
    def q = { v -> "'" + v.toString().replace("'", "'\\''") + "'" }
    """
    compute_alignment_stats.py \\
        --proteome ${proteome} \\
        --query_gff ${target_gff} --subject_gff ${reference_gff} \\
        --query_name ${target_name} --subject_name ${reference_name} \\
        --query_source ${q(query_source)} --query_species ${q(query_species)} \\
        --subject_source ${q(subject_source)} --subject_species ${q(subject_species)} \\
        --proteome_source ${q(proteome_source)} --proteome_species ${q(proteome_species)} \\
        --out ${target_name}.${reference_name}.stats.json
    """
}

workflow BUILD_SYNTENY {
    take:
    target_gff              // tuple(name, path gff), seqids already renamed (RENAME_GFF)
    reference_gff
    target_chrom_sizes      // tuple(name, path chrom.sizes) -- hits on sequences not listed are dropped
    reference_chrom_sizes
    proteome                // path -- the same proteome fasta MINIPROT_ALIGN aligned against both genomes
    min_identity            // '' auto-tunes; otherwise an explicit 0-1 floor
    max_gap                 // max gene-rank step between consecutive chain members
    min_block               // '' auto-tunes; minimum block size (distinct loci) for links.tsv
    input_sources           // tuple(query source, query species, subject source, subject species,
                            //       proteome source, proteome species) -- for the stats panel

    main:
    // one EXTRACT_HITS call over both genomes (a DSL2 process can't be
    // invoked twice in one scope), split back by the role name each tuple
    // carries
    hits = EXTRACT_HITS(target_gff.join(target_chrom_sizes).mix(reference_gff.join(reference_chrom_sizes))).hits
    hits_by_role = hits.branch {
        target: it[0] == 'target'
        reference: it[0] == 'reference'
    }

    chained = CHAIN_CROSS(hits_by_role.target, hits_by_role.reference, min_identity, max_gap, min_block)
    stats   = COMPUTE_ALIGNMENT_STATS(target_gff, reference_gff, proteome, input_sources).stats

    // both genomes are always self-chained (homeologs); the page decides
    // whether to show them (its Show self-links switch)
    homeolog_slider = CHAIN_SELF(hits, min_identity, max_gap).slider

    emit:
    cross_slider    = chained.slider
    cross_links     = chained.links
    homeolog_slider = homeolog_slider
    stats           = stats
    hits            = hits
}
