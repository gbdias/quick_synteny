// Synteny blocks straight from the two miniprot GFFs: every protein's hit
// positions in the two genomes are the anchors. bin/extract_hits.py makes a
// hit table per genome, and bin/chain.js (via bin/chain_blocks.mjs) chains
// them -- the same chainer the page runs client-side. Contract: chain.js's
// header.
//
// On an arm64 host CHAIN_* need the native arm64 node image: node under
// x86-64 emulation miscompiles chain.js (see chain_blocks.mjs).
//
// The page doesn't read any links file. links.tsv is the result at the
// run's parameters; *.slider_*links.tsv are chained with --min_block 5:
// extraction doesn't depend on the minimum block size (chain.js's CHAINING
// RULES, rule 6), so any higher threshold's blocks are exactly that file's
// blocks with score >= it.

process EXTRACT_HITS {
    tag "${name}"
    label 'process_low'
    label 'env_python'
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
    // exit 3: chain_blocks.mjs caught the node runtime miscomputing (see its
    // --no-maglev comment); a fresh process computes it correctly
    errorStrategy { task.exitStatus == 3 ? 'retry' : 'terminate' }
    maxRetries 2
    label 'env_nodejs'
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(target_name), path(target_hits)
    tuple val(reference_name), path(reference_hits)
    val min_identity  // '' takes chain.js's default (0.7); otherwise an explicit 0-1 floor
    val max_gap
    val min_block     // '' takes chain.js's default (15); the slider file always uses 5 (see top of file)

    output:
    tuple val(target_name), val(reference_name), path("${target_name}.${reference_name}.slider_links.tsv")
    tuple val(target_name), val(reference_name), path("${target_name}.${reference_name}.links.tsv")

    script:
    def idFlag = min_identity != '' ? "--min_identity ${min_identity}" : ''
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
    // exit 3: chain_blocks.mjs caught the node runtime miscomputing (see its
    // --no-maglev comment); a fresh process computes it correctly
    errorStrategy { task.exitStatus == 3 ? 'retry' : 'terminate' }
    maxRetries 2
    label 'env_nodejs'
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(name), path(hits)
    val min_identity
    val max_gap

    output:
    tuple val(name), path("${name}.slider_homeolog_links.tsv")

    script:
    def idFlag = min_identity != '' ? "--min_identity ${min_identity}" : ''
    """
    chain_blocks.mjs --query_hits ${hits} --self \\
        ${idFlag} --max_gap ${max_gap} --min_block 5 \\
        --out ${name}.slider_homeolog_links.tsv
    """
}

// The page's stats panel: proteome size plus each genome's aligned-protein
// count and mean identity, from the miniprot GFFs.
process COMPUTE_ALIGNMENT_STATS {
    tag "${target_name} vs ${reference_name}"
    label 'process_low'
    label 'env_python'
    publishDir "${params.outdir}/synteny", mode: 'copy'

    input:
    tuple val(target_name), path(target_gff)
    tuple val(reference_name), path(reference_gff)
    path proteome
    // see main.nf's input_sources
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
    proteome                // path -- the proteome aligned against both genomes
    min_identity            // '' takes chain.js's default (0.7); otherwise an explicit 0-1 floor
    max_gap                 // max rank step between consecutive chain members
    min_block               // '' takes chain.js's default (15); minimum block size (distinct loci) for links.tsv
    input_sources           // tuple(query source, query species, subject source, subject species,
                            //       proteome source, proteome species) -- for the stats panel

    main:
    hits = EXTRACT_HITS(target_gff.join(target_chrom_sizes).mix(reference_gff.join(reference_chrom_sizes))).hits
    hits_by_role = hits.branch {
        target: it[0] == 'target'
        reference: it[0] == 'reference'
    }

    // CHAIN_* only publish their links files
    CHAIN_CROSS(hits_by_role.target, hits_by_role.reference, min_identity, max_gap, min_block)
    CHAIN_SELF(hits, min_identity, max_gap)
    stats = COMPUTE_ALIGNMENT_STATS(target_gff, reference_gff, proteome, input_sources).stats

    emit:
    stats = stats
    hits  = hits
}
