// Chromosome-size table used to size pyGenomeViz tracks and to filter out
// sequences shorter than --min_seq_size. The jcvi-layout chromEnd-marker
// trick from the original bash script is gone: pyGenomeViz takes chromosome
// lengths directly, it doesn't need a fake BED feature to infer track extent.

process CHROM_SIZES {
    tag "${name}"
    label 'process_low'
    container 'quay.io/biocontainers/ucsc-fasize:482--h0b57e2e_0'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    tuple val(name), path(genome_fasta)
    val min_seq_size

    output:
    tuple val(name), path("${name}.chrom.sizes"), emit: chrom_sizes

    script:
    """
    faSize -detailed ${genome_fasta} | awk -v s=${min_seq_size} '\$2 >= s' | sort -k2,2nr > ${name}.chrom.sizes
    """
}
