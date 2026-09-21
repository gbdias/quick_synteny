// Chromosome-size table used to size pyGenomeViz tracks and to filter out
// sequences shorter than --min_seq_size. The jcvi-layout chromEnd-marker
// trick from the original bash script is gone: pyGenomeViz takes chromosome
// lengths directly, it doesn't need a fake BED feature to infer track extent.
//
// Emitted in the FASTA's own (natural) sequence order -- deliberately NOT
// sorted by size here (an earlier version did `sort -k2,2nr`). The plot
// script needs both orderings (its "Order by size" switch, default on,
// toggles between them for the ring and dotplot alike -- see
// plot_synteny_interactive.py's Dataset.__init__ and SYN.dpNaturalOrder),
// and size order is trivial to derive from natural order in Python, but the
// reverse -- recovering natural order once it's been sorted away -- is not
// possible at all. Sorting exactly once, downstream, in the one place that
// actually needs both, is the only version of this that can't lose data.

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
    faSize -detailed ${genome_fasta} | awk -v s=${min_seq_size} '\$2 >= s' > ${name}.chrom.sizes
    """
}

// Assembly gaps (runs of N's) per genome -- scanned from the same renamed
// FASTA CHROM_SIZES uses, so gap coordinates land on the same chr1/chr2/...
// names as chrom.sizes/links.tsv rather than the original download's raw
// accessions. Plain Python (find_assembly_gaps.py), not a bioinformatics
// tool -- this is just a streaming character scan, nothing that needs a
// specialized container.
process FIND_ASSEMBLY_GAPS {
    tag "${name}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    tuple val(name), path(genome_fasta)
    val min_asm_gap

    output:
    tuple val(name), path("${name}.gaps.tsv"), emit: gaps

    script:
    """
    find_assembly_gaps.py --fasta ${genome_fasta} --min_gap ${min_asm_gap} --out ${name}.gaps.tsv
    """
}
