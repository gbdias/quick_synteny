// Single full read of the genome: renames sequence IDs to plot-friendly
// labels (chr1, chr2, ...) and, in the same pass, computes the chromosome
// size table and the assembly-gap table that used to each cost a separate
// full read (faSize, then a per-character Python gap scan -- see
// rename_sequences.py's docstring). No renamed FASTA is written at all:
// MINIPROT_INDEX indexes the original genome_fasta, and RENAME_GFF renames
// the seqid column of MINIPROT_ALIGN's GFF afterwards using the lookup this
// process emits.
process RENAME_SEQUENCES {
    tag "${name}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    tuple val(name), path(genome_fasta)
    val min_seq_size
    val min_asm_gap

    output:
    tuple val(name), path("${name}.rename_lookup.tsv"), emit: lookup
    tuple val(name), path("${name}.chrom.sizes"),        emit: chrom_sizes
    tuple val(name), path("${name}.gaps.tsv"),           emit: gaps

    script:
    """
    rename_sequences.py \\
        --fasta ${genome_fasta} \\
        --out_lookup ${name}.rename_lookup.tsv \\
        --out_sizes ${name}.chrom.sizes --min_seq_size ${min_seq_size} \\
        --out_gaps ${name}.gaps.tsv --min_gap ${min_asm_gap}
    """
}
