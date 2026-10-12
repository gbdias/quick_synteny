// One read of the genome: plot-friendly names (chr1, chr2, ...), chromosome
// sizes and assembly gaps. No renamed FASTA is written: miniprot indexes
// the original, and RENAME_GFF applies this lookup to its GFF.
process RENAME_SEQUENCES {
    tag "${name}"
    label 'process_low'
    label 'env_python'
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
