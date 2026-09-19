process RENAME_SEQUENCES {
    tag "${name}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: '*.rename_lookup.tsv'

    input:
    tuple val(name), path(genome_fasta)

    output:
    tuple val(name), path("${name}.renamed.fa"), emit: fasta
    path "${name}.rename_lookup.tsv",            emit: lookup

    script:
    """
    rename_sequences.py \\
        --fasta ${genome_fasta} \\
        --out_fasta ${name}.renamed.fa \\
        --out_lookup ${name}.rename_lookup.tsv
    """
}
