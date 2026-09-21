// Generic (genome, proteome) -> gff aligner, called twice from main.nf
// (target, comparison) rather than as two separate modules.
//
// --outs/-N are loosened from miniprot's stock defaults (0.99 / 30) so
// secondary alignments -- a protein's OTHER hits, at a lower score than its
// best one -- are actually retained in the output. build_synteny_blocks.nf
// uses exactly these secondary hits as its raw synteny/homeolog anchors,
// so this alignment doubles as the "ortholog search" a separate aligner
// would otherwise be needed for.
//
// -M (miniprot_m, optional) trades sensitivity for RAM: it samples 1/2^M of
// genomic k-mers when building the index, so it belongs on the -d line, not
// the alignment line -- the sampling is baked into the index at build time,
// not redone per query. Unset by default (miniprot's own default applies);
// see main.nf --miniprot_m and benchmark/miniprot_m_sweep/ for the measured
// RAM-vs-concordance tradeoff this buys on RAM-constrained hardware.

process MINIPROT_ALIGN {
    tag "${name}"
    label 'process_high'
    container 'quay.io/biocontainers/miniprot:0.18--h577a1d6_0'

    input:
    tuple val(name), path(genome_fasta), path(proteome_faa)
    val miniprot_m   // '' keeps miniprot's own default; otherwise an int k-mer sampling exponent

    output:
    tuple val(name), path("${name}.gff"), emit: gff

    script:
    def mFlag = miniprot_m ? "-M ${miniprot_m}" : ''
    """
    miniprot -t ${task.cpus} ${mFlag} -d ${name}.mpi ${genome_fasta}
    miniprot -t ${task.cpus} -I -N 5 --outs=0.7 --gff ${name}.mpi ${proteome_faa} > ${name}.gff
    """
}
