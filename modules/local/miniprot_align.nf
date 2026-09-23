// Generic (genome, proteome) -> gff aligner, called twice from main.nf
// (target, comparison) rather than as two separate modules.
//
// --outs/-N are loosened from miniprot's stock defaults (0.99 / 30) so
// secondary alignments -- a protein's OTHER hits, at a lower score than its
// best one -- are actually retained in the output. build_synteny.nf uses
// exactly these secondary hits as its raw synteny/homeolog anchors,
// so this alignment doubles as the "ortholog search" a separate aligner
// would otherwise be needed for.
//
// A single miniprot invocation, indexing the genome on the fly instead of
// with a separate `-d` step: the old two-call form (`-d` to build a .mpi
// index file, then align against it) wrote that index to disk only to read
// it straight back for the one alignment that ever uses it -- on axolotl,
// 72 GB written and 100 GB read for no benefit. All indexing options
// (including -M) still apply here; miniprot accepts them on the alignment
// line and builds the index in memory before aligning.
//
// Runs on the ORIGINAL genome fasta, not a renamed copy (RENAME_SEQUENCES
// no longer produces one -- see that module and rename_sequences.py), so
// the GFF this emits still carries the genome's original sequence names;
// RENAME_GFF swaps those for the renamed IDs right after.
//
// -M (miniprot_m, optional) trades sensitivity for RAM: it samples 1/2^M of
// genomic k-mers when building the index. Unset by default (miniprot's own
// default applies); see main.nf --miniprot_m and benchmark/miniprot_m_sweep/
// for the measured RAM-vs-concordance tradeoff this buys on RAM-constrained
// hardware.

process MINIPROT_ALIGN {
    tag "${name}"
    label 'process_high'
    container 'quay.io/biocontainers/miniprot:0.18--h577a1d6_0'

    input:
    tuple val(name), path(genome_fasta), path(proteome_faa)
    val miniprot_m   // '' keeps miniprot's own default; otherwise an int k-mer sampling exponent

    output:
    tuple val(name), path("${name}.raw.gff"), emit: gff

    script:
    def mFlag = miniprot_m ? "-M ${miniprot_m}" : ''
    """
    miniprot -t ${task.cpus} ${mFlag} -I -N 5 --outs=0.7 --gff ${genome_fasta} ${proteome_faa} > ${name}.raw.gff
    """
}

// Renames the aligner's GFF seqid column (original genome IDs -> chr1/chr2/
// ... etc) using the same lookup RENAME_SEQUENCES produced from the same
// genome -- see rename_gff.py. Plain Python, not a bioinformatics tool, so
// it runs in the generic python container like RENAME_SEQUENCES itself.
process RENAME_GFF {
    tag "${name}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'

    input:
    tuple val(name), path(raw_gff), path(lookup)

    output:
    tuple val(name), path("${name}.gff"), emit: gff

    script:
    """
    rename_gff.py --gff ${raw_gff} --lookup ${lookup} --out ${name}.gff
    """
}
