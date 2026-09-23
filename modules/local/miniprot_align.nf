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

// ---------------------------------------------------------------------------
// Chunked path (params.miniprot_chunk_gb, opt-in): SPLIT_GENOME ->
// MINIPROT_ALIGN_CHUNK (one task per chunk) -> MERGE_MINIPROT_GFF, wired in
// main.nf to emit the same tuple(name, gff) shape MINIPROT_ALIGN does, so
// RENAME_GFF and everything after it is unchanged either way.
//
// The point is peak RAM: miniprot's index is ~10 GB per Gb of genome, so
// aligning against N smaller chunks instead of one whole-genome index caps
// peak RSS per task at roughly (chunk size / N) instead of the whole
// genome. Per-chunk outputs are a
// superset of the whole-genome output (same --outs/-N filters, applied to
// a smaller index), so MERGE_MINIPROT_GFF re-applies both filters globally
// to recover a result equivalent to one whole-genome run -- see that
// script's docstring for why this is exact, and its two gotchas (an
// explicit -G computed from the WHOLE genome's length, since -I would
// otherwise derive max intron size from a chunk's own, smaller length;
// and per-chunk k-mer statistics that can shift sensitivity slightly,
// which is what the plan's acceptance checks measure).
// ---------------------------------------------------------------------------

process SPLIT_GENOME {
    tag "${name}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'

    input:
    tuple val(name), path(genome_fasta)
    val target_chunk_bp

    output:
    tuple val(name), path("chunk_*.fa"), path("chunks.tsv"), path("total_length.txt"), emit: chunks

    script:
    """
    split_genome.py --fasta ${genome_fasta} --chunk_bp ${target_chunk_bp} --outdir .
    """
}

// cpus/memory/errorStrategy/maxRetries are set directly on the process
// (not via a label) because the memory formula needs `chunk_bp`, a process
// input -- a dynamic directive in nextflow.config can only see task.* and
// params.*, not a process's own named inputs, so this can't be split across
// a label (whose memory Nextflow config would otherwise win over anything
// set here). cpus is lower than MINIPROT_ALIGN's process_high default
// because chunking already spreads the alignment across N parallel tasks
// instead of one -- the original single-job benchmark ran 32 threads at
// only 392% CPU (APPLICATION_NOTE_PLAN.md Sec 2), so a handful of threads
// per chunk is already proportionate.
process MINIPROT_ALIGN_CHUNK {
    tag "${name}:${chunk_fasta.baseName}"
    container 'quay.io/biocontainers/miniprot:0.18--h577a1d6_0'
    cpus 4
    memory { "${Math.ceil((chunk_bp / 1e9 * params.miniprot_gb_per_gb + 4) * Math.pow(1.5, task.attempt - 1))} GB" }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(name), path(chunk_fasta), val(chunk_bp), val(total_length), path(proteome_faa)
    val miniprot_m

    output:
    tuple val(name), path("${chunk_fasta.baseName}.raw.gff"), emit: gff

    script:
    def mFlag = miniprot_m ? "-M ${miniprot_m}" : ''
    // -I derives max intron size from the length of the INPUT genome; a
    // chunk is shorter than the whole genome, so -G is computed here from
    // the whole genome's length instead, replicating miniprot's own
    // formula (mp_mapopt_set_max_intron, confirmed empirically against -I's
    // own logged value: it's a ceiling, not a round or floor).
    def G = Math.min(Math.max(Math.ceil(3.6 * Math.sqrt(total_length as double)) as long, 10000L), 300000L)
    """
    miniprot -t ${task.cpus} ${mFlag} -G ${G} -N 5 --outs=0.7 --gff ${chunk_fasta} ${proteome_faa} > ${chunk_fasta.baseName}.raw.gff
    """
}

process MERGE_MINIPROT_GFF {
    tag "${name}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'

    input:
    tuple val(name), path(chunk_gffs)

    output:
    tuple val(name), path("${name}.raw.gff"), emit: gff

    script:
    """
    merge_miniprot_gff.py --gff ${chunk_gffs} --out ${name}.raw.gff
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
