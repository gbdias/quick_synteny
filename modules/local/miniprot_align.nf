// Generic (genome, proteome) -> gff aligner, called twice from main.nf
// (target, reference) rather than as two separate modules: MINIPROT_INDEX
// indexes the genome, then MINIPROT_ALIGN aligns the proteome against it.
//
// --outs/-N are loosened from miniprot's stock defaults (0.99 / 30) so
// secondary alignments -- a protein's OTHER hits, at a lower score than its
// best one -- are actually retained in the output. build_synteny.nf uses
// exactly these secondary hits as its raw synteny/homeolog anchors,
// so this alignment doubles as the "ortholog search" a separate aligner
// would otherwise be needed for.
//
// Indexing is its own task so it can run on one thread while the alignment
// runs on many. miniprot's index builder gives each thread a sequence strand
// to decode into a buffer of its own and collect k-mers from (build_worker
// in miniprot's index.c), so the build's peak RAM grows by a decoded
// chromosome strand plus its k-mers per thread -- 1 thread instead of 8 cut
// it by 36% on TAIR10 -- while the index itself is byte-identical at any
// thread count. The price is a disk round trip for the index (72 GB on
// axolotl) that indexing on the fly inside the alignment call would avoid.
//
// MINIPROT_INDEX reads the ORIGINAL genome fasta, not a renamed copy
// (RENAME_SEQUENCES no longer produces one -- see that module and
// rename_sequences.py), so the GFF MINIPROT_ALIGN emits still carries the
// genome's original sequence names; RENAME_GFF swaps those for the renamed
// IDs right after.
//
// -M (miniprot_m, optional) trades sensitivity for RAM: it samples 1/2^M of
// genomic k-mers when building the index, so it belongs on the -d line --
// the sampling is baked into the index, and miniprot ignores -M when
// aligning against a prebuilt one. Unset by default (miniprot's own default
// applies); see main.nf --miniprot_m.

// process_high supplies memory and time; cpus is pinned to 1 in
// nextflow.config.
process MINIPROT_INDEX {
    tag "${name}"
    label 'process_high'
    label 'env_miniprot'

    input:
    tuple val(name), path(genome_fasta)
    val miniprot_m   // '' keeps miniprot's own default; otherwise an int k-mer sampling exponent

    output:
    tuple val(name), path("${name}.mpi"), emit: index

    script:
    def mFlag = miniprot_m ? "-M ${miniprot_m}" : ''
    """
    miniprot -t ${task.cpus} ${mFlag} -d ${name}.mpi ${genome_fasta}
    """
}

process MINIPROT_ALIGN {
    tag "${name}"
    label 'process_high'
    label 'env_miniprot'

    input:
    tuple val(name), path(genome_index), path(proteome_faa)

    output:
    tuple val(name), path("${name}.raw.gff"), emit: gff

    script:
    """
    miniprot -t ${task.cpus} -I -N 5 --outs=0.7 --gff ${genome_index} ${proteome_faa} > ${name}.raw.gff
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
// and per-chunk k-mer statistics that can shift sensitivity slightly).
// ---------------------------------------------------------------------------

process SPLIT_GENOME {
    tag "${name}"
    label 'process_low'
    label 'env_python'

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
// instead of one -- a single whole-genome job given 32 threads used only
// ~4 cores' worth, so a handful of threads per chunk is already proportionate.
process MINIPROT_ALIGN_CHUNK {
    tag "${name}:${chunk_fasta.baseName}"
    label 'env_miniprot'
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
    label 'env_python'

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
    label 'env_python'

    input:
    tuple val(name), path(raw_gff), path(lookup)

    output:
    tuple val(name), path("${name}.gff"), emit: gff

    script:
    """
    rename_gff.py --gff ${raw_gff} --lookup ${lookup} --out ${name}.gff
    """
}
