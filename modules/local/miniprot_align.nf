// MINIPROT_INDEX indexes a genome, then MINIPROT_ALIGN aligns the proteome
// against it.
//
// -N/--outs are loosened from miniprot's defaults (30 / 0.99) so secondary
// hits are kept: build_synteny.nf uses them as its synteny/homeolog anchors.
//
// Indexing is its own task so it can run on one thread: each index-builder
// thread decodes a whole sequence strand into its own buffer (build_worker
// in miniprot's index.c), so peak RAM grows with threads -- 1 instead of 8
// cut it by 36% on TAIR10 -- while the index is byte-identical. The price
// is a disk round trip for the index (72 GB on axolotl).
//
// The GFF carries the genome's original sequence names; RENAME_GFF swaps
// in the renamed IDs.
//
// -M and -L belong on the -d line: they decide what goes into the index,
// and miniprot ignores them when aligning against a prebuilt one. The
// pipeline's defaults (-M 2 -L 45, against miniprot's 1 and 30) cut the
// index's peak RAM by about half; on rice and on two beetles 80 My apart,
// large blocks were unchanged and small ones (homeologs most) lost a few
// anchors.

// '' keeps miniprot's own default; 0 is a real value (-M 0 indexes every k-mer)
def indexFlags(miniprot_m, miniprot_l) {
    [miniprot_m != '' ? "-M ${miniprot_m}" : '', miniprot_l != '' ? "-L ${miniprot_l}" : ''].findAll().join(' ')
}

// process_high sets memory and time; cpus is pinned to 1 in nextflow.config
process MINIPROT_INDEX {
    tag "${name}"
    label 'process_high'
    label 'env_miniprot'

    input:
    tuple val(name), path(genome_fasta)
    val miniprot_m   // k-mer sampling exponent (-M), or ''
    val miniprot_l   // shortest ORF to index in aa (-L), or ''

    output:
    tuple val(name), path("${name}.mpi"), emit: index

    script:
    """
    miniprot -t ${task.cpus} ${indexFlags(miniprot_m, miniprot_l)} -d ${name}.mpi ${genome_fasta}
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
// Chunked path (--miniprot_chunk_gb): SPLIT_GENOME -> MINIPROT_ALIGN_CHUNK
// (one task per chunk) -> MERGE_MINIPROT_GFF, emitting the same
// tuple(name, gff) as MINIPROT_ALIGN.
//
// The point is peak RAM: miniprot's index is ~3.5 GB per Gb of genome at
// -M 2 -L 45 (~7.5 at miniprot's own -M 1 -L 30), so each task holds one
// chunk's index. MERGE_MINIPROT_GFF re-applies --outs/-N across chunks, and
// -G comes from the whole genome's length.
//
// Close to, but NOT identical to, a whole-genome run, and can't be without
// changing miniprot: it decides per index (1) which k-mers are too
// repetitive to seed from (mp_cal_max_occ, which -c can only lower) and
// (2) which candidates survive -p/-N before DP. So a chunk can gain or lose
// seeds at a locus and keep candidates a whole-genome run pruned. On
// C. elegans vs C. portoensis (C. brenneri proteome, 3-6 chunks): the best
// hit matched for ~98% of proteins and ~97% of blocks; ~7-10% of proteins
// differed somewhere in their hit set (weak hits, rank 4-6 secondaries,
// tandem arrays, score ties). More chunks drift further and cost more CPU
// (every chunk aligns the whole proteome), so use the largest chunk that
// fits. Loosening -N/-p/--outs per chunk moves further from whole-genome.
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

// RAM (GB) to request per Gb of chunk, on top of a fixed 4 GB for the
// proteome's side (a few huge isoforms can need several GB each):
// --miniprot_gb_per_gb, or else by -M. miniprot's peak per Gb of chunk,
// 4 threads building the chunk's index then aligning: ~5.5 GB at -M 2 -L 45
// (0.1-0.5 Gb chunks of rice, O. australiensis and two beetles; smaller
// chunks add a fixed ~0.5 GB) and ~12.5 at -M 1 -L 30, which 11 + 4 still
// covers up to ~2.5 Gb chunks. -M 0 is unmeasured: twice -M 1's.
def chunkGbPerGb() {
    if (params.miniprot_gb_per_gb != null) return params.miniprot_gb_per_gb
    def m = params.miniprot_m != null ? params.miniprot_m as int : 1
    m <= 0 ? 22 : m == 1 ? 11 : 6
}

// No label: the memory formula needs `chunk_bp`, a process input that a
// config directive can't see, and a label's config memory would override
// it. 4 cpus: the chunks already run in parallel, and a whole-genome job
// given 32 threads used only ~4 cores' worth.
process MINIPROT_ALIGN_CHUNK {
    tag "${name}:${chunk_fasta.baseName}"
    label 'env_miniprot'
    cpus 4
    memory { "${Math.ceil((chunk_bp / 1e9 * chunkGbPerGb() + 4) * Math.pow(1.5, task.attempt - 1))} GB" }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(name), path(chunk_fasta), val(chunk_bp), val(total_length), path(proteome_faa)
    val miniprot_m
    val miniprot_l

    output:
    tuple val(name), path("${chunk_fasta.baseName}.raw.gff"), emit: gff

    script:
    // -I would size max introns from the chunk's length: replicate its
    // formula (mp_mapopt_set_max_intron, a ceiling) on the whole genome's
    def G = Math.min(Math.max(Math.ceil(3.6 * Math.sqrt(total_length as double)) as long, 10000L), 300000L)
    """
    miniprot -t ${task.cpus} ${indexFlags(miniprot_m, miniprot_l)} -G ${G} -N 5 --outs=0.7 --gff ${chunk_fasta} ${proteome_faa} > ${chunk_fasta.baseName}.raw.gff
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

// original sequence IDs -> RENAME_SEQUENCES' names (chr1, chr2, ...)
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
