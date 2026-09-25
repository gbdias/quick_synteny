nextflow.enable.dsl = 2

include { RESOLVE_TAXONOMY }        from './modules/local/resolve_taxonomy.nf'
include { FIND_REFERENCE_ASSEMBLY } from './modules/local/find_reference_assembly.nf'
include { FIND_PROTEOME_ASSEMBLY }  from './modules/local/find_proteome_assembly.nf'
include { DOWNLOAD_GENOME }         from './modules/local/download_assembly.nf'
include { DOWNLOAD_PROTEIN }        from './modules/local/download_assembly.nf'
include { RENAME_SEQUENCES }        from './modules/local/rename_sequences.nf'
include { MINIPROT_INDEX; MINIPROT_ALIGN; RENAME_GFF; SPLIT_GENOME; MINIPROT_ALIGN_CHUNK; MERGE_MINIPROT_GFF } from './modules/local/miniprot_align.nf'
include { BUILD_SYNTENY }           from './modules/local/build_synteny.nf'
include { PYGENOMEVIZ_PLOT }        from './modules/local/pygenomeviz_plot.nf'

def RANKS() {
    return ['species', 'genus', 'family', 'order', 'class', 'phylum']
}

def helpMessage() {
    log.info """
    quick_synteny — taxonomy-guided synteny plotting

    Usage:
      nextflow run main.nf --taxid <TAXID> --assembly <target.fa> [options]

    Required:
      --taxid <int>            NCBI taxid of the target organism. Drives reference/
                                proteome discovery. Not required if both --reference
                                and --proteome are given.
      --assembly <path>        FASTA of the target genome to place on the synteny plot.

    Manual overrides (skip NCBI discovery entirely):
      --reference <path>       Reference genome FASTA to use instead of auto-discovery.
      --proteome <path>        Proteome FASTA to use instead of auto-discovery.

    Discovery behavior:
      --max_rank <rank>        How far to climb the taxonomy ladder before giving up.
                                One of: ${RANKS().join(',')} (default: order).
      --min_seq_size <int>     Minimum sequence length (bp) to include in the synteny
                                plot (default: 500000; pass 0 to disable filtering).
      --exclude_target         Never pick a reference genome or proteome source of the
                                same species as the target, even a chromosome-level one.
                                Off by default: a same-species result is kept -- it's a
                                legitimate outcome (e.g. a second, independently submitted
                                assembly of the same organism, or the only chromosome-level
                                resource available for a sparsely-sequenced group) -- as
                                long as it's chromosome-level or better; only a lower-
                                quality same-species candidate is dropped either way.
      --min_asm_gap <int>       Minimum run of N's (bp) counted as an assembly gap, shown
                                 by the plot's "Show gaps" switch (default: 100).

    Synteny chaining (the interactive page can re-chain with any values; these
    set its starting point and the published links.tsv):
      --min_identity <float>   Per-hit identity (miniprot Positive=, 0-1) below which a
                                hit is ignored. Default: auto-tuned from the weaker
                                genome's mean best-hit identity, clamped to [0.3, 0.9],
                                so a divergent pair isn't forced through a threshold
                                tuned for close relatives. Pass a value to override.
      --max_gap <int>          Max number of anchors skipped between consecutive anchors of
                                a synteny block, on either genome (default: 25). Counted
                                in anchors, not bp, so it means the same thing in an
                                anchor-dense and an anchor-sparse genome.
      --min_block <int>        Min anchors (distinct loci on both genomes) per block in
                                links.tsv. Default: auto -- 15 for close relatives (weaker
                                mean identity >= 0.8), 5 otherwise.

    Resource tuning:
      --miniprot_m <int>       Miniprot k-mer sampling exponent: samples 1/2^INT
                                of genomic k-mers when building its index. Higher
                                values cut peak alignment RAM at some cost to
                                sensitivity -- useful on RAM-constrained hardware
                                (a laptop, a small VM). Default: unset (miniprot's
                                own default applies).
      --miniprot_chunk_gb <float>  Align the proteome against genome chunks of about
                                this many Gb each, in parallel, instead of one
                                whole-genome index -- caps miniprot's peak RAM
                                (~10 GB per Gb of genome) without changing results
                                (equivalent to a single whole-genome run; see
                                modules/local/miniprot_align.nf). Default: unset
                                (one whole-genome alignment, as above).
      --miniprot_gb_per_gb <float>  RAM (GB) requested per Gb of CHUNK when
                                --miniprot_chunk_gb is set, plus a fixed 4 GB
                                overhead; multiplied by 1.5 per retry attempt
                                on a 137/140 (OOM) exit. Default: 11.

    Output:
      --outdir <path>          Output directory (default: results).

    Execution (-profile slurm only):
      --slurm_queue <name>            SLURM queue/partition (default: normal).
      --singularity_cache_dir <path>  Shared Apptainer/Singularity image cache dir.

    Profiles:
      -profile standard         Docker + local executor (default for a laptop/CI).
      -profile slurm             Apptainer + SLURM executor (HPC).
      -profile test              Docker + local executor with tiny resource caps.
    """.stripIndent()
}

def validateParams() {
    if (!RANKS().contains(params.max_rank)) {
        exit 1, "ERROR: --max_rank must be one of ${RANKS().join(',')}, got '${params.max_rank}'"
    }

    if (!params.assembly) {
        exit 1, "ERROR: --assembly <target genome fasta> is required. Run with --help for usage."
    }
    if (!file(params.assembly).exists()) {
        exit 1, "ERROR: --assembly file not found: ${params.assembly}"
    }

    def skipDiscovery = params.reference && params.proteome
    if (!params.taxid && !skipDiscovery) {
        exit 1, "ERROR: --taxid is required unless both --reference and --proteome are given. Run with --help for usage."
    }

    if (params.reference && !file(params.reference).exists()) {
        exit 1, "ERROR: --reference file not found: ${params.reference}"
    }
    if (params.proteome && !file(params.proteome).exists()) {
        exit 1, "ERROR: --proteome file not found: ${params.proteome}"
    }
}

def readSelection(selectionFile) {
    return new groovy.json.JsonSlurper().parse(selectionFile)
}

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    validateParams()

    log.info "quick_synteny: assembly=${params.assembly} taxid=${params.taxid ?: '(skipped, using overrides)'} outdir=${params.outdir}"

    def skipDiscovery = params.reference && params.proteome

    // target is always user-supplied (no auto-discovery for it), so its
    // display name -- shown as a plot subtitle so a viewer can tell which
    // physical input file "target"/"reference" refer to -- is just its
    // filename, known immediately
    target_display_name = Channel.value(file(params.assembly).name)
    // species names for the page's alignment summary -- known only for what
    // NCBI told us about: the target's from --taxid's lineage, a discovered
    // reference genome's/proteome's from its selection. Empty otherwise.
    target_species = Channel.value('')
    reference_species = Channel.value('')
    proteome_species = Channel.value('')
    // proteome discovery prefers the discovered reference species' own
    // proteins when it has an annotated assembly (discovery otherwise ranks
    // by assembly quality, not relatedness); a user-supplied reference
    // genome has no known species, so it keeps the plain ranking
    proteome_prefer_taxid = Channel.value('')

    // ---- reference + proteome accession discovery (or manual override) ----
    reference_fasta = null
    proteome_fasta   = null
    // reference genome's display name is its filename when given manually, or
    // (since a downloaded reference genome is always staged as the generic
    // "genome.fna", which tells a viewer nothing) the NCBI accession it was
    // discovered from
    reference_display_name = null
    // proteome's display name for the stats panel: likewise its filename when
    // given manually, or the accession it was discovered from (its species
    // goes alongside, in proteome_species)
    proteome_display_name = null

    if (params.reference) {
        reference_fasta = Channel.value(file(params.reference))
        reference_display_name = Channel.value(file(params.reference).name)
    }
    if (params.proteome) {
        proteome_fasta = Channel.value(file(params.proteome))
        proteome_display_name = Channel.value(file(params.proteome).name)
    }

    if (!skipDiscovery) {
        lineage_ch = RESOLVE_TAXONOMY(Channel.value(params.taxid))
        target_species = lineage_ch.map { f ->
            def row = f.readLines().collect { it.split('\t') }.find { it[0] == 'species' }
            row && row.size() > 2 ? row[2] : ''
        }

        if (!params.reference) {
            reference_selection = FIND_REFERENCE_ASSEMBLY(lineage_ch, params.max_rank, params.exclude_target).selection
            reference_selection.map { f ->
                def sel = readSelection(f)
                def sameSpeciesNote = sel.same_species_as_target ? ' [SAME SPECIES AS TARGET]' : ''
                log.info "quick_synteny: reference accession=${sel.accession} (${sel.organism_name}) rank=${sel.rank} (${sel.name})${sameSpeciesNote} from ${sel.candidate_count} candidate(s)"
                sel
            }.view()
            reference_fasta = DOWNLOAD_GENOME(reference_selection.map { readSelection(it).accession })
            reference_display_name = reference_selection.map { readSelection(it).accession }
            reference_species = reference_selection.map { readSelection(it).organism_name ?: '' }
            proteome_prefer_taxid = reference_selection.map { readSelection(it).organism_taxid ?: '' }
        }

        if (!params.proteome) {
            prot_selection = FIND_PROTEOME_ASSEMBLY(lineage_ch, params.max_rank, params.exclude_target,
                                                    proteome_prefer_taxid).selection
            prot_selection.map { f ->
                def sel = readSelection(f)
                def sameSpeciesNote = sel.same_species_as_target ? ' [SAME SPECIES AS TARGET]' : ''
                log.info "quick_synteny: proteome accession=${sel.accession} (${sel.organism_name}) rank=${sel.rank} (${sel.name}) annotated=${sel.annotated}${sel.preferred_species ? ' [REFERENCE SPECIES]' : ''}${sameSpeciesNote} from ${sel.candidate_count} candidate(s)"
                sel
            }.view()
            proteome_fasta = DOWNLOAD_PROTEIN(prot_selection.map { readSelection(it).accession })
            proteome_display_name = prot_selection.map { readSelection(it).accession }
            proteome_species = prot_selection.map { readSelection(it).organism_name ?: '' }
        }
    }

    // Each per-genome step below (rename+sizes+gaps, align, rename the
    // alignment output) is a single process definition invoked ONCE on a
    // channel carrying both the target and the reference genome tagged by
    // role name ('target'/'reference') -- Nextflow DSL2 does not allow
    // calling the same process twice in one workflow scope, so
    // target+reference are mixed into one channel and split back apart
    // with .branch{} wherever a downstream step needs them as two separate
    // arguments.

    // ---- rename sequences (both genomes) -- the only full read of each
    // genome; chrom sizes and assembly gaps come out of the same pass
    // instead of two more full reads (faSize, then a separate gap scan) --
    // see rename_sequences.nf/.py ----
    genomes_in = Channel.value('target').combine(Channel.value(file(params.assembly)))
        .mix(Channel.value('reference').combine(reference_fasta))

    renamed = RENAME_SEQUENCES(genomes_in, params.min_seq_size, params.min_asm_gap)

    chrom_sizes_by_role = renamed.chrom_sizes.branch {
        target: it[0] == 'target'
        reference: it[0] == 'reference'
    }
    target_chrom_sizes     = chrom_sizes_by_role.target
    reference_chrom_sizes = chrom_sizes_by_role.reference

    gaps_by_role = renamed.gaps.branch {
        target: it[0] == 'target'
        reference: it[0] == 'reference'
    }
    target_gaps     = gaps_by_role.target
    reference_gaps = gaps_by_role.reference

    // ---- align proteome against both ORIGINAL genomes (also doubles as
    // the raw synteny/homeolog anchor source -- see build_synteny.nf), then
    // rename each resulting GFF's seqid column using the lookup from the
    // matching RENAME_SEQUENCES call (join()'d by role name) ----
    // params.miniprot_m defaults to null (miniprot's own default); normalized
    // to '' here for the same reason as min_identity below -- see that comment
    def miniprot_m = params.miniprot_m ?: ''

    if (params.miniprot_chunk_gb) {
        // Chunked path: split each genome into ~miniprot_chunk_gb-sized
        // pieces, align each in its own (smaller-index, lower-RAM) miniprot
        // task, then merge back into a GFF equivalent to a whole-genome run
        // -- see modules/local/miniprot_align.nf.
        def target_chunk_bp = Math.round(params.miniprot_chunk_gb * 1_000_000_000)

        split_ch = SPLIT_GENOME(genomes_in, target_chunk_bp)

        // one row per (name, chunk fasta), carrying that chunk's own length
        // (dynamic memory directive) and the WHOLE genome's length (the -G
        // formula, computed from the genome, not the chunk -- see
        // MINIPROT_ALIGN_CHUNK)
        chunk_rows = split_ch.chunks.flatMap { name, chunk_fastas, chunks_tsv, total_length_file ->
            def total_length = total_length_file.text.trim() as long
            def length_by_file = [:]
            chunks_tsv.readLines().each { line ->
                def (fname, len) = line.split('\t')
                length_by_file[fname] = len as long
            }
            def fastas = chunk_fastas instanceof List ? chunk_fastas : [chunk_fastas]
            fastas.collect { fa -> [name, fa, length_by_file[fa.name], total_length] }
        }
        n_chunks_by_name = split_ch.chunks.map { name, chunk_fastas, chunks_tsv, total_length_file ->
            def fastas = chunk_fastas instanceof List ? chunk_fastas : [chunk_fastas]
            [name, fastas.size()]
        }

        align_chunk_in = chunk_rows.combine(proteome_fasta)
        raw_gff_chunks = MINIPROT_ALIGN_CHUNK(align_chunk_in, miniprot_m).gff

        // groupTuple() with no `size:` waits for the (finite, known-size)
        // upstream channel to close before emitting each group, so this
        // already gathers every chunk regardless of target/reference
        // having different chunk counts; the join+check below is purely a
        // fail-fast guard against a silently-missing chunk task.
        grouped_gff = raw_gff_chunks.groupTuple().join(n_chunks_by_name).map { name, gffs, n_chunks ->
            if (gffs.size() != n_chunks) {
                error "quick_synteny: ${name} produced ${gffs.size()} chunk GFF(s), expected ${n_chunks} -- a MINIPROT_ALIGN_CHUNK task must be missing"
            }
            [name, gffs]
        }
        raw_gff = MERGE_MINIPROT_GFF(grouped_gff).gff
    } else {
        genome_index = MINIPROT_INDEX(genomes_in, miniprot_m).index
        align_in = genome_index.combine(proteome_fasta)
        raw_gff = MINIPROT_ALIGN(align_in).gff
    }

    gff = RENAME_GFF(raw_gff.join(renamed.lookup)).gff
    gff_by_role = gff.branch {
        target: it[0] == 'target'
        reference: it[0] == 'reference'
    }
    target_gff     = gff_by_role.target
    reference_gff = gff_by_role.reference

    // ---- synteny blocks + final plot ----
    // params.min_identity defaults to null (auto-tune); Groovy's null is not
    // safe to pass through a process `val` input the same way every call site
    // expects (falsy-but-interpolatable), so it's normalized to '' here once
    def min_identity = params.min_identity ?: ''
    def min_block = params.min_block ?: ''
    // what the page's alignment summary names each input by -- (target
    // source, target species, reference source, reference species,
    // proteome source, proteome species); combine() flattens into one tuple
    input_sources = target_display_name
        .combine(target_species)
        .combine(reference_display_name)
        .combine(reference_species)
        .combine(proteome_display_name)
        .combine(proteome_species)
    synteny = BUILD_SYNTENY(target_gff, reference_gff, target_chrom_sizes, reference_chrom_sizes,
                             proteome_fasta, min_identity, params.max_gap, min_block,
                             input_sources)

    PYGENOMEVIZ_PLOT(
        synteny.hits,
        target_chrom_sizes, reference_chrom_sizes,
        min_identity, params.max_gap, min_block,
        synteny.stats,
        target_display_name, reference_display_name,
        target_gaps, reference_gaps,
    )
}
