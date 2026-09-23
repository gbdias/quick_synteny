nextflow.enable.dsl = 2

include { RESOLVE_TAXONOMY }        from './modules/local/resolve_taxonomy.nf'
include { FIND_COMPARISON_ASSEMBLY } from './modules/local/find_comparison_assembly.nf'
include { FIND_PROTEOME_ASSEMBLY }  from './modules/local/find_proteome_assembly.nf'
include { DOWNLOAD_GENOME }         from './modules/local/download_assembly.nf'
include { DOWNLOAD_PROTEIN }        from './modules/local/download_assembly.nf'
include { RENAME_SEQUENCES }        from './modules/local/rename_sequences.nf'
include { MINIPROT_ALIGN; RENAME_GFF } from './modules/local/miniprot_align.nf'
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
      --taxid <int>            NCBI taxid of the target organism. Drives comparison/
                                proteome discovery. Not required if both --comparison
                                and --proteome are given.
      --assembly <path>        FASTA of the target genome to place on the synteny plot.

    Manual overrides (skip NCBI discovery entirely):
      --comparison <path>      Comparison genome FASTA to use instead of auto-discovery.
      --proteome <path>        Proteome FASTA to use instead of auto-discovery.

    Discovery behavior:
      --max_rank <rank>        How far to climb the taxonomy ladder before giving up.
                                One of: ${RANKS().join(',')} (default: order).
      --min_seq_size <int>     Minimum sequence length (bp) to include in the synteny
                                plot (default: 500000; pass 0 to disable filtering).
      --exclude_target         Never pick a comparison genome or proteome source of the
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
      --max_gap <int>          Max number of genes skipped between consecutive genes of
                                a synteny block, on either genome (default: 25). Counted
                                in genes, not bp, so it means the same thing in a
                                gene-dense and a gene-sparse genome.
      --min_block <int>        Min genes (distinct loci on both genomes) per block in
                                links.tsv. Default: auto -- 15 for close relatives (weaker
                                mean identity >= 0.8), 5 otherwise.

    Resource tuning:
      --miniprot_m <int>       Miniprot k-mer sampling exponent: samples 1/2^INT
                                of genomic k-mers when building its index. Higher
                                values cut peak alignment RAM at some cost to
                                sensitivity -- useful on RAM-constrained hardware
                                (a laptop, a small VM). Default: unset (miniprot's
                                own default applies). See benchmark/miniprot_m_sweep/
                                for the measured RAM-vs-concordance tradeoff.

    Polyploid support:
      --show_homeologs <mode>  Also self-compare one or both genomes' own proteome
                                hits to find and draw homeologous chromosome pairs
                                (e.g. for an allopolyploid genome). One of: target,
                                comparison, both. The synteny/homeolog detection
                                itself needs no ploidy ratio -- it's built directly
                                from miniprot's own multi-hit alignments, so it
                                picks up whatever multiplicity is actually in the
                                data. Default: both. Pass '' to turn it off entirely.

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

    def skipDiscovery = params.comparison && params.proteome
    if (!params.taxid && !skipDiscovery) {
        exit 1, "ERROR: --taxid is required unless both --comparison and --proteome are given. Run with --help for usage."
    }

    if (params.comparison && !file(params.comparison).exists()) {
        exit 1, "ERROR: --comparison file not found: ${params.comparison}"
    }
    if (params.proteome && !file(params.proteome).exists()) {
        exit 1, "ERROR: --proteome file not found: ${params.proteome}"
    }

    if (params.show_homeologs && !['target', 'comparison', 'both'].contains(params.show_homeologs)) {
        exit 1, "ERROR: --show_homeologs must be one of target,comparison,both, got '${params.show_homeologs}'"
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

    def skipDiscovery = params.comparison && params.proteome

    // target is always user-supplied (no auto-discovery for it), so its
    // display name -- shown as a plot subtitle so a viewer can tell which
    // physical input file "target"/"comparison" refer to -- is just its
    // filename, known immediately
    target_display_name = Channel.value(file(params.assembly).name)

    // ---- comparison + proteome accession discovery (or manual override) ----
    comparison_fasta = null
    proteome_fasta   = null
    // comparison genome's display name is its filename when given manually, or
    // (since a downloaded comparison genome is always staged as the generic
    // "genome.fna", which tells a viewer nothing) the NCBI accession it was
    // discovered from
    comparison_display_name = null
    // proteome's origin for the stats panel: its filename when given manually,
    // or (species name is more useful than an accession here -- the proteome
    // is background context, not one of the two genomes being compared) the
    // species name it was discovered from
    proteome_display_name = null

    if (params.comparison) {
        comparison_fasta = Channel.value(file(params.comparison))
        comparison_display_name = Channel.value(file(params.comparison).name)
    }
    if (params.proteome) {
        proteome_fasta = Channel.value(file(params.proteome))
        proteome_display_name = Channel.value(file(params.proteome).name)
    }

    if (!skipDiscovery) {
        lineage_ch = RESOLVE_TAXONOMY(Channel.value(params.taxid))

        if (!params.comparison) {
            comparison_selection = FIND_COMPARISON_ASSEMBLY(lineage_ch, params.max_rank, params.exclude_target).selection
            comparison_selection.map { f ->
                def sel = readSelection(f)
                def sameSpeciesNote = sel.same_species_as_target ? ' [SAME SPECIES AS TARGET]' : ''
                log.info "quick_synteny: comparison accession=${sel.accession} rank=${sel.rank} (${sel.name})${sameSpeciesNote} from ${sel.candidate_count} candidate(s)"
                sel
            }.view()
            comparison_fasta = DOWNLOAD_GENOME(comparison_selection.map { readSelection(it).accession })
            comparison_display_name = comparison_selection.map { readSelection(it).accession }
        }

        if (!params.proteome) {
            prot_selection = FIND_PROTEOME_ASSEMBLY(lineage_ch, params.max_rank, params.exclude_target).selection
            prot_selection.map { f ->
                def sel = readSelection(f)
                def sameSpeciesNote = sel.same_species_as_target ? ' [SAME SPECIES AS TARGET]' : ''
                log.info "quick_synteny: proteome accession=${sel.accession} rank=${sel.rank} (${sel.name}) annotated=${sel.annotated}${sameSpeciesNote} from ${sel.candidate_count} candidate(s)"
                sel
            }.view()
            proteome_fasta = DOWNLOAD_PROTEIN(prot_selection.map { readSelection(it).accession })
            proteome_display_name = prot_selection.map { readSelection(it).name }
        }
    }

    // Each per-genome step below (rename+sizes+gaps, align, rename the
    // alignment output) is a single process definition invoked ONCE on a
    // channel carrying both the target and the comparison genome tagged by
    // role name ('target'/'comparison') -- Nextflow DSL2 does not allow
    // calling the same process twice in one workflow scope, so
    // target+comparison are mixed into one channel and split back apart
    // with .branch{} wherever a downstream step needs them as two separate
    // arguments.

    // ---- rename sequences (both genomes) -- the only full read of each
    // genome; chrom sizes and assembly gaps come out of the same pass
    // instead of two more full reads (faSize, then a separate gap scan) --
    // see rename_sequences.nf/.py ----
    genomes_in = Channel.value('target').combine(Channel.value(file(params.assembly)))
        .mix(Channel.value('comparison').combine(comparison_fasta))

    renamed = RENAME_SEQUENCES(genomes_in, params.min_seq_size, params.min_asm_gap)

    chrom_sizes_by_role = renamed.chrom_sizes.branch {
        target: it[0] == 'target'
        comparison: it[0] == 'comparison'
    }
    target_chrom_sizes     = chrom_sizes_by_role.target
    comparison_chrom_sizes = chrom_sizes_by_role.comparison

    gaps_by_role = renamed.gaps.branch {
        target: it[0] == 'target'
        comparison: it[0] == 'comparison'
    }
    target_gaps     = gaps_by_role.target
    comparison_gaps = gaps_by_role.comparison

    // ---- align proteome against both ORIGINAL genomes (also doubles as
    // the raw synteny/homeolog anchor source -- see build_synteny.nf), then
    // rename each resulting GFF's seqid column using the lookup from the
    // matching RENAME_SEQUENCES call (join()'d by role name) ----
    align_in = genomes_in.combine(proteome_fasta)
    // params.miniprot_m defaults to null (miniprot's own default); normalized
    // to '' here for the same reason as min_identity below -- see that comment
    def miniprot_m = params.miniprot_m ?: ''
    raw_gff = MINIPROT_ALIGN(align_in, miniprot_m).gff
    gff = RENAME_GFF(raw_gff.join(renamed.lookup)).gff
    gff_by_role = gff.branch {
        target: it[0] == 'target'
        comparison: it[0] == 'comparison'
    }
    target_gff     = gff_by_role.target
    comparison_gff = gff_by_role.comparison

    // ---- synteny blocks + final plot ----
    // params.min_identity defaults to null (auto-tune); Groovy's null is not
    // safe to pass through a process `val` input the same way every call site
    // expects (falsy-but-interpolatable), so it's normalized to '' here once
    def min_identity = params.min_identity ?: ''
    def min_block = params.min_block ?: ''
    synteny = BUILD_SYNTENY(target_gff, comparison_gff, target_chrom_sizes, comparison_chrom_sizes,
                             proteome_fasta, params.show_homeologs, min_identity, params.max_gap, min_block,
                             proteome_display_name)

    PYGENOMEVIZ_PLOT(
        synteny.hits,
        target_chrom_sizes, comparison_chrom_sizes,
        min_identity, params.max_gap, min_block,
        synteny.stats,
        target_display_name, comparison_display_name,
        target_gaps, comparison_gaps,
    )
}
