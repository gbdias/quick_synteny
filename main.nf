include { validateParameters } from 'plugin/nf-schema'

include { RESOLVE_TAXONOMY }        from './modules/local/resolve_taxonomy.nf'
include { RESOLVE_TAXONOMY as RESOLVE_REFERENCE_TAXONOMY } from './modules/local/resolve_taxonomy.nf'
include { FIND_REFERENCE_ASSEMBLY; FETCH_ACCESSION_SUMMARY; DESCRIBE_REFERENCE_ACCESSION } from './modules/local/find_reference_assembly.nf'
include { FIND_PROTEOME_ASSEMBLY }  from './modules/local/find_proteome_assembly.nf'
include { DOWNLOAD_GENOME }         from './modules/local/download_assembly.nf'
include { DOWNLOAD_PROTEIN }        from './modules/local/download_assembly.nf'
include { RENAME_SEQUENCES }        from './modules/local/rename_sequences.nf'
include { MINIPROT_INDEX; MINIPROT_ALIGN; RENAME_GFF; SPLIT_GENOME; MINIPROT_ALIGN_CHUNK; MERGE_MINIPROT_GFF } from './modules/local/miniprot_align.nf'
include { BUILD_SYNTENY }           from './modules/local/build_synteny.nf'
include { PLOT_SYNTENY }            from './modules/local/plot_synteny.nf'

// chunk size for a bare --miniprot_chunk_gb
def DEFAULT_CHUNK_GB() {
    return 1.0
}

// --miniprot_chunk_gb -> chunk size in Gb, or null for the whole-genome
// path. A bare flag arrives as boolean true.
def chunkGb() {
    def v = params.miniprot_chunk_gb
    if (v == null || v == false) return null
    if (v == true) return DEFAULT_CHUNK_GB()
    return v as double
}

// The reference is --reference, --reference_accession, --reference_taxid,
// or else discovered outward from --taxid; the proteome is --proteome, or
// else discovered from --taxid's relatives. --taxid is needed whenever
// something is discovered from it.
def referenceGiven() {
    return params.reference || params.reference_accession || params.reference_taxid
}

def needsTargetTaxid() {
    return !referenceGiven() || !params.proteome
}

// What nextflow_schema.json can't say: which parameters go together.
// Everything else is checked by validateParameters() against the schema.
def validateParams() {
    def references = ['reference', 'reference_accession', 'reference_taxid'].findAll { params[it] }
    if (references.size() > 1) {
        exit 1, "ERROR: give only one of --reference, --reference_accession and --reference_taxid, got ${references.collect { '--' + it }.join(' and ')}."
    }
    if (params.taxid == null && needsTargetTaxid()) {
        exit 1, "ERROR: --taxid is required unless --proteome is given along with a reference genome (--reference, --reference_accession or --reference_taxid). Run with --help for usage."
    }
}

def readSelection(selectionFile) {
    return new groovy.json.JsonSlurper().parse(selectionFile)
}

workflow {
    // nf-schema prints --help/--helpFull but doesn't stop the run: without
    // the exit, validation would then fail on the missing --assembly
    if (params.help || params.containsKey('helpFull') && params.helpFull) {
        exit 0
    }
    validateParameters()
    validateParams()

    log.info "quick_synteny: assembly=${params.assembly} taxid=${params.taxid ?: '(skipped, using overrides)'} outdir=${params.outdir}"

    // what the page names each input by: its filename if user-supplied, else
    // the NCBI accession it came from (a download is always "genome.fna")
    target_display_name = Channel.value(file(params.assembly).name)
    // species names, known only for what NCBI told us about; '' otherwise
    target_species = Channel.value('')
    reference_species = Channel.value('')
    proteome_species = Channel.value('')
    // proteome discovery prefers a discovered reference species' own
    // proteins; a user-supplied reference has no known species
    proteome_prefer_taxid = Channel.value('')
    proteome_prefer_rank_taxid = Channel.value('')

    reference_fasta = null
    proteome_fasta   = null
    reference_display_name = null
    proteome_display_name = null

    if (params.reference) {
        reference_fasta = Channel.value(file(params.reference))
        reference_display_name = Channel.value(file(params.reference).name)
    }
    if (params.proteome) {
        proteome_fasta = Channel.value(file(params.proteome))
        proteome_display_name = Channel.value(file(params.proteome).name)
    }

    lineage_ch = null
    if (needsTargetTaxid()) {
        lineage_ch = RESOLVE_TAXONOMY(Channel.value(params.taxid), '')
        target_species = lineage_ch.map { f ->
            def row = f.readLines().collect { it.split('\t') }.find { it[0] == 'species' }
            row && row.size() > 2 ? row[2] : ''
        }
    }

    // ---- the reference, unless --reference gave a file ----
    reference_selection = null
    if (params.reference_accession) {
        summary = FETCH_ACCESSION_SUMMARY(Channel.value(params.reference_accession))
        reference_selection = DESCRIBE_REFERENCE_ACCESSION(Channel.value(params.reference_accession), summary).selection
    } else if (params.reference_taxid) {
        reference_lineage = RESOLVE_REFERENCE_TAXONOMY(Channel.value(params.reference_taxid), 'reference_')
        reference_selection = FIND_REFERENCE_ASSEMBLY(reference_lineage, params.max_rank, false, true).selection
    } else if (!params.reference) {
        reference_selection = FIND_REFERENCE_ASSEMBLY(lineage_ch, params.max_rank, params.exclude_target, false).selection
    }
    if (reference_selection) {
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
        proteome_prefer_rank_taxid = reference_selection.map { readSelection(it).taxid ?: '' }
    }

    // ---- the proteome, unless --proteome gave a file ----
    if (!params.proteome) {
        prot_selection = FIND_PROTEOME_ASSEMBLY(lineage_ch, params.max_rank, params.exclude_target,
                                                proteome_prefer_taxid, proteome_prefer_rank_taxid).selection
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

    // DSL2 can't call a process twice in one scope, so each per-genome step
    // runs once on a channel of both genomes tagged by role
    // ('target'/'reference'), split back apart with branch{} where needed.

    // ---- rename sequences, plus chrom sizes and gaps, in one read ----
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

    lookup_by_role = renamed.lookup.branch {
        target: it[0] == 'target'
        reference: it[0] == 'reference'
    }

    // ---- align the proteome against both original genomes, then rename
    // each GFF's seqids with RENAME_SEQUENCES' lookup ----
    // null params are passed to processes as '' (see min_identity below)
    def miniprot_m = params.miniprot_m != null ? params.miniprot_m : ''

    def chunk_gb = chunkGb()
    if (chunk_gb) {
        // align ~chunk_gb pieces separately, then merge -- close to, not
        // identical to, a whole-genome run; see miniprot_align.nf
        def target_chunk_bp = Math.round(chunk_gb * 1_000_000_000)

        split_ch = SPLIT_GENOME(genomes_in, target_chunk_bp)

        // one row per chunk, with its own length (for memory) and the whole
        // genome's (for -G, see MINIPROT_ALIGN_CHUNK)
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

        // groupTuple() without `size:` waits for the channel to close, so it
        // already gathers every chunk; the count check only guards against a
        // silently missing chunk task
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
    // null (chain.js's own defaults) doesn't pass through a `val` input
    // cleanly, so it becomes ''. A --min_identity 0 stays 0: the processes
    // test for '', not falsiness.
    def min_identity = params.min_identity != null ? params.min_identity : ''
    def min_block = params.min_block ?: ''
    // (source, species) for target, reference and proteome, as one tuple
    input_sources = target_display_name
        .combine(target_species)
        .combine(reference_display_name)
        .combine(reference_species)
        .combine(proteome_display_name)
        .combine(proteome_species)
    synteny = BUILD_SYNTENY(target_gff, reference_gff, target_chrom_sizes, reference_chrom_sizes,
                             proteome_fasta, min_identity, params.max_gap, min_block,
                             input_sources)

    PLOT_SYNTENY(
        synteny.hits,
        target_chrom_sizes, reference_chrom_sizes,
        min_identity, params.max_gap, min_block,
        synteny.stats,
        target_display_name, reference_display_name,
        target_gaps, reference_gaps,
        lookup_by_role.target, lookup_by_role.reference,
    )
}
