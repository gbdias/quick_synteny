nextflow.enable.dsl = 2

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
include { PYGENOMEVIZ_PLOT }        from './modules/local/pygenomeviz_plot.nf'

// chunk size used when --miniprot_chunk_gb is given without a value
def DEFAULT_CHUNK_GB() {
    return 1.0
}

// --miniprot_chunk_gb -> chunk size in Gb, or null for the whole-genome path.
// A bare flag arrives as boolean true; false/null/unset mean whole-genome.
def chunkGb() {
    def v = params.miniprot_chunk_gb
    if (v == null || v == false) return null
    if (v == true) return DEFAULT_CHUNK_GB()
    return v as double
}

// The reference genome comes from one of four places: --reference (a
// file), --reference_accession (one exact NCBI assembly), --reference_taxid
// (the best assembly of a taxon you choose), or else discovery outward
// from --taxid. The proteome is --proteome, or else discovered from
// --taxid's relatives. So --taxid -- the target's organism -- is needed
// whenever something is discovered from it.
def referenceGiven() {
    return params.reference || params.reference_accession || params.reference_taxid
}

def needsTargetTaxid() {
    return !referenceGiven() || !params.proteome
}

// What nextflow_schema.json can't say: which parameters go together.
// Everything else about each parameter -- types, ranges, allowed values,
// files that must exist, unknown parameters -- is checked by
// validateParameters() against the schema, which is also where --help
// comes from.
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
    // --help/--helpFull are nf-schema's, built from nextflow_schema.json (see
    // nextflow.config). It prints them as the run starts, but doesn't stop
    // this workflow: without the exit, validation would then fail on the
    // missing --assembly
    if (params.help || params.containsKey('helpFull') && params.helpFull) {
        exit 0
    }
    validateParameters()
    validateParams()

    log.info "quick_synteny: assembly=${params.assembly} taxid=${params.taxid ?: '(skipped, using overrides)'} outdir=${params.outdir}"


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
    proteome_prefer_rank_taxid = Channel.value('')

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

    lineage_ch = null
    if (needsTargetTaxid()) {
        lineage_ch = RESOLVE_TAXONOMY(Channel.value(params.taxid), '')
        target_species = lineage_ch.map { f ->
            def row = f.readLines().collect { it.split('\t') }.find { it[0] == 'species' }
            row && row.size() > 2 ? row[2] : ''
        }
    }

    // ---- the reference: an exact assembly, the best of a chosen taxon, or
    // discovered outward from the target (--reference, a file, is above) ----
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

    // ---- the proteome: discovered from the target's relatives, preferring
    // the reference species' own proteins when NCBI has them annotated ----
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

    lookup_by_role = renamed.lookup.branch {
        target: it[0] == 'target'
        reference: it[0] == 'reference'
    }

    // ---- align proteome against both ORIGINAL genomes (also doubles as
    // the raw synteny/homeolog anchor source -- see build_synteny.nf), then
    // rename each resulting GFF's seqid column using the lookup from the
    // matching RENAME_SEQUENCES call (join()'d by role name) ----
    // params.miniprot_m defaults to null (miniprot's own default); normalized
    // to '' here for the same reason as min_identity below -- see that comment
    def miniprot_m = params.miniprot_m != null ? params.miniprot_m : ''

    def chunk_gb = chunkGb()
    if (chunk_gb) {
        // Chunked path: split each genome into ~chunk_gb-sized pieces, align
        // each in its own (smaller-index, lower-RAM) miniprot task, then
        // merge back into one GFF -- close to, but not identical to, a
        // whole-genome run; see modules/local/miniprot_align.nf.
        def target_chunk_bp = Math.round(chunk_gb * 1_000_000_000)

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
    // params.min_identity defaults to null (chain.js's own 0.7) and
    // params.min_block to null (auto-tune); Groovy's null is not safe to pass
    // through a process `val` input the same way every call site expects
    // (falsy-but-interpolatable), so they're normalized to '' here once. A
    // --min_identity 0 stays 0 (the call sites test for '', not falsiness).
    def min_identity = params.min_identity != null ? params.min_identity : ''
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
        lookup_by_role.target, lookup_by_role.reference,
    )
}
