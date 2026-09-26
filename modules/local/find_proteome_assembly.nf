// Same ladder-walk as find_reference_assembly.nf, filtered for annotation
// availability instead of assembly level -- a scaffold-level annotated
// assembly is fine here since only its proteins are used. Unlike that
// query, this one has no --assembly-level filter of its own, so
// find_closest_assembly.py's chromosome-level gate on same-species
// candidates (see that script) is the only thing enforcing it here. The
// preferred (reference) species is also queried on its own, since it may
// join the lineage above the taxon where the climb stops.

process QUERY_GENOME_LADDER_PROTEOME {
    tag "max_rank=${max_rank}"
    label 'process_low'
    label 'env_ncbi_datasets'

    input:
    path lineage
    val max_rank
    val prefer_taxid

    output:
    path 'ladder', emit: ladder

    script:
    """
    query_genome_ladder.sh ${lineage} ${max_rank} '${prefer_taxid}' --annotated
    """
}

process SELECT_PROTEOME_ASSEMBLY {
    label 'process_low'
    label 'env_python'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'proteome_*'

    input:
    path lineage
    path ladder
    val max_rank
    val exclude_target  // true: drop every same-species candidate outright (see find_closest_assembly.py)
    val prefer_taxid    // the reference genome's species taxid ('' if user-supplied): its own
                        // annotation wins over the quality ranking when it has one
    val prefer_rank_taxid // the lineage taxon that species joins at (the reference's selection taxid)

    output:
    path 'proteome_selection.json', emit: selection
    path 'proteome_candidates.tsv'
    path 'proteome_search_log.tsv'

    script:
    def excludeFlag = exclude_target ? '--exclude_target' : ''
    def preferFlag = prefer_taxid ? "--prefer_taxid ${prefer_taxid} --prefer_rank_taxid ${prefer_rank_taxid}" : ''
    """
    find_closest_assembly.py --lineage ${lineage} --max_rank ${max_rank} --ladder_dir ${ladder} --outprefix proteome ${excludeFlag} ${preferFlag}
    """
}

workflow FIND_PROTEOME_ASSEMBLY {
    take:
    lineage
    max_rank
    exclude_target
    prefer_taxid      // species taxid whose own annotation to prefer ('' for none)
    prefer_rank_taxid // lineage taxid that species joins the target's lineage at

    main:
    ladder = QUERY_GENOME_LADDER_PROTEOME(lineage, max_rank, prefer_taxid)
    sel    = SELECT_PROTEOME_ASSEMBLY(lineage, ladder.ladder, max_rank, exclude_target, prefer_taxid, prefer_rank_taxid)

    emit:
    selection = sel.selection
}
