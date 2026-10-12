// Climbs the lineage ladder (species -> params.max_rank, every taxon in
// between -- see parse_lineage.py) for a chromosome/complete-level assembly.
// Query and selection are separate processes because the datasets image
// has no Python.
//
// --assembly-level here and --require_chromosome_level in the selection are
// deliberately redundant: the selection script enforces it itself rather
// than relying on the query's filter.
//
// `chosen` is --reference_taxid: the lineage is the chosen taxon's, and
// only its first row -- that taxon -- is searched.

process QUERY_GENOME_LADDER_REFERENCE {
    tag "${chosen ? 'chosen taxon' : "max_rank=${max_rank}"}"
    label 'process_low'
    label 'env_ncbi_datasets'

    input:
    path lineage
    val max_rank
    val chosen

    output:
    path 'ladder', emit: ladder

    script:
    // a chosen taxon: its own row alone, so the ladder stops there (phylum
    // only keeps the one-row ladder from being cut short by a major rank)
    if (chosen)
        """
        head -n 1 ${lineage} > chosen_lineage.tsv
        query_genome_ladder.sh chosen_lineage.tsv phylum '' --assembly-level chromosome,complete
        """
    else
        """
        query_genome_ladder.sh ${lineage} ${max_rank} '' --assembly-level chromosome,complete
        """
}

process SELECT_REFERENCE_ASSEMBLY {
    label 'process_low'
    label 'env_python'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'reference_*'

    input:
    path lineage
    path ladder
    val max_rank
    val exclude_target  // true: drop every same-species candidate outright (see find_closest_assembly.py)
    val chosen          // --reference_taxid: the lineage is the chosen taxon's (see above)

    output:
    path 'reference_selection.json', emit: selection
    path 'reference_candidates.tsv'
    path 'reference_search_log.tsv'

    script:
    def excludeFlag = exclude_target ? '--exclude_target' : ''
    if (chosen)
        """
        head -n 1 ${lineage} > chosen_lineage.tsv
        find_closest_assembly.py --lineage chosen_lineage.tsv --max_rank phylum --ladder_dir ${ladder} --outprefix reference \\
            --require_chromosome_level --chosen
        """
    else
        """
        find_closest_assembly.py --lineage ${lineage} --max_rank ${max_rank} --ladder_dir ${ladder} --outprefix reference \\
            --require_chromosome_level ${excludeFlag}
        """
}

// --reference_accession: no search. Its NCBI summary becomes a
// reference_selection.json of the same shape the search writes, so
// everything downstream treats it alike. Any assembly level: it's the
// user's pick.
process FETCH_ACCESSION_SUMMARY {
    tag "${accession}"
    label 'process_low'
    label 'env_ncbi_datasets'

    input:
    val accession

    output:
    path 'summary.jsonl'

    script:
    """
    datasets summary genome accession ${accession} --as-json-lines > summary.jsonl
    """
}

process DESCRIBE_REFERENCE_ACCESSION {
    tag "${accession}"
    label 'process_low'
    label 'env_python'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'reference_*'

    input:
    val accession
    path summary

    output:
    path 'reference_selection.json', emit: selection

    script:
    """
    describe_assembly.py --summary ${summary} --accession ${accession} --out reference_selection.json
    """
}

workflow FIND_REFERENCE_ASSEMBLY {
    take:
    lineage
    max_rank
    exclude_target
    chosen          // true: lineage is --reference_taxid's, searched at that taxon only

    main:
    ladder = QUERY_GENOME_LADDER_REFERENCE(lineage, max_rank, chosen)
    sel    = SELECT_REFERENCE_ASSEMBLY(lineage, ladder.ladder, max_rank, exclude_target, chosen)

    emit:
    selection = sel.selection
}
