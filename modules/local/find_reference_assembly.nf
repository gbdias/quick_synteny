// Climbs the lineage ladder (species -> params.max_rank, through every
// intermediate taxon NCBI has -- see parse_lineage.py) looking for a
// chromosome/complete-level assembly. Split into a ladder-query process
// (datasets container, one query per taxon) and a selection process (plain
// Python container) for the same reason as resolve_taxonomy.nf: the datasets
// CLI image has no Python to do the JSON selection logic in-process.
//
// query_genome_ladder.sh stops at the first taxon with another species'
// assembly and caps each query's size -- see that script for how the cap
// avoids keeping an arbitrary subset of a large taxon.
//
// --assembly-level chromosome,complete here is belt-and-suspenders with
// find_closest_assembly.py's own --require_chromosome_level flag (always
// passed below): a low-quality genome makes a poor synteny comparison
// regardless of species, so this is enforced explicitly in the selection
// script itself, not left to only ever be true as an incidental side effect
// of this query's own filter. Proteome discovery passes neither -- it only
// needs a candidate's protein sequences, so a scaffold-level annotated
// assembly is fine there (see find_proteome_assembly.nf).

process QUERY_GENOME_LADDER_REFERENCE {
    tag "max_rank=${max_rank}"
    label 'process_low'
    label 'env_ncbi_datasets'

    input:
    path lineage
    val max_rank

    output:
    path 'ladder', emit: ladder

    script:
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

    output:
    path 'reference_selection.json', emit: selection
    path 'reference_candidates.tsv'
    path 'reference_search_log.tsv'

    script:
    def excludeFlag = exclude_target ? '--exclude_target' : ''
    """
    find_closest_assembly.py --lineage ${lineage} --max_rank ${max_rank} --ladder_dir ${ladder} --outprefix reference \\
        --require_chromosome_level ${excludeFlag}
    """
}

workflow FIND_REFERENCE_ASSEMBLY {
    take:
    lineage
    max_rank
    exclude_target

    main:
    ladder = QUERY_GENOME_LADDER_REFERENCE(lineage, max_rank)
    sel    = SELECT_REFERENCE_ASSEMBLY(lineage, ladder.ladder, max_rank, exclude_target)

    emit:
    selection = sel.selection
}
