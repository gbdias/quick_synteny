// Climbs the lineage ladder (species -> params.max_rank) looking for a
// chromosome/complete-level assembly. Split into a ladder-query process
// (datasets container, one query per rank) and a selection process (plain
// Python container) for the same reason as resolve_taxonomy.nf: the datasets
// CLI image has no Python to do the JSON selection logic in-process.
//
// --limit 200: for a heavily-sampled species (e.g. E. coli, human) this
// query can otherwise return tens of thousands of records, which both takes
// a long time and can OOM-kill the downstream Python selection step. 200
// candidates is more than enough to pick a good one -- this tool doesn't
// need the single best assembly out of every one ever submitted.
//
// --assembly-level chromosome,complete here is belt-and-suspenders with
// find_closest_assembly.py's own --require_chromosome_level flag (always
// passed below): a low-quality genome makes a poor synteny comparison
// regardless of species, so this is enforced explicitly in the selection
// script itself, not left to only ever be true as an incidental side effect
// of this query's own filter. Proteome discovery passes neither -- it only
// needs a candidate's protein sequences, so a scaffold-level annotated
// assembly is fine there (see find_proteome_assembly.nf).

process QUERY_GENOME_LADDER_COMPARISON {
    tag "max_rank=${max_rank}"
    label 'process_low'
    container 'quay.io/staphb/ncbi-datasets:18.35.0'

    input:
    path lineage
    val max_rank

    output:
    path '*.jsonl', emit: jsonl

    script:
    """
    RANKS=(species genus family order class phylum)
    for r in "\${RANKS[@]}"; do
        taxid=\$(awk -F'\\t' -v rank="\$r" '\$1==rank{print \$2}' ${lineage})
        if [ -n "\$taxid" ]; then
            datasets summary genome taxon "\$taxid" --assembly-level chromosome,complete --as-json-lines --limit 200 > "\${r}.jsonl" || true
        fi
        if [ "\$r" == "${max_rank}" ]; then
            break
        fi
    done
    """
}

process SELECT_COMPARISON_ASSEMBLY {
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'comparison_*'

    input:
    path lineage
    path jsonl_files
    val max_rank
    val exclude_target  // true: drop every same-species candidate outright (see find_closest_assembly.py)

    output:
    path 'comparison_selection.json', emit: selection
    path 'comparison_candidates.tsv'
    path 'comparison_search_log.tsv'

    script:
    def excludeFlag = exclude_target ? '--exclude_target' : ''
    """
    find_closest_assembly.py --lineage ${lineage} --max_rank ${max_rank} --jsonl_dir . --outprefix comparison \\
        --require_chromosome_level ${excludeFlag}
    """
}

workflow FIND_COMPARISON_ASSEMBLY {
    take:
    lineage
    max_rank
    exclude_target

    main:
    ladder = QUERY_GENOME_LADDER_COMPARISON(lineage, max_rank)
    sel    = SELECT_COMPARISON_ASSEMBLY(lineage, ladder.jsonl, max_rank, exclude_target)

    emit:
    selection = sel.selection
}
