// Same ladder-walk as find_reference_assembly.nf, filtered for annotation
// availability instead of assembly level -- a scaffold-level annotated
// assembly is fine here since only its proteins are used. Unlike that
// query, this one has no --assembly-level filter of its own, so
// find_closest_assembly.py's chromosome-level gate on same-species
// candidates (see that script) is the only thing enforcing it here.

process QUERY_GENOME_LADDER_PROTEOME {
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
            datasets summary genome taxon "\$taxid" --annotated --as-json-lines --limit 200 > "\${r}.jsonl" || true
        fi
        if [ "\$r" == "${max_rank}" ]; then
            break
        fi
    done
    """
}

process SELECT_PROTEOME_ASSEMBLY {
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'proteome_*'

    input:
    path lineage
    path jsonl_files
    val max_rank
    val exclude_target  // true: drop every same-species candidate outright (see find_closest_assembly.py)
    val prefer_taxid    // the reference genome's species taxid ('' if user-supplied): its own
                        // annotation wins over the quality ranking when it has one

    output:
    path 'proteome_selection.json', emit: selection
    path 'proteome_candidates.tsv'
    path 'proteome_search_log.tsv'

    script:
    def excludeFlag = exclude_target ? '--exclude_target' : ''
    def preferFlag = prefer_taxid ? "--prefer_taxid ${prefer_taxid}" : ''
    """
    find_closest_assembly.py --lineage ${lineage} --max_rank ${max_rank} --jsonl_dir . --outprefix proteome ${excludeFlag} ${preferFlag}
    """
}

workflow FIND_PROTEOME_ASSEMBLY {
    take:
    lineage
    max_rank
    exclude_target
    prefer_taxid     // species taxid whose own annotation to prefer ('' for none)

    main:
    ladder = QUERY_GENOME_LADDER_PROTEOME(lineage, max_rank)
    sel    = SELECT_PROTEOME_ASSEMBLY(lineage, ladder.jsonl, max_rank, exclude_target, prefer_taxid)

    emit:
    selection = sel.selection
}
