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
//
// `chosen` is --reference_taxid: the lineage is then the chosen taxon's,
// and only its first row -- that taxon itself -- is searched, for its best
// chromosome-level assembly, instead of climbing outward from the target.

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
    def excludeFlag = exclude_target && !chosen ? '--exclude_target' : ''
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

// --reference_accession: one exact assembly, no search. Its summary from
// NCBI becomes a reference_selection.json of the same shape the search
// writes, so everything downstream (the download, the proteome's preferred
// species, the page's labels) treats it alike. Any assembly level: it's the
// user's pick. Fetched and described in two processes, like the search
// above: the datasets image has no Python.
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
