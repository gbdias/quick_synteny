// Ranks are queried outward from species; see bin/parse_lineage.py for the
// exact ladder order (species outward to phylum, every intermediate taxon
// included -- hence --parents, which adds each ancestor's own report).
//
// Split into two processes because the datasets image (envs/ncbi_datasets.yml)
// has no Python: fetch the raw taxonomy JSON in the datasets container, parse
// it in the Python one (envs/python.yml).

process FETCH_TAXONOMY_JSON {
    tag "${taxid}"
    label 'process_low'
    label 'env_ncbi_datasets'

    input:
    val taxid

    output:
    path 'taxonomy.json'

    script:
    """
    datasets summary taxonomy taxon ${taxid} --parents > taxonomy.json
    """
}

process PARSE_LINEAGE {
    tag "${taxonomy_json}"
    label 'process_low'
    label 'env_python'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path taxonomy_json
    val prefix      // '' for the target's lineage, 'reference_' for --reference_taxid's

    output:
    path "${prefix}lineage.tsv"

    script:
    """
    parse_lineage.py ${taxonomy_json} > ${prefix}lineage.tsv
    """
}

// called twice with --reference_taxid (main.nf includes it a second time as
// RESOLVE_REFERENCE_TAXONOMY): `prefix` keeps the two published lineages apart
workflow RESOLVE_TAXONOMY {
    take:
    taxid
    prefix

    main:
    taxonomy_json = FETCH_TAXONOMY_JSON(taxid)
    lineage       = PARSE_LINEAGE(taxonomy_json, prefix)

    emit:
    lineage
}
