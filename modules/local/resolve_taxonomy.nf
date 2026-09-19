// Ranks are queried outward from species; see bin/parse_lineage.py for the
// exact ladder order (species,genus,family,order,class,phylum).
//
// Split into two processes because the ncbi-datasets-cli biocontainers image
// (quay.io/staphb/ncbi-datasets, used in place of quay.io/biocontainers/
// ncbi-datasets-cli -- that image's CA bundle is broken, see README) has no
// Python: fetch the raw taxonomy JSON in the datasets container, parse it in
// a plain Python container.

process FETCH_TAXONOMY_JSON {
    tag "${taxid}"
    label 'process_low'
    container 'quay.io/staphb/ncbi-datasets:18.35.0'

    input:
    val taxid

    output:
    path 'taxonomy.json'

    script:
    """
    datasets summary taxonomy taxon ${taxid} > taxonomy.json
    """
}

process PARSE_LINEAGE {
    tag "${taxonomy_json}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.13.7'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path taxonomy_json

    output:
    path 'lineage.tsv'

    script:
    """
    parse_lineage.py ${taxonomy_json} > lineage.tsv
    """
}

workflow RESOLVE_TAXONOMY {
    take:
    taxid

    main:
    taxonomy_json = FETCH_TAXONOMY_JSON(taxid)
    lineage       = PARSE_LINEAGE(taxonomy_json)

    emit:
    lineage
}
