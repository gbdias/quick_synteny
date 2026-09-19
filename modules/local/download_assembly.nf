// The comparison genome only ever needs its genome FASTA (miniprot aligns
// the proteome against it; its own GFF3/protein are never used) and
// proteome only ever needs its protein FASTA -- see the plan's "download
// scope trimmed" note. Because each accession's role only ever needs ONE
// distinct file type, a same-accession comparison-genome+proteome pair each
// fetch different content (genome vs. protein), so no separate dedup branch
// is needed: there is no redundant download to avoid.

process DOWNLOAD_GENOME {
    tag "${accession}"
    label 'process_low'
    container 'quay.io/staphb/ncbi-datasets:18.35.0'

    input:
    val accession

    output:
    path 'genome.fna', emit: fasta

    script:
    """
    datasets download genome accession ${accession} --include genome --no-progressbar --filename dl.zip
    unzip -q dl.zip
    mv ncbi_dataset/data/${accession}/*_genomic.fna genome.fna
    """
}

process DOWNLOAD_PROTEIN {
    tag "${accession}"
    label 'process_low'
    container 'quay.io/staphb/ncbi-datasets:18.35.0'

    input:
    val accession

    output:
    path 'protein.faa', emit: fasta

    script:
    """
    datasets download genome accession ${accession} --include protein --no-progressbar --filename dl.zip
    unzip -q dl.zip
    mv ncbi_dataset/data/${accession}/protein.faa protein.faa
    """
}
