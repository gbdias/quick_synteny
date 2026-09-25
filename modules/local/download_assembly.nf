// The reference genome only ever needs its genome FASTA (miniprot aligns
// the proteome against it; its own GFF3/protein are never used) and
// proteome only ever needs its protein FASTA. Because each accession's role only ever needs ONE
// distinct file type, a same-accession reference-genome+proteome pair each
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
    # dehydrated: the zip holds only a manifest; rehydrate then fetches the
    # sequence files as separate requests. NCBI drops long single-stream
    # downloads now and then (HTTP/2 stream resets), and a rerun of
    # rehydrate only fetches what is still missing, so retry that step
    datasets download genome accession ${accession} --include genome --dehydrated --no-progressbar --filename dl.zip
    unzip -q dl.zip
    n=0
    until datasets rehydrate --directory . --no-progressbar; do
        n=\$((n + 1))
        if [ \$n -ge 3 ]; then exit 1; fi
        sleep \$((n * 60))
    done
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
    # dehydrated: the zip holds only a manifest; rehydrate then fetches the
    # sequence files as separate requests. NCBI drops long single-stream
    # downloads now and then (HTTP/2 stream resets), and a rerun of
    # rehydrate only fetches what is still missing, so retry that step
    datasets download genome accession ${accession} --include protein --dehydrated --no-progressbar --filename dl.zip
    unzip -q dl.zip
    n=0
    until datasets rehydrate --directory . --no-progressbar; do
        n=\$((n + 1))
        if [ \$n -ge 3 ]; then exit 1; fi
        sleep \$((n * 60))
    done
    mv ncbi_dataset/data/${accession}/protein.faa protein.faa
    """
}
