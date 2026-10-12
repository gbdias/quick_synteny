// The reference needs only its genome FASTA and the proteome only its
// protein FASTA, so a reference and proteome from the same accession never
// download the same file.
//
// Both download dehydrated, then rehydrate fetches the sequence files. NCBI
// drops downloads now and then (HTTP/2 stream resets), so each step is
// retried in the script itself, not with errorStrategy: that would resubmit
// the whole task (a new SLURM job), and conf/slurm.config's errorStrategy
// only retries OOM kills.

// `retry CMD...`: at most 4 tries, waiting 1, 2 then 3 minutes in between
def retryFunction() {
    return '''
    retry() {
        n=0
        until "$@"; do
            status=$?
            n=$((n + 1))
            if [ $n -gt 3 ]; then return $status; fi
            echo "quick_synteny: '$*' failed (exit $status), retry $n of 3 in $((n * 60)) s" >&2
            sleep $((n * 60))
        done
    }
    '''.trim()   // keeps its inner indentation, so it lines up with the script around it
}

process DOWNLOAD_GENOME {
    tag "${accession}"
    label 'process_low'
    label 'env_ncbi_datasets'

    input:
    val accession

    output:
    path 'genome.fna', emit: fasta

    script:
    """
    ${retryFunction()}
    # a failed download can leave a partial zip behind
    download() {
        rm -f dl.zip
        datasets download genome accession ${accession} --include genome --dehydrated --no-progressbar --filename dl.zip
    }
    retry download
    unzip -q dl.zip
    # a rerun of rehydrate only fetches what is still missing
    retry datasets rehydrate --directory . --no-progressbar
    mv ncbi_dataset/data/${accession}/*_genomic.fna genome.fna
    """
}

process DOWNLOAD_PROTEIN {
    tag "${accession}"
    label 'process_low'
    label 'env_ncbi_datasets'

    input:
    val accession

    output:
    path 'protein.faa', emit: fasta

    script:
    """
    ${retryFunction()}
    # a failed download can leave a partial zip behind
    download() {
        rm -f dl.zip
        datasets download genome accession ${accession} --include protein --dehydrated --no-progressbar --filename dl.zip
    }
    retry download
    unzip -q dl.zip
    # a rerun of rehydrate only fetches what is still missing
    retry datasets rehydrate --directory . --no-progressbar
    mv ncbi_dataset/data/${accession}/protein.faa protein.faa
    """
}
