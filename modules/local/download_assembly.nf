// The reference genome only ever needs its genome FASTA (miniprot aligns
// the proteome against it; its own GFF3/protein are never used) and
// proteome only ever needs its protein FASTA. Because each accession's role only ever needs ONE
// distinct file type, a same-accession reference-genome+proteome pair each
// fetch different content (genome vs. protein), so no separate dedup branch
// is needed: there is no redundant download to avoid.
//
// Both download dehydrated: the zip holds only a manifest, and rehydrate
// then fetches the sequence files as separate requests. NCBI drops
// downloads now and then (HTTP/2 stream resets, "stream error: stream ID 5;
// INTERNAL_ERROR; received from peer"), the dehydrated zip as well as the
// files rehydrate fetches, so each step is retried in the script itself.
// Not with errorStrategy: that would resubmit the whole task (a new SLURM
// job), and conf/slurm.config's own errorStrategy only retries OOM kills.

// `retry CMD...`: runs CMD until it succeeds, at most 4 times, waiting 1, 2
// then 3 minutes in between; the last failure's exit status is the task's
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
