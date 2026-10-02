#!/usr/bin/env bash
# Parameter validation: nextflow_schema.json (nf-schema's validateParameters)
# plus main.nf's validateParams(), and nf-schema's --help. Each case runs
# `nextflow run main.nf` with a set of parameters and checks the exit status
# and the message. Valid cases use -preview, which runs main.nf (and so the
# validation) without running any process: no data, Docker or conda needed.
#
#   tests/params.sh            # uses `nextflow` on PATH
#   NEXTFLOW=/path/to/nextflow tests/params.sh

set -u
NEXTFLOW=${NEXTFLOW:-nextflow}
MAIN=$(cd "$(dirname "$0")/.." && pwd)/main.nf
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
printf '>chr1\nACGTACGTACGT\n' > genome.fa
printf '>p1\nMKVLAAGIV\n' > proteome.faa
export NXF_DISABLE_CHECK_LATEST=true

pass=0
fail=0

run() {
    "$NEXTFLOW" run "$MAIN" "$@" -ansi-log false > out.txt 2>&1
    status=$?
    # nf-schema colours its list even without a terminal
    sed 's/\x1b\[[0-9;]*m//g' out.txt > clean.txt
}

report() {
    if [ "$1" = ok ]; then
        pass=$((pass + 1)); echo "ok    $2"
    else
        fail=$((fail + 1)); echo "FAIL  $2"; sed 's/^/      | /' clean.txt | tail -15
    fi
}

# fails CASE PATTERN ARGS...: exits non-zero, and the output matches PATTERN
fails() {
    local name=$1 pattern=$2; shift 2
    run "$@"
    if [ $status -ne 0 ] && grep -qE -- "$pattern" clean.txt; then report ok "$name"; else report fail "$name (exit $status)"; fi
}

# passes CASE ARGS...: validation passes (run with -preview)
passes() {
    local name=$1; shift
    run "$@" -preview
    if [ $status -eq 0 ]; then report ok "$name"; else report fail "$name (exit $status)"; fi
}

echo "== --help"
run --help
if [ $status -eq 0 ] && grep -q -- '--max_rank' clean.txt && grep -q 'taxonomy-guided synteny plotting' clean.txt \
    && ! grep -q ERROR clean.txt; then report ok "--help lists the parameters, exits 0"; else report fail "--help"; fi
if ! grep -qE -- '--(min_identity|max_gap|min_block)\b' clean.txt; then report ok "--help hides the chaining parameters"; else report fail "--help hides the chaining parameters"; fi
run --help max_rank
if [ $status -eq 0 ] && grep -q 'help_text' clean.txt; then report ok "--help max_rank shows its details"; else report fail "--help max_rank"; fi

echo "== invalid parameters"
fails "no --assembly"                    'Missing required parameter\(s\): assembly'
fails "--assembly that doesn't exist"    "--assembly \(nope.fa\): the file or directory 'nope.fa' does not exist" --assembly nope.fa --taxid 4932
fails "--taxid 000000 (arrives as 0)"    '--taxid .*must be an NCBI taxid' --assembly genome.fa --taxid 000000
fails "--taxid abc"                      '--taxid \(abc\).*must be an NCBI taxid' --assembly genome.fa --taxid abc
fails "bare --taxid"                     '--taxid \(true\).*must be an NCBI taxid' --assembly genome.fa --taxid
fails "no --taxid and no overrides"      '--taxid is required unless --proteome is given along with a reference genome' --assembly genome.fa
fails "--reference without --proteome"   '--taxid is required unless --proteome is given along with a reference genome' --assembly genome.fa --reference genome.fa
fails "--reference_taxid, no --proteome"  '--taxid is required unless --proteome is given' --assembly genome.fa --reference_taxid 4932
fails "--reference_taxid abc"            '--reference_taxid \(abc\).*must be an NCBI taxid' --assembly genome.fa --taxid 4932 --reference_taxid abc
fails "--reference_accession not one"    '--reference_accession \(GCA_1\).*must be an NCBI assembly accession' --assembly genome.fa --taxid 4932 --reference_accession GCA_1
fails "two references at once"           'give only one of --reference, --reference_accession and --reference_taxid, got --reference and --reference_taxid' --assembly genome.fa --taxid 4932 --reference genome.fa --reference_taxid 4932
fails "--max_rank kingdom"               '--max_rank \(kingdom\): Expected any of \[species, genus, family, order, class, phylum\]' --assembly genome.fa --taxid 4932 --max_rank kingdom
fails "--min_identity 1.5"               '--min_identity \(1.5\): 1.5 is greater than 1' --assembly genome.fa --taxid 4932 --min_identity 1.5
fails "--min_block 1"                    '--min_block \(1\): 1 is less than 2' --assembly genome.fa --taxid 4932 --min_block 1
fails "--miniprot_chunk_gb 0"            '--miniprot_chunk_gb .*must be a positive number of Gb' --assembly genome.fa --taxid 4932 --miniprot_chunk_gb 0
fails "--miniprot_chunk_gb -1"           '--miniprot_chunk_gb \(-1\).*must be a positive number of Gb' --assembly genome.fa --taxid 4932 --miniprot_chunk_gb -1
fails "misspelled --max-rank"            '--maxRank: genus' --assembly genome.fa --taxid 4932 --max-rank genus
fails "unknown --colour"                 '--colour: red' --assembly genome.fa --taxid 4932 --colour red

echo "== valid parameters"
passes "--taxid"                                --assembly genome.fa --taxid 4932
passes "--reference and --proteome, no taxid"   --assembly genome.fa --reference genome.fa --proteome proteome.faa
passes "--reference, proteome from NCBI"        --assembly genome.fa --reference genome.fa --taxid 4932
passes "--reference_taxid"                      --assembly genome.fa --taxid 4932 --reference_taxid 27291
passes "--reference_accession"                  --assembly genome.fa --taxid 4932 --reference_accession GCA_056824455.1
passes "a chosen reference and --proteome, no taxid" --assembly genome.fa --reference_accession GCF_000146045.2 --proteome proteome.faa
passes "every discovery and chaining option"    --assembly genome.fa --taxid 4932 --max_rank genus --exclude_target \
    --min_seq_size 0 --min_asm_gap 10 --min_identity 0.5 --max_gap 10 --min_block 3
passes "bare --miniprot_chunk_gb"               --assembly genome.fa --taxid 4932 --miniprot_chunk_gb
passes "--miniprot_chunk_gb 0.5"                --assembly genome.fa --taxid 4932 --miniprot_chunk_gb 0.5
passes "--miniprot_m 0"                         --assembly genome.fa --taxid 4932 --miniprot_m 0

echo
echo "$pass passed, $fail failed"
[ $fail -eq 0 ]
